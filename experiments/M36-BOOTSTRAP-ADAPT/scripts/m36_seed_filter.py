"""S2 seed construction from the miniprot alignments (protocol 22 §2): keep complete protein-to-genome models with
  coverage >= 0.95 (aligned query span / protein length), Identity >= 0.70, no frameshift and no in-frame stop (miniprot fs/st tags),
  a stop codon closing the aligned CDS (miniprot `stop_codon` feature; already inside the terminal exon in miniprot 0.18, otherwise appended),
  and a chain that passes the frozen `manifest.chain_checks` (ATG start, terminal stop, no internal stop, length % 3, GT/GC-AG introns, ACGT only,
  terminal exons >= 3 bp) with <= 128 exons;  overlapping candidates on the same strand: keep the highest identity (ties: coverage, CDS bp, id).
Only label-set A contigs are present in the miniprot input; the script re-asserts that no evaluation-set B contig appears.
Output (data/m36/seeds/): seeds_S2_<species>.gff3, seeds_S2_<species>.json (seed rows for m36_seed_labels), seed_stats_S2_<species>.json
usage: m36_seed_filter.py <species> [raw.gff.gz]
"""
import gzip
import json
import sys
import time
from collections import Counter, defaultdict

import numpy as np

import m36_common as W
import manifest as MF                       # frozen: chain_checks, read_allowed


def attrs(s):
    return dict(x.split('=', 1) for x in s.rstrip(';').split(';') if '=' in x)


def parse_miniprot_gff(path):
    """yield dict(paf=..., mrna=..., cds=[(start1, end1, phase)], stop=[(start1,end1)]) per alignment"""
    cur = None
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('##PAF'):
                if cur is not None:
                    yield cur
                p = line.rstrip('\n').split('\t')
                tags = {t[:2]: t[5:] for t in p[13:]}
                cur = dict(query=p[1], qlen=int(p[2]), qs=int(p[3]), qe=int(p[4]), ctg=p[6], ts=int(p[8]), te=int(p[9]),
                           fs=int(tags.get('fs', 0)), st=int(tags.get('st', 0)), mrna=None, cds=[], stop=[])
            elif line[0] == '#':
                continue
            else:
                p = line.rstrip('\n').split('\t')
                if len(p) < 9 or cur is None:
                    continue
                if p[2] == 'mRNA':
                    a = attrs(p[8])
                    cur['mrna'] = dict(ctg=p[0], start=int(p[3]), end=int(p[4]), strand=p[6], id=a.get('ID'), identity=float(a.get('Identity', 'nan')),
                                       positive=float(a.get('Positive', 'nan')), frameshift=int(a.get('Frameshift', cur['fs'])), stopcodon=int(a.get('StopCodon', cur['st'])))
                elif p[2] == 'CDS':
                    cur['cds'].append((int(p[3]), int(p[4]), p[7]))
                elif p[2] == 'stop_codon':
                    cur['stop'].append((int(p[3]), int(p[4])))
    if cur is not None:
        yield cur


def chain_from_alignment(al):
    """-> (CDS [(a0, b0, phase)] ascending, has_stop) with the stop codon merged into the terminal exon when contiguous."""
    strand = al['mrna']['strand']
    ex = sorted((s - 1, e) for s, e, _ in al['cds'])
    has_stop = False
    if al['stop']:
        stops = sorted((s - 1, e) for s, e in al['stop'])
        if len(stops) == 1 and stops[0][1] - stops[0][0] == 3:
            s0, e0 = stops[0]
            # miniprot >= 0.12 already includes the stop codon in the terminal CDS exon (stop interval inside the exon); older layouts put it just outside
            if strand == '+' and (e0 == ex[-1][1] and s0 >= ex[-1][0]):
                has_stop = True
            elif strand == '+' and s0 == ex[-1][1]:
                ex[-1] = (ex[-1][0], e0); has_stop = True
            elif strand == '-' and (s0 == ex[0][0] and e0 <= ex[0][1]):
                has_stop = True
            elif strand == '-' and e0 == ex[0][0]:
                ex[0] = (s0, ex[0][1]); has_stop = True
        if not has_stop:                                           # split stop codon or odd feature: keep the pieces, chain_checks will reject
            ex = sorted(set(ex) | set(stops))
    order = ex if strand == '+' else ex[::-1]
    phases, off = {}, 0
    for a, b in order:
        phases[(a, b)] = str((3 - off % 3) % 3)
        off += b - a
    return [(a, b, phases[(a, b)]) for a, b in ex], has_stop


def main():
    W.pin_and_nice()
    sp = sys.argv[1]
    raw = sys.argv[2] if len(sys.argv) > 2 else W.DATA / 'miniprot_raw' / f'{sp}.gff.gz'
    t0 = time.time()
    lens = W.a_lengths(sp)
    seqs = MF.read_allowed(W.M1 / sp / 'genome.fa', lens)
    split_of = {c: ('val' if c in W.TARGETS[sp]['val'] else 'train') for c in lens}
    funnel = Counter()
    fail_reasons = Counter()
    cands = []
    for al in parse_miniprot_gff(raw):
        funnel['alignments'] += 1
        W.assert_not_B(sp, [al['ctg']])
        if al['mrna'] is None or not al['cds'] or al['ctg'] not in lens:
            funnel['no_mrna_or_unknown_contig'] += 1
            continue
        cov = (al['qe'] - al['qs']) / al['qlen']
        if cov < W.COVERAGE_MIN:
            funnel['fail_coverage'] += 1
            continue
        funnel['pass_coverage'] += 1
        if not (al['mrna']['identity'] >= W.IDENTITY_MIN):
            funnel['fail_identity'] += 1
            continue
        funnel['pass_identity'] += 1
        if al['mrna']['frameshift'] or al['mrna']['stopcodon']:
            funnel['fail_frameshift_or_internal_stop'] += 1
            continue
        funnel['pass_no_fs_stop'] += 1
        cds, has_stop = chain_from_alignment(al)
        if not has_stop:
            funnel['fail_no_stop_codon_feature'] += 1
            continue
        funnel['pass_stop_codon'] += 1
        if len(cds) > W.MAX_EXONS:
            funnel['fail_exon_cap'] += 1
            continue
        strand = al['mrna']['strand']
        failures, motifs = MF.chain_checks(dict(strand=strand, CDS=cds), seqs[al['ctg']])
        if failures:
            funnel['fail_chain_checks'] += 1
            for r in failures:
                fail_reasons[r] += 1
            continue
        funnel['pass_chain_checks'] += 1
        species, pid = al['query'].split('|', 1)
        cands.append(dict(seqid=al['ctg'], strand=strand, CDS=cds, span=[cds[0][0], cds[-1][1]], source_protein=pid, source_species=species,
                          identity=al['mrna']['identity'], positive=al['mrna']['positive'], coverage=round(cov, 4), n_exons=len(cds),
                          cds_bp=sum(b - a for a, b, _ in cds), splice_motifs=motifs, miniprot_id=al['mrna']['id']))
    # same-locus de-duplication: greedy by identity (then coverage, CDS bp), per contig & strand
    cands.sort(key=lambda r: (-r['identity'], -r['coverage'], -r['cds_bp'], r['seqid'], r['span'][0], r['source_protein']))
    accepted = defaultdict(list)       # (seqid, strand) -> [(a, b)]
    seeds = []
    for r in cands:
        key = (r['seqid'], r['strand'])
        a, b = r['span']
        if any(x < b and y > a for x, y in accepted[key]):
            funnel['dropped_same_locus_lower_identity'] += 1
            continue
        accepted[key].append((a, b))
        seeds.append(r)
    funnel['seeds'] = len(seeds)
    seeds.sort(key=lambda r: (r['seqid'], r['span'][0], r['strand']))
    for i, r in enumerate(seeds):
        r['id'] = f'{sp[:3]}S2_{i + 1:06d}'
        r['gene_id'] = r['id']
    # cross-strand overlaps (allowed in the per-strand labels; reported only)
    xs = 0
    for c in lens:
        plus = sorted(r['span'] for r in seeds if r['seqid'] == c and r['strand'] == '+')
        minus = [r['span'] for r in seeds if r['seqid'] == c and r['strand'] == '-']
        starts = [a for a, _ in plus]
        import bisect
        for a, b in minus:
            i = bisect.bisect_left(starts, b)
            if any(x < b and y > a for x, y in plus[max(0, i - 50):i]):
                xs += 1
    outd = W.DATA / 'seeds'
    outd.mkdir(parents=True, exist_ok=True)
    with open(outd / f'seeds_S2_{sp}.gff3', 'w') as f:
        f.write('##gff-version 3\n')
        for c, n in lens.items():
            f.write(f'##sequence-region {c} 1 {n}\n')
        for r in seeds:
            a, b = r['span']
            common = f"source_protein={r['source_protein']};source_species={r['source_species']};identity={r['identity']};coverage={r['coverage']};positive={r['positive']};n_exons={r['n_exons']}"
            f.write(f"{r['seqid']}\tm36_miniprot\tgene\t{a + 1}\t{b}\t.\t{r['strand']}\t.\tID={r['id']};{common}\n")
            f.write(f"{r['seqid']}\tm36_miniprot\tmRNA\t{a + 1}\t{b}\t.\t{r['strand']}\t.\tID={r['id']}.m;Parent={r['id']};{common}\n")
            for j, (x, y, p) in enumerate(r['CDS']):
                f.write(f"{r['seqid']}\tm36_miniprot\tCDS\t{x + 1}\t{y}\t.\t{r['strand']}\t{p}\tID={r['id']}.cds{j};Parent={r['id']}.m\n")
    json.dump(seeds, open(outd / f'seeds_S2_{sp}.json', 'w'))
    by_contig = {c: dict(split=split_of[c], Mb=round(n / 1e6, 2), seeds=sum(r['seqid'] == c for r in seeds)) for c, n in lens.items()}
    strata = lambda rows: {'1': sum(r['n_exons'] == 1 for r in rows), '2-5': sum(2 <= r['n_exons'] <= 5 for r in rows), '>=6': sum(r['n_exons'] >= 6 for r in rows)}
    ident = np.array([r['identity'] for r in seeds]) if seeds else np.zeros(1)
    stats = dict(species=sp, raw=str(raw), thresholds=dict(coverage_min=W.COVERAGE_MIN, identity_min=W.IDENTITY_MIN, max_exons=W.MAX_EXONS),
                 A_contigs=list(lens), funnel=dict(funnel), chain_check_failures=dict(fail_reasons), by_contig=by_contig,
                 by_split={s: sum(split_of[r['seqid']] == s for r in seeds) for s in ('train', 'val')},
                 exon_strata=strata(seeds), exon_strata_train=strata([r for r in seeds if split_of[r['seqid']] == 'train']),
                 by_source_species=dict(Counter(r['source_species'] for r in seeds)),
                 identity_quantiles={q: float(np.quantile(ident, q)) for q in (0.05, 0.25, 0.5, 0.75, 0.95)},
                 seeds_per_Mb_train=len([r for r in seeds if split_of[r['seqid']] == 'train']) / (sum(n for c, n in lens.items() if split_of[c] == 'train') / 1e6),
                 cross_strand_overlapping_seeds=xs, splice_motifs=dict(sum((Counter(r['splice_motifs']) for r in seeds), Counter())),
                 seconds=time.time() - t0)
    json.dump(stats, open(outd / f'seed_stats_S2_{sp}.json', 'w'), indent=1)
    print(json.dumps({k: stats[k] for k in ('species', 'funnel', 'by_split', 'exon_strata', 'seeds_per_Mb_train', 'seconds')}), flush=True)


if __name__ == '__main__':
    main()
