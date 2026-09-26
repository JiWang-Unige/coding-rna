#!/usr/bin/env python3
"""Fixed-prediction reference representation control, not independent truth."""
import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import partial_semantics as P
S = P.S
SEQID, LENGTH = 'NC_003074.8', 23459830
SPECIES = 'arabidopsis_thaliana'


def fasta_target(path, target):
    chunks, found = [], False
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if line.startswith('>'):
                if found:
                    break
                found = line[1:].split()[0] == target
            elif found:
                chunks.append(line.strip())
    if not found:
        raise ValueError(f'Missing FASTA target {target}')
    return ''.join(chunks).upper()


def normalize_araport(line):
    if line.startswith('#'):
        return None
    fields = line.rstrip('\n').split('\t')
    if len(fields) != 9:
        raise ValueError('Expected GFF3 row')
    if fields[0] != 'Chr3':
        return None
    fields[0] = SEQID
    if fields[2] == 'gene':
        attrs = S.E.parse_attrs(fields[8])
        biotype = attrs.get('locus_type')
        if biotype is None:
            raise ValueError('Araport gene lacks locus_type')
        fields[8] += ';gene_biotype=' + biotype
    return '\t'.join(fields) + '\n'


def pick_primary(choices):
    return [min(ts, key=lambda t: (-sum(b-a for a,b,_ in t['CDS']), t['id']))
            for ts in choices.values() if ts]


def source_views(path, sequence):
    coding = S.E.parse_annotation(path, {SEQID: LENGTH}, protein_coding_only=True)
    raw_g, raw_t, raw_c, tags = P.raw_features(open(path), SEQID)
    parent = S.E.primary_transcripts(coding)
    eligible = defaultdict(list)
    paired = {}
    for tx in coding['transcripts'].values():
        if not tx['CDS']:
            continue
        gene = tx['gene_id']
        gene_attrs = raw_g[gene]
        locus = gene_attrs.get('locus_tag', gene)
        cp = any(S.E.is_partial(a) for a in raw_c[tx['id']])
        semantic = not cp and P.compatible(tx, S.V.transcript_sequence(tx, {SEQID: sequence}).upper())
        parent_ok = not coding['genes'][gene]['partial'] and not tx['partial']
        if semantic:
            eligible[gene].append(tx)
        identity = (locus, S.key(tx))
        row = paired.setdefault(identity, {'parent_eligible': False, 'CDS_assessable': False})
        row['parent_eligible'] |= parent_ok
        row['CDS_assessable'] |= semantic
    cds = pick_primary(eligible)
    parent_genes = {t['gene_id'] for t in parent}
    cds_genes = {t['gene_id'] for t in cds}
    return {'parent_filter': parent, 'CDS_assessable': cds}, coding, paired, {
        'gene_count': len(coding['genes']), 'raw_feature_row_flags': tags,
        'parent_primary': len(parent), 'CDS_assessable_primary': len(cds),
        'new_genes_under_CDS_policy': len(cds_genes-parent_genes),
        'old_genes_without_compatible_CDS': len(parent_genes-cds_genes),
        'same_gene_primary_chain_changes': sum(S.key(t) != S.key({x['gene_id']:x for x in cds}[t['gene_id']]) for t in parent if t['gene_id'] in cds_genes)}


def background_mask(txs):
    mask = np.zeros(LENGTH, bool)
    for t in txs:
        parts = t['exon'] or t['CDS']
        mask[parts[0][0]:parts[-1][1]] = True
    return mask


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', required=True)
    out = Path(parser.parse_args().output_dir)
    out.mkdir(parents=True, exist_ok=True)
    if (out/'reference_pair.json').exists():
        raise FileExistsError('Do not overwrite completed results')
    data = S.ROOT/'data/m26_reference_pair_20260926'
    sequence = S.V.read_fasta(S.ROOT/'data/m1_screen'/SPECIES/'genome.fa')[SEQID].upper()
    assert len(sequence) == LENGTH
    assert fasta_target(data/'TAIR10_chr_all.fas.gz', 'Chr3') == sequence, 'FASTA mismatch: stop coordinate comparison'
    normalized = out/'araport_chr3.normalized.gff3'
    # This dated TAIR file contains byte 0x91 in curator_summary; Latin-1 preserves bytes.
    with gzip.open(data/'Araport11.20240701.gff.gz', 'rt', encoding='latin-1') as inp, normalized.open('x', encoding='utf-8') as handle:
        for line in inp:
            row = normalize_araport(line)
            if row:
                handle.write(row)
    paths = {'RefSeq': S.ROOT/'data/m1_screen'/SPECIES/'reference.gff3', 'Araport11_20240701': normalized}
    previous = json.loads((S.ROOT/'outputs/M26-SAME-SCOPE-MECHANISM-R1/summary.json').read_text())
    predictions = {m: S.E.primary_transcripts(S.E.parse_annotation(S.ROOT/t.format(species=SPECIES), {SEQID:LENGTH}), complete_only=False) for m,t in S.METHODS.items()}
    result = {'status':'development_reference_representation_control_shared_annotation_ancestry',
              'species':SPECIES, 'seqid':SEQID, 'length':LENGTH, 'FASTA_exact_uppercase_match':True,
              'setaria_access':False, 'sources':{}}
    paired_sources = {}
    for source, path in paths.items():
        views, coding, paired, stats = source_views(path, sequence)
        paired_sources[source] = paired
        if source == 'RefSeq':
            assert len(views['parent_filter']) == 4151
        all_coding = np.zeros(LENGTH, bool)
        for g in coding['genes'].values():
            all_coding[g['start']:g['end']] = True
        stats['fixed_all_coding_gene_background_bp'] = int((~all_coding).sum())
        stats['policies'] = {}
        for policy, refs in views.items():
            background = ~background_mask(refs)
            ref_keys = {S.key(t) for t in refs}
            per_method = {}
            for method, txs in predictions.items():
                keys = {S.key(t) for t in txs}
                span, cds = P.coverage_masks(txs, LENGTH)
                metrics = S.prf(len(keys & ref_keys), len(keys), len(ref_keys))
                counts = {'span_fp_bp':int((span & background).sum()), 'background_bp':int(background.sum()),
                          'CDS_coverage_on_background_bp':int((cds & background).sum())}
                fpr = counts['span_fp_bp']/counts['background_bp']
                if source == 'RefSeq' and policy == 'parent_filter':
                    old = previous['methods'][method]['per_species'][SPECIES]
                    assert metrics == old['exact_chain']
                    assert abs(fpr-old['frozen_metrics']['intergenic_FPR']) < 1e-12
                per_method[method] = {'chain':metrics, 'policy_span_background_counts':counts,
                                      'policy_span_FPR':fpr,
                                      'all_coding_background_span_bp':int((span & ~all_coding).sum()),
                                      'all_coding_background_CDS_bp':int((cds & ~all_coding).sum())}
            stats['policies'][policy] = {'reference_unique_chains':len(ref_keys), 'methods':per_method}
        result['sources'][source] = stats
    n, a = paired_sources.values()
    shared = n.keys() & a.keys()
    same_cds = {k for k in shared if n[k]['CDS_assessable'] and a[k]['CDS_assessable']}
    discordant = {k for k in same_cds if not n[k]['parent_eligible'] and a[k]['parent_eligible']}
    result['shared_chain_control'] = {
        'same_locus_and_CDS_chain_records':len(shared),
        'CDS_assessable_in_both':len(same_cds),
        'CDS_assessable_parent_eligibility_RefSeq_Araport':dict(Counter(f'{int(n[k]["parent_eligible"])},{int(a[k]["parent_eligible"])}' for k in same_cds)),
        'RefSeq_excluded_Araport_retained_shared_chains':len(discordant),
        'RefSeq_excluded_Araport_retained_shared_loci':len({k[0] for k in discordant}),
        'methods':{}}
    for method, txs in predictions.items():
        keys = {S.key(t) for t in txs}
        matched = {k for k in discordant if k[1] in keys}
        result['shared_chain_control']['methods'][method] = {'exact_discordant_chains':len(matched), 'exact_discordant_loci':len({k[0] for k in matched})}
    with (out/'shared_chain_eligibility.tsv').open('x') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['locus','strand','CDS_0based_halfopen','RefSeq_parent','Araport_parent','CDS_assessable_both'])
        for locus, key in sorted(shared):
            writer.writerow([locus,key[1],';'.join(f'{x}-{y}' for x,y in key[2]),n[(locus,key)]['parent_eligible'],a[(locus,key)]['parent_eligible'],(locus,key) in same_cds])
    with (out/'reference_pair.json').open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
