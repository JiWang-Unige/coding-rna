"""Protein-homology audit of the FULL panel proteomes (12 species, all NCBI protein records, isoforms included) against the M34 training proteomes
(T10 = 8 new species + Arabidopsis + rice; T2 = Arabidopsis + rice).  Same search as M34 (`m34_homology_audit`: MMseqs2 -s 4.0 -e 1e-5 --max-seqs 300,
query coverage >= 0.7), but every protein instead of a 1,500 sample, and the whole best-hit identity distribution instead of three thresholds.
Purely descriptive (nothing is filtered).  CPU only.
usage: panel_homology.py [threads]
"""
import gzip
import json
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime

import panel_common as P
import m34_common as C
from m34_homology_audit import faa, MMSEQS
from m34_schedule import LEGACY, NEW

OUTD = P.HOMOLOGY
BINS = [('no_hit', None), ('<30', (0, 30)), ('30-40', (30, 40)), ('40-50', (40, 50)), ('50-70', (50, 70)), ('70-90', (70, 90)), ('>=90', (90, 101))]


def binned(best):
    """best: best-hit identity in PERCENT (qcov >= 0.7) or None -> bin label."""
    if best is None:
        return 'no_hit'
    for lab, (lo, hi) in BINS[1:]:
        if lo <= best < hi:
            return lab
    raise ValueError(best)


def summarize(n, best):
    """n proteins, best = list of best identities (percent) for the proteins WITH a hit -> distribution dict."""
    hist = Counter(binned(b) for b in best)
    hist['no_hit'] = n - len(best)
    srt = sorted(best)
    return dict(n_proteins=n, n_with_hit=len(best), frac={lab: hist[lab] / n for lab, _ in BINS},
                frac_ge={str(t): sum(b >= t for b in best) / n for t in (30, 50, 70, 90)},
                median_best_identity_with_hit=(srt[len(srt) // 2] if srt else None), mean_best_identity_with_hit=(sum(srt) / len(srt) if srt else None))


def build_inputs():
    OUTD.mkdir(parents=True, exist_ok=True)
    rows = C.new_species()
    train = C.DATA / 'homology' / 'train_T10.faa'            # M34's target set (10 species, '<key>|' prefixed ids); rebuilt here if absent
    if not train.exists():
        train = OUTD / 'train_T10.faa'
        train.unlink(missing_ok=True)
        for key in NEW:
            r = rows[key]
            faa(key, C.raw_dir(r) / f"{r['accession']}_{r['assembly_name']}_protein.faa.gz", train)
        faa('arabidopsis_thaliana', C.RAW / 'legacy_proteins/GCF_000001735.4_TAIR10.1_protein.faa.gz', train)
        faa('oryza_sativa', C.RAW / 'legacy_proteins/GCF_034140825.1_ASM3414082v1_protein.faa.gz', train)
    query = OUTD / 'panel_full.faa'
    query.unlink(missing_ok=True)
    for key, row in P.species().items():
        faa(key, P.raw_dir(row) / f"{P.stem(row)}_protein.faa.gz", query)
    return query, train


def main(threads=6):
    query, train = build_inputs()
    res = OUTD / 'best_hits_full.tsv'
    gz = res.with_suffix('.tsv.gz')                    # the 6 GB hit table is kept gzipped after the run (disk budget)
    if not res.exists() and not gz.exists():
        subprocess.run([MMSEQS, 'easy-search', str(query), str(train), str(res), str(OUTD / 'tmp'), '--threads', str(threads), '-s', '4.0', '-e', '1e-5',
                        '--max-seqs', '300', '--format-output', 'query,target,pident,qcov,tcov,evalue,bits'], check=True)
        subprocess.run(['rm', '-rf', str(OUTD / 'tmp')])
    nq, ntarget = Counter(), Counter()
    for path, cnt in ((query, nq), (train, ntarget)):
        with open(path) as f:
            for line in f:
                if line[0] == '>':
                    cnt[line[1:].split('|')[0]] += 1
    t2, t10 = set(LEGACY), set(LEGACY + NEW)
    best10, best2, best_sp = {}, {}, defaultdict(dict)          # query -> best pident (percent); query -> {train species: best pident}
    with (gzip.open(gz, 'rt') if not res.exists() else open(res)) as f:
        for line in f:
            q, t, pid, qc, tc, ev, bits = line.rstrip('\n').split('\t')
            if float(qc) < 0.7:
                continue
            pid, ts = float(pid), t.split('|')[0]
            if ts in t10 and pid > best10.get(q, -1):
                best10[q] = pid
            if ts in t2 and pid > best2.get(q, -1):
                best2[q] = pid
            if pid > best_sp[q].get(ts, -1):
                best_sp[q][ts] = pid
    out = {}
    for key in P.species():
        qs10 = [v for q, v in best10.items() if q.split('|')[0] == key]
        qs2 = [v for q, v in best2.items() if q.split('|')[0] == key]
        nearest = Counter(max(d, key=d.get) for q, d in best_sp.items() if q.split('|')[0] == key and d)
        out[key] = dict(T10=summarize(nq[key], qs10), T2=summarize(nq[key], qs2),
                        nearest_training_species=dict(nearest.most_common()),
                        per_training_species_frac_ge70={sp: sum(d.get(sp, 0) >= 70 for q, d in best_sp.items() if q.split('|')[0] == key) / nq[key] for sp in sorted(t10)})
    meta = dict(search='mmseqs easy-search -s 4.0 -e 1e-5 --max-seqs 300; best hit with query coverage >= 0.7; identity = mmseqs pident over the aligned region',
                queries={k: nq[k] for k in P.species()}, targets=dict(ntarget), target_file=str(train), threads=threads,
                computed=datetime.now().astimezone().isoformat(timespec='seconds'), note='all NCBI protein records (isoforms included); descriptive only, nothing filtered')
    json.dump(dict(meta=meta, species=out), open(OUTD / 'homology_audit_full.json', 'w'), indent=1)
    with open(OUTD / 'homology_table.tsv', 'w') as f:
        cols = ['species', 'n_proteins', 'T10_no_hit', 'T10_<30', 'T10_30-40', 'T10_40-50', 'T10_50-70', 'T10_70-90', 'T10_>=90', 'T10_ge50', 'T10_ge70', 'T10_ge90', 'T10_median_id',
                'T2_no_hit', 'T2_ge50', 'T2_ge70', 'T2_ge90', 'T2_median_id', 'nearest_training_species_top3']
        f.write('\t'.join(cols) + '\n')
        for key, d in out.items():
            a, b = d['T10'], d['T2']
            row = [key, a['n_proteins']] + [f"{a['frac'][lab]:.3f}" for lab, _ in BINS] + [f"{a['frac_ge'][t]:.3f}" for t in ('50', '70', '90')] + [a['median_best_identity_with_hit']]
            row += [f"{b['frac']['no_hit']:.3f}"] + [f"{b['frac_ge'][t]:.3f}" for t in ('50', '70', '90')] + [b['median_best_identity_with_hit']]
            row += [';'.join(f'{s}:{n}' for s, n in list(d['nearest_training_species'].items())[:3])]
            f.write('\t'.join(map(str, row)) + '\n')
    for key, d in out.items():
        print(key, 'n', d['T10']['n_proteins'], 'T10 ge70', round(d['T10']['frac_ge']['70'], 3), 'no_hit', round(d['T10']['frac']['no_hit'], 3), 'T2 ge70', round(d['T2']['frac_ge']['70'], 3), flush=True)
    print('PANEL_HOMOLOGY_DONE', flush=True)


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
