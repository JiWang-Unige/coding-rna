"""Report 37 part B: panel long-gene failure stratification, T2 vs Carbon-A unfiltered (descriptive; plan fixed before computing).
Error categories: the M29 R4 structural rules (m29_r4_error_types: Index, ov_bp, classify_local; same priority order, no confidence_filtered category).
usage: long_gene_strata.py <c3_data_root> <out.json>      (same inputs as c3_calib.py)
"""
import glob
import json
import pickle
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'M30-HOLDOUT-ADAPT' / 'scripts'))
import m30_common as M30                                 # noqa: E402,F401
import m29_common as M29                                 # noqa: E402
import m29_r4_error_types as R4                          # noqa: E402

CORE10 = ['dryas_octopetala', 'eschscholzia_californica', 'akebia_trifoliata', 'hakea_prostrata', 'zizania_latifolia', 'themeda_triandra', 'juncus_effusus',
          'rhynchospora_breviuscula', 'acorus_gramineus', 'amborella_trichopoda']
MIN_BP = R4.MIN_BP


def ex_bin(n):
    return '1' if n == 1 else '2' if n == 2 else '3-5' if n <= 5 else '6-9' if n <= 9 else '10-19' if n <= 19 else '20+'


def len_bin(bp):
    kb = bp / 1000
    return '<2kb' if kb < 2 else '2-5kb' if kb < 5 else '5-10kb' if kb < 10 else '10-20kb' if kb < 20 else '20-50kb' if kb < 50 else '>=50kb'


EXB = ['1', '2', '3-5', '6-9', '10-19', '20+']
LENB = ['<2kb', '2-5kb', '5-10kb', '10-20kb', '20-50kb', '>=50kb']


def classify(r, k, Er, pred_idx, preds, ref_index, iso_keys):
    if k in preds:
        return 'hit'
    if len(Er) < 2:
        return 'single_exon_miss'
    strand = r['strand']
    lo, hi = Er[0][0], Er[-1][1]
    cands = [(R4.ov_bp(c['exons'], Er), c) for c in pred_idx.query(strand, lo, hi)]
    cands = [(o, c) for o, c in cands if o > 0]
    if any(c['key'] in iso_keys and c['key'] != k for o, c in cands):
        return 'other_isoform'
    if not cands:
        other = '-' if strand == '+' else '+'
        oc = [c for c in pred_idx.query(other, lo, hi) if R4.ov_bp(c['exons'], Er) > 0]
        return 'wrong_strand_only' if oc else 'absent_no_call'
    cands.sort(key=lambda x: -x[0])
    best = cands[0][1]
    others = [x for x in ref_index.query(strand, best['exons'][0][0], best['exons'][-1][1]) if x['id'] != r['id'] and R4.ov_bp(x['exons'], best['exons']) >= MIN_BP]
    if others:
        return 'merge'
    if len([c for o, c in cands if o >= MIN_BP]) >= 2:
        return 'split'
    return R4.classify_local(Er, best['exons'], strand)[0]


def main():
    root, outp = Path(sys.argv[1]), sys.argv[2]
    post_root = root / 'M34-NTV2-SCALEUP' / 'panel' / 'post_T2'
    rows = []                                            # (species, n_exons, span_bp, cat_T2, cat_CA)
    for sp in CORE10:
        seqs = sorted({Path(f).name.rsplit('_', 1)[0] for f in glob.glob(str(post_root / sp / 'tau1' / '*.json'))})
        for s in seqs:
            ref = pickle.load(open(root / 'PANEL-P6' / '_refcache' / f'ref_{s}.pkl', 'rb'))
            M29.SEQ2SPECIES[s], M29.LENGTH[s] = sp, ref['length']
            t2 = {}
            for st in '+-':
                for g in json.load(open(post_root / sp / 'tau1' / f'{s}_{st}.json'))['genes']:
                    ex = [(int(a), int(b)) for a, b in g['exons']]
                    if sum(b - a for a, b in ex) >= 96:
                        k = (s, st, tuple(ex))
                        t2.setdefault(k, dict(key=k, strand=st, exons=ex))
            ca = {}
            for tx in M29.load_txs(str(root / 'PANEL-P6' / 'carbonA' / sp / s / 'prediction.gff'), s):
                ex = [(a, b) for a, b, _ in tx['CDS']]
                if sum(b - a for a, b in ex) >= 96:
                    k = M29.key(tx)
                    ca.setdefault(k, dict(key=k, strand=tx['strand'], exons=ex))
            idx = {'T2': R4.Index(list(t2.values())), 'CA': R4.Index(list(ca.values()))}
            iso = defaultdict(set)
            for t in ref['isoforms']:
                iso[t['gene_id']].add(M29.key(t))
            ref_index = R4.Index([dict(strand=t['strand'], exons=[(a, b) for a, b, _ in t['CDS']], id=t['id']) for t in ref['refs']])
            for r in ref['refs']:
                k = M29.key(r)
                Er = [(a, b) for a, b, _ in r['CDS']]
                c1 = classify(r, k, Er, idx['T2'], t2, ref_index, iso[r['gene_id']])
                c2 = classify(r, k, Er, idx['CA'], ca, ref_index, iso[r['gene_id']])
                rows.append((sp, len(Er), Er[-1][1] - Er[0][0], c1, c2))
        print(sp, 'done', flush=True)
    res = dict(n_ref=len(rows), species={})

    def block(sub):
        n = len(sub)
        out = dict(n=n)
        for nm, fn, bins in (('by_exons', lambda r: ex_bin(r[1]), EXB), ('by_span', lambda r: len_bin(r[2]), LENB)):
            tab = {}
            for b in bins:
                x = [r for r in sub if fn(r) == b]
                h1 = sum(r[3] == 'hit' for r in x); h2 = sum(r[4] == 'hit' for r in x)
                tab[b] = dict(n=len(x), recall_T2=100 * h1 / max(len(x), 1), recall_CA=100 * h2 / max(len(x), 1),
                              contrib_to_recall_gap=100 * (h2 - h1) / max(n, 1))
            out[nm] = tab
        out['recall_T2'] = 100 * sum(r[3] == 'hit' for r in sub) / max(n, 1)
        out['recall_CA'] = 100 * sum(r[4] == 'hit' for r in sub) / max(n, 1)
        gap = [r for r in sub if r[3] != 'hit' and r[4] == 'hit']
        rev = [r for r in sub if r[4] != 'hit' and r[3] == 'hit']
        out['T2_miss_CA_hit'] = dict(n=len(gap), categories=dict(Counter(r[3] for r in gap).most_common()),
                                     by_exons={b: dict(Counter(r[3] for r in gap if ex_bin(r[1]) == b).most_common()) for b in EXB})
        out['CA_miss_T2_hit'] = dict(n=len(rev), categories=dict(Counter(r[4] for r in rev).most_common()))
        out['all_T2_misses'] = dict(Counter(r[3] for r in sub if r[3] != 'hit').most_common())
        out['all_CA_misses'] = dict(Counter(r[4] for r in sub if r[4] != 'hit').most_common())
        return out

    for sp in CORE10:
        res['species'][sp] = block([r for r in rows if r[0] == sp])
    res['pooled'] = block(rows)
    json.dump(res, open(outp, 'w'), indent=1)
    for sp in CORE10 + ['pooled']:
        b = res['species'].get(sp, res['pooled'])
        g = b['T2_miss_CA_hit']
        top = ', '.join(f'{k} {100 * v / max(g["n"], 1):.0f}%' for k, v in list(g['categories'].items())[:5])
        print(f"{sp:26s} n {b['n']:6d} recall T2 {b['recall_T2']:5.1f} CA {b['recall_CA']:5.1f} | gap genes {g['n']:5d} (reverse {b['CA_miss_T2_hit']['n']:5d}): {top}")
    P = res['pooled']
    for nm in ('by_exons', 'by_span'):
        print(nm, {b: (v['n'], round(v['recall_T2'], 1), round(v['recall_CA'], 1), round(v['contrib_to_recall_gap'], 2)) for b, v in P[nm].items()})


if __name__ == '__main__':
    main()
