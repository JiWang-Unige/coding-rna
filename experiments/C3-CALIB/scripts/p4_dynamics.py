"""Report 37 part A: P4 training-dynamics decomposition (descriptive; plan fixed before computing).
usage: p4_dynamics.py <data_root> <out.json>
  data_root/<arm>/metrics.jsonl, selection.json, val_<step>/<contig>/prediction/full_grammar.gff3   arms: pc2_T10_k50U[_s2], pc2_T10_k50F[_s2], pc2_T10
  data_root/_cache/ref_<contig>.pkl   (the M34 ruler objects of the 8 new-species validation segments)
"""
import json
import pickle
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'M30-HOLDOUT-ADAPT' / 'scripts'))
import m30_common as M30                                 # noqa: E402,F401
import m29_common as M29                                 # noqa: E402
import m30_score as SC                                   # noqa: E402

STEPS = (768, 1536, 2304, 3072)
SEEDS = {'s1': ('pc2_T10_k50U', 'pc2_T10_k50F'), 's2': ('pc2_T10_k50U_s2', 'pc2_T10_k50F_s2')}
BLK, DRAWS, SEED = 1_000_000, 2000, 0


def main():
    root, outp = Path(sys.argv[1]), sys.argv[2]
    refs, species = {}, {}
    for p in sorted((root / '_cache').glob('ref_*__val8M.pkl')):
        o = pickle.load(open(p, 'rb'))
        c = o['seq']
        M29.SEQ2SPECIES[c], M29.LENGTH[c] = o['species'], o['length']
        refs[c] = {M29.key(t) for t in o['refs']}
        species[c] = o['species']
    contigs = sorted(refs)
    nref = sum(len(v) for v in refs.values())
    chains = {}
    arms = sorted({a for v in SEEDS.values() for a in v} | {'pc2_T10'})
    for a in arms:
        for s in STEPS:
            chains[a, s] = {c: set(SC.load_gff_chains(str(root / a / f'val_{s}' / c / 'prediction' / 'full_grammar.gff3'), c)) for c in contigs}
    # sanity against metrics.jsonl (own_species)
    check = []
    for a in arms:
        for line in open(root / a / 'metrics.jsonl'):
            m = json.loads(line)
            if m.get('tag') != 'own_species' or m['step'] not in STEPS:
                continue
            for reg, v in m['per_region'].items():
                c = reg.split(':')[0]
                ch = chains[a, m['step']][c]
                got = (len(ch & refs[c]), len(ch), len(refs[c]))
                check.append(got == (v['tp'], v['pred'], v['ref']))
    sel = {a: json.load(open(root / a / 'selection.json'))['best_step'] for a in arms}

    def stats(a, s, cs=contigs):
        tp = sum(len(chains[a, s][c] & refs[c]) for c in cs)
        pr = sum(len(chains[a, s][c]) for c in cs)
        rf = sum(len(refs[c]) for c in cs)
        return dict(tp=tp, pred=pr, ref=rf, F1=200 * tp / (pr + rf), precision=100 * tp / max(pr, 1), recall=100 * tp / rf)

    # blocks: (contig, first CDS start // 1 Mb) of references and of every scored chain
    bl = set()
    for c in contigs:
        bl.update((c, k[2][0][0] // BLK) for k in refs[c])
        for (a, s), d in chains.items():
            bl.update((c, k[2][0][0] // BLK) for k in d[c])
    bl = sorted(bl)
    idx = {b: i for i, b in enumerate(bl)}
    rfv = np.zeros(len(bl))
    for c in contigs:
        for k in refs[c]:
            rfv[idx[c, k[2][0][0] // BLK]] += 1

    def vec(a, s):
        tp, pr = np.zeros(len(bl)), np.zeros(len(bl))
        for c in contigs:
            for k in chains[a, s][c]:
                i = idx[c, k[2][0][0] // BLK]
                pr[i] += 1
                tp[i] += k in refs[c]
        return tp, pr

    vecs = {k: vec(*k) for k in chains}
    rng = np.random.default_rng(SEED)
    draws = [rng.integers(0, len(bl), len(bl)) for _ in range(DRAWS)]
    f1 = lambda v, ii: 200 * v[0][ii].sum() / (v[1][ii].sum() + rfv[ii].sum())

    def boot(u, f):
        d = np.array([f1(vecs[u], ii) - f1(vecs[f], ii) for ii in draws])
        return [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))]

    res = dict(n_ref=nref, contigs={c: species[c] for c in contigs}, sanity_metrics_match=f'{sum(check)}/{len(check)}', selected=sel, seeds={})
    for sd, (u, f) in SEEDS.items():
        rows = {}
        for s in STEPS:
            U, F, T = stats(u, s), stats(f, s), stats('pc2_T10', s)
            d = U['F1'] - F['F1']
            tp_term = 200 * (U['tp'] - F['tp']) / (F['pred'] + F['ref'])
            pos = sum(stats(u, s, [c])['F1'] > stats(f, s, [c])['F1'] for c in contigs)
            rows[s] = dict(U=U, F=F, T10=T, diff=d, diff_CI95=boot((u, s), (f, s)), tp_term=tp_term, pred_term=d - tp_term,
                           species_positive=pos, U_minus_T10=U['F1'] - T['F1'], F_minus_T10=F['F1'] - T['F1'])
        grid = {f'{su}|{sf}': stats(u, su)['F1'] - stats(f, sf)['F1'] for su in STEPS for sf in STEPS}
        res['seeds'][sd] = dict(arms=(u, f), per_step=rows, selected=(sel[u], sel[f]), selected_diff=grid[f'{sel[u]}|{sel[f]}'],
                                grid=grid, grid_positive=sum(v > 0 for v in grid.values()), grid_min=min(grid.values()), grid_max=max(grid.values()))
    json.dump(res, open(outp, 'w'), indent=1)
    print('sanity', res['sanity_metrics_match'], 'selected', sel)
    for sd, r in res['seeds'].items():
        print(f'== {sd} {r["arms"]} selected {r["selected"]} diff {r["selected_diff"]:+.2f}; grid positive {r["grid_positive"]}/16 [{r["grid_min"]:+.2f}, {r["grid_max"]:+.2f}]')
        for s, x in r['per_step'].items():
            print(f"{s:5d} U {x['U']['F1']:.2f} (tp {x['U']['tp']} pred {x['U']['pred']} R {x['U']['recall']:.1f} P {x['U']['precision']:.1f}) "
                  f"F {x['F']['F1']:.2f} (tp {x['F']['tp']} pred {x['F']['pred']} R {x['F']['recall']:.1f} P {x['F']['precision']:.1f}) "
                  f"diff {x['diff']:+.2f} {[round(v, 2) for v in x['diff_CI95']]} tp_term {x['tp_term']:+.2f} pred_term {x['pred_term']:+.2f} "
                  f"pos {x['species_positive']}/8 T10 {x['T10']['F1']:.2f}")


if __name__ == '__main__':
    main()
