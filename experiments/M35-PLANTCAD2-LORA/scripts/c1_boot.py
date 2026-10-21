"""M35 block bootstrap (1-Mb blocks, 2000 reps, seed 0) of pooled / per-species exact-chain F1 on DEV and of PAIRED differences of PlantCAD2+LoRA (M35) to
PlantCAD2 frozen (M32), NT-v2 A frozen / NT-v2 B + LoRA (M31, the pre-declared comparator), GENERanno fresh-joint and Carbon-A filtered (same blocks,
same resampled indices).  Usage: c1_boot.py <m35_dev_dir> <out.json>"""
import json, sys, os
import numpy as np
import m35_common as M
import a1_common as K
sys.path.insert(0, str(K.CA_SRC))
import a1_eval as E              # noqa: E402
import a1_compare as AC          # noqa: E402

run, outp = sys.argv[1], sys.argv[2]
O = os.path.expanduser('~/coding-rna/outputs')
seqs = list(M.DEV_SEQS)
SPECIES = {'NC_003074.8': 'Arabidopsis', 'NC_089041.1': 'rice'}
sets = {'PC2_LoRA_M35': {s: E.load_chains(run, s) for s in seqs},
        'PC2_frozen_M32': {s: E.load_chains(O + '/M32-C1-PLANTCAD2/runA/dev_2048', s) for s in seqs},
        'NTv2_A_frozen_M31': {s: E.load_chains(O + '/M31-C1-NTV2/dev_final', s) for s in seqs},
        'NTv2_B_LoRA_M31': {s: E.load_chains(O + '/M31-C1-NTV2/devB_final', s) for s in seqs},
        'GENERanno_fresh_joint': {s: AC.gff_chains(AC.DEVGFF['B_freshjoint'][s], s) for s in seqs},
        'CarbonA_filtered': {s: AC.carbon_dev(s, True) for s in seqs}}
BLK = 1_000_000
blocks = {}
for s in seqs:
    for k in AC.refs.keys[s]:
        blocks.setdefault((s, k[2][0][0] // BLK), None)
    for name, d in sets.items():
        for k in d[s]:
            blocks.setdefault((s, k[2][0][0] // BLK), None)
bl = sorted(blocks); idx = {b: i for i, b in enumerate(bl)}
bl_seq = np.array([seqs.index(b[0]) for b in bl])


def counts(name):
    tp = np.zeros(len(bl)); pr = np.zeros(len(bl)); rf = np.zeros(len(bl))
    for s in seqs:
        for k in sets[name][s]:
            i = idx[(s, k[2][0][0] // BLK)]; pr[i] += 1; tp[i] += k in AC.refs.keys[s]
        for k in AC.refs.keys[s]:
            rf[idx[(s, k[2][0][0] // BLK)]] += 1
    return tp, pr, rf


C = {n: counts(n) for n in sets}
rng = np.random.default_rng(0)


def f1(c, ii, species=None):
    tp, p, r = (x[ii] for x in c)
    if species is not None:
        m = bl_seq[ii] == species
        tp, p, r = tp[m], p[m], r[m]
    return 200 * tp.sum() / max(p.sum() + r.sum(), 1)


full = np.arange(len(bl))
scopes = {'pooled': None, **{SPECIES[s]: i for i, s in enumerate(seqs)}}
res = {n: {sc: {'F1': f1(c, full, sp)} for sc, sp in scopes.items()} for n, c in C.items()}
draw = {n: {sc: [] for sc in scopes} for n in sets}
for _ in range(2000):
    ii = rng.integers(0, len(bl), len(bl))
    for n, c in C.items():
        for sc, sp in scopes.items():
            draw[n][sc].append(f1(c, ii, sp))
ME = 'PC2_LoRA_M35'
for n in sets:
    for sc in scopes:
        res[n][sc]['CI95'] = [float(np.percentile(draw[n][sc], 2.5)), float(np.percentile(draw[n][sc], 97.5))]
        if n != ME:
            d = np.array(draw[ME][sc]) - np.array(draw[n][sc])
            res[n][sc]['diff_M35_minus_this_point'] = res[ME][sc]['F1'] - res[n][sc]['F1']
            res[n][sc]['diff_M35_minus_this_CI95'] = [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))]
res['blocks'] = len(bl)
json.dump(res, open(outp, 'w'), indent=1)
print(json.dumps({n: {sc: {k: (round(v, 3) if isinstance(v, float) else [round(x, 3) for x in v]) for k, v in d.items()} for sc, d in r.items()} for n, r in res.items() if n != 'blocks'}))
