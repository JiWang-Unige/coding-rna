"""M32 block bootstrap (1-Mb blocks, 2000 reps, seed 0) of pooled exact-chain F1 on DEV and of paired differences of PlantCAD2-C1 to
NT-v2 A / NT-v2 B (M31 predictions), GENERanno fresh-joint and Carbon-A filtered (same blocks).  Usage: c1_boot.py <m32_dev_dir> <out.json>"""
import json, sys, os
import numpy as np
import m32_common as M
import a1_common as K
sys.path.insert(0, str(K.CA_SRC))
import a1_eval as E              # noqa: E402
import a1_compare as AC          # noqa: E402

run, outp = sys.argv[1], sys.argv[2]
M31 = os.path.expanduser('~/coding-rna/outputs/M31-C1-NTV2')
seqs = list(M.DEV_SEQS)
sets = {'PlantCAD2_C1': {s: E.load_chains(run, s) for s in seqs},
        'NTv2_A_frozen': {s: E.load_chains(M31 + '/dev_final', s) for s in seqs},
        'NTv2_B_LoRA': {s: E.load_chains(M31 + '/devB_final', s) for s in seqs},
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
f1 = lambda tp, p, r: 200 * tp.sum() / (p.sum() + r.sum())
res = {n: {'F1': f1(*c)} for n, c in C.items()}
draws = {n: [] for n in sets}; diff = {n: [] for n in sets if n != 'PlantCAD2_C1'}
for _ in range(2000):
    ii = rng.integers(0, len(bl), len(bl))
    f = {n: f1(c[0][ii], c[1][ii], c[2][ii]) for n, c in C.items()}
    for n in sets:
        draws[n].append(f[n])
    for n in diff:
        diff[n].append(f['PlantCAD2_C1'] - f[n])
for n in sets:
    res[n]['CI95'] = [float(np.percentile(draws[n], 2.5)), float(np.percentile(draws[n], 97.5))]
for n in diff:
    res[n]['diff_PC2_minus_this_CI95'] = [float(np.percentile(diff[n], 2.5)), float(np.percentile(diff[n], 97.5))]
    res[n]['diff_point'] = res['PlantCAD2_C1']['F1'] - res[n]['F1']
res['blocks'] = len(bl)
json.dump(res, open(outp, 'w'), indent=1)
print(json.dumps(res))
