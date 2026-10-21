"""Block bootstrap (1-Mb genomic blocks, resampled with replacement, 2000 reps) of pooled exact-chain F1 on the DEV chromosomes, and of the
paired difference to the GENERanno fresh-joint predictions (same blocks).  Usage: c1_boot.py <c1_run_dir> <out.json>"""
import json, sys
import numpy as np
import m31_common as M
import a1_common as K
sys.path.insert(0, str(K.CA_SRC))
import a1_eval as E              # noqa: E402
import a1_compare as AC          # noqa: E402

run, outp = sys.argv[1], sys.argv[2]
seqs = list(M.DEV_SEQS)
BLK = 1_000_000
sets = {'NTv2_C1': {s: E.load_chains(run, s) for s in seqs},
        'GENERanno_fresh_joint': {s: AC.gff_chains(AC.DEVGFF['B_freshjoint'][s], s) for s in seqs},
        'CarbonA_filtered': {s: AC.carbon_dev(s, True) for s in seqs}}
blocks = {}
for s in seqs:
    for k in AC.refs.keys[s]:
        blocks.setdefault((s, k[2][0][0] // BLK), None)
    for name, d in sets.items():
        for k in d[s]:
            blocks.setdefault((s, k[2][0][0] // BLK), None)
bl = sorted(blocks)
idx = {b: i for i, b in enumerate(bl)}
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
draws = {n: [] for n in sets}; diff = {n: [] for n in sets if n != 'NTv2_C1'}
for _ in range(2000):
    ii = rng.integers(0, len(bl), len(bl))
    f = {n: f1(c[0][ii], c[1][ii], c[2][ii]) for n, c in C.items()}
    for n in sets:
        draws[n].append(f[n])
    for n in diff:
        diff[n].append(f['NTv2_C1'] - f[n])
for n in sets:
    res[n]['CI95'] = [float(np.percentile(draws[n], 2.5)), float(np.percentile(draws[n], 97.5))]
for n in diff:
    res[n]['diff_NTv2_minus_this_CI95'] = [float(np.percentile(diff[n], 2.5)), float(np.percentile(diff[n], 97.5))]
    res[n]['diff_point'] = res['NTv2_C1']['F1'] - res[n]['F1']
res['blocks'] = len(bl)
json.dump(res, open(outp, 'w'), indent=1)
print(json.dumps(res))
