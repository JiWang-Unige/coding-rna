"""Paired 1-Mb block bootstrap (2000 reps, seed 0; same blocks as M31 `c1_boot.py`) of pooled DEV exact-chain F1: CI of every run, and of the paired differences
of the first run to all the others (incl. the existing GENERanno fresh-joint / Carbon-A-filtered DEV predictions).
Usage: m34_boot_pair.py <out.json> NAME=run_dir [NAME=run_dir ...]      (first = reference run, e.g. T10=.../T10/devM34_T10  T2=.../T2/devM34_T2  M31B=.../devB_final)"""
import json, sys
import numpy as np
import m34_train_common as T       # sets sys.path (M31 scripts, M29 helpers)
import m31_common as M
import a1_common as K
sys.path.insert(0, str(K.CA_SRC))
import a1_eval as E                # noqa: E402
import a1_compare as AC            # noqa: E402

outp, specs = sys.argv[1], [a.split('=', 1) for a in sys.argv[2:]]
seqs = list(M.DEV_SEQS)
BLK = 1_000_000
sets = {n: {s: E.load_chains(p, s) for s in seqs} for n, p in specs}
sets['GENERanno_fresh_joint'] = {s: AC.gff_chains(AC.DEVGFF['B_freshjoint'][s], s) for s in seqs}
sets['CarbonA_filtered'] = {s: AC.carbon_dev(s, True) for s in seqs}
blocks = {}
for s in seqs:
    for k in AC.refs.keys[s]:
        blocks.setdefault((s, k[2][0][0] // BLK), None)
    for d in sets.values():
        for k in d[s]:
            blocks.setdefault((s, k[2][0][0] // BLK), None)
bl = sorted(blocks)
idx = {b: i for i, b in enumerate(bl)}


def counts(name):
    tp, pr, rf = np.zeros(len(bl)), np.zeros(len(bl)), np.zeros(len(bl))
    for s in seqs:
        for k in sets[name][s]:
            i = idx[(s, k[2][0][0] // BLK)]; pr[i] += 1; tp[i] += k in AC.refs.keys[s]
        for k in AC.refs.keys[s]:
            rf[idx[(s, k[2][0][0] // BLK)]] += 1
    return tp, pr, rf


C = {n: counts(n) for n in sets}
f1 = lambda tp, p, r: 200 * tp.sum() / (p.sum() + r.sum())
ref = specs[0][0]
rng = np.random.default_rng(0)
draws = {n: [] for n in sets}
diff = {n: [] for n in sets if n != ref}
for _ in range(2000):
    ii = rng.integers(0, len(bl), len(bl))
    f = {n: f1(c[0][ii], c[1][ii], c[2][ii]) for n, c in C.items()}
    for n in sets:
        draws[n].append(f[n])
    for n in diff:
        diff[n].append(f[ref] - f[n])
res = {n: {'F1': f1(*c), 'CI95': [float(np.percentile(draws[n], 2.5)), float(np.percentile(draws[n], 97.5))]} for n, c in C.items()}
for n in diff:
    res[n]['diff_ref_minus_this_point'] = res[ref]['F1'] - res[n]['F1']
    res[n]['diff_ref_minus_this_CI95'] = [float(np.percentile(diff[n], 2.5)), float(np.percentile(diff[n], 97.5))]
res['reference_run'] = ref
res['blocks'] = len(bl)
json.dump(res, open(outp, 'w'), indent=1)
print(json.dumps(res))
