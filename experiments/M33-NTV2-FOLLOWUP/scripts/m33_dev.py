"""M33 step 1b: DEV report of the tau*/thr* posterior filter chosen on the internal validation (m33_select.py).  DEV is never used for selection.
Reports unfiltered / filtered pooled + per-species + strata exact-chain F1, CDS-base F1 and background burden (same functions as M31), AUROC / AP / ECE with
reliability bins for tau* (and for any other tau present, descriptively), 1-Mb block bootstrap (as M31 `c1_boot.py`) of filtered F1 and of the differences to
unfiltered B, fresh-joint (GENERanno) and Carbon-A filtered.
Usage: m33_dev.py <dev_run_dir> <dev_post_dir> <selection.json> <out_json>
"""
import json, sys
import numpy as np
import m33_common as C
import a1_eval as E
import a1_compare as AC
from m29_r3_posterior_filter import curve

run, post, selp, outp = sys.argv[1:5]
sel = json.load(open(selp))['selection']
tau, thr = sel['tau'], sel['thr']
seqs = list(C.DEV)
chains = {s: E.load_chains(run, s) for s in seqs}
refs = AC.refs


def score(ch):
    r = E.summarize(ch, refs, seqs)
    r['base'] = AC.base_metrics(ch, seqs)
    return r


res = dict(selection=sel, unfiltered=score(chains))
lp = C.lp_by_key(post, tau, seqs)
filt = {s: {k: v for k, v in chains[s].items() if lp[k] >= thr} for s in seqs}
res['filtered'] = score(filt)
res['filtered']['kept_fraction'] = sum(len(v) for v in filt.values()) / sum(len(v) for v in chains.values())
# calibration / ranking at every tau that has DEV posteriors (tau* is the pre-selected one; others are descriptive)
res['calibration'] = {}
for t in C.TAU_GRID:
    try:
        l2 = C.lp_by_key(post, t, seqs)
    except FileNotFoundError:
        continue
    rows = [(l2[k], int(k in refs.keys[s]), s) for s in seqs for k in chains[s]]
    sc = np.array([r[0] for r in rows]); lb = np.array([r[1] for r in rows]); sq = np.array([r[2] for r in rows])
    cal = C.calibration(np.exp(sc), lb)
    cal['per_species'] = {s: {k: v for k, v in C.calibration(np.exp(sc[sq == s]), lb[sq == s]).items() if k != 'bins'} for s in seqs}
    nref = sum(len(refs.keys[s]) for s in seqs)
    tt, kept, tp, f1 = curve(sc, lb, nref)
    cal['oracle_best_F1_on_DEV_not_for_selection'] = 100 * float(f1.max())
    cal['mean_P'] = float(np.exp(sc).mean())
    res['calibration'][f'{t:g}'] = cal
    print(f'DEV tau={t:g} ECE={cal["ECE"]:.4f} AUROC={cal["AUROC"]:.4f} AP={cal["AP"]:.4f}', flush=True)
# 1-Mb block bootstrap, same blocks and seed as M31 c1_boot.py
BLK = 1_000_000
sets = {'B_unfiltered': chains, 'B_filtered': filt,
        'GENERanno_fresh_joint': {s: AC.gff_chains(AC.DEVGFF['B_freshjoint'][s], s) for s in seqs},
        'CarbonA_filtered': {s: AC.carbon_dev(s, True) for s in seqs}}
blocks = {}
for s in seqs:
    for k in refs.keys[s]:
        blocks.setdefault((s, k[2][0][0] // BLK), None)
    for d in sets.values():
        for k in d[s]:
            blocks.setdefault((s, k[2][0][0] // BLK), None)
bl = sorted(blocks)
idx = {b: i for i, b in enumerate(bl)}


def counts(d):
    tp = np.zeros(len(bl)); pr = np.zeros(len(bl)); rf = np.zeros(len(bl))
    for s in seqs:
        for k in d[s]:
            i = idx[(s, k[2][0][0] // BLK)]; pr[i] += 1; tp[i] += k in refs.keys[s]
        for k in refs.keys[s]:
            rf[idx[(s, k[2][0][0] // BLK)]] += 1
    return tp, pr, rf


Cn = {n: counts(d) for n, d in sets.items()}
f1f = lambda tp, p, r: 200 * tp.sum() / (p.sum() + r.sum())
rng = np.random.default_rng(0)
draws = {n: [] for n in sets}; diffs = {n: [] for n in ('B_unfiltered', 'GENERanno_fresh_joint', 'CarbonA_filtered')}
for _ in range(2000):
    ii = rng.integers(0, len(bl), len(bl))
    f = {n: f1f(c[0][ii], c[1][ii], c[2][ii]) for n, c in Cn.items()}
    for n in sets:
        draws[n].append(f[n])
    for n in diffs:
        diffs[n].append(f['B_filtered'] - f[n])
res['bootstrap'] = {n: dict(F1=f1f(*Cn[n]), CI95=[float(np.percentile(draws[n], 2.5)), float(np.percentile(draws[n], 97.5))]) for n in sets}
for n in diffs:
    res['bootstrap'][n]['diff_Bfiltered_minus_this'] = dict(point=f1f(*Cn['B_filtered']) - f1f(*Cn[n]), CI95=[float(np.percentile(diffs[n], 2.5)), float(np.percentile(diffs[n], 97.5))])
res['bootstrap']['blocks'] = len(bl)
res['reference_numbers_existing'] = dict(Helixer=62.851, Tiberius=71.805, ANNEVO=75.703, B_unfiltered_M31=63.61)
json.dump(res, open(outp, 'w'), indent=1)
u, f = res['unfiltered'], res['filtered']
print('DEV unfiltered', round(u['pooled']['F1'], 3), {s: round(d['F1'], 2) for s, d in u['by_species'].items()}, {k: round(v, 2) for k, v in u['base'].items()})
print('DEV filtered  ', round(f['pooled']['F1'], 3), {s: round(d['F1'], 2) for s, d in f['by_species'].items()}, {k: round(v, 2) for k, v in f['base'].items()}, 'kept', round(f['kept_fraction'], 3))
print('strata unf', {a: round(b['F1'], 2) for a, b in u['by_strata'].items()}, 'filt', {a: round(b['F1'], 2) for a, b in f['by_strata'].items()})
print(json.dumps(res['bootstrap']))
