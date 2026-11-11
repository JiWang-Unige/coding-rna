"""M41 pre-registered DEV comparison (protocol docs/review_20261016/33 section 4): paired 1-Mb block bootstrap of the exact-chain F1 of the M41 arm vs T2,
pooled AND per exon-count stratum (1 / 2-5 / >=6), with the T2 seed replicates and Carbon-A (filtered) as descriptive references, and the verdict.
Same blocks / draws / statistic as `M34-NTV2-SCALEUP/scripts/m34_boot_pair.py` (blocks = (chromosome, first-CDS-start // 1 Mb) over the reference and every
run, 2,000 draws of rng.integers with default_rng(0), F1 = 200 tp / (pred + ref)); a stratum's F1 uses the stratum's tp / pred / ref of the resampled blocks
(prediction stratum = its own exon count, reference stratum = the reference's, as c1_eval / a1_eval).  M34 / M29 files are imported, not edited.
Usage (ont, generanno env, from M34-NTV2-SCALEUP/scripts):  python m41_compare.py <out.json> T2_ILS=<devM41 dir> T2=<pc2_T2/devM34_T2> T2_s2=<...> T2_s3=<...>
(first = the M41 arm, second = the paired reference T2; the rest are descriptive)
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'M34-NTV2-SCALEUP' / 'scripts'))
import m34_train_common as T      # noqa: E402,F401  (sys.path: M31 scripts, M29 helpers)
import m31_common as M            # noqa: E402
import a1_common as K             # noqa: E402
sys.path.insert(0, str(K.CA_SRC))
import a1_eval as E               # noqa: E402
import a1_compare as AC           # noqa: E402

BLK = 1_000_000
STRATA = ('1', '2-5', '>=6')
# pre-registered thresholds (protocol 33 section 4.2): seed SDs of T2 (seeds 0/1/2 = T2, T2_s2, T2_s3)
SD_POOLED, SD_GE6, SD_SETB = 1.35, 1.29, 1.29


def strat(n):
    return '1' if n == 1 else '2-5' if n <= 5 else '>=6'


def verdict(d_ge6, ci_ge6, d_pooled):
    """Pre-registered decision rule on DEV (protocol 33 section 4.2).  Set-B guard is applied by the caller."""
    if d_ge6 >= 2 * SD_GE6 and ci_ge6[0] > 0 and d_pooled > -SD_POOLED:
        return 'supported'
    if d_ge6 <= -2 * SD_GE6 or d_pooled <= -2 * SD_POOLED:
        return 'harmful'
    if d_ge6 >= SD_GE6 and ci_ge6[0] > 0:
        return 'grey_zone'
    return 'no_detectable_effect'


def main():
    outp, specs = sys.argv[1], [a.split('=', 1) for a in sys.argv[2:]]
    seqs = list(M.DEV_SEQS)
    sets = {n: {s: E.load_chains(p, s) for s in seqs} for n, p in specs}
    sets['CarbonA_filtered'] = {s: {k: dict(nex=len(k[2])) for k in AC.carbon_dev(s, True)} for s in seqs}
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
        c = {z: np.zeros((3, len(bl))) for z in STRATA}           # tp, pred, ref per block
        for s in seqs:
            for k, v in sets[name][s].items():
                z = strat(v['nex'])
                i = idx[(s, k[2][0][0] // BLK)]
                c[z][1, i] += 1
                c[z][0, i] += k in AC.refs.keys[s]
            for k in AC.refs.keys[s]:
                c[strat(AC.refs.nex[k])][2, idx[(s, k[2][0][0] // BLK)]] += 1
        c['pooled'] = sum(c[z] for z in STRATA)
        return c

    C = {n: counts(n) for n in sets}
    f1 = lambda m, ii=slice(None): 200 * m[0, ii].sum() / max(1., m[1, ii].sum() + m[2, ii].sum())
    scopes = ('pooled',) + STRATA
    arm, ref = specs[0][0], specs[1][0]
    point = {n: {z: f1(C[n][z]) for z in scopes} for n in sets}
    rng = np.random.default_rng(0)
    diffs = {z: [] for z in scopes}
    draws = {n: {z: [] for z in scopes} for n in sets}
    for _ in range(2000):
        ii = rng.integers(0, len(bl), len(bl))
        for n in sets:
            for z in scopes:
                draws[n][z].append(f1(C[n][z], ii))
        for z in scopes:
            diffs[z].append(draws[arm][z][-1] - draws[ref][z][-1])
    ci = lambda v: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
    res = dict(arm=arm, paired_reference=ref, blocks=len(bl),
               runs={n: {z: dict(F1=point[n][z], CI95=ci(draws[n][z])) for z in scopes} for n in sets},
               paired_diff={z: dict(point=point[arm][z] - point[ref][z], CI95=ci(diffs[z])) for z in scopes})
    seeds = [n for n, _ in specs[1:]]
    for z in scopes:
        v = np.array([point[n][z] for n in seeds])
        res.setdefault('t2_seed_summary', {})[z] = dict(runs=seeds, mean=float(v.mean()), sd=float(v.std(ddof=1)) if len(v) > 1 else None,
                                                         z_arm=float((point[arm][z] - v.mean()) / v.std(ddof=1)) if len(v) > 1 else None)
    res['carbonA_gap_ge6'] = dict(T2=point['CarbonA_filtered']['>=6'] - point[ref]['>=6'], arm=point['CarbonA_filtered']['>=6'] - point[arm]['>=6'])
    pd = res['paired_diff']
    res['dev_verdict'] = verdict(pd['>=6']['point'], pd['>=6']['CI95'], pd['pooled']['point'])
    res['rule'] = ('supported: d(>=6) >= 2 x 1.29 and paired CI(>=6) lower > 0 and d(pooled) > -1.35;  harmful: d(>=6) <= -2.58 or d(pooled) <= -2.70;  '
                   'grey_zone: 1.29 <= d(>=6) < 2.58 with CI lower > 0;  otherwise no_detectable_effect.  Set-B guard (protocol 33 4.2) applied in the report.')
    Path(outp).write_text(json.dumps(res, indent=1))
    print(json.dumps(dict(paired_diff=pd, dev_verdict=res['dev_verdict'])))


if __name__ == '__main__':
    main()
