"""M34 P2b scoring: zero-shot exact-chain F1 of two models on the M30 time-held-out evaluation set B (Dryas octopetala NC_141367.1/69/71/73, Eschscholzia californica
NW_028238309.1/461.1/483.1) and the pre-registered paired 1-Mb block bootstrap of the T10 - T2 difference (protocol: docs/review_20261016/19_M34_scaleup_plan.md section 3.1, bullet P2b).

Ruler: `m30_score.Evaluator` (M26 CDS-assessable primaries, exact CDS chain) with `SC.load_gff_chains` (unique exact chains, CDS >= 96 bp, NO posterior filter), i.e. the M33
`m33_evalB.py` counting.  Bootstrap statistic: the DEV one of `c1_boot.py` / `m34_boot_pair.py` / M35 `c1_boot.py` -- blocks (contig, first-CDS-start // 1 Mb) taken from the
reference chains and from every given run's chains, 2,000 draws of `rng.integers(0, n_blocks, n_blocks)` with `default_rng(0)`, F1 = 200 tp / (pred + ref) of the resampled
blocks; the same resampled indices serve every run and every scope (pooled, per species: the species F1 uses the resampled blocks that lie in that species), 95% CI = 2.5/97.5
percentiles.  Criterion P2b: pooled diff >= +2.0 and pooled CI lower bound > 0, and both species' point differences > 0.

Usage: m34_p2b_score.py --new NAME=DIR[,DIR...] --old NAME=DIR[,DIR...] [--extra NAME=DIR[,DIR...] ...] --out RESULT.json [--label P2b] [--z-arm NAME]
  DIR = a run directory with <contig>/prediction/full_grammar.gff3 (m34_infer layout); several DIRs are allowed when the contigs are spread over several run directories
  (M33 ran one directory per species).  Blocks are the union of the reference and ALL given runs.  Diff = F1(new) - F1(old) (and new - extra, descriptive).
Run with the generanno env (the M30 ruler's environment), e.g.  taskset -c 25-31 nice -n 10 ~/envs/generanno/bin/python m34_p2b_score.py ...
"""
import argparse
import glob
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1]
_M30 = str(EXP / 'M30-HOLDOUT-ADAPT' / 'scripts')
if _M30 not in sys.path:
    sys.path.insert(0, _M30)
import m30_common as M30                      # noqa: E402
import m30_score as SC                        # noqa: E402  (registers the M30 contigs in the M29 helpers)

BLK = 1_000_000
DRAWS, SEED = 2000, 0
DIFF_MIN, SPECIES_MIN = 2.0, 0.0              # pre-registered: pooled diff >= +2.0, CI lower bound > 0, both species' point diffs > 0
SPECIES = {'dryas_octopetala': ['NC_141367.1', 'NC_141369.1', 'NC_141371.1', 'NC_141373.1'],
           'eschscholzia_californica': ['NW_028238309.1', 'NW_028238461.1', 'NW_028238483.1']}
CONTIGS = [c for cs in SPECIES.values() for c in cs]
for _sp, _cs in SPECIES.items():              # the registered contig lists must be exactly the M30 evaluation set B
    assert sorted(_cs) == sorted(M30.SPECIES[_sp]['B']), (_sp, _cs, M30.SPECIES[_sp]['B'])
SCOPES = ['pooled'] + list(SPECIES)


# ---- loading / counting ------------------------------------------------------------------------------------------------------------------------
def gff_path(run, contig):
    return Path(run) / contig / 'prediction' / 'full_grammar.gff3'


def load_run(dirs, contigs=CONTIGS):
    """{contig: {key: dict(nex, exons, strand)}} -- each contig must be found in exactly one of the run directories."""
    out = {}
    for c in contigs:
        hit = [d for d in dirs if gff_path(d, c).exists()]
        assert len(hit) == 1, f'{c}: GFF3 found in {len(hit)} of {list(dirs)}'
        out[c] = SC.load_gff_chains(str(gff_path(hit[0], c)), c)
    return out


def provenance(dirs):
    """model / step / throughput recorded by m34_infer / m33_infer in the DONE_*.json of each run directory."""
    rows = []
    for d in dirs:
        for f in sorted(glob.glob(str(Path(d) / 'DONE_*.json'))):
            j = json.load(open(f))
            rows.append({k: j.get(k) for k in ('model', 'step', 'backend', 'regions', 'wall_seconds', 'encode_seconds', 'completed')} | {'file': f})
    return rows


def block_of(contig, key):
    return (contig, key[2][0][0] // BLK)          # key = (seqid, strand, ((a, b), ...)): first CDS interval start


def f1(tp, pred, ref):
    return 200 * tp / max(pred + ref, 1)


class Scorer:
    """One Evaluator for the seven B contigs + the block count vectors of every added run."""

    def __init__(self, contigs=CONTIGS):
        self.contigs = list(contigs)
        self.ev = SC.Evaluator(self.contigs)
        self.chains, self.dirs = {}, {}

    def add(self, name, dirs):
        dirs = list(dirs)
        self.chains[name] = load_run(dirs, self.contigs)
        self.dirs[name] = dirs
        self.ev.add(name, self.chains[name])

    def species_counts(self, name):
        per = self.ev.methods[name]['per']
        out = {}
        for sp, cs in SPECIES.items():
            out[sp] = {k: sum(per[c][k] for c in cs if c in per) for k in ('tp', 'pred', 'ref')}
        out['pooled'] = {k: sum(v[k] for v in out.values()) for k in ('tp', 'pred', 'ref')}
        return out

    def table(self, name):
        """per-scope tp/pred/ref/F1/P/R + per contig, from the M30 Evaluator counting"""
        sc = self.species_counts(name)
        per = self.ev.methods[name]['per']
        return dict({s: dict(v, F1=f1(v['tp'], v['pred'], v['ref']), precision=100 * v['tp'] / max(v['pred'], 1), recall=100 * v['tp'] / max(v['ref'], 1)) for s, v in sc.items()},
                    per_contig={c: dict(p, F1=f1(p['tp'], p['pred'], p['ref'])) for c, p in per.items()})

    # ---- blocks ----------------------------------------------------------------------------------------------------------------------------------
    def blocks(self):
        bl = set()
        for c in self.contigs:
            bl.update(block_of(c, k) for k in self.ev.refkeys[c])
        for ch in self.chains.values():
            for c in self.contigs:
                bl.update(block_of(c, k) for k in ch[c])
        return sorted(bl)

    def block_counts(self, bl):
        idx = {b: i for i, b in enumerate(bl)}
        rf = np.zeros(len(bl))
        for c in self.contigs:
            for k in self.ev.refkeys[c]:
                rf[idx[block_of(c, k)]] += 1
        C = {}
        for n, ch in self.chains.items():
            tp, pr = np.zeros(len(bl)), np.zeros(len(bl))
            for c in self.contigs:
                for k in ch[c]:
                    i = idx[block_of(c, k)]
                    pr[i] += 1
                    tp[i] += k in self.ev.refkeys[c]
            C[n] = (tp, pr, rf)
            t = self.species_counts(n)['pooled']                  # the block counts must add up to the Evaluator counts
            assert (int(tp.sum()), int(pr.sum()), int(rf.sum())) == (t['tp'], t['pred'], t['ref']), (n, (tp.sum(), pr.sum(), rf.sum()), t)
        return C


def bootstrap(C, bl, ref_name, draws=DRAWS, seed=SEED):
    """C: {run: (tp, pred, ref) per block}.  -> {run: {scope: {F1, CI95}}}, {run: {scope: {point, CI95, frac_positive}}} for ref_name - run (every run != ref_name)."""
    bl_sp = np.array([next(i for i, (sp, cs) in enumerate(SPECIES.items()) if b[0] in cs) for b in bl])
    scopes = {'pooled': None, **{sp: i for i, sp in enumerate(SPECIES)}}

    def stat(c, ii, species=None):
        tp, p, r = (x[ii] for x in c)
        if species is not None:
            m = bl_sp[ii] == species
            tp, p, r = tp[m], p[m], r[m]
        return 200 * tp.sum() / max(p.sum() + r.sum(), 1)

    full = np.arange(len(bl))
    point = {n: {sc: stat(c, full, sp) for sc, sp in scopes.items()} for n, c in C.items()}
    rng = np.random.default_rng(seed)
    draw = {n: {sc: [] for sc in scopes} for n in C}
    for _ in range(draws):
        ii = rng.integers(0, len(bl), len(bl))                       # same blocks / same resample for every run and scope
        for n, c in C.items():
            for sc, sp in scopes.items():
                draw[n][sc].append(stat(c, ii, sp))
    ci = lambda v: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
    runs = {n: {sc: dict(F1=float(point[n][sc]), CI95=ci(draw[n][sc])) for sc in scopes} for n in C}
    diffs = {}
    for n in C:
        if n == ref_name:
            continue
        diffs[n] = {}
        for sc in scopes:
            d = np.array(draw[ref_name][sc]) - np.array(draw[n][sc])
            diffs[n][sc] = dict(point=float(point[ref_name][sc] - point[n][sc]), CI95=ci(d), frac_positive=float((d > 0).mean()))
    return runs, diffs


def verdict(d):
    """d = diffs[old] (T10 - T2): the pre-registered P2b criterion"""
    checks = dict(pooled_diff_ge_2=bool(d['pooled']['point'] >= DIFF_MIN), pooled_ci_lower_gt_0=bool(d['pooled']['CI95'][0] > 0),
                  both_species_point_positive=bool(all(d[sp]['point'] > SPECIES_MIN for sp in SPECIES)))
    return dict(criterion=f'pooled diff >= +{DIFF_MIN} and pooled CI lower bound > 0 and both species point differences > 0', checks=checks, met=bool(all(checks.values())))


def run_scoring(new, old, extra=(), draws=DRAWS, seed=SEED, label='P2b', z_arm=None):
    """new / old: (name, [dirs]); extra: [(name, [dirs])].  -> result dict"""
    t0 = time.time()
    S = Scorer()
    for n, dirs in [new, old, *extra]:
        S.add(n, dirs)
    bl = S.blocks()
    C = S.block_counts(bl)
    runs, diffs = bootstrap(C, bl, new[0], draws, seed)
    res = dict(label=label, protocol='docs/review_20261016/19_M34_scaleup_plan.md section 3.1 P2b; 13_M30_protocol.md section 2 (evaluation set B)',
               ruler='m30_score.Evaluator + load_gff_chains (exact CDS chain, CDS>=96 bp, no posterior filter)',
               bootstrap=dict(block_bp=BLK, blocks=len(bl), draws=draws, seed=seed, statistic='200 tp/(pred+ref) of resampled 1-Mb blocks; same indices for all runs/scopes; 95% CI = 2.5/97.5 percentiles'),
               contigs=SPECIES, new=new[0], old=old[0], runs={})
    for n, dirs in [new, old, *extra]:
        res['runs'][n] = dict(dirs=dirs, provenance=provenance(dirs), table=S.table(n), bootstrap_F1={sc: runs[n][sc] for sc in SCOPES})
    res['paired_new_minus'] = diffs
    res['verdict'] = dict(verdict(diffs[old[0]]), name=f'{new[0]} - {old[0]}', point={sc: diffs[old[0]][sc]['point'] for sc in SCOPES},
                          CI95={sc: diffs[old[0]][sc]['CI95'] for sc in SCOPES})
    if z_arm:
        res['m36_zero_shot_arm_Z'] = dict(name=z_arm, run_dirs=res['runs'][z_arm]['dirs'], note='zero-shot B-set output of the M36 base model, reusable as the M36 Z arm')
    res['score_seconds'] = round(time.time() - t0, 1)
    return res


def parse_run(s):
    name, _, dirs = s.partition('=')
    assert name and dirs, f'expected NAME=DIR[,DIR...], got {s!r}'
    return name, dirs.split(',')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--new', required=True)
    ap.add_argument('--old', required=True)
    ap.add_argument('--extra', action='append', default=[])
    ap.add_argument('--out', required=True)
    ap.add_argument('--label', default='P2b')
    ap.add_argument('--z-arm', default=None)
    args = ap.parse_args()
    res = run_scoring(parse_run(args.new), parse_run(args.old), [parse_run(e) for e in args.extra], label=args.label, z_arm=args.z_arm)
    tmp = args.out + '.tmp'
    json.dump(res, open(tmp, 'w'), indent=1)
    Path(tmp).replace(args.out)
    for n, r in res['runs'].items():
        t = r['table']
        print(f"{n}: " + ' | '.join(f"{sc} F1 {t[sc]['F1']:.2f} (tp {t[sc]['tp']} pred {t[sc]['pred']} ref {t[sc]['ref']}) CI [{r['bootstrap_F1'][sc]['CI95'][0]:.2f},{r['bootstrap_F1'][sc]['CI95'][1]:.2f}]" for sc in SCOPES))
    for n, d in res['paired_new_minus'].items():
        print(f"{res['new']} - {n}: " + ' | '.join(f"{sc} {d[sc]['point']:+.2f} [{d[sc]['CI95'][0]:+.2f},{d[sc]['CI95'][1]:+.2f}]" for sc in SCOPES))
    v = res['verdict']
    print(json.dumps(dict(label=res['label'], verdict=v['name'], met=v['met'], checks=v['checks'], blocks=res['bootstrap']['blocks'])))


if __name__ == '__main__':
    main()
