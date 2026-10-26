"""M36 final scoring on evaluation set B (protocol 22 §4-§5).  The ONLY place that reads the B reference annotation; one-shot (sentinels, protocol 22 §8).
Ruler: `m30_score.Evaluator` (M26 CDS-assessable primaries, exact CDS chain) with `SC.load_gff_chains` (unique exact chains, CDS >= 96 bp, NO posterior filter) -- the M33
`m33_evalB.py` / M34 P2b counting.  Per species (a model is adapted and evaluated on ONE species): Z = the P2b T10 zero-shot run, the trained arms S1 / S2 / S3 / F1 at step 512,
Carbon-A filtered / unfiltered read from the M30 outputs (unchanged).
Statistic (user decision m36-ambiguities, 2026-10-26, before any M36 output; protocol 22 section 9): the section 5 verdict uses the M30 / M33 gene-unit paired bootstrap
(`Evaluator.bootstrap`: reference genes + 100-kb orphan bins, 2,000 draws, seed 0) -- the implementation protocol 22 section 4 names.  The paired 1-Mb block bootstrap of the M34 P2b /
DEV `c1_boot.py` (block = (contig, first-CDS-start // 1 Mb), union of the reference and every run's chains, 2,000 draws of `default_rng(0).integers(0, n, n)`) is reported next to it
(`verdict_block_bootstrap`); if the two give different labels the paper reports the result as sensitive to the statistic.
closure ratio = (F1_arm - F1_Z) / (F1_CarbonA_filt - F1_Z)  (per species; F1_CarbonA_filt = M30 value measured on the same B contigs, Dryas 73.98 / Esch 73.32).
Verdict (protocol 22 §5; see `verdict`):  GO-strong: S2 or S3 within 2.0 points of Carbon-A filtered in BOTH species with a paired CI that does not exclude parity;  GO: S1 or S3 with closure
>= 30% and a paired (arm - Z) CI lower bound > 0 in BOTH species;  NO-GO: some species in which no arm closes 30%;  C2 retracted: F1 - S1 >= -2.0 points.
usage: m36_score.py main  [--out RESULT.json] [--root OUT]          one-shot scoring of the unfiltered exact-chain F1 (verdict)
       m36_score.py post  [--post-root DIR] [--out RESULT.json]     one-shot descriptive scoring of the posterior-filtered chains (frozen thr, tau) + calibration
       m36_score.py carbonA                                         reproduce the M30 Carbon-A B numbers (sanity; reads Carbon-A predictions only)
Run with the generanno env:  taskset -c 18-23 nice -n 10 ~/envs/generanno/bin/python m36_score.py main
"""
import argparse
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

import m36_common as W
import m36_sentinel as SN
import m30_common as M30
import m30_score as SC
import m34_p2b_score as PS

BLK, DRAWS, SEED = PS.BLK, PS.DRAWS, PS.SEED
B = {sp: list(M30.SPECIES[sp]['B']) for sp in W.SPECIES}
CARBON_A_DIR = M30.OUT / 'carbonA_raw'
CARBON_A_FILT_EXPECTED = {'dryas_octopetala': 73.98, 'eschscholzia_californica': 73.32}        # M30 / M33 results (results_<species>.json, 2 decimals)
CLOSURE_GO, STRONG_MARGIN, C2_MARGIN = 0.30, 2.0, -2.0                                           # protocol 22 §5 (fixed before any M36 output)
S_ARMS = ('S1', 'S2', 'S3')
ALL_ARMS = ('S1', 'S2', 'S3', 'F1')


def ci(v):
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]


# ---- runs ---------------------------------------------------------------------------------------------------------------------------------
def carbon_chains(species, filtered=True):
    return {c: SC.load_gff_chains(str(CARBON_A_DIR / species / c / 'prediction.gff'), c, cds96=False, conf_filter='carbon' if filtered else None) for c in B[species]}


def gff_chains(run_dir, species):
    return {c: SC.load_gff_chains(str(PS.gff_path(run_dir, c)), c) for c in B[species]}


class SpeciesScorer:
    """One Evaluator for the B contigs of one species + per-block count vectors of every added run."""

    def __init__(self, species):
        self.species, self.contigs = species, list(B[species])
        self.ev = SC.Evaluator(self.contigs)
        self.chains = {}

    def add(self, name, chains):
        self.chains[name] = chains
        return self.ev.add(name, chains)

    def blocks(self):
        bl = set()
        for c in self.contigs:
            bl.update(PS.block_of(c, k) for k in self.ev.refkeys[c])
        for ch in self.chains.values():
            for c in self.contigs:
                bl.update(PS.block_of(c, k) for k in ch[c])
        return sorted(bl)

    def block_counts(self, bl):
        idx = {b: i for i, b in enumerate(bl)}
        rf = np.zeros(len(bl))
        for c in self.contigs:
            for k in self.ev.refkeys[c]:
                rf[idx[PS.block_of(c, k)]] += 1
        out = {}
        for n, ch in self.chains.items():
            tp, pr = np.zeros(len(bl)), np.zeros(len(bl))
            for c in self.contigs:
                for k in ch[c]:
                    i = idx[PS.block_of(c, k)]
                    pr[i] += 1
                    tp[i] += k in self.ev.refkeys[c]
            t = self.ev.table(n)['pooled']
            assert (int(tp.sum()), int(pr.sum()), int(rf.sum())) == (t['tp'], t['pred'], t['ref']), (n, (tp.sum(), pr.sum(), rf.sum()), t)
            out[n] = (tp, pr, rf)
        return out


def bootstrap_draws(counts, n_blocks, draws=DRAWS, seed=SEED):
    """counts: {run: (tp, pred, ref) per block}.  -> (point F1 {run}, draws {run: array})  -- the resample loop of m34_p2b_score.bootstrap (one `integers(0, n, n)` call per draw)."""
    full = np.arange(n_blocks)
    stat = lambda c, ii: 200 * c[0][ii].sum() / max(c[1][ii].sum() + c[2][ii].sum(), 1)
    point = {n: float(stat(c, full)) for n, c in counts.items()}
    rng = np.random.default_rng(seed)
    out = {n: np.empty(draws) for n in counts}
    for d in range(draws):
        ii = rng.integers(0, n_blocks, n_blocks)
        for n, c in counts.items():
            out[n][d] = stat(c, ii)
    return point, out


def pair(point, draws, a, b):
    d = draws[a] - draws[b]
    return dict(point=point[a] - point[b], CI95=ci(d), frac_positive=float((d > 0).mean()))


def closure(f_arm, f_z, f_ca):
    den = f_ca - f_z
    return (f_arm - f_z) / den if den > 0 else float('nan')


# ---- verdict (protocol 22 §5), pure function of the per-species numbers ---------------------------------------------------------------------
def verdict(S, f1_minus_s1=None, species=W.SPECIES):
    """S[sp][arm] = dict(diff_vs_Z=dict(point, CI95), diff_vs_CA=dict(point, CI95), closure=float).  f1_minus_s1[sp] = dict(point, CI95) (only species with an F1 control).
    -> dict of flags and the label.  Conventions (the §5 text leaves them open; recorded in the result):
      * 'CI does not exclude parity' (GO-strong) = the upper bound of the paired (arm - Carbon-A filtered) CI is >= 0;
      * GO: one arm of {S1, S3} meets (closure >= 30% and paired (arm - Z) CI lower bound > 0) in every species;
      * NO-GO: there is a species in which NO arm of {S1, S2, S3} reaches closure 30% (also reported: the stricter reading restricted to {S1, S3});
      * C2 retracted if F1 - S1 >= -2.0 points in at least one species with an F1 control."""
    def has(sp, arm):
        return sp in S and arm in S[sp]

    strong = {arm: {sp: (S[sp][arm]['diff_vs_CA']['point'] >= -STRONG_MARGIN and S[sp][arm]['diff_vs_CA']['CI95'][1] >= 0) if has(sp, arm) else None for sp in species}
              for arm in ('S2', 'S3')}
    go_arm = {arm: {sp: (S[sp][arm]['closure'] >= CLOSURE_GO and S[sp][arm]['diff_vs_Z']['CI95'][0] > 0) if has(sp, arm) else None for sp in species} for arm in ('S1', 'S3')}
    go_strong = any(all(v is True for v in per.values()) for per in strong.values())
    go = any(all(v is True for v in per.values()) for per in go_arm.values())
    best_closure = {sp: max((S[sp][a]['closure'] for a in S_ARMS if has(sp, a)), default=float('nan')) for sp in species}
    s13_closure = {sp: max((S[sp][a]['closure'] for a in ('S1', 'S3') if has(sp, a)), default=float('nan')) for sp in species}
    no_go = any(not (v >= CLOSURE_GO) for v in best_closure.values())                  # NaN (missing) counts as not reaching 30%
    no_go_s13 = any(not (v >= CLOSURE_GO) for v in s13_closure.values())
    c2 = {sp: (f1_minus_s1[sp]['point'] >= C2_MARGIN) for sp in (f1_minus_s1 or {})}
    label = 'GO-strong' if go_strong else 'GO' if go else 'NO-GO' if no_go else 'INCONCLUSIVE (between GO and NO-GO; not covered by the protocol 22 section 5 table)'
    return dict(label=label, go_strong=go_strong, go=go, no_go=no_go, no_go_reading_S1_S3_only=no_go_s13, c2_retracted=bool(c2) and any(c2.values()),
                c2_by_species=c2, go_strong_by_arm=strong, go_by_arm=go_arm, best_closure_S1_S2_S3=best_closure, best_closure_S1_S3=s13_closure,
                rules=dict(closure_go=CLOSURE_GO, strong_margin_points=STRONG_MARGIN, c2_margin_points=C2_MARGIN,
                           go_strong='S2 or S3: (arm - CarbonA_filt) point >= -2.0 and paired CI upper bound >= 0, both species',
                           go='S1 or S3: closure >= 0.30 and (arm - Z) paired CI lower bound > 0, both species',
                           no_go='a species where max closure over S1/S2/S3 < 0.30', c2='F1 - S1 point >= -2.0 (species with an F1 control)'))


# ---- inputs --------------------------------------------------------------------------------------------------------------------------------
def m30_reference_rows(species):
    """protocol 22 §4 'reference rows': M30 adaptation with real labels (k = 25% U, k = 25% F, k = 100%) on the same B contigs -- read from the M30 results, unchanged"""
    try:
        t = json.load(open(M30.OUT / f'results_{species}_final.json'))['tables']
        return {k: t[f'{species}_{k}']['pooled']['F1'] for k in ('k25_U', 'k25_F', 'k100_FULL')}
    except Exception as e:
        return f'not available: {e!r}'


def trained_arms(root):
    """[(arm, species, key, model_path, run_dir)] of the trained arms whose B inference finished"""
    out = []
    for arm, sp in sorted({(a, s) for a in ALL_ARMS for s in W.SPECIES}):
        key = W.arm_key(arm, sp)
        run = root / 'evalB' / key
        if list(run.glob('DONE_*.json')):
            out.append((arm, sp, key, root / 'train' / key / f"model_{W.TRAIN['main_step']:06d}.pt", run))
    return out


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()


def provenance_of(run_dir):
    return PS.provenance([str(run_dir)])


def sentinel_identity(root):
    ident = dict(models={}, runs={}, z=str(W.Z_DIR))
    for arm, sp, key, model, run in trained_arms(root):
        ident['models'][key] = md5(model)
        ident['runs'][key] = [(r['model'], r['regions']) for r in provenance_of(run)]
    ident['z_runs'] = [(r['model'], r['regions']) for r in provenance_of(W.Z_DIR)]
    return ident


def seed_coverage(root):
    """descriptive (protocol 22 §2 'covered reference-gene fraction'): per arm, seeds on the A training contigs vs the CDS-assessable reference primaries of those contigs"""
    import m36_train_inputs as TI
    out = {}
    M29 = None
    for sp in W.SPECIES:
        arms = [a for a in ALL_ARMS if (root / 'train' / W.arm_key(a, sp) / 'm36_arm.json').exists()]
        if not arms:
            continue
        M29 = M29 or M30.register()
        train = set(W.TARGETS[sp]['train'])
        refkeys = set()
        for c in W.TARGETS[sp]['train']:
            refkeys |= {M29.key(t) for t in M29.build_reference(c)['refs']}
        for arm in arms:
            inputs = TI.load_arm_inputs(arm, sp)
            keys = {(r['seqid'], r['strand'], tuple((a, b) for a, b, _ in r['CDS'])) for r in inputs['rows'] if r['seqid'] in train}
            hit = len(keys & refkeys)
            out[W.arm_key(arm, sp)] = dict(seeds_train_contigs=len(keys), reference_primaries_train_contigs=len(refkeys), exact_chain_hits=hit,
                                           seed_precision=hit / max(len(keys), 1), reference_recall=hit / max(len(refkeys), 1))
    return out


# ---- main scoring ----------------------------------------------------------------------------------------------------------------------------
def score_main(root, out_path, draws=DRAWS, seed=SEED, with_seed_coverage=True):
    t0 = time.time()
    arms = trained_arms(root)
    res = dict(label='M36', protocol='docs/review_20261016/22_M36_bootstrap_adaptation_protocol.md sections 4-5; 26_M36_cpu_preparation.md section 7',
               ruler='m30_score.Evaluator + load_gff_chains (exact CDS chain, CDS >= 96 bp, no posterior filter); Carbon-A: cds96=False, conf_filter=carbon',
               bootstrap=dict(kind='M30 / M33 gene-unit paired bootstrap (Evaluator.bootstrap, reference genes + 100-kb orphan bins), per species -- drives `verdict`',
                              draws=draws, seed=seed, decision='m36-ambiguities 2026-10-26 (protocol 22 section 9)',
                              alt='paired 1-Mb block bootstrap (m34_p2b_score / DEV c1_boot statistic) -- reported as `verdict_block_bootstrap`', block_bp=BLK),
               z_run=str(W.Z_DIR), species={}, arms_scored=[a[2] for a in arms])
    summary, f1_minus_s1, alt_summary = {}, {}, {}
    for sp in W.SPECIES:
        S = SpeciesScorer(sp)
        S.add('Z', gff_chains(W.Z_DIR, sp))
        S.add('CarbonA_filt', carbon_chains(sp, True))
        S.add('CarbonA_unf', carbon_chains(sp, False))
        mine = [a for a in arms if a[1] == sp]
        for arm, _, key, model, run in mine:
            S.add(arm, gff_chains(run, sp))
        bl = S.blocks()
        counts = S.block_counts(bl)
        point, dr = bootstrap_draws(counts, len(bl), draws, seed)
        names = list(counts)
        pairs = {f'{a} - {b}': pair(point, dr, a, b) for a in names for b in names if a != b}
        tabs = {n: S.ev.table(n) for n in names}
        # gene-unit (M30 / M33) bootstrap of the same differences, as a robustness check
        gp = [(a, 'Z') for a, *_ in mine] + [(a, 'CarbonA_filt') for a, *_ in mine] + ([('F1', 'S1')] if {'F1', 'S1'} <= {a[0] for a in mine} else [])
        g_ci, g_diffs = S.ev.bootstrap(gp, B=draws, seed=seed) if gp else ({}, {})
        f_ca, f_z = point['CarbonA_filt'], point['Z']
        arm_rows, alt_rows = {}, {}
        for arm, _, key, model, run in mine:
            d_z, d_ca = pairs[f'{arm} - Z'], pairs[f'{arm} - CarbonA_filt']
            den = f_ca - f_z
            clo = closure(point[arm], f_z, f_ca)
            arm_rows[arm] = dict(key=key, F1=point[arm], F1_CI95=ci(dr[arm]), diff_vs_Z=d_z, diff_vs_CA=d_ca, closure=clo,
                                 closure_CI95=[x / den for x in d_z['CI95']] if den > 0 else [float('nan')] * 2, closure_denominator=den,
                                 table=tabs[arm]['pooled'], model=str(model), provenance=provenance_of(run))
            g_z, g_ca = g_diffs.get(f'{arm} - Z'), g_diffs.get(f'{arm} - CarbonA_filt')
            alt_rows[arm] = dict(diff_vs_Z=dict(point=g_z['point'], CI95=[g_z['lo'], g_z['hi']]), diff_vs_CA=dict(point=g_ca['point'], CI95=[g_ca['lo'], g_ca['hi']]), closure=clo)
        summary[sp], alt_summary[sp] = arm_rows, alt_rows
        if 'F1 - S1' in pairs and 'F1' in point and 'S1' in point:
            f1_minus_s1[sp] = pairs['F1 - S1']
        res['species'][sp] = dict(B_contigs=B[sp], blocks=len(bl), reference_chains=int(sum(len(S.ev.refkeys[c]) for c in S.contigs)), F1_Z=f_z, F1_CarbonA_filt=f_ca,
                                  F1_CarbonA_unf=point['CarbonA_unf'], CarbonA_filt_expected_M30=CARBON_A_FILT_EXPECTED[sp],
                                  CarbonA_filt_matches_M30=round(f_ca, 2) == CARBON_A_FILT_EXPECTED[sp], gap_CarbonA_filt_minus_Z=f_ca - f_z,
                                  F1_with_CI95={n: dict(F1=point[n], CI95=ci(dr[n])) for n in names}, pairs=pairs, tables=tabs, arms=arm_rows,
                                  m30_reference_rows=m30_reference_rows(sp),
                                  gene_unit_bootstrap=dict(F1_CI95={n: [v['lo'], v['hi']] for n, v in g_ci.items()}, diffs=g_diffs))
    res['verdict_block_bootstrap'] = verdict(summary, f1_minus_s1)
    alt_f1_minus_s1 = {sp: dict(point=res['species'][sp]['gene_unit_bootstrap']['diffs']['F1 - S1']['point'],
                                CI95=[res['species'][sp]['gene_unit_bootstrap']['diffs']['F1 - S1']['lo'], res['species'][sp]['gene_unit_bootstrap']['diffs']['F1 - S1']['hi']])
                       for sp in f1_minus_s1}
    res['verdict'] = verdict(alt_summary, alt_f1_minus_s1)                     # primary (gene-unit, M30 / M33 implementation)
    res['verdict_statistic_sensitive'] = res['verdict']['label'] != res['verdict_block_bootstrap']['label']
    res['n_trained_arms'] = len(arms)
    if with_seed_coverage:
        try:
            res['seed_coverage_descriptive'] = seed_coverage(root)
        except Exception as e:                                                  # descriptive only; never blocks the verdict
            res['seed_coverage_descriptive'] = f'not computed: {e!r}'
    res['score_seconds'] = round(time.time() - t0, 1)
    tmp = str(out_path) + '.tmp'
    json.dump(res, open(tmp, 'w'), indent=1)
    Path(tmp).replace(out_path)
    return res


def lp_by_key(post_dir, contigs, tau=W.S1_TAU):
    out = {}
    for c in contigs:
        for st in '+-':
            for g in json.load(open(Path(post_dir) / f'tau{tau:g}' / f'{c}_{st}.json'))['genes']:
                out[(c, st, tuple(tuple(e) for e in g['exons']))] = g['log_p']
    return out


def score_post(root, post_root, out_path, base_result):
    """descriptive: F1 of each run after the frozen posterior filter (log P >= thr*), closure on those values, AUROC / ECE of the posterior on B"""
    sys.path.insert(0, str(W.EXPS / 'M33-NTV2-FOLLOWUP' / 'scripts'))
    import m33_common as C33
    res = dict(label='M36 posterior-filtered (descriptive)', thr=W.S1_LOGP_THR, tau=W.S1_TAU, species={})
    for sp in W.SPECIES:
        S = SpeciesScorer(sp)
        runs = {'Z': (W.Z_DIR, f'Z_{sp}')} | {a: (r, key) for a, s, key, m, r in trained_arms(root) if s == sp}
        raw = {}
        for n, (run, pname) in runs.items():
            ch = gff_chains(run, sp)
            lp = lp_by_key(Path(post_root) / pname, B[sp])
            missing = sum(1 for c in B[sp] for k in ch[c] if k not in lp)
            raw[n] = ch
            S.add(n, ch)
            S.add(n + '_post', {c: {k: v for k, v in ch[c].items() if lp.get(k, -1e9) >= W.S1_LOGP_THR} for c in B[sp]})
            rows = [(lp[k], int(k in S.ev.refkeys[c])) for c in B[sp] for k in ch[c] if k in lp]
            cal = C33.calibration(np.exp([r[0] for r in rows]), np.array([r[1] for r in rows]))
            res['species'].setdefault(sp, {}).setdefault('calibration', {})[n] = {k: v for k, v in cal.items() if k != 'bins'} | dict(chains_without_posterior=missing)
        S.add('CarbonA_filt', carbon_chains(sp, True))
        bl = S.blocks()
        counts = S.block_counts(bl)
        point, dr = bootstrap_draws(counts, len(bl))
        den = point['CarbonA_filt'] - point['Z_post'] if 'Z_post' in point else float('nan')
        res['species'][sp].update(F1={n: dict(F1=point[n], CI95=ci(dr[n])) for n in point},
                                  closure_post={n: closure(point[n + '_post'], point['Z_post'], point['CarbonA_filt']) for n in runs if n != 'Z'},
                                  diff_post_vs_Z_post={n: pair(point, dr, n + '_post', 'Z_post') for n in runs if n != 'Z'})
    tmp = str(out_path) + '.tmp'
    json.dump(res, open(tmp, 'w'), indent=1)
    Path(tmp).replace(out_path)
    return res


def carbon_a_numbers():
    """F1 / counts of Carbon-A filtered and unfiltered on the B contigs of both species (sanity: must equal the M30 / M33 recorded values)"""
    out = {}
    for sp in W.SPECIES:
        S = SpeciesScorer(sp)
        for n, filt in (('CarbonA_filt', True), ('CarbonA_unf', False)):
            out.setdefault(sp, {})[n] = S.add(n, carbon_chains(sp, filt))
    return out


def summary_lines(res):
    L = []
    for sp, r in res['species'].items():
        L.append(f"{sp}: Z {r['F1_Z']:.2f}  CarbonA_filt {r['F1_CarbonA_filt']:.2f}  gap {r['gap_CarbonA_filt_minus_Z']:.2f}  blocks {r['blocks']}")
        for a, x in r['arms'].items():
            L.append(f"   {a:2s} F1 {x['F1']:.2f}  -Z {x['diff_vs_Z']['point']:+.2f} [{x['diff_vs_Z']['CI95'][0]:+.2f},{x['diff_vs_Z']['CI95'][1]:+.2f}]  closure {100 * x['closure']:.1f}%  "
                     f"-CarbonA {x['diff_vs_CA']['point']:+.2f} [{x['diff_vs_CA']['CI95'][0]:+.2f},{x['diff_vs_CA']['CI95'][1]:+.2f}]")
    v = res['verdict']
    L.append(f"verdict {v['label']}  go_strong {v['go_strong']}  go {v['go']}  no_go {v['no_go']}  c2_retracted {v['c2_retracted']} {v['c2_by_species']}")
    return L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=('main', 'post', 'carbonA'))
    ap.add_argument('--root', default=str(W.OUT))
    ap.add_argument('--out', default=None)
    ap.add_argument('--post-root', default=None)
    a = ap.parse_args()
    W.K.repoint()
    root = Path(a.root)
    if a.cmd == 'carbonA':
        print(json.dumps(carbon_a_numbers(), indent=1))
        return 0
    results = root / 'results'
    results.mkdir(parents=True, exist_ok=True)
    if a.cmd == 'main':
        out = Path(a.out or results / 'm36_result.json')
        sent = results / 'SCORE_MAIN_STARTED.json'
        ident = sentinel_identity(root)
        expected = {W.arm_key(arm, sp) for arm, sp in W.order()}
        if not expected <= set(ident['models']):
            print(f"trained arms with a finished B inference missing: {sorted(expected - set(ident['models']))}", file=sys.stderr)
            return 2
        rc = SN.start(sent, ident)
        if rc != 0:
            print({10: 'main scoring already finished (B is scored once)', 11: 'refused: inputs differ from the first scoring start'}[rc], file=sys.stderr)
            return rc
        res = score_main(root, out)
        SN.done(sent)
        print('\n'.join(summary_lines(res)))
        return 0
    if a.cmd == 'post':
        base = results / 'm36_result.json'
        assert (results / 'SCORE_MAIN_STARTED.json.done').exists(), 'main scoring has not finished'
        out = Path(a.out or results / 'm36_result_posterior_filtered.json')
        sent = results / 'SCORE_POST_STARTED.json'
        rc = SN.start(sent, sentinel_identity(root))
        if rc != 0:
            print('post scoring: already finished or refused', rc, file=sys.stderr)
            return rc
        res = score_post(root, Path(a.post_root or root / 'postB'), out, json.load(open(base)))
        SN.done(sent)
        print(json.dumps({sp: r['closure_post'] for sp, r in res['species'].items()}))
        return 0


if __name__ == '__main__':
    sys.exit(main())
