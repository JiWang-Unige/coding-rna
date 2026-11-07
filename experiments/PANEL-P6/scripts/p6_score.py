"""PANEL-P6 scorer: the P6 baselines and our T2 on the core-10 independent panel, the C1 table (report 28 §6) and the seed rule (report 28 §7).
Protocol: docs/review_20261016/31_P6_panel_baselines.md §2 (written before any panel score).

Nothing of the metric is re-implemented:
  converter   M38 `vip_score.convert` (stdlib; validated on Carbon-A / GENERanno / Vipsania GFF dialects in M38, here also on Helixer / Tiberius / ANNEVO) writes every tool's
              GFF into <run>/<contig>/prediction/full_grammar.gff3 (the m34_infer layout of our own runs); Carbon-A filtered = genes with confidence >= 0.1 or NA (M29 `carbon_conf`).
  ruler       P2b / M30: `m30_score.Evaluator` (M26 CDS-assessable primary references of the frozen common_ruler) + `m30_score.load_gff_chains` (unique exact CDS chains of the
              primary transcript of each predicted gene, CDS >= 96 bp rule; the rule-off value is reported as well).  For the panel sequences the references are built by the
              same code (M29 `build_reference`, copied with the path pointed at data/panel/prepared/<key>; cache outputs/PANEL-P6/_refcache).
  bootstrap   the P2b statistic generalised to many species: blocks (contig, first CDS start // 1 Mb) of the references and of every scored run, 2,000 joint draws
              `rng.integers(0, n, n)` with default_rng(0); each species' F1 from the resampled blocks of that species; per draw also the panel median / mean of the species F1.
Sub-commands (generanno env, CPU):
  convert   --tool T --key K                 native outputs of one tool/species -> runs/<tool>/<key>  (tools: Helixer Tiberius ANNEVO CarbonA Vipsania_ft Vipsania_pre)
  validate                                   re-score published DEV / set B values through the converter and this code path (touches only DEV and set-B sequences)
  panel     [--allow-incomplete NO]          the panel table, C1 verdict and seed-rule trigger; refuses to run before every baseline is complete for all 10 species
"""
import argparse
import hashlib
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import p6_common as P                                  # noqa: E402
EXP = P.EXP
sys.path.insert(0, str(EXP / 'M30-HOLDOUT-ADAPT' / 'scripts'))
import m30_score as SC                                 # noqa: E402  (registers the M30 contigs in the M29 helpers; puts the Carbon-A pipeline copy on sys.path)
from m29_r3_posterior_filter import carbon_conf        # noqa: E402
M29 = SC.M29
sys.path.insert(0, str(EXP / 'M34-NTV2-SCALEUP' / 'scripts'))
import m34_panel_common as PC                          # noqa: E402
PC.register(M29)
sys.path.insert(0, str(EXP / 'M38-VIPSANIA' / 'scripts'))
import vip_score as V                                  # noqa: E402  (converter only)

BLK, DRAWS, SEED = 1_000_000, 2000, 0
THETA, SEED_RULE_MARGIN, PARTIAL_MAX = 6.0, 2.5, 10.0
CRITERION_TOOLS = ['Helixer', 'Tiberius', 'Vipsania_ft']          # T2 median must be >= each of these medians (report 28 §6)
CARBON_CRIT = 'CarbonA_filt'                                       # and median(Carbon-A filtered) - median(T2) <= THETA
TOOLS = ['Helixer', 'Tiberius', 'ANNEVO', 'CarbonA_unf', 'CarbonA_filt', 'Vipsania_ft', 'Vipsania_pre']
RUN_ORDER = ['T2'] + TOOLS
TIBERIUS_OUT_OF_SCOPE = ['amborella_trichopoda']                   # Mesangiospermae model (report 10 §1.3): Amborellales outside its range -> sensitivity row only
REFCACHE = P.OUT / '_refcache'


# ======================================================================================================================================
# references of the panel sequences (M29 build_reference with the panel path)
# ======================================================================================================================================
_orig_build = M29.build_reference


def build_reference_panel(seq):
    key = PC._REG['seq2key'].get(seq)
    if key is None:
        return _orig_build(seq)
    p = REFCACHE / f'ref_{seq}.pkl'
    if p.exists():
        return pickle.loads(p.read_bytes())
    C, length, path = M29.C, M29.LENGTH[seq], P.PREP / key
    sequences = C.M.read_allowed(path / 'genome.fa', {seq: length})
    coding = C.S.E.parse_annotation(path / 'reference.gff3', {seq: length}, protein_coding_only=True)
    annotation = C.S.E.parse_annotation(path / 'reference.gff3', {seq: length})
    with (path / 'reference.gff3').open() as handle:
        _, _, raw_cds, _ = C.P.raw_features(handle, seq)
    refs, isoforms = C.select_reference(coding, raw_cds, sequences)
    genes = [(g['start'], g['end']) for g in annotation['genes'].values()]
    obj = dict(seq=seq, species=key, length=length, refs=refs, isoforms=isoforms, gene_spans=genes, annotation=annotation)
    REFCACHE.mkdir(parents=True, exist_ok=True)
    p.write_bytes(pickle.dumps(obj, protocol=4))
    return obj


M29.build_reference = build_reference_panel


# ======================================================================================================================================
# conversion of native outputs
# ======================================================================================================================================
def gff_path(run_dir, contig):
    return Path(run_dir) / contig / 'prediction' / 'full_grammar.gff3'


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def native_sources(tool, key):
    """-> [(gff, [contigs read from it])] and the completion sentinel of one tool / species"""
    contigs = list(P.lengths(key))
    if tool in ('Helixer', 'Tiberius', 'ANNEVO'):
        d = P.NATIVE / tool / key
        return [(d / 'prediction.gff3', contigs)], d / 'DONE'
    if tool in ('CarbonA_unf', 'CarbonA_filt'):
        d = P.CARBON / key
        return [(d / c / 'prediction.gff', [c]) for c in contigs], d / 'DONE'
    if tool in ('Vipsania_ft', 'Vipsania_pre'):
        d = P.VIP / 'runs' / key / tool.split('_')[1]
        return [(d / 'vipsania.gff3', contigs)], d / f"{tool.split('_')[1]}_DONE"
    raise KeyError(tool)


def convert_tool(tool, key, out_root=None):
    out = Path(out_root or P.RUNS) / tool / key
    srcs, done = native_sources(tool, key)
    assert done.exists(), f'{tool} {key}: {done} missing (run not complete)'
    rep = dict(tool=tool, key=key, sources=[], contigs={})
    for gff, contigs in srcs:
        gf = None
        if tool == 'CarbonA_filt':
            conf = carbon_conf(str(gff))
            gf = (lambda g, conf=conf: conf[g] is None or conf[g] >= 0.1)
        rep['contigs'].update(V.convert(str(gff), out, contigs, gene_filter=gf))
        rep['sources'].append(dict(gff=str(gff), sha256=sha256(gff)))
    rep['time'] = time.strftime('%Y-%m-%dT%H:%M:%S')
    json.dump(rep, open(out / 'CONVERTED.json', 'w'), indent=1)
    return out


# ======================================================================================================================================
# scoring
# ======================================================================================================================================
def f1(tp, pred, ref):
    return 200 * tp / max(pred + ref, 1)


def strat(n):
    return '1' if n == 1 else '2-5' if n <= 5 else '>=6'


class Panel:
    """One Evaluator over every contig of every group (group = species) + the chains of every added run."""

    def __init__(self, groups):
        self.groups = {g: list(cs) for g, cs in groups.items()}
        self.contigs = [c for cs in self.groups.values() for c in cs]
        self.group_of = {c: g for g, cs in self.groups.items() for c in cs}
        self.ev = SC.Evaluator(self.contigs)
        self.chains, self.dirs = {}, {}

    def _dir(self, dirs, c):
        return dirs[self.group_of[c]] if isinstance(dirs, dict) else dirs

    def load(self, dirs, cds96=True):
        out = {}
        for c in self.contigs:
            p = gff_path(self._dir(dirs, c), c)
            assert p.exists(), p
            out[c] = SC.load_gff_chains(str(p), c, cds96=cds96)
        return out

    def add(self, name, dirs, cds96=True):
        ch = self.load(dirs, cds96)
        self.chains[name], self.dirs[name] = ch, dirs
        self.ev.add(name, ch)

    def group_counts(self, chains):
        out = {}
        for g, cs in self.groups.items():
            tp = sum(sum(1 for k in chains[c] if k in self.ev.refkeys[c]) for c in cs)
            pred = sum(len(chains[c]) for c in cs)
            ref = sum(len(self.ev.refkeys[c]) for c in cs)
            out[g] = dict(tp=tp, pred=pred, ref=ref, F1=f1(tp, pred, ref), precision=100 * tp / max(pred, 1), recall=100 * tp / max(ref, 1))
        t = {k: sum(v[k] for v in out.values()) for k in ('tp', 'pred', 'ref')}
        out['pooled'] = dict(t, F1=f1(t['tp'], t['pred'], t['ref']), precision=100 * t['tp'] / max(t['pred'], 1), recall=100 * t['tp'] / max(t['ref'], 1))
        return out

    def strata(self, name):
        out = {}
        for g, cs in self.groups.items():
            s = {k: [0, 0, 0] for k in ('1', '2-5', '>=6')}
            for c in cs:
                for k, v in self.chains[name][c].items():
                    s[strat(v['nex'])][1] += 1
                    s[strat(v['nex'])][0] += k in self.ev.refkeys[c]
                for k in self.ev.refkeys[c]:
                    s[strat(self.ev.nex[k])][2] += 1
            out[g] = {k: dict(tp=a, pred=b, ref=r, F1=f1(a, b, r), recall=100 * a / max(r, 1), precision=100 * a / max(b, 1)) for k, (a, b, r) in s.items()}
        return out

    def per_contig(self, name):
        return {c: dict(p, F1=f1(p['tp'], p['pred'], p['ref'])) for c, p in self.ev.methods[name]['per'].items()}

    # ---- joint block bootstrap (P2b statistic, many groups) ---------------------------------------------------------------------------------------
    def bootstrap(self, names, ref_name, draws=DRAWS, seed=SEED):
        bl = set()
        for c in self.contigs:
            bl.update((c, k[2][0][0] // BLK) for k in self.ev.refkeys[c])
        for n in names:
            for c in self.contigs:
                bl.update((c, k[2][0][0] // BLK) for k in self.chains[n][c])
        bl = sorted(bl)
        idx = {b: i for i, b in enumerate(bl)}
        G = list(self.groups)
        gi = np.array([G.index(self.group_of[b[0]]) for b in bl])
        rf = np.zeros(len(bl))
        for c in self.contigs:
            for k in self.ev.refkeys[c]:
                rf[idx[(c, k[2][0][0] // BLK)]] += 1
        C = {}
        for n in names:
            tp, pr = np.zeros(len(bl)), np.zeros(len(bl))
            for c in self.contigs:
                for k in self.chains[n][c]:
                    i = idx[(c, k[2][0][0] // BLK)]
                    pr[i] += 1
                    tp[i] += k in self.ev.refkeys[c]
            C[n] = (tp, pr)

        def scores(n, ii):
            """-> (pooled F1, per-group F1 array) of the blocks ii"""
            tp, pr = C[n][0][ii], C[n][1][ii]
            r, g = rf[ii], gi[ii]
            tg, pg, rg = (np.bincount(g, weights=x, minlength=len(G)) for x in (tp, pr, r))
            return 200 * tp.sum() / max(pr.sum() + r.sum(), 1), 200 * tg / np.maximum(pg + rg, 1)

        full = np.arange(len(bl))
        point = {n: scores(n, full) for n in names}
        rng = np.random.default_rng(seed)
        D = {n: dict(pooled=np.empty(draws), groups=np.empty((draws, len(G)))) for n in names}
        for d in range(draws):
            ii = rng.integers(0, len(bl), len(bl))
            for n in names:
                D[n]['pooled'][d], D[n]['groups'][d] = scores(n, ii)
        ci = lambda v: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
        res = dict(blocks=len(bl), draws=draws, seed=seed, groups=G, runs={}, diffs={})
        for n in names:
            pm, pg = point[n]
            res['runs'][n] = dict(pooled=dict(F1=float(pm), CI95=ci(D[n]['pooled'])),
                                  **{g: dict(F1=float(pg[j]), CI95=ci(D[n]['groups'][:, j])) for j, g in enumerate(G)},
                                  panel_median=dict(F1=float(np.median(pg)), CI95=ci(np.median(D[n]['groups'], axis=1))),
                                  panel_mean=dict(F1=float(np.mean(pg)), CI95=ci(np.mean(D[n]['groups'], axis=1))))
        for n in names:
            if n == ref_name:
                continue
            dd = {}
            for j, g in enumerate(['pooled'] + G):
                a = D[ref_name]['pooled'] if g == 'pooled' else D[ref_name]['groups'][:, j - 1]
                b = D[n]['pooled'] if g == 'pooled' else D[n]['groups'][:, j - 1]
                pa = point[ref_name][0] if g == 'pooled' else point[ref_name][1][j - 1]
                pb = point[n][0] if g == 'pooled' else point[n][1][j - 1]
                dd[g] = dict(point=float(pa - pb), CI95=ci(a - b), frac_positive=float(((a - b) > 0).mean()))
            for stat, fn in (('panel_median', np.median), ('panel_mean', np.mean)):
                a, b = fn(D[ref_name]['groups'], axis=1), fn(D[n]['groups'], axis=1)
                dd[stat] = dict(point=float(fn(point[ref_name][1]) - fn(point[n][1])), CI95=ci(a - b), frac_positive=float(((a - b) > 0).mean()))
            res['diffs'][n] = dd
        return res


# ======================================================================================================================================
# C1 verdict (report 28 §6, outline §0) and seed rule (report 28 §7) -- pure functions of the species F1
# ======================================================================================================================================
def c1_verdict(med):
    """med: {run: panel median F1}.  -> dict with the conditions, the verdict and the seed-rule trigger"""
    t2 = med['T2']
    lead = {t: t2 - med[t] for t in CRITERION_TOOLS}
    gap = med[CARBON_CRIT] - t2
    cond_tools = all(v >= 0 for v in lead.values())
    cond_gap = gap <= THETA
    if gap > PARTIAL_MAX:
        verdict = '不成立'
    elif cond_tools and cond_gap:
        verdict = '成立'
    else:
        verdict = '部分成立'
    trig1 = min(lead.values()) < SEED_RULE_MARGIN
    trig2 = abs(gap - THETA) <= SEED_RULE_MARGIN
    return dict(T2_median=t2, lead_T2_minus=lead, gap_CarbonA_filt_minus_T2=gap, cond_T2_ge_Helixer_Tiberius_Vipsania=cond_tools, cond_gap_le_theta=cond_gap,
                theta=THETA, verdict=verdict, note_unevaluated='outline §0 "不成立" also covers "T2 loses its advantage over NT-v2 + LoRA on the panel": NT-v2 + LoRA was not run on the panel, so that clause is not evaluated',
                seed_rule=dict(trigger_min_lead_lt_2p5=trig1, trigger_gap_within_theta_pm_2p5=trig2, triggered=bool(trig1 or trig2), min_lead=min(lead.values()),
                               rule='report 28 §7: min lead over Helixer / Tiberius / Vipsania-ft < 2.5, or |gap to Carbon-A filtered - 6.0| <= 2.5'))


def species_tests(F, tools):
    """paired species-level tests of T2 - tool (report 10 §4.1): Wilcoxon signed-rank (two-sided) + sign test, Holm over the tools"""
    from scipy.stats import binomtest, wilcoxon
    out = {}
    for t in tools:
        d = np.array([F['T2'][g] - F[t][g] for g in F['T2']])
        nz = d[d != 0]
        w = wilcoxon(nz) if len(nz) else None
        s = binomtest(int((nz > 0).sum()), len(nz)) if len(nz) else None
        out[t] = dict(n=len(d), wins=int((d > 0).sum()), losses=int((d < 0).sum()), ties=int((d == 0).sum()), median_diff=float(np.median(d)),
                      wilcoxon_p=float(w.pvalue) if w else None, sign_p=float(s.pvalue) if s else None)
    ps = sorted((v['wilcoxon_p'], t) for t, v in out.items() if v['wilcoxon_p'] is not None)
    m, run = len(ps), 0.0
    for i, (p, t) in enumerate(ps):
        run = max(run, min(1.0, (m - i) * p))
        out[t]['wilcoxon_p_holm'] = run
    return out


# ======================================================================================================================================
# validate: published DEV / set B values through the converter and this code path
# ======================================================================================================================================
def cmd_validate(args):
    out = P.OUT / 'validation'
    rep, ok_all = dict(checks={}), True
    VR = out / 'runs'

    def check(name, got, want):
        nonlocal ok_all
        good = (abs(got - want) <= 1e-9) if isinstance(want, float) else got == want
        ok_all &= bool(good)
        rep['checks'][name] = dict(got=got, expected=want, ok=bool(good))
        print(f"[{'OK ' if good else 'FAIL'}] {name}: got {got} expected {want}", flush=True)

    M34 = P.HOME / 'coding-rna' / 'outputs' / 'M34-NTV2-SCALEUP'
    M38 = P.HOME / 'coding-rna' / 'outputs' / 'M38-VIPSANIA'
    # ---- A. DEV (frozen M29 references of NC_003074.8 / NC_089041.1) ------------------------------------------------------------------------------
    dev = Panel({'arabidopsis_thaliana': ['NC_003074.8'], 'oryza_sativa': ['NC_089041.1']})
    pub = {'ANNEVO': ((5644, 7183, 7728), 75.70250150895312, {'arabidopsis_thaliana': (3941, 4898, 5437), 'oryza_sativa': (1703, 2285, 2291)}),
           'Tiberius': ((5441, 7427, 7728), 71.80468492246783, {'arabidopsis_thaliana': (3865, 4984, 5437), 'oryza_sativa': (1576, 2443, 2291)}),
           'Helixer': ((5112, 8539, 7728), 62.85117108255979, {'arabidopsis_thaliana': (3673, 5453, 5437), 'oryza_sativa': (1439, 3086, 2291)})}
    nat = P.OUT / 'validation' / 'M28-NATIVE-BASELINES'
    rep['dev_rule_on'] = {}
    for tool, (tot, fpub, per) in pub.items():
        d = VR / 'DEV' / tool
        for sp, c in (('arabidopsis_thaliana', 'NC_003074.8'), ('oryza_sativa', 'NC_089041.1')):
            V.convert(str(nat / tool / sp / 'prediction.gff3'), d, [c])
        g = dev.group_counts(dev.load(d, cds96=False))
        check(f'DEV {tool} rule off: pooled tp/pred/ref == published', (g['pooled']['tp'], g['pooled']['pred'], g['pooled']['ref']), tot)
        check(f'DEV {tool} rule off: pooled F1 == published', g['pooled']['F1'], fpub)
        for sp, t in per.items():
            check(f'DEV {tool} rule off: {sp} tp/pred/ref == published', (g[sp]['tp'], g[sp]['pred'], g[sp]['ref']), t)
        g96 = dev.group_counts(dev.load(d, cds96=True))
        rep['dev_rule_on'][tool] = dict(pooled=g96['pooled'], rule_effect=g96['pooled']['F1'] - g['pooled']['F1'])
        print(f'   DEV {tool} with the CDS >= 96 rule: {g96["pooled"]["F1"]:.4f} (rule effect {g96["pooled"]["F1"] - g["pooled"]["F1"]:+.4f})')
    cag = {'NC_003074.8': M29.OUTS / 'M29-CARBON-A-GAP' / 'G0' / 'ara_pred' / 'prediction.gff', 'NC_089041.1': M29.OUTS / 'M29-CARBON-A-GAP' / 'G0' / 'rice_pred' / 'prediction.gff'}
    for name, filt, tot in (('CarbonA_unf', False, (5730, 7503, 7728)), ('CarbonA_filt', True, (5687, 6966, 7728))):
        d = VR / 'DEV' / name
        for c, gff in cag.items():
            conf = carbon_conf(str(gff)) if filt else None
            V.convert(str(gff), d, [c], gene_filter=(lambda g, conf=conf: conf[g] is None or conf[g] >= 0.1) if filt else None)
        g = dev.group_counts(dev.load(d, cds96=False))
        check(f'DEV {name} rule off: pooled tp/pred/ref == published', (g['pooled']['tp'], g['pooled']['pred'], g['pooled']['ref']), tot)
        g96 = dev.group_counts(dev.load(d, cds96=True))
        check(f'DEV {name}: rule on == rule off (pipeline has its own >= 96 bp floor)', g96['pooled']['tp'] * 10**6 + g96['pooled']['pred'], g['pooled']['tp'] * 10**6 + g['pooled']['pred'])
    bpub = json.load(open(M34 / 'boot_pair_T10_T2_DEV.json'))
    g = dev.group_counts(dev.load(M34 / 'pc2_T2' / 'devM34_T2', cds96=True))
    check('DEV T2 (P5 model, rule on) pooled F1 == published 73.88', g['pooled']['F1'], bpub['T2']['F1'])
    vdev = json.load(open(M38 / 'results' / 'dev' / 'vip_dev_result.json'))
    d = VR / 'DEV' / 'Vipsania_ft'
    for sp, c in (('arabidopsis_thaliana', 'NC_003074.8'), ('oryza_sativa', 'NC_089041.1')):
        V.convert(str(M38 / 'runs' / sp / 'ft' / 'vipsania.gff3'), d, [c])
    g = dev.group_counts(dev.load(d, cds96=True))
    check('DEV Vipsania ft (rule on) pooled F1 == M38 69.47', g['pooled']['F1'], vdev['variants']['VIP_ft']['pooled']['F1'])
    # ---- B. set B with the PANEL references (data/panel/prepared) --------------------------------------------------------------------------------
    B = {'dryas_octopetala': ['NC_141367.1', 'NC_141369.1', 'NC_141371.1', 'NC_141373.1'],
         'eschscholzia_californica': ['NW_028238309.1', 'NW_028238461.1', 'NW_028238483.1']}
    for c in [c for cs in B.values() for c in cs]:                     # panel references == the M30 (data/m1_screen) references used by P2b
        a = {M29.key(t) for t in build_reference_panel(c)['refs']}
        b = {M29.key(t) for t in _orig_build(c)['refs']}
        check(f'set B {c}: panel reference keys == M30 reference keys ({len(a)})', a == b, True)
    sb = Panel(B)
    p2b = json.load(open(M34 / 'p2b' / 'p2b_result.json'))
    sb.add('T10', M34 / 'p2b' / 'T10')
    sb.add('T2', M34 / 'p2b' / 'T2')
    for n in ('T10', 'T2'):
        g = sb.group_counts(sb.chains[n])
        for sc in ('pooled', 'dryas_octopetala', 'eschscholzia_californica'):
            want = p2b['runs'][n]['table'][sc]
            check(f'set B {n} [{sc}] tp/pred/ref == P2b', (g[sc]['tp'], g[sc]['pred'], g[sc]['ref']), (want['tp'], want['pred'], want['ref']))
            check(f'set B {n} [{sc}] F1 == P2b', g[sc]['F1'], want['F1'])
    bs = sb.bootstrap(['T10', 'T2'], 'T10')
    check('set B bootstrap blocks == P2b', bs['blocks'], p2b['bootstrap']['blocks'])
    for n in ('T10', 'T2'):
        for sc in ('pooled', 'dryas_octopetala', 'eschscholzia_californica'):
            check(f'set B bootstrap {n} [{sc}] CI95 == P2b', bs['runs'][n][sc]['CI95'], p2b['runs'][n]['bootstrap_F1'][sc]['CI95'])
    for sc in ('pooled', 'dryas_octopetala', 'eschscholzia_californica'):
        check(f'set B T10 - T2 [{sc}] CI95 == P2b', bs['diffs']['T2'][sc]['CI95'], p2b['paired_new_minus']['T2'][sc]['CI95'])
    m33 = {'CarbonA_unf': {'dryas_octopetala': (6885, 9532, 9508), 'eschscholzia_californica': (8683, 11838, 12367)},        # M30/M33 set-B tables
           'CarbonA_filt': {'dryas_octopetala': (6807, 8895, 9508), 'eschscholzia_californica': (8584, 11048, 12367)}}      # (experiments/M33-NTV2-FOLLOWUP/results)
    for name, key, filt in (('CarbonA_unf', 'CarbonA_unf', False), ('CarbonA_filt', 'CarbonA_filt', True)):
        d = VR / 'B' / name
        for sp, cs in B.items():
            for c in cs:
                gff = P.HOME / 'coding-rna' / 'outputs' / 'M30-HOLDOUT-ADAPT' / 'carbonA_raw' / sp / c / 'prediction.gff'
                conf = carbon_conf(str(gff)) if filt else None
                V.convert(str(gff), d, [c], gene_filter=(lambda g, conf=conf: conf[g] is None or conf[g] >= 0.1) if filt else None)
        g = sb.group_counts(sb.load(d, cds96=True))
        for sp in B:
            check(f'set B {name} {sp} tp/pred/ref == M30/M33', (g[sp]['tp'], g[sp]['pred'], g[sp]['ref']), m33[key][sp])
    vsb = json.load(open(M38 / 'results' / 'setb' / 'vip_setb_result.json'))
    d = VR / 'B' / 'Vipsania_pre'
    for sp, cs in B.items():
        V.convert(str(M38 / 'runs' / sp / 'pre' / 'vipsania.gff3'), d, cs)
    g = sb.group_counts(sb.load(d, cds96=True))
    for sc in ('pooled', 'dryas_octopetala', 'eschscholzia_californica'):
        want = vsb['variants']['VIP_pre']['table'][sc]
        check(f'set B Vipsania pre [{sc}] tp/pred/ref == M38', (g[sc]['tp'], g[sc]['pred'], g[sc]['ref']), (want['tp'], want['pred'], want['ref']))
    # ---- C. verdict logic on synthetic medians ---------------------------------------------------------------------------------------------------
    base = dict(T2=70.0, Helixer=60.0, Tiberius=65.0, Vipsania_ft=66.0, CarbonA_filt=72.0)
    check('C1 synthetic: gap 2, leads >= 4 -> 成立, no trigger', (c1_verdict(base)['verdict'], c1_verdict(base)['seed_rule']['triggered']), ('成立', False))
    x = dict(base, CarbonA_filt=77.0)
    check('C1 synthetic: gap 7 -> 部分成立, trigger (gap within 6 +- 2.5)', (c1_verdict(x)['verdict'], c1_verdict(x)['seed_rule']['triggered']), ('部分成立', True))
    x = dict(base, Tiberius=70.5)
    check('C1 synthetic: below Tiberius -> 部分成立, trigger (lead < 2.5)', (c1_verdict(x)['verdict'], c1_verdict(x)['seed_rule']['triggered']), ('部分成立', True))
    x = dict(base, CarbonA_filt=80.5)
    check('C1 synthetic: gap 10.5 -> 不成立', c1_verdict(x)['verdict'], '不成立')
    x = dict(base, CarbonA_filt=76.0, Vipsania_ft=67.6)
    check('C1 synthetic: gap 6.0, lead 2.4 -> 成立, trigger', (c1_verdict(x)['verdict'], c1_verdict(x)['seed_rule']['triggered']), ('成立', True))
    rep.update(ok=bool(ok_all), time=time.strftime('%Y-%m-%dT%H:%M:%S'))
    out.mkdir(parents=True, exist_ok=True)
    json.dump(rep, open(out / 'validation.json', 'w'), indent=1, default=str)
    print('VALIDATION', 'PASSED' if ok_all else 'FAILED', f'({sum(c["ok"] for c in rep["checks"].values())}/{len(rep["checks"])} checks)')
    return 0 if ok_all else 1


# ======================================================================================================================================
# panel
# ======================================================================================================================================
def completeness():
    """-> {tool: {key: bool}} of the native outputs, and T2 (P5)"""
    st = {'T2': {k: (P.T2_PANEL.parent / 'PANEL_INFER_DONE').exists() and all(gff_path(P.T2_PANEL / k, c).exists() for c in P.lengths(k)) for k in P.CORE10}}
    for t in TOOLS:
        st[t] = {k: native_sources(t, k)[1].exists() for k in P.CORE10}
    return st


def cmd_status(args):
    st = completeness()
    for t, d in st.items():
        print(f"{t:14s} {sum(d.values())}/{len(d)}  missing: {[k for k, v in d.items() if not v]}")
    return 0


def cmd_refs(args):
    """build (cache) the references of the 119 panel sequences and compare the per-species CDS-assessable primary counts with report 25 (reference_stats.tsv);
    reads only the reference annotation, scores nothing"""
    want = {'dryas_octopetala': 21467, 'eschscholzia_californica': 26580, 'akebia_trifoliata': 23926, 'hakea_prostrata': 23486, 'zizania_latifolia': 35671,
            'themeda_triandra': 41205, 'juncus_effusus': 17779, 'rhynchospora_breviuscula': 19976, 'acorus_gramineus': 21365, 'amborella_trichopoda': 16467}
    out, ok = {}, True
    for k in P.CORE10:
        t0 = time.time()
        n = sum(len(build_reference_panel(c)['refs']) for c in P.lengths(k))
        out[k] = dict(refs=n, report25=want[k], equal=n == want[k], seconds=round(time.time() - t0, 1))
        ok &= n == want[k]
        print(f"{k}: {n} CDS-assessable primary references (report 25: {want[k]}) {'OK' if n == want[k] else 'DIFF'} {out[k]['seconds']} s", flush=True)
    out['total'] = sum(v['refs'] for v in out.values())
    P.RESULTS.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(P.RESULTS / 'panel_reference_counts.json', 'w'), indent=1)
    print('REFS', 'EQUAL' if ok else 'DIFFER', out['total'])
    return 0 if ok else 1


def cmd_convert(args):
    o = convert_tool(args.tool, args.key)
    print(json.dumps(json.load(open(o / 'CONVERTED.json'))['contigs'])[:2000])
    return 0


def cmd_panel(args):
    P.RESULTS.mkdir(parents=True, exist_ok=True)
    res_path = P.RESULTS / 'p6_panel_result.json'
    if res_path.exists() and not args.overwrite:
        raise SystemExit(f'{res_path} exists: the panel is scored once (protocol §2); --overwrite only to redo a failed scoring')
    st = completeness()
    missing = {t: [k for k, v in d.items() if not v] for t, d in st.items() if not all(d.values())}
    if missing:
        raise SystemExit(f'not every run is complete, the panel is not scored before all baselines are done (protocol §2): {missing}')
    t0 = time.time()
    for t in TOOLS:
        for k in P.CORE10:
            if not (P.RUNS / t / k / 'CONVERTED.json').exists():
                convert_tool(t, k)
    groups = {k: list(P.lengths(k)) for k in P.CORE10}
    pn = Panel(groups)
    dirs = {'T2': {k: P.T2_PANEL / k for k in P.CORE10}, **{t: {k: P.RUNS / t / k for k in P.CORE10} for t in TOOLS}}
    for n in RUN_ORDER:
        pn.add(n, dirs[n], cds96=True)
        print(f'loaded {n} {time.time() - t0:.0f} s', flush=True)
    tables = {n: pn.group_counts(pn.chains[n]) for n in RUN_ORDER}
    rule_off = {n: pn.group_counts(pn.load(dirs[n], cds96=False)) for n in RUN_ORDER}
    strata = {n: pn.strata(n) for n in RUN_ORDER}
    F = {n: {k: tables[n][k]['F1'] for k in P.CORE10} for n in RUN_ORDER}
    med = {n: float(np.median(list(F[n].values()))) for n in RUN_ORDER}
    boot = pn.bootstrap(RUN_ORDER, 'T2')
    verdict = c1_verdict(med)
    sp = P.species()
    sens = {}
    inscope = [k for k in P.CORE10 if k not in TIBERIUS_OUT_OF_SCOPE]
    sens['Tiberius_in_scope_9'] = dict(species=inscope, T2_median=float(np.median([F['T2'][k] for k in inscope])), Tiberius_median=float(np.median([F['Tiberius'][k] for k in inscope])),
                                       verdict_with_this_Tiberius_median=c1_verdict(dict(med, Tiberius=float(np.median([F['Tiberius'][k] for k in inscope]))))['verdict'])
    for label, keep in (('TH', lambda k: sp[k]['tier'] == 'TH'), ('SH', lambda k: sp[k]['tier'] == 'SH'), ('AH', lambda k: sp[k]['tier'] == 'AH'),
                        ('Poaceae', lambda k: sp[k]['family'] == 'Poaceae'), ('non-Poaceae', lambda k: sp[k]['family'] != 'Poaceae'),
                        ('RefSeq', lambda k: sp[k]['source'] == 'RefSeq'), ('GenBank', lambda k: sp[k]['source'] != 'RefSeq')):
        ks = [k for k in P.CORE10 if keep(k)]
        if ks:
            sens[f'stratum_{label}'] = dict(species=ks, medians={n: float(np.median([F[n][k] for k in ks])) for n in RUN_ORDER})
    res = dict(label='PANEL-P6 core-10 panel', protocol='docs/review_20261016/31_P6_panel_baselines.md §2; C1 = report 28 §6; seed rule = report 28 §7',
               ruler='m30_score.Evaluator + load_gff_chains (exact CDS chain of the primary transcript, CDS >= 96 bp rule, no posterior filter); panel references = frozen common_ruler on data/panel/prepared',
               species=P.CORE10, contigs=groups, runs={n: {k: str(v) for k, v in dirs[n].items()} for n in RUN_ORDER},
               tables=tables, rule_off=rule_off, strata=strata, per_contig={n: pn.per_contig(n) for n in RUN_ORDER},
               species_F1=F, panel_median=med, panel_mean={n: float(np.mean(list(F[n].values()))) for n in RUN_ORDER},
               panel_min={n: float(min(F[n].values())) for n in RUN_ORDER}, bootstrap=boot, C1=verdict, sensitivity=sens,
               species_tests=species_tests(F, TOOLS), score_seconds=round(time.time() - t0, 1), time=time.strftime('%Y-%m-%dT%H:%M:%S'))
    json.dump(res, open(res_path, 'w'), indent=1, default=str)
    write_tsv(res, P.RESULTS / 'p6_panel_table.tsv')
    print(json.dumps(dict(panel_median=med, C1=verdict), indent=1, ensure_ascii=False))
    return 0


def write_tsv(res, path):
    runs = RUN_ORDER
    rows = ['species\t' + '\t'.join(runs)]
    for k in res['species']:
        rows.append(k + '\t' + '\t'.join(f"{res['species_F1'][n][k]:.2f}" for n in runs))
    for stat in ('panel_median', 'panel_mean', 'panel_min'):
        rows.append(stat + '\t' + '\t'.join(f"{res[stat][n]:.2f}" for n in runs))
    open(path, 'w').write('\n'.join(rows) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('convert'); p.add_argument('--tool', required=True, choices=TOOLS); p.add_argument('--key', required=True, choices=P.CORE10); p.set_defaults(fn=cmd_convert)
    p = sub.add_parser('validate'); p.set_defaults(fn=cmd_validate)
    p = sub.add_parser('status'); p.set_defaults(fn=cmd_status)
    p = sub.add_parser('refs'); p.set_defaults(fn=cmd_refs)
    p = sub.add_parser('panel'); p.add_argument('--overwrite', action='store_true'); p.set_defaults(fn=cmd_panel)
    a = ap.parse_args()
    sys.exit(a.fn(a) or 0)


if __name__ == '__main__':
    main()
