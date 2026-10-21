"""M33 step 2 scoring: zero-shot model B on the M30 evaluation-set B chromosomes of one species, side by side with the M30 arms, same ruler,
same Evaluator (`m30_score.Evaluator`: M26 CDS-assessable primaries, exact CDS chain, CDS>=96 bp, gene-unit paired bootstrap, 2000 draws, seed 0).
Arms:  CarbonA_unf, CarbonA_filt, Gate, Gate_post (read from the M30 outputs, unchanged)  +  NTv2B (unfiltered), NTv2B_post (posterior filter, tau*/thr* from
the internal validation, `selection.json`).  Evaluation set B is used for scoring only: no selection of tau / threshold / arm happens here.
Usage: m33_evalB.py <species> <run_dir of m33_infer on the B contigs> <post_dir> <selection.json> [<out_json>]
"""
import json, sys
import numpy as np
import m33_common as C
import m30_common as M30
import m30_score as SC
from m29_r3_posterior_filter import curve

species, run, post, selp = sys.argv[1:5]
outp = sys.argv[5] if len(sys.argv) > 5 else str(C.OUT / f'results_{species}.json')
sel = json.load(open(selp))['selection']
tau, thr = sel['tau'], sel['thr']
B = M30.SPECIES[species]['B']
ev = SC.Evaluator(B)
m30thr = json.load(open(M30.OUT / 'gate_posterior_threshold.json'))
res = {}


def gate_lp(c):
    lp = {}
    for st in '+-':
        for g in json.load(open(M30.OUT / 'posterior_gate_tau0.03_B' / f'{c}_{st}.json'))['genes']:
            lp[(c, st, tuple(tuple(e) for e in g['exons']))] = g['log_p']
    return lp


arms = {}
arms['CarbonA_unf'] = {c: SC.load_gff_chains(str(M30.OUT / 'carbonA_raw' / species / c / 'prediction.gff'), c, cds96=False) for c in B}
arms['CarbonA_filt'] = {c: SC.load_gff_chains(str(M30.OUT / 'carbonA_raw' / species / c / 'prediction.gff'), c, cds96=False, conf_filter='carbon') for c in B}
arms['Gate'] = {c: SC.load_gff_chains(str(M30.OUT / 'evalB' / 'gate' / c / 'prediction' / 'full_grammar.gff3'), c) for c in B}
arms['Gate_post'] = {}
for c in B:
    lp = gate_lp(c)
    arms['Gate_post'][c] = {k: v for k, v in arms['Gate'][c].items() if lp[k] >= m30thr['threshold']}
arms['NTv2B'] = {c: SC.load_gff_chains(f'{run}/{c}/prediction/full_grammar.gff3', c) for c in B}
lpB = C.lp_by_key(post, tau, B)
arms['NTv2B_post'] = {c: {k: v for k, v in arms['NTv2B'][c].items() if lpB[k] >= thr} for c in B}
for n, ch in arms.items():
    res[n] = ev.add(n, ch)
    print(n, res[n], flush=True)
names = list(ev.methods)
pairs = [('NTv2B_post', 'NTv2B'), ('NTv2B', 'CarbonA_unf'), ('NTv2B', 'CarbonA_filt'), ('NTv2B_post', 'CarbonA_filt'), ('NTv2B_post', 'CarbonA_unf'),
         ('NTv2B_post', 'Gate_post'), ('NTv2B_post', 'Gate'), ('NTv2B', 'Gate'), ('Gate_post', 'CarbonA_filt'), ('CarbonA_filt', 'CarbonA_unf')]
cis, diffs = ev.bootstrap(pairs)
out = dict(species=species, B=B, selection=sel, tables={n: ev.table(n) for n in names}, bootstrap_F1=cis, bootstrap_diffs=diffs)
# calibration of model B's posterior on this new species (descriptive: nothing here feeds back into tau / thr)
rows = [(lpB[k], int(k in ev.refkeys[c])) for c in B for k in arms['NTv2B'][c]]
sc = np.array([r[0] for r in rows]); lb = np.array([r[1] for r in rows])
out['posterior_on_B'] = dict(C.calibration(np.exp(sc), lb), kept_fraction=float(np.mean([lpB[k] >= thr for c in B for k in arms['NTv2B'][c]])), mean_P=float(np.exp(sc).mean()))
nref = sum(len(ev.refkeys[c]) for c in B)
tt, kept, tp, f1 = curve(sc, lb, nref)
out['posterior_on_B']['oracle_best_F1_not_for_selection'] = 100 * float(f1.max())
out['posterior_on_B']['oracle_best_thr_not_for_selection'] = float(tt[int(np.argmax(f1))])
json.dump(out, open(outp, 'w'), indent=1)
for n in names:
    t = out['tables'][n]
    sts = ', '.join(k + ':' + format(v['F1'], '.1f') for k, v in t['strata'].items())
    print(f"{n}: F1 {t['pooled']['F1']:.2f} CI [{cis[n]['lo']:.2f},{cis[n]['hi']:.2f}] P/R {100*t['pooled']['tp']/max(t['pooled']['pred'],1):.1f}/{100*t['pooled']['tp']/t['pooled']['ref']:.1f} strata {{{sts}}} base {t['base']['base_F1']:.2f} bg {t['base']['bg_pred_kb']:.0f}kb {t['base']['bg_burden_pct']:.2f}%")
for k, v in diffs.items():
    print(f"  {k}: {v['point']:+.2f} [{v['lo']:+.2f},{v['hi']:+.2f}]")
print('posterior_on_B', json.dumps({k: v for k, v in out['posterior_on_B'].items() if k != 'bins'}))
