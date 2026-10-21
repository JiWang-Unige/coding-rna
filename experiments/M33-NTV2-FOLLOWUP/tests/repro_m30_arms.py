"""Sanity: re-score the four M30 arms with the M33 loader code path; must equal M30's published values (Dryas 73.98/72.32/71.97/74.09, Esch 73.32/71.75/71.38/72.52).
Usage: repro_m30_arms.py <species>"""
import json, sys
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[1] / 'scripts'))
import m33_common as C
import m30_common as M30
import m30_score as SC
species = sys.argv[1]
B = M30.SPECIES[species]['B']
ev = SC.Evaluator(B)
thr = json.load(open(M30.OUT / 'gate_posterior_threshold.json'))['threshold']
arms = {
    'CarbonA_unf': {c: SC.load_gff_chains(str(M30.OUT / 'carbonA_raw' / species / c / 'prediction.gff'), c, cds96=False) for c in B},
    'CarbonA_filt': {c: SC.load_gff_chains(str(M30.OUT / 'carbonA_raw' / species / c / 'prediction.gff'), c, cds96=False, conf_filter='carbon') for c in B},
    'Gate': {c: SC.load_gff_chains(str(M30.OUT / 'evalB' / 'gate' / c / 'prediction' / 'full_grammar.gff3'), c) for c in B}}
arms['Gate_post'] = {}
for c in B:
    lp = {}
    for st in '+-':
        for g in json.load(open(M30.OUT / 'posterior_gate_tau0.03_B' / f'{c}_{st}.json'))['genes']:
            lp[(c, st, tuple(tuple(e) for e in g['exons']))] = g['log_p']
    arms['Gate_post'][c] = {k: v for k, v in arms['Gate'][c].items() if lp[k] >= thr}
for n, ch in arms.items():
    r = ev.add(n, ch)
    print(species, n, round(r['F1'], 2), r, flush=True)
