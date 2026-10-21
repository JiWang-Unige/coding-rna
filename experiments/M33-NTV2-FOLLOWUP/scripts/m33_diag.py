"""Descriptive failure-mode split (not used for any selection): for every reference primary gene of a contig set, is it (a) hit exactly, (b) overlapped by >=1 predicted
chain on the same strand but never exactly (structure / junction error), (c) not overlapped by any predicted chain (missed).  Same gene units as the gene-unit bootstrap.
Usage: m33_diag.py <out_json>   (model B on DEV, Dryas B, Eschscholzia B; Carbon-A filtered and gate+posterior on the two new species)"""
import json, sys
import numpy as np
import m33_common as C
import m30_common as M30
import m30_score as SC
import a1_eval as E
out = {}


def split(ev, name, contigs):
    tp, pred = ev.arrays(name)
    n = ev.n_gene_units
    res = {}
    nexs = np.array([ev.nex[u[2]] for u in ev.units[:n]])
    strat = np.array([SC.strat(x) for x in nexs])
    for lab, m in (('all', np.ones(n, bool)), ('1', strat == '1'), ('2-5', strat == '2-5'), ('>=6', strat == '>=6')):
        t = tp[:n][m] > 0; o = (pred[:n][m] > 0) & ~t; z = pred[:n][m] == 0
        res[lab] = dict(n=int(m.sum()), exact=float(t.mean()), overlap_not_exact=float(o.mean()), missed=float(z.mean()))
    return res


for sp in M30.SPECIES:
    B = M30.SPECIES[sp]['B']
    ev = SC.Evaluator(B)
    thrG = json.load(open(M30.OUT / 'gate_posterior_threshold.json'))['threshold']
    ev.add('CarbonA_filt', {c: SC.load_gff_chains(str(M30.OUT / 'carbonA_raw' / sp / c / 'prediction.gff'), c, cds96=False, conf_filter='carbon') for c in B})
    ev.add('NTv2B', {c: SC.load_gff_chains(str(C.OUT / ('evalB_dryas' if sp.startswith('dryas') else 'evalB_esch') / c / 'prediction' / 'full_grammar.gff3'), c) for c in B})
    out[sp] = {n: split(ev, n, B) for n in ('CarbonA_filt', 'NTv2B')}
    print(sp, json.dumps(out[sp]), flush=True)
seqs = list(C.DEV)
ev = SC.Evaluator(seqs)
ev.add('NTv2B', {s: E.load_chains(C.OUT / 'dev', s) for s in seqs})
out['DEV'] = {'NTv2B': split(ev, 'NTv2B', seqs)}
print('DEV', json.dumps(out['DEV']), flush=True)
json.dump(out, open(sys.argv[1], 'w'), indent=1)
