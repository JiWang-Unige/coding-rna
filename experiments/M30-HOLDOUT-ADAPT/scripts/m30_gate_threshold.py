"""Posterior threshold for 'gate + temperature-calibrated posterior' (protocol §5-3): chosen on the six Arabidopsis/rice TRAIN chromosomes' gate chains
(tau=0.03 posterior), maximising pooled exact-chain F1; never sees the new species.  -> outputs/M30/gate_posterior_threshold.json"""
import json, sys
import numpy as np
import m30_common as M
import a1_common as K
sys.path.insert(0, str(K.CA_SRC))
import a1_eval as E                       # noqa: E402
from m29_r3_posterior_filter import curve, auroc_ap   # noqa: E402

refs = E.Refs()
rows, nref = [], 0
for s in E.TRAIN:
    ch = E.load_chains(K.A1 / 'gate_ownpath', s)
    lp = {}
    for st in '+-':
        for g in json.load(open(M.OUT / 'posterior_gate_tau0.03_TRAIN' / f'{s}_{st}.json'))['genes']:
            lp[(s, st, tuple(tuple(e) for e in g['exons']))] = g['log_p']
    for k in ch:
        rows.append((lp[k], int(k in refs.keys[s])))
    nref += len(refs.keys[s])
sc = np.array([a for a, _ in rows]); lb = np.array([b for _, b in rows])
t, kept, tp, f1 = curve(sc, lb, nref)
i = int(np.argmax(f1))
auc, ap = auroc_ap(sc, lb)
unf = 100 * 2 * lb.sum() / (len(lb) + nref)
out = dict(tau=0.03, threshold=float(t[i]), kept=int(kept[i]), tp=int(tp[i]), F1_at_threshold=100 * float(f1[i]), F1_unfiltered=float(unf), AUROC=auc, AP=ap, n_chains=len(rows), n_ref=nref)
json.dump(out, open(M.OUT / 'gate_posterior_threshold.json', 'w'), indent=1)
print(out)
