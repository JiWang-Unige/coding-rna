"""Summarise gene-posterior diagnostic: AUROC / AP / ECE (10 bins) of P = exp(log p) per tau, pooled over region jobs (and per species)."""
import glob, json, sys
import numpy as np
sys.path.insert(0, '.')
import a2_common as C2
import a1_common as K
sys.path.insert(0, str(K.CA_SRC))
from m29_r3_posterior_filter import auroc_ap, reliability

rows = {}
for f in sorted(glob.glob(str(C2.A2 / 'posterior_diag' / '*.json'))):
    d = json.load(open(f))
    sp = 'ara' if d['seq'].startswith('NC_0030') else 'rice'
    for tau, v in d['taus'].items():
        rows.setdefault(tau, []).extend((x, l, sp) for x, l in zip(v['log_p'], d['labels']))
out = {}
for tau in sorted(rows, key=float, reverse=True):
    lp = np.array([r[0] for r in rows[tau]]); y = np.array([r[1] for r in rows[tau]]); sp = np.array([r[2] for r in rows[tau]])
    p = np.exp(lp)
    auc, ap = auroc_ap(lp, y)
    rel, ece = reliability(p, y)
    out[tau] = dict(n=len(y), precision=float(y.mean()), mean_P=float(p.mean()), AUROC=auc, AP=ap, ECE=ece, frac_P_eq_1=float((lp > -1e-6).mean()),
                    AUROC_ara=auroc_ap(lp[sp == 'ara'], y[sp == 'ara'])[0], AUROC_rice=auroc_ap(lp[sp == 'rice'], y[sp == 'rice'])[0],
                    ECE_ara=reliability(p[sp == 'ara'], y[sp == 'ara'])[1], ECE_rice=reliability(p[sp == 'rice'], y[sp == 'rice'])[1], reliability=rel)
    print(f"tau {tau}: n={len(y)} prec={y.mean():.3f} meanP={p.mean():.3f} AUROC={auc:.3f} AP={ap:.3f} ECE={ece:.3f} (ara {out[tau]['ECE_ara']:.3f}, rice {out[tau]['ECE_rice']:.3f}) frac(P>=1-1e-6)={out[tau]['frac_P_eq_1']:.2f}")
json.dump(out, open(C2.A2 / 'posterior_diag_summary.json', 'w'), indent=1)
