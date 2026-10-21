"""M33 step 1a: choose tau and the posterior filter threshold on the TRAIN INTERNAL-VALIDATION chromosomes only (NC_003071.7, NC_089044.1).
Rule (declared in protocol.md before any posterior of model B existed):
  tau*  = argmin over TAU_GRID of the pooled ECE (10 equal-width bins) of P(g) = exp(log_p) against the exact-chain hit label
  thr*  = log_p value maximising pooled exact-chain F1 of 'keep chains with log_p >= thr' on those chromosomes at tau*   (same `curve` as M30)
Also reported (descriptive, not used for selection): AUROC / AP / ECE for every tau, per species; cross-species transfer of thr*; best F1 per tau.
Usage: m33_select.py <val_run_dir> <post_dir> <out_json>
"""
import json, sys
import numpy as np
import m33_common as C
import a1_eval as E
from m29_r3_posterior_filter import curve

run, post, outp = sys.argv[1], sys.argv[2], sys.argv[3]
seqs = list(C.HELDOUT)
refs = E.Refs()
chains = {s: E.load_chains(run, s) for s in seqs}
nref = {s: len(refs.keys[s]) for s in seqs}
NREF = sum(nref.values())
unf = E.summarize(chains, refs, seqs)['pooled']
res = dict(val_run=str(run), seqs=seqs, n_ref=nref, unfiltered=unf, taus={})
for tau in C.TAU_GRID:
    try:
        lp = C.lp_by_key(post, tau, seqs)
    except FileNotFoundError:
        continue
    rows = []                                                 # (log_p, hit, seq)
    for s in seqs:
        for k in chains[s]:
            rows.append((lp[k], int(k in refs.keys[s]), s))
    sc = np.array([r[0] for r in rows]); lb = np.array([r[1] for r in rows]); sq = np.array([r[2] for r in rows])
    t, kept, tp, f1 = curve(sc, lb, NREF)
    i = int(np.argmax(f1))
    cal = C.calibration(np.exp(sc), lb)
    per = {}
    for s in seqs:
        m = sq == s
        c = C.calibration(np.exp(sc[m]), lb[m])
        per[s] = dict(AUROC=c['AUROC'], AP=c['AP'], ECE=c['ECE'], n=c['n'], hits=c['hits'])
    d = dict(thr=float(t[i]), kept=int(kept[i]), tp=int(tp[i]), F1_filtered=100 * float(f1[i]), F1_unfiltered=unf['F1'],
             mean_P=float(np.exp(sc).mean()), frac_P_ge_0_999=float((np.exp(sc) >= 0.999).mean()), per_species=per, **cal)
    # cross-species transfer of the pooled threshold's twin: threshold fitted on one species, applied to the other
    xs = {}
    for s in seqs:
        o = [x for x in seqs if x != s][0]
        ms, mo = sq == s, sq == o
        ts, ks, tps, fs = curve(sc[ms], lb[ms], nref[s])
        thr_s = float(ts[int(np.argmax(fs))])
        keep = mo & (sc >= thr_s)
        xs[f'thr_from_{s}_applied_to_{o}'] = dict(thr=thr_s, F1=100 * 2 * float(lb[keep].sum()) / (keep.sum() + nref[o]),
                                                  F1_unfiltered=100 * 2 * float(lb[mo].sum()) / (mo.sum() + nref[o]))
    d['cross_species'] = xs
    res['taus'][f'{tau:g}'] = d
    print(f"tau={tau:g} ECE={cal['ECE']:.4f} AUROC={cal['AUROC']:.4f} AP={cal['AP']:.4f} thr={d['thr']:.4f} F1 {d['F1_unfiltered']:.2f}->{d['F1_filtered']:.2f} kept {d['kept']}/{cal['n']} meanP {d['mean_P']:.3f}", flush=True)
ok = {k: v for k, v in res['taus'].items() if len(v['bins']) > 0}
tau_star = min(ok, key=lambda k: ok[k]['ECE'])
grid = [f'{x:g}' for x in C.TAU_GRID if f'{x:g}' in ok]
res['selection'] = dict(rule='argmin pooled ECE over TAU_GRID on internal-validation chains', tau=float(tau_star), thr=ok[tau_star]['thr'],
                        at_grid_edge=(tau_star in (grid[0], grid[-1])), grid_evaluated=grid, ECE=ok[tau_star]['ECE'], AUROC=ok[tau_star]['AUROC'],
                        F1_filtered_internal=ok[tau_star]['F1_filtered'], F1_unfiltered_internal=ok[tau_star]['F1_unfiltered'])
json.dump(res, open(outp, 'w'), indent=1)
print('SELECTED', json.dumps(res['selection']))
