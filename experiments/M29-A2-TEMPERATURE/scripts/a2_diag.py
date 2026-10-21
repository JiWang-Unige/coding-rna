"""A2 step 1 (CPU only): saturation diagnostic of the CRF target at inverse temperature tau, on the GATE model (head = 0).
E_tau = tau * (offset + head).  For N stratified TRAIN windows with observations and each tau:
  p(y_obs) = exp(-loss), gradient norm wrt the head logits, fraction of windows with (numerically) zero gradient,
  gradient mass on windows that the offset already decodes correctly (p~1) -> 'solved windows still give gradient'.
Usage: a2_diag.py N taus(comma) workers
"""
import json, sys, time
import numpy as np
import a2_common as C2
import a1_common as K
import a1_model as MD, a1_offset as O, a1_dp as DP, a1_train as T
import prepare_full_grammar_fresh_draws as PP, fit_full_grammar_objective_contrast as C

def main():
    N, taus, workers = int(sys.argv[1]), [float(x) for x in sys.argv[2].split(',')], int(sys.argv[3])
    plan = [json.loads(x) for x in (K.OUTS / 'M28-FULL-GRAMMAR-FRESH-DRAWS-PREPARATION/training_plan.jsonl').read_text().splitlines()]
    windows, strata, catalogs, original = PP.load_population()
    sequences, _ = C.P.read_train_DNA()
    tracks = {s: O.load_tracks(s) for s in {w['seqid'] for w in windows.values()}}
    rng = np.random.default_rng(0)
    cand = [r for r in plan[:4608] if r['fine_positions'] + r['coarse_positions'] > 0]
    seen, pick = set(), []
    for i in rng.permutation(len(cand)):
        r = cand[i]
        if r['window_id'] in seen: continue
        seen.add(r['window_id']); pick.append(r)
        if len(pick) == N: break
    jobs, meta = [], []
    for r in pick:
        w = windows[r['window_id']]
        chromosome, geometry, context, offset, targets, counts = T.example(w, sequences, catalogs)
        a, b = geometry['genomic_context']
        off = O.context_offset(tracks[w['seqid']], w['strand'], a, b, len(chromosome))
        tg = {k: np.asarray(targets[k]) for k in ('seg', 'region')}
        seg = tg['seg'][:w['valid_bases']]
        meta.append(dict(window=w['id'], species=w['species'], seqid=w['seqid'], known=counts[0] + counts[1], gene_roles=int((seg >= 1).sum())))
        jobs.append((context, off, offset, offset + w['valid_bases'], tg))
    pool = DP.DPPool(workers, C2.CPUS)
    out = dict(N=len(jobs), taus={})
    for tau in taus:
        t0 = time.time()
        futs = [pool.submit(c, (tau * lg).astype(np.float32), cs, ce, tg) for c, lg, cs, ce, tg in jobs]
        rows = []
        for m, f in zip(meta, futs):
            r = f.get()
            g = tau * r['grad'] / max(1, m['known'])             # gradient of the per-observation loss wrt the (unscaled) head logits
            rows.append(dict(m, loss=r['loss'], p=float(np.exp(-r['loss'])), gnorm=float(np.linalg.norm(g)), gmax=float(np.abs(r['grad']).max() * tau)))
        L = np.array([r['loss'] for r in rows]); P = np.array([r['p'] for r in rows]); G = np.array([r['gnorm'] for r in rows])
        has_gene = np.array([r['gene_roles'] > 0 for r in rows])
        out['taus'][str(tau)] = {
            'mean_p_label': float(P.mean()), 'median_p_label': float(np.median(P)), 'frac_p_gt_0.99': float((P > .99).mean()),
            'frac_p_lt_0.01': float((P < .01).mean()), 'frac_zero_grad': float((G < 1e-12).mean()),
            'mean_gnorm': float(G.mean()), 'median_gnorm': float(np.median(G)),
            'mean_loss_per_obs': float(np.mean([r['loss'] / max(1, r['known']) for r in rows])),
            'gene_windows': int(has_gene.sum()), 'frac_zero_grad_gene_windows': float((G[has_gene] < 1e-12).mean()),
            'frac_zero_grad_genefree': float((G[~has_gene] < 1e-12).mean()) if (~has_gene).any() else None,
            'gnorm_share_from_windows_p_gt_0.5': float(G[P > .5].sum() / max(G.sum(), 1e-30)), 'seconds': time.time() - t0, 'rows': rows}
        s = out['taus'][str(tau)]
        print(f"tau {tau}: mean p(y)={s['mean_p_label']:.3f} med={s['median_p_label']:.3f} P>.99:{s['frac_p_gt_0.99']:.2f} P<.01:{s['frac_p_lt_0.01']:.2f} zero-grad:{s['frac_zero_grad']:.2f} (genes {s['frac_zero_grad_gene_windows']:.2f}, genefree {s['frac_zero_grad_genefree']}) mean|g|={s['mean_gnorm']:.4g} share_g_from_p>.5:{s['gnorm_share_from_windows_p_gt_0.5']:.2f} ({s['seconds']:.0f}s)", flush=True)
    pool.close()
    C2.A2.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(C2.A2 / 'diag_windows.json', 'w'))

if __name__ == '__main__':
    main()
