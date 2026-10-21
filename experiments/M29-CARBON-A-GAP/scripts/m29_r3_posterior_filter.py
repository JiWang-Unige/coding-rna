"""R3: posterior / confidence filtering curves on the six original TRAIN contigs (CPU, read-only).

Predictors compared on the same ruler (unique exact CDS chain = hit iff it equals a primary reference key; denominator 22,688):
  * EM_pool posterior  (log P(g included), whole-contig gene-set partition function; summary.json of ENERGY-POOL-POSTERIOR-TRAIN)
  * EM_pool internal gain (GFF/summary 'gain')
  * Carbon-A confidence (GFF gene attribute `confidence`; NA = edge-truncated gene, kept by the 0.1-or-NA operating point)
Outputs  outputs/M29-CARBON-A-GAP/R3/{summary.json, curve_*.tsv, reliability_*.tsv, figure png}
"""
import json
import math
import sys
from collections import defaultdict
import numpy as np
import m29_common as M

OUT = M.M29 / 'R3'
OUT.mkdir(parents=True, exist_ok=True)
POSTERIOR = M.OUTS / 'M28-FULL-GRAMMAR-ENERGY-POOL-POSTERIOR-TRAIN'


def carbon_conf(path):
    res = {}
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        r = line.rstrip('\n').split('\t')
        if r[2] != 'gene':
            continue
        a = dict(x.split('=', 1) for x in r[8].split(';'))
        res[a['ID']] = None if a['confidence'] == 'NA' else float(a['confidence'])
    return res


def auroc_ap(score, label):
    """AUROC (tie-aware, rank based) and AP over a prediction pool."""
    s, y = np.asarray(score, float), np.asarray(label, int)
    order = np.argsort(-s, kind='stable')
    s, y = s[order], y[order]
    # AP with whole-tie groups
    ap, tp, n = 0.0, 0, 0
    npos = int(y.sum())
    i = 0
    while i < len(s):
        j = i
        while j < len(s) and s[j] == s[i]:
            j += 1
        tp_new = tp + int(y[i:j].sum())
        n_new = n + (j - i)
        # tie group: precision at group end weighted by positives gained
        ap += (tp_new - tp) * (tp_new / n_new)
        tp, n = tp_new, n_new
        i = j
    ap /= npos
    from scipy.stats import rankdata
    r = rankdata(s)  # average ranks ascending
    pos = r[y == 1].sum()
    nn = int(y.sum()); nneg = len(y) - nn
    auc = (pos - nn * (nn + 1) / 2) / (nn * nneg)
    return float(auc), float(ap)


def curve(score, label, n_ref, thresholds=None):
    """F1 for 'keep predictions with score >= t' for every distinct t (ties kept together)."""
    s, y = np.asarray(score, float), np.asarray(label, int)
    order = np.argsort(-s, kind='stable')
    s, y = s[order], y[order]
    ctp = np.cumsum(y)
    ends = np.r_[np.where(np.diff(s) != 0)[0], len(s) - 1]  # last index of each tie group
    kept = ends + 1
    tp = ctp[ends]
    f1 = 2 * tp / (kept + n_ref)
    return s[ends], kept, tp, f1


def best_threshold_pick(rows_train, n_ref_train):
    """argmax-F1 threshold on a training set of (score,label)."""
    sc = np.array([r[0] for r in rows_train]); lb = np.array([r[1] for r in rows_train])
    t, kept, tp, f1 = curve(sc, lb, n_ref_train)
    k = int(np.argmax(f1))
    return float(t[k]), float(f1[k])


def apply_threshold(rows, thr, n_ref):
    tp = sum(l for s, l in rows if s >= thr)
    kept = sum(1 for s, l in rows if s >= thr)
    return tp, kept, 100 * 2 * tp / (kept + n_ref)


def reliability(prob, label, bins=10):
    p, y = np.asarray(prob, float), np.asarray(label, float)
    edges = np.linspace(0, 1, bins + 1)
    idx = np.minimum((p * bins).astype(int), bins - 1)
    out, ece = [], 0.0
    for b in range(bins):
        m = idx == b
        if m.sum():
            out.append(dict(bin=f'{edges[b]:.1f}-{edges[b+1]:.1f}', n=int(m.sum()), mean_pred=float(p[m].mean()), hit_rate=float(y[m].mean())))
            ece += m.sum() / len(p) * abs(p[m].mean() - y[m].mean())
    return out, float(ece)


def main():
    pools = {'posterior': [], 'gain': [], 'carbonA_conf': [], 'carbonA_conf_NA1': []}
    meta = []  # (seq, species)
    per = defaultdict(lambda: defaultdict(list))   # tool -> seq -> [(score,label)]
    nref = {}
    for sp, seq, L in M.SCOPE:
        ref = M.build_reference(seq)
        rk = {M.key(t) for t in ref['refs']}
        nref[seq] = len(rk)
        # --- EM pool posterior
        post = json.loads((POSTERIOR / seq / 'posterior/summary.json').read_text())
        txs = {t['id']: t for t in M.load_txs(M.gff_path('EM_pool', seq), seq)}
        seen = {}
        for g in post['genes']:
            tx = txs[g['id']]
            k = M.key(tx)
            lab = int(k in rk)
            # unique chains: keep first record
            if k in seen:
                continue
            seen[k] = 1
            lp = g['log_probability']
            per['posterior'][seq].append((math.exp(lp), lab))
            per['gain'][seq].append((g['gain'], lab))
        assert len(seen) == len(post['genes']), (seq, len(seen), len(post['genes']))
        # --- Carbon-A confidence
        gff = M.gff_path('CarbonA_unfiltered', seq)
        conf = carbon_conf(gff)
        ctx = M.load_txs(gff, seq)
        seen = {}
        for tx in ctx:
            k = M.key(tx)
            if k in seen:
                continue
            seen[k] = 1
            c = conf[tx['gene_id']]
            lab = int(k in rk)
            per['carbonA_conf'][seq].append((-1.0 if c is None else c, lab, c is None))
            per['carbonA_conf_NA1'][seq].append((2.0 if c is None else c, lab))
        print(seq, 'posterior pool', len(post['genes']), 'carbonA unique', len(seen), flush=True)
    seqs = [s for _, s, _ in M.SCOPE]
    species_of = M.SEQ2SPECIES
    R = {'n_ref_primary': sum(nref.values()), 'predictors': {}}

    def pooled(tool):
        rows = []
        for s in seqs:
            rows += [(r[0], r[1]) for r in per[tool][s]]
        return rows

    # ---- ranking quality (prediction-pool AUROC / AP, primary-exact label)
    for tool in ('posterior', 'gain', 'carbonA_conf_NA1'):
        rows = pooled(tool)
        sc = [r[0] for r in rows]; lb = [r[1] for r in rows]
        auc, ap = auroc_ap(sc, lb)
        sp_stats = {}
        for sp in ('arabidopsis_thaliana', 'oryza_sativa'):
            rr = [(r[0], r[1]) for s in seqs if species_of[s] == sp for r in per[tool][s]]
            a2, p2 = auroc_ap([x[0] for x in rr], [x[1] for x in rr])
            sp_stats[sp] = dict(auroc=a2, ap=p2, n=len(rr), hits=sum(x[1] for x in rr))
        R['predictors'][tool] = dict(n_pred=len(rows), hits=int(sum(lb)), precision_all=sum(lb) / len(lb),
                                     auroc=auc, ap=ap, per_species=sp_stats,
                                     F1_unfiltered=100 * 2 * sum(lb) / (len(rows) + R['n_ref_primary']))
    # Carbon-A: excluding NA genes (not rankable by confidence)
    rr = [(r[0], r[1]) for s in seqs for r in per['carbonA_conf'][s] if not r[2]]
    a, p = auroc_ap([x[0] for x in rr], [x[1] for x in rr])
    R['predictors']['carbonA_conf_excluding_NA'] = dict(n_pred=len(rr), hits=sum(x[1] for x in rr), auroc=a, ap=p,
                                                        n_NA=sum(1 for s in seqs for r in per['carbonA_conf'][s] if r[2]),
                                                        NA_hits=sum(r[1] for s in seqs for r in per['carbonA_conf'][s] if r[2]))

    # ---- F1 vs retained fraction / threshold curves
    nref_all = R['n_ref_primary']
    curves = {}
    for tool in ('posterior', 'gain', 'carbonA_conf_NA1'):
        rows = pooled(tool)
        t, kept, tp, f1 = curve([r[0] for r in rows], [r[1] for r in rows], nref_all)
        curves[tool] = (t, kept, tp, f1)
        with open(OUT / f'curve_{tool}.tsv', 'w') as f:
            f.write('threshold\tkept\ttp\tF1_pct\tprecision_pct\trecall_pct\tkept_fraction\n')
            n = len(rows)
            for a, b, c, d in zip(t, kept, tp, f1):
                f.write(f'{a:.8g}\t{b}\t{c}\t{100*d:.4f}\t{100*c/b:.4f}\t{100*c/nref_all:.4f}\t{b/n:.5f}\n')
        k = int(np.argmax(f1))
        R['predictors'][tool]['best_posthoc'] = dict(threshold=float(t[k]), kept=int(kept[k]), tp=int(tp[k]), F1=100 * float(f1[k]))
        # F1 at fixed retained fractions
        fr = {}
        n = len(rows)
        for q in (0.5, 0.75, 0.8, 0.85, 0.9, 0.93, 0.95, 0.97, 1.0):
            # first tie-group end with kept >= q n (tie groups kept together, as in the protocol)
            i = int(np.searchsorted(kept, q * n))
            i = min(i, len(kept) - 1)
            fr[str(q)] = dict(kept=int(kept[i]), tp=int(tp[i]), F1=100 * float(f1[i]), thr=float(t[i]))
        R['predictors'][tool]['F1_at_retained_fraction'] = fr

    # Carbon-A operating point actually used by the tool (conf>=0.1 or NA): reproduce
    rows = [(r[0], r[1]) for s in seqs for r in per['carbonA_conf_NA1'][s]]
    tp, kept, f1 = apply_threshold(rows, 0.1, nref_all)
    R['carbonA_conf_ge_0.1_or_NA'] = dict(tp=tp, kept=kept, F1=f1)
    # equal-retention head-to-head: posterior kept = same #predictions as Carbon-A filter
    t, kp, tp_, f1_ = curves['posterior']
    i = int(np.searchsorted(kp, kept))
    R['posterior_at_same_kept_as_CarbonA_filter'] = dict(kept=int(kp[i]), tp=int(tp_[i]), F1=100 * float(f1_[i]))

    # ---- threshold transfer (realizable gain): fit on one set of contigs, apply on the held-out ones
    transfer = {}
    for tool in ('posterior', 'carbonA_conf_NA1'):
        res = {}
        # leave-one-contig-out
        tot_tp = tot_kept = 0
        thrs = {}
        for held in seqs:
            train_rows = [(r[0], r[1]) for s in seqs if s != held for r in per[tool][s]]
            n_tr = sum(nref[s] for s in seqs if s != held)
            thr, f1_tr = best_threshold_pick(train_rows, n_tr)
            test_rows = [(r[0], r[1]) for r in per[tool][held]]
            tp, kept_, _ = apply_threshold(test_rows, thr, nref[held])
            tot_tp += tp; tot_kept += kept_; thrs[held] = thr
        res['leave_one_contig_out'] = dict(thresholds=thrs, tp=tot_tp, kept=tot_kept, F1=100 * 2 * tot_tp / (tot_kept + nref_all))
        # cross-species
        for src, dst in (('arabidopsis_thaliana', 'oryza_sativa'), ('oryza_sativa', 'arabidopsis_thaliana')):
            train_rows = [(r[0], r[1]) for s in seqs if species_of[s] == src for r in per[tool][s]]
            n_tr = sum(nref[s] for s in seqs if species_of[s] == src)
            thr, f1_tr = best_threshold_pick(train_rows, n_tr)
            dst_seqs = [s for s in seqs if species_of[s] == dst]
            test_rows = [(r[0], r[1]) for s in dst_seqs for r in per[tool][s]]
            n_te = sum(nref[s] for s in dst_seqs)
            tp, kept_, f1v = apply_threshold(test_rows, thr, n_te)
            base_tp = sum(l for _, l in test_rows)
            res[f'fit_{src}_apply_{dst}'] = dict(threshold=thr, F1_on_fit=100 * f1_tr, tp=tp, kept=kept_, F1_applied=f1v,
                                                 F1_unfiltered_dst=100 * 2 * base_tp / (len(test_rows) + n_te))
        transfer[tool] = res
    R['threshold_transfer'] = transfer

    # ---- per-species/contig unfiltered vs best post-hoc, for the record
    # ---- reliability (calibration)
    rel = {}
    rows = pooled('posterior')
    bins, ece = reliability([r[0] for r in rows], [r[1] for r in rows])
    rel['posterior_pooled'] = dict(bins=bins, ECE=ece)
    rows = [(r[0], r[1]) for s in seqs for r in per['carbonA_conf'][s] if not r[2]]
    bins, ece = reliability([r[0] for r in rows], [r[1] for r in rows])
    rel['carbonA_conf_pooled_excluding_NA'] = dict(bins=bins, ECE=ece)
    for sp in ('arabidopsis_thaliana', 'oryza_sativa'):
        rr = [(r[0], r[1]) for s in seqs if species_of[s] == sp for r in per['posterior'][s]]
        rel[f'posterior_{sp}'] = dict(zip(('bins', 'ECE'), reliability([x[0] for x in rr], [x[1] for x in rr])))
        rr = [(r[0], r[1]) for s in seqs if species_of[s] == sp for r in per['carbonA_conf'][s] if not r[2]]
        rel[f'carbonA_{sp}'] = dict(zip(('bins', 'ECE'), reliability([x[0] for x in rr], [x[1] for x in rr])))
    R['reliability'] = rel
    with open(OUT / 'summary.json', 'w') as f:
        json.dump(R, f, indent=1)
    # ---- figure
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
        n = R['predictors']['posterior']['n_pred']
        for tool, lab in (('posterior', 'E/M pool posterior'), ('gain', 'E/M pool gain'), ('carbonA_conf_NA1', 'Carbon-A confidence (NA kept first)')):
            t, kp_, tp_, f1_ = curves[tool]
            ax[0].plot(kp_ / R['predictors'][tool]['n_pred'], 100 * f1_, label=lab)
        ax[0].axhline(R['predictors']['posterior']['F1_unfiltered'], ls=':', c='gray')
        ax[0].set_xlabel('retained fraction'); ax[0].set_ylabel('exact-chain F1 (%) pooled'); ax[0].legend(fontsize=7); ax[0].set_title('F1 vs retained fraction')
        for tool, lab in (('posterior', 'E/M pool posterior'), ('carbonA_conf_NA1', 'Carbon-A confidence')):
            t, kp_, tp_, f1_ = curves[tool]
            ax[1].plot(kp_, 100 * tp_ / kp_, label=lab)
        ax[1].set_xlabel('# predictions retained'); ax[1].set_ylabel('precision (%)'); ax[1].legend(fontsize=7); ax[1].set_title('precision vs retained #')
        for key, lab in (('posterior_pooled', 'E/M pool posterior'), ('carbonA_conf_pooled_excluding_NA', 'Carbon-A confidence')):
            b = rel[key]['bins']
            ax[2].plot([x['mean_pred'] for x in b], [x['hit_rate'] for x in b], 'o-', label=f"{lab} (ECE {rel[key]['ECE']:.3f})")
        ax[2].plot([0, 1], [0, 1], 'k:'); ax[2].set_xlabel('mean predicted'); ax[2].set_ylabel('primary-exact hit rate'); ax[2].legend(fontsize=7); ax[2].set_title('reliability')
        fig.tight_layout(); fig.savefig(OUT / 'r3_filter_curves.png', dpi=130)
    except Exception as e:
        print('figure failed', e)
    print(json.dumps({k: R[k] for k in R if k != 'reliability'}, indent=1)[:6000])


if __name__ == '__main__':
    main()
