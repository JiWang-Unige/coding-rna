"""C3 final calibration analysis (report 36; plan fixed before computing).  T2 gene-chain posteriors (tau = 1) on the core-10 panel, scored against the panel
references of report 31 (cached ruler objects), with Carbon-A's own confidence and the decode gain as comparators.  CPU only, no model.
usage: c3_calib.py <data_root> <out.json>
  data_root/M34-NTV2-SCALEUP/panel/post_T2/<species>/tau1/<contig>_<strand>.json
  data_root/PANEL-P6/_refcache/ref_<contig>.pkl
  data_root/PANEL-P6/carbonA/<species>/<contig>/prediction.gff
Sanity: per-species unfiltered exact-chain F1 of T2 and Carbon-A must equal experiments/PANEL-P6/results/p6_panel_table.tsv.
"""
import glob
import json
import pickle
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1]
sys.path.insert(0, str(EXP / 'M30-HOLDOUT-ADAPT' / 'scripts'))
import m30_common as M30                                 # noqa: E402  (a1_common: pipeline / M29 paths)
import m29_common as M29                                 # noqa: E402
from m29_r3_posterior_filter import reliability, auroc_ap, carbon_conf   # noqa: E402

CORE10 = ['dryas_octopetala', 'eschscholzia_californica', 'akebia_trifoliata', 'hakea_prostrata', 'zizania_latifolia', 'themeda_triandra', 'juncus_effusus',
          'rhynchospora_breviuscula', 'acorus_gramineus', 'amborella_trichopoda']
THR = -1.400270636077039                                  # M33 frozen thr* (log_p)
COVER = (1.0, .9, .8, .7, .5)
BLK, DRAWS, SEED = 1_000_000, 2000, 0


def key_of(seq, strand, exons):
    return (seq, strand, tuple((int(a), int(b)) for a, b in exons))


def main():
    root, outp = Path(sys.argv[1]), sys.argv[2]
    post_root = root / 'M34-NTV2-SCALEUP' / 'panel' / 'post_T2'
    rows = []                                            # (species, seq, key, log_p, gain, hit)
    carbon = []                                          # (species, seq, key, conf or None, hit)
    nref = {}
    refkeys = {}
    for sp in CORE10:
        files = sorted(glob.glob(str(post_root / sp / 'tau1' / '*.json')))
        seqs = sorted({Path(f).name.rsplit('_', 1)[0] for f in files})
        nref[sp] = 0
        for s in seqs:
            ref = pickle.load(open(root / 'PANEL-P6' / '_refcache' / f'ref_{s}.pkl', 'rb'))
            M29.SEQ2SPECIES[s], M29.LENGTH[s] = sp, ref['length']
            rk = {M29.key(t) for t in ref['refs']}
            refkeys[s] = rk
            nref[sp] += len(rk)
            seen = {}
            for st in '+-':
                d = json.load(open(post_root / sp / 'tau1' / f'{s}_{st}.json'))
                for g in d['genes']:
                    ex = [(int(a), int(b)) for a, b in g['exons']]
                    if sum(b - a for a, b in ex) < 96:
                        continue
                    k = key_of(s, st, ex)
                    if k not in seen:
                        seen[k] = (float(g['log_p']), float(g['gain']))
            for k, (lp, gn) in seen.items():
                rows.append((sp, s, k, lp, gn, int(k in rk)))
            gff = root / 'PANEL-P6' / 'carbonA' / sp / s / 'prediction.gff'
            conf = carbon_conf(str(gff))
            cseen = {}
            for tx in M29.load_txs(str(gff), s):
                ex = [(a, b) for a, b, _ in tx['CDS']]
                if sum(b - a for a, b in ex) < 96:
                    continue
                cseen.setdefault(M29.key(tx), conf[tx['gene_id']])
            for k, c in cseen.items():
                carbon.append((sp, s, k, c, int(k in rk)))
    res = dict(n_ref=nref, species=CORE10, per_species={}, pooled={})

    def f1(tp, pred, ref):
        return 200 * tp / max(pred + ref, 1)

    def summarize(sub, csub, ref):
        lp = np.array([r[3] for r in sub]); gn = np.array([r[4] for r in sub]); y = np.array([r[5] for r in sub])
        p = np.exp(lp)
        bins, ece = reliability(p, y)
        au, ap = auroc_ap(lp, y)
        aug, apg = auroc_ap(gn, y)
        hi = p >= .9
        out = dict(n=len(y), hits=int(y.sum()), F1=f1(y.sum(), len(y), ref), precision=100 * y.mean(), recall=100 * y.sum() / ref,
                   ECE=ece, Brier=float(np.mean((p - y) ** 2)), mean_P=float(p.mean()), AUROC=au, AP=ap, AUROC_gain=aug, AP_gain=apg, bins=bins,
                   hiconf=dict(n=int(hi.sum()), mean_P=float(p[hi].mean()) if hi.any() else None, hit_rate=float(y[hi].mean()) if hi.any() else None))
        order = np.argsort(-lp, kind='stable')
        cov = {}
        for c in COVER:
            m = int(round(c * len(y)))
            t = y[order[:m]].sum()
            cov[f'{int(c * 100)}'] = dict(kept=m, precision=100 * t / max(m, 1), recall=100 * t / ref, F1=f1(t, m, ref))
        out['coverage'] = cov
        keep = lp >= THR
        out['frozen_thr'] = dict(kept_frac=float(keep.mean()), F1=f1(y[keep].sum(), keep.sum(), ref), dF1=f1(y[keep].sum(), keep.sum(), ref) - out['F1'])
        cy = np.array([r[4] for r in csub]); cc = [r[3] for r in csub]
        known = np.array([c is not None for c in cc])
        cs = np.array([c if c is not None else -1.0 for c in cc])
        cout = dict(n=len(cy), hits=int(cy.sum()), F1=f1(cy.sum(), len(cy), ref), precision=100 * cy.mean(), conf_NA=int((~known).sum()))
        filt = (~known) | (cs >= .1)
        cout['filt_rule'] = dict(kept_frac=float(filt.mean()), precision=100 * cy[filt].mean(), F1=f1(cy[filt].sum(), filt.sum(), ref))
        if known.sum():
            cout['AUROC_conf_known'] = auroc_ap(cs[known], cy[known])[0]
            cb, cece = reliability(np.clip(cs[known], 0, 1), cy[known])
            cout['ECE_conf_known'] = cece
        corder = np.argsort(-cs, kind='stable')                     # NA ranked last
        for c in (.7, .5):
            m = int(round(c * len(cy)))
            cout[f'precision_at_{int(c * 100)}'] = 100 * cy[corder[:m]].sum() / max(m, 1)
        out['carbonA'] = cout
        return out

    by_sp = defaultdict(list); cby = defaultdict(list)
    for r in rows:
        by_sp[r[0]].append(r)
    for r in carbon:
        cby[r[0]].append(r)
    for sp in CORE10:
        res['per_species'][sp] = summarize(by_sp[sp], cby[sp], nref[sp])
    res['pooled'] = summarize(rows, carbon, sum(nref.values()))
    # ---- 1-Mb block bootstrap (pooled) ------------------------------------------------------------------------------------------------------
    blk = lambda s, k: (s, k[2][0][0] // BLK)
    bl = sorted({blk(r[1], r[2]) for r in rows} | {blk(r[1], r[2]) for r in carbon})
    idx = {b: i for i, b in enumerate(bl)}
    rb = np.array([idx[blk(r[1], r[2])] for r in rows]); cb_ = np.array([idx[blk(r[1], r[2])] for r in carbon])
    lp = np.array([r[3] for r in rows]); gn = np.array([r[4] for r in rows]); y = np.array([r[5] for r in rows])
    cy = np.array([r[4] for r in carbon]); cs = np.array([r[3] if r[3] is not None else -1.0 for r in carbon])
    by_r = [np.flatnonzero(rb == i) for i in range(len(bl))]; by_c = [np.flatnonzero(cb_ == i) for i in range(len(bl))]
    rng = np.random.default_rng(SEED)
    ece_d, dau_d, dprec_d = [], [], []

    def prec_at(score, lab, c):
        m = int(round(c * len(lab)))
        o = np.argsort(-score, kind='stable')
        return 100 * lab[o[:m]].sum() / max(m, 1)

    for _ in range(DRAWS):
        ii = rng.integers(0, len(bl), len(bl))
        r_ = np.concatenate([by_r[i] for i in ii]); c_ = np.concatenate([by_c[i] for i in ii])
        ece_d.append(reliability(np.exp(lp[r_]), y[r_])[1])
        dau_d.append(auroc_ap(lp[r_], y[r_])[0] - auroc_ap(gn[r_], y[r_])[0])
        dprec_d.append(prec_at(lp[r_], y[r_], .7) - prec_at(cs[c_], cy[c_], .7))
    ci = lambda v: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
    P = res['pooled']
    res['bootstrap'] = dict(blocks=len(bl), draws=DRAWS, ECE_CI95=ci(ece_d), AUROC_logp_minus_gain=P['AUROC'] - P['AUROC_gain'], AUROC_diff_CI95=ci(dau_d),
                            prec70_T2_minus_CarbonA=P['coverage']['70']['precision'] - P['carbonA']['precision_at_70'], prec70_diff_CI95=ci(dprec_d))
    json.dump(res, open(outp, 'w'), indent=1)
    print(f"{'species':26s} {'n':>6s} {'F1':>6s} {'caF1':>6s} {'ECE':>6s} {'Brier':>6s} {'AUROC':>6s} {'AUgain':>6s} {'hiP':>5s} {'hiHit':>5s} {'P@70':>6s} {'caP@70':>6s} {'dF1thr':>6s}")
    for sp in CORE10 + ['pooled']:
        r = res['per_species'].get(sp, res['pooled'])
        h = r['hiconf']
        print(f"{sp:26s} {r['n']:6d} {r['F1']:6.2f} {r['carbonA']['F1']:6.2f} {r['ECE']:6.3f} {r['Brier']:6.3f} {r['AUROC']:6.3f} {r['AUROC_gain']:6.3f} "
              f"{(h['mean_P'] or 0):5.2f} {(h['hit_rate'] or 0):5.2f} {r['coverage']['70']['precision']:6.2f} {r['carbonA']['precision_at_70']:6.2f} {r['frozen_thr']['dF1']:+6.2f}")
    print(json.dumps(res['bootstrap']))


if __name__ == '__main__':
    main()
