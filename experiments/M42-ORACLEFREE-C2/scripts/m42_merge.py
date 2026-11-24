"""M42 post-hoc merge analysis (protocol 35 §8.4; descriptive, rules fixed before running): OLD-release chains on E + model chains that
  G0: share no CDS base with any OLD CDS (either strand), or
  G1: G0 and the chain's CDS span lies entirely inside the E evidence-unknown region (miniprot hits cov >= .5, id >= .4, span +-2 kb; gate-0 V1 rule on E).
Scored with the same ruler / classes / block bootstrap as m42_score.py.
usage: m42_merge.py <species_dir> --old DIR --run NAME=DIR [--run ...] --out RESULT.json
"""
import argparse
import bisect
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.append(str(HERE))
import m42_score as S                                    # noqa: E402  (puts the M30 scripts on sys.path)
import m42_common as Z                                   # noqa: E402

EXPAND = 2000


def merged_intervals(iv):
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return ([x for x, _ in out], [y for _, y in out])


def inside(iv, a, b):
    s, e = iv
    i = bisect.bisect_right(s, a) - 1
    return i >= 0 and e[i] >= b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('species_dir')
    ap.add_argument('--old', required=True)
    ap.add_argument('--run', action='append', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    D = Path(args.species_dir)
    cfg = json.load(open(D / 'm42_config.json'))
    E = json.load(open(D / 'gate0' / 'split.json'))['E']
    import m30_score as SC
    Z.register(SC.M29)
    ev_ = SC.Evaluator(E)
    load = lambda d: {c: SC.load_gff_chains(str(Path(d) / c / 'prediction' / 'full_grammar.gff3'), c) for c in E}
    old = load(args.old)
    runs = {n: load(d) for n, d in (x.split('=', 1) for x in args.run)}
    oc, ochains = S.old_cds(D / cfg['old_gff'], E)
    plen = S.protein_lengths(cfg['queries'])
    raw = defaultdict(list)
    keep = set(E)
    for line in S.opn(D / cfg['miniprot_gff']):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'mRNA' or f[0] not in keep:
            continue
        a = S.attrs(f[8])
        t = a.get('Target', '').split()
        if len(t) < 3 or t[0] not in plen:
            continue
        if (int(t[2]) - int(t[1]) + 1) / plen[t[0]] >= S.EV_COV and float(a.get('Identity', 0)) >= S.EV_ID:
            raw[f[0]].append((max(0, int(f[3]) - 1 - EXPAND), int(f[4]) + EXPAND))
    unk = {c: merged_intervals(raw.get(c, [])) for c in E}
    unk_bp = sum(sum(y - x for x, y in zip(*unk[c])) for c in E)
    ev = S.evidence(D / cfg['miniprot_gff'], plen, E)
    cls = {}
    for c in E:
        for k in ev_.refkeys[c]:
            ex = k[2]
            if not any(S.overlaps(oc.get(c), a, b) for a, b in ex):
                cls[k] = 'novel_ev' if S.ev_hit(ev, c, ex[0][0], ex[-1][1]) else 'novel_noev'
            elif (c, k[1], tuple(ex)) in ochains:
                cls[k] = 'identical'
            else:
                cls[k] = 'revised'
    counts = defaultdict(int)
    for v in cls.values():
        counts[v] += 1
    n_novel = counts['novel_ev'] + counts['novel_noev']
    sets = {'OLD': old}
    added_info = {}
    for n, ch in runs.items():
        for g in ('G0', 'G1'):
            m, add_n, add_tp = {}, 0, 0
            for c in E:
                m[c] = dict(old[c])
                for k, v in ch[c].items():
                    if k in m[c] or any(S.overlaps(oc.get(c), a, b) for a, b in v['exons']):
                        continue
                    if g == 'G1' and not inside(unk[c], v['exons'][0][0], v['exons'][-1][1]):
                        continue
                    m[c][k] = v
                    add_n += 1
                    add_tp += k in ev_.refkeys[c]
            name = f'{n}+{g}'
            sets[name] = m
            added_info[name] = dict(added=add_n, added_tp=add_tp, added_precision=100 * add_tp / max(add_n, 1))
    bl = set()
    for c in E:
        bl.update(S.block(c, k) for k in ev_.refkeys[c])
        for ch in sets.values():
            bl.update(S.block(c, k) for k in ch[c])
    bl = sorted(bl)
    idx = {b: i for i, b in enumerate(bl)}
    nb = len(bl)
    rf, nrf = np.zeros(nb), np.zeros(nb)
    for c in E:
        for k in ev_.refkeys[c]:
            rf[idx[S.block(c, k)]] += 1
            nrf[idx[S.block(c, k)]] += cls[k].startswith('novel')
    res = dict(E_unknown_region_bp=unk_bp, ref_classes=dict(counts), novel=n_novel, per_set={}, blocks=nb)
    vec = {}
    for n, ch in sets.items():
        t = ev_.add(n, ch)
        tp, pr, ntp = np.zeros(nb), np.zeros(nb), np.zeros(nb)
        rec = defaultdict(int)
        for c in E:
            for k in ch[c]:
                i = idx[S.block(c, k)]
                pr[i] += 1
                if k in ev_.refkeys[c]:
                    tp[i] += 1
                    rec[cls[k]] += 1
                    ntp[i] += cls[k].startswith('novel')
        vec[n] = (tp, pr, ntp)
        nov = rec['novel_ev'] + rec['novel_noev']
        res['per_set'][n] = dict(F1=t['F1'], precision=100 * t['tp'] / max(t['pred'], 1), recall=100 * t['tp'] / max(t['ref'], 1), pred=t['pred'], tp=t['tp'],
                                 novel_recall=100 * nov / max(n_novel, 1), recall_by_class={k: 100 * rec[k] / max(counts[k], 1) for k in counts},
                                 **added_info.get(n, {}))
    rng = np.random.default_rng(S.SEED)
    f1 = lambda v, ii: 200 * v[0][ii].sum() / max(v[1][ii].sum() + rf[ii].sum(), 1)
    nr = lambda v, ii: 100 * v[2][ii].sum() / max(nrf[ii].sum(), 1)
    pairs = [('U+G1', 'OLD'), ('U+G1', 'Z+G1'), ('U+G1', 'F+G1'), ('U+G0', 'OLD'), ('Z+G1', 'OLD')]
    pairs = [p for p in pairs if p[0] in vec and p[1] in vec]
    dr = {p: (np.empty(S.DRAWS), np.empty(S.DRAWS)) for p in pairs}
    for d in range(S.DRAWS):
        ii = rng.integers(0, nb, nb)
        for p in pairs:
            dr[p][0][d] = f1(vec[p[0]], ii) - f1(vec[p[1]], ii)
            dr[p][1][d] = nr(vec[p[0]], ii) - nr(vec[p[1]], ii)
    full = np.arange(nb)
    ci = lambda x: [float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))]
    res['paired'] = {f'{a} - {b}': dict(F1=f1(vec[a], full) - f1(vec[b], full), F1_CI95=ci(dr[(a, b)][0]),
                                        novel_recall=nr(vec[a], full) - nr(vec[b], full), novel_recall_CI95=ci(dr[(a, b)][1])) for a, b in pairs}
    json.dump(res, open(args.out, 'w'), indent=1)
    print('E unknown region bp', unk_bp, 'novel', n_novel, dict(counts))
    for n, r in res['per_set'].items():
        print(f"{n:8s} F1 {r['F1']:.2f} P {r['precision']:.2f} R {r['recall']:.2f} novel {r['novel_recall']:.2f} ident {r['recall_by_class']['identical']:.1f} "
              f"nov_ev {r['recall_by_class']['novel_ev']:.1f} nov_noev {r['recall_by_class']['novel_noev']:.1f} added {r.get('added', '-')} add_prec {r.get('added_precision', float('nan')):.1f}")
    for k, v in res['paired'].items():
        print(k, round(v['F1'], 2), [round(x, 2) for x in v['F1_CI95']], round(v['novel_recall'], 2), [round(x, 2) for x in v['novel_recall_CI95']])


if __name__ == '__main__':
    main()
