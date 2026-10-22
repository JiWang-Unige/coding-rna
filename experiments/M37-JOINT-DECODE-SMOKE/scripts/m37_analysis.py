"""M37 analysis (protocol §4).  DEV is never read.
 4.1 error independence: for every CDS-assessable primary reference of A (NC_089044.1, 1,601) that is fully mappable to B (all CDS exon boundaries
     mapped, one strand sign, co-linear), A-alone exact hit  x  exact hit of a B-alone chain projected into A coordinates -> 2x2 table, headroom
     = P(B-projected hit | A miss).  Symmetric direction: B references (CM000135.2 CDS-assessable primaries), A chains projected into B.
 4.2 potential-averaged joint decoding of A (lambda 0 / 0.25 / 0.5): project ruler exact-chain F1 vs A alone, 1-Mb block paired bootstrap (2,000
     draws, seed 0, as M31 c1_boot / M33 m33_dev), exon-count strata, CDS-base F1, background burden.
 4.3 descriptive: SNP / indel density (from the cs tags) in the gene spans of corrected / lost / unchanged references; mappable reference share.
Usage: m37_analysis.py [--out OUT/results.json]
"""
import argparse, json
from collections import defaultdict
from pathlib import Path
import numpy as np
import m37_common as C

M = C.register()
A, B = C.A_CONTIG, C.B_CONTIG


def load_run_chains(seq, chain_dir):
    by = {st: C.load_chain_json(Path(chain_dir) / f'{seq}_{st}.json')[0] for st in '+-'}
    return C.chain_dict(seq, by)


def project_exons(exons, strand, x2y, x2y_s, seq_y):
    """exact projection of a CDS chain through the base map; None unless every exon boundary maps with one strand sign and the result is co-linear"""
    sign = None; out = []
    for s, e in exons:
        ys, ye = int(x2y[s]), int(x2y[e - 1])
        if ys < 0 or ye < 0:
            return None
        sg = int(x2y_s[s])
        if sg != int(x2y_s[e - 1]) or (sign is not None and sg != sign):
            return None
        sign = sg
        if sg > 0:
            if ye < ys:
                return None
            out.append((ys, ye + 1))
        else:
            if ys < ye:
                return None
            out.append((ye, ys + 1))
    out.sort()
    if any(a[1] > b[0] for a, b in zip(out, out[1:])):
        return None
    st = strand if sign > 0 else ('-' if strand == '+' else '+')
    return (seq_y, st, tuple(out))


def all_bases_mapped(exons, x2y, x2y_s):
    sg = None
    for s, e in exons:
        seg = x2y[s:e]
        if (seg < 0).any():
            return False
        ss = np.unique(x2y_s[s:e])
        if len(ss) != 1 or (sg is not None and ss[0] != sg):
            return False
        sg = ss[0]
    return True


def two_by_two(ref_keys, hitA, hitB, mappable):
    n = {'both': 0, 'A_only': 0, 'B_only': 0, 'neither': 0}
    for k in ref_keys:
        if not mappable[k]:
            continue
        a, b = k in hitA, k in hitB
        n['both' if a and b else 'A_only' if a else 'B_only' if b else 'neither'] += 1
    N = sum(n.values())
    a_miss = n['B_only'] + n['neither']; b_miss = n['A_only'] + n['neither']
    a_rate = (n['both'] + n['A_only']) / N; b_rate = (n['both'] + n['B_only']) / N
    exp_B_only = a_miss * b_rate                     # expected count of (A miss, B hit) under independence
    oddsr = (n['both'] * n['neither']) / max(1, n['A_only'] * n['B_only'])
    phi = (n['both'] * n['neither'] - n['A_only'] * n['B_only']) / np.sqrt(max(1e-9, (n['both'] + n['A_only']) * (n['B_only'] + n['neither']) * (n['both'] + n['B_only']) * (n['A_only'] + n['neither'])))
    return dict(n_mappable=N, n_total=len(ref_keys), mappable_frac=N / len(ref_keys), table=n,
                A_recall=a_rate, B_projected_recall=b_rate, union_recall=(N - n['neither']) / N,
                headroom=n['B_only'] / a_miss if a_miss else None, headroom_other_direction=n['A_only'] / b_miss if b_miss else None,
                B_only_expected_if_independent=exp_B_only, B_only_ratio_obs_over_indep=n['B_only'] / exp_B_only if exp_B_only else None,
                odds_ratio=oddsr, phi=float(phi), A_miss=a_miss, B_miss=b_miss)


def score(chains, keys, nex_ref):
    tp = sum(1 for k in chains if k in keys)
    st = {s: [0, 0, 0] for s in ('1', '2-5', '>=6')}
    for k, v in chains.items():
        st[C.strat(v['nex'])][1] += 1
        if k in keys:
            st[C.strat(v['nex'])][0] += 1
    for k in keys:
        st[C.strat(nex_ref[k])][2] += 1
    return dict(tp=tp, pred=len(chains), ref=len(keys), F1=C.f1(tp, len(chains), len(keys)), precision=100 * tp / max(1, len(chains)), recall=100 * tp / len(keys),
                strata={s: dict(tp=v[0], pred=v[1], ref=v[2], F1=C.f1(*v)) for s, v in st.items()})


def base_metrics(chains, seq, ref):
    L = M.LENGTH[seq]
    cds = {'+': np.zeros(L, bool), '-': np.zeros(L, bool)}
    for t in list(ref['refs']) + list(ref['isoforms']):
        for a, b, _ in t['CDS']:
            cds[t['strand']][a:b] = True
    gene = np.zeros(L, bool)
    for a, b in ref['gene_spans']:
        gene[max(0, a - 1):b] = True
    pred = {'+': np.zeros(L, bool), '-': np.zeros(L, bool)}
    for v in chains.values():
        for a, b in v['exons']:
            pred[v['strand']][a:b] = True
    tp = pp = rp = 0
    for st in '+-':
        tp += int((pred[st] & cds[st]).sum()); pp += int(pred[st].sum()); rp += int(cds[st].sum())
    anyp = pred['+'] | pred['-']
    bg = int((anyp & ~gene).sum()); bgsize = int((~gene).sum())
    return dict(base_P=100 * tp / pp, base_R=100 * tp / rp, base_F1=100 * 2 * tp / (pp + rp), bg_pred_kb=bg / 1000, bg_burden_pct=100 * bg / bgsize)


def block_bootstrap(sets, keys, baseline, n_draw=2000, seed=0):
    blocks = set(k[2][0][0] // C.BLK for k in keys)
    for d in sets.values():
        blocks |= set(k[2][0][0] // C.BLK for k in d)
    bl = sorted(blocks); idx = {b: i for i, b in enumerate(bl)}
    rf = np.zeros(len(bl))
    for k in keys:
        rf[idx[k[2][0][0] // C.BLK]] += 1
    cnt = {}
    for n, d in sets.items():
        tp = np.zeros(len(bl)); pr = np.zeros(len(bl))
        for k in d:
            i = idx[k[2][0][0] // C.BLK]; pr[i] += 1; tp[i] += k in keys
        cnt[n] = (tp, pr)
    f1f = lambda tp, p, r: 200 * tp.sum() / (p.sum() + r.sum())
    rng = np.random.default_rng(seed)
    draws = {n: [] for n in sets}
    for _ in range(n_draw):
        ii = rng.integers(0, len(bl), len(bl))
        for n, (tp, pr) in cnt.items():
            draws[n].append(f1f(tp[ii], pr[ii], rf[ii]))
    draws = {n: np.array(v) for n, v in draws.items()}
    out = {n: dict(F1=f1f(cnt[n][0], cnt[n][1], rf), CI95=[float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]) for n, v in draws.items()}
    for n in sets:
        if n != baseline:
            d = draws[n] - draws[baseline]
            out[n]['diff_vs_' + baseline] = dict(point=out[n]['F1'] - out[baseline]['F1'], CI95=[float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))],
                                                 frac_positive=float((d > 0).mean()))
    out['blocks'] = len(bl); out['draws'] = n_draw; out['seed'] = seed
    return out


def density(ref_t, a_state, a_ins):
    s, e = ref_t['CDS'][0][0], ref_t['CDS'][-1][1]
    kb = (e - s) / 1000
    seg = a_state[s:e]
    return dict(span_kb=kb, mismatch_per_kb=float((seg == 2).sum() / kb), indel_per_kb=float(((seg == 3).sum() + a_ins[s:e].sum()) / kb),
                unmapped_frac=float((seg == 0).mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=str(C.OUT / 'results.json'))
    ap.add_argument('--joint', default=str(C.OUT / 'joint'))
    a = ap.parse_args()
    mp_ = np.load(C.OUT / 'align' / 'map.npz')
    a2b, a2b_s, b2a, b2a_s = mp_['a2b'], mp_['a2b_strand'], mp_['b2a'], mp_['b2a_strand']
    refA, keysA, nexA = C.ref_keys(A)
    refB, keysB, nexB = C.ref_keys(B)
    chA = load_run_chains(A, C.M33OUT / 'val' / 'chains')
    chB = load_run_chains(B, C.OUT / 'inferB' / 'chains')
    # cross-check: the M33 val GFF (the file M33 scored) gives the same unique chain set as the chains json
    gff_keys = set()
    for tx in M.load_txs(C.M33OUT / 'val' / A / 'prediction' / 'full_grammar.gff3', A):
        if sum(b - a for a, b, _ in tx['CDS']) >= 96:
            gff_keys.add(M.key(tx))
    assert gff_keys == set(chA), (len(gff_keys), len(chA), len(gff_keys ^ set(chA)))
    res = dict(A=A, B=B, n_refA=len(keysA), n_refB=len(keysB), n_chainsA=len(chA), n_chainsB=len(chB), map_stats=json.load(open(C.OUT / 'align' / 'map.json')))
    # ---- 4.1 error independence (A references; B chains projected into A) ----
    projB = {}
    for k, v in chB.items():
        pk = project_exons(v['exons'], v['strand'], b2a, b2a_s, A)
        if pk is not None:
            projB[pk] = k
    mappableA = {k: project_exons(k[2], k[1], a2b, a2b_s, B) is not None for k in keysA}
    allbaseA = {k: all_bases_mapped(k[2], a2b, a2b_s) for k in keysA}
    res['A_alone'] = score(chA, keysA, nexA); res['A_alone']['base'] = base_metrics(chA, A, refA)
    res['B_alone_on_B_refs'] = score(chB, keysB, nexB); res['B_alone_on_B_refs']['base'] = base_metrics(chB, B, refB)
    res['independence_A'] = two_by_two(keysA, set(chA), set(projB), mappableA)
    res['independence_A']['n_B_chains_projected'] = len(projB)
    res['independence_A']['n_refs_all_bases_mapped'] = int(sum(allbaseA.values()))
    res['independence_A']['two_by_two_all_bases_mapped'] = two_by_two(keysA, set(chA), set(projB), allbaseA)['table']
    byst = {}
    for s in ('1', '2-5', '>=6'):
        sub = [k for k in keysA if C.strat(nexA[k]) == s]
        byst[s] = two_by_two(sub, set(chA), set(projB), mappableA) if sub else None
    res['independence_A']['by_strata'] = byst
    # symmetric direction: B references, A chains projected into B
    projA = {}
    for k, v in chA.items():
        pk = project_exons(v['exons'], v['strand'], a2b, a2b_s, B)
        if pk is not None:
            projA[pk] = k
    mappableB = {k: project_exons(k[2], k[1], b2a, b2a_s, A) is not None for k in keysB}
    res['independence_B'] = two_by_two(keysB, set(chB), set(projA), mappableB)
    res['independence_B']['n_A_chains_projected'] = len(projA)
    # ---- 4.2 joint decoding ----
    sets = {'A_alone': chA}
    joint = {}
    for d in sorted(Path(a.joint).glob('lambda*')):
        lam = d.name[len('lambda'):]
        if not all((d / 'chains' / f'{A}_{st}.json').exists() for st in '+-'):
            continue
        ch = load_run_chains(A, d / 'chains')
        joint[lam] = ch; sets['lambda' + lam] = ch
        r = score(ch, keysA, nexA); r['base'] = base_metrics(ch, A, refA)
        r['identical_to_A_alone'] = set(ch) == set(chA)
        r['n_chains_in_mapped_bases_frac'] = float(np.mean([(a2b[v['exons'][0][0]] >= 0) for v in ch.values()]))
        res['joint_' + lam] = r
    res['bootstrap'] = block_bootstrap(sets, keysA, 'A_alone')
    # ---- 4.3 descriptive: which references change, and the local divergence ----
    a_state, a_ins = mp_['a_state'], mp_['a_ins']
    tx_of = {M.key(t): t for t in refA['refs']}
    desc = {}
    for lam, ch in joint.items():
        if lam == '0':
            continue
        cat = defaultdict(list)
        for k in keysA:
            c = ('hit' if k in chA else 'miss') + '->' + ('hit' if k in ch else 'miss')
            cat[c].append(k)
        d = {}
        for c, ks in cat.items():
            dens = [density(tx_of[k], a_state, a_ins) for k in ks]
            d[c] = dict(n=len(ks), mappable=int(sum(mappableA[k] for k in ks)), B_projected_hit=int(sum(k in projB for k in ks)),
                        mismatch_per_kb_median=float(np.median([x['mismatch_per_kb'] for x in dens])) if dens else None,
                        indel_per_kb_median=float(np.median([x['indel_per_kb'] for x in dens])) if dens else None,
                        unmapped_frac_mean=float(np.mean([x['unmapped_frac'] for x in dens])) if dens else None,
                        strata={s: sum(1 for k in ks if C.strat(nexA[k]) == s) for s in ('1', '2-5', '>=6')})
        new_fp = [k for k in ch if k not in keysA and k not in chA]
        gone_fp = [k for k in chA if k not in keysA and k not in ch]
        d['false_chains'] = dict(new=len(new_fp), removed=len(gone_fp), new_in_mapped_frac=float(np.mean([a2b[k[2][0][0]] >= 0 for k in new_fp])) if new_fp else None)
        desc[lam] = d
    res['transitions'] = desc
    all_dens = [density(tx_of[k], a_state, a_ins) for k in keysA]
    res['reference_divergence'] = dict(mismatch_per_kb_median=float(np.median([x['mismatch_per_kb'] for x in all_dens])),
                                       indel_per_kb_median=float(np.median([x['indel_per_kb'] for x in all_dens])),
                                       unmapped_frac_mean=float(np.mean([x['unmapped_frac'] for x in all_dens])),
                                       mappable_refs_frac=float(np.mean(list(mappableA.values()))), all_bases_mapped_frac=float(np.mean(list(allbaseA.values()))))
    json.dump(res, open(a.out, 'w'), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k not in ('transitions',)}, indent=1))
    print(json.dumps(res['transitions'], indent=1))


if __name__ == '__main__':
    main()
