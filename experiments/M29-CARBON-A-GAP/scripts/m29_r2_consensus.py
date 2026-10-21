"""R2: tool-consensus upper bound and >=k-tool voting F1 on the six original TRAIN contigs (CPU)."""
import json
from collections import Counter, defaultdict
from itertools import combinations
import m29_common as M
from m29_load_preds import tool_chains

OUT = M.M29 / 'R2'
OUT.mkdir(parents=True, exist_ok=True)
TOOLS = ['CarbonA_unfiltered', 'ANNEVO', 'Tiberius', 'Helixer', 'B_freshjoint', 'E_lora', 'M_denseqv', 'EM_pool']
CORE4 = ['CarbonA', 'ANNEVO', 'Tiberius', 'E_lora']


def overlap(a, b):
    return a[0] < b[1] and b[0] < a[1]


def greedy_nonoverlap(cands, rank):
    """cands: list of chain dicts; keep best-ranked chains with no same-strand span overlap."""
    kept = []
    by_strand = defaultdict(list)
    for c in sorted(cands, key=rank):
        span = (c['exons'][0][0], c['exons'][-1][1])
        if any(overlap(span, s) for s in by_strand[c['strand']]):
            continue
        by_strand[c['strand']].append(span)
        kept.append(c)
    return kept


def main():
    data = {}      # seq -> tool -> {key: rec}
    refs = {}
    for sp, seq, L in M.SCOPE:
        r = M.build_reference(seq)
        refs[seq] = {M.key(t) for t in r['refs']}
        data[seq] = {t: tool_chains(t, seq) for t in TOOLS}
        # filtered Carbon-A = confidence>=0.1 or NA
        data[seq]['CarbonA'] = {k: v for k, v in data[seq]['CarbonA_unfiltered'].items() if v['na'] or v['conf'] >= 0.1}
    nref = sum(len(v) for v in refs.values())
    species = M.SEQ2SPECIES
    R = {'n_ref': nref, 'single_tool': {}, 'union_upper_bound': {}, 'vote': {}, 'vote_nonoverlap': {}, 'hybrid': {}}

    def score(selected_per_seq, label=None):
        tp = pred = 0
        per_sp = defaultdict(lambda: [0, 0, 0])
        for sp, seq, L in M.SCOPE:
            s = selected_per_seq[seq]
            t = sum(1 for k in s if k in refs[seq])
            tp += t; pred += len(s)
            per_sp[sp][0] += t; per_sp[sp][1] += len(s); per_sp[sp][2] += len(refs[seq])
        out = dict(tp=tp, pred=pred, precision=tp / pred if pred else 0, recall=tp / nref, F1=100 * 2 * tp / (pred + nref))
        out['per_species_F1'] = {s: 100 * 2 * v[0] / (v[1] + v[2]) for s, v in per_sp.items()}
        return out

    all_tools = ['CarbonA', 'CarbonA_unfiltered'] + [t for t in TOOLS if t != 'CarbonA_unfiltered']
    for t in all_tools:
        R['single_tool'][t] = score({seq: set(data[seq][t]) for _, seq, _ in M.SCOPE})

    # ---- union upper bound (recall ceiling): reference hit by >=1 tool of a set
    def hit_sets(toolset):
        res = {}
        for seq in refs:
            hit = set()
            for t in toolset:
                hit |= set(data[seq][t]) & refs[seq]
            res[seq] = hit
        return res
    sets = {
        'core4=CarbonA_filtered+ANNEVO+Tiberius+E': ['CarbonA', 'ANNEVO', 'Tiberius', 'E_lora'],
        'core4_unfilteredCA': ['CarbonA_unfiltered', 'ANNEVO', 'Tiberius', 'E_lora'],
        'core3_noE=CarbonA+ANNEVO+Tiberius': ['CarbonA_unfiltered', 'ANNEVO', 'Tiberius'],
        'ANNEVO+Tiberius+E': ['ANNEVO', 'Tiberius', 'E_lora'],
        'ours=B+E+M+EMpool': ['B_freshjoint', 'E_lora', 'M_denseqv', 'EM_pool'],
        'CarbonA_unfiltered+E': ['CarbonA_unfiltered', 'E_lora'],
        'CarbonA_unfiltered+ANNEVO': ['CarbonA_unfiltered', 'ANNEVO'],
        'all8': TOOLS,
    }
    for name, ts in sets.items():
        h = hit_sets(ts)
        n = sum(len(v) for v in h.values())
        sp = {}
        for s in ('arabidopsis_thaliana', 'oryza_sativa'):
            nn = sum(len(h[q]) for q in refs if species[q] == s); rr = sum(len(refs[q]) for q in refs if species[q] == s)
            sp[s] = dict(hits=nn, ref=rr, recall=nn / rr)
        R['union_upper_bound'][name] = dict(tools=ts, union_hits=n, recall_ceiling=n / nref,
                                            F1_if_union_hits_and_zero_FP=100 * 2 * n / (n + nref), per_species=sp)
    # reference never hit by any tool
    h = hit_sets(TOOLS)
    miss_all = sum(len(refs[s] - h[s]) for s in refs)
    R['never_hit_by_any_of_8'] = miss_all
    # per-ref hit-count histogram on core4 (CarbonA unfiltered, ANNEVO, Tiberius, E)
    hist = Counter()
    for seq in refs:
        for k in refs[seq]:
            hist[sum(k in data[seq][t] for t in ('CarbonA_unfiltered', 'ANNEVO', 'Tiberius', 'E_lora'))] += 1
    R['ref_hit_count_hist_core4_unfilteredCA'] = dict(sorted(hist.items()))

    # ---- voting: chains predicted by >=k tools of a set (unique-chain votes)
    vote_sets = {
        'core4(CA_filt,ANNEVO,Tib,E)': ['CarbonA', 'ANNEVO', 'Tiberius', 'E_lora'],
        'core4(CA_unfilt,ANNEVO,Tib,E)': ['CarbonA_unfiltered', 'ANNEVO', 'Tiberius', 'E_lora'],
        'core3_external(CA_unfilt,ANNEVO,Tib)': ['CarbonA_unfiltered', 'ANNEVO', 'Tiberius'],
        'core3_external_filt(CA_filt,ANNEVO,Tib)': ['CarbonA', 'ANNEVO', 'Tiberius'],
        'core5(+Helixer)': ['CarbonA_unfiltered', 'ANNEVO', 'Tiberius', 'E_lora', 'Helixer'],
        'ANNEVO,Tib,E': ['ANNEVO', 'Tiberius', 'E_lora'],
        'CA_unfilt,E,ANNEVO': ['CarbonA_unfiltered', 'E_lora', 'ANNEVO'],
        'ours(B,E,M)': ['B_freshjoint', 'E_lora', 'M_denseqv'],
    }
    for name, ts in vote_sets.items():
        votes = {seq: Counter() for seq in refs}
        recs = {seq: {} for seq in refs}
        for seq in refs:
            for t in ts:
                for k, v in data[seq][t].items():
                    votes[seq][k] += 1
                    recs[seq].setdefault(k, v)
        res = {}
        res_no = {}
        for kmin in range(1, len(ts) + 1):
            sel = {seq: {k for k, c in votes[seq].items() if c >= kmin} for seq in refs}
            res[f'>={kmin}'] = score(sel)
            # non-overlapping version: rank = (-votes, CarbonA-membership first, -chain length)
            sel2 = {}
            for seq in refs:
                cands = [recs[seq][k] for k in sel[seq]]
                ca = data[seq]['CarbonA_unfiltered']
                kept = greedy_nonoverlap(cands, lambda c: (-votes[seq][c['key']], 0 if c['key'] in ca else 1,
                                                           -(c['exons'][-1][1] - c['exons'][0][0]), c['key']))
                sel2[seq] = {c['key'] for c in kept}
            res_no[f'>={kmin}'] = score(sel2)
        R['vote'][name] = res
        R['vote_nonoverlap'][name] = res_no
    # ---- hybrid: Carbon-A filtered union chains supported by >=2 of {ANNEVO,Tiberius,E} not in Carbon-A, non-overlapping w/ CA chains
    for name, base in (('CarbonA_filtered', 'CarbonA'), ('CarbonA_unfiltered', 'CarbonA_unfiltered')):
        for need in (2, 3):
            sel = {}
            for seq in refs:
                ca = data[seq][base]
                extra = Counter()
                recs = {}
                for t in ('ANNEVO', 'Tiberius', 'E_lora'):
                    for k, v in data[seq][t].items():
                        if k not in ca:
                            extra[k] += 1; recs[k] = v
                cand = [recs[k] for k, c in extra.items() if c >= need]
                spans = defaultdict(list)
                for v in ca.values():
                    spans[v['strand']].append((v['exons'][0][0], v['exons'][-1][1]))
                add = []
                for c in sorted(cand, key=lambda c: (-extra[c['key']], c['key'])):
                    sp_ = (c['exons'][0][0], c['exons'][-1][1])
                    if any(overlap(sp_, s) for s in spans[c['strand']]):
                        continue
                    spans[c['strand']].append(sp_); add.append(c['key'])
                sel[seq] = set(ca) | set(add)
            R['hybrid'][f'{name}+nonoverlap_chains_in>={need}_of(ANNEVO,Tiberius,E)'] = score(sel)
    # ---- precision of "all agree" hits: chains predicted by CA(unfilt) and E and ANNEVO etc.
    with open(OUT / 'summary.json', 'w') as f:
        json.dump(R, f, indent=1)
    def line(n, v): print(f"{n:55s} tp={v['tp']:6d} pred={v['pred']:6d} P={100*v['precision']:.2f} R={100*v['recall']:.2f} F1={v['F1']:.3f}")
    for t, v in R['single_tool'].items(): line(t, v)
    print('--- union upper bound')
    for n, v in R['union_upper_bound'].items(): print(f"{n:50s} hits={v['union_hits']} recall_ceiling={100*v['recall_ceiling']:.2f}% F1(0FP)={v['F1_if_union_hits_and_zero_FP']:.2f}")
    print('never hit by any of 8:', miss_all, 'hist', R['ref_hit_count_hist_core4_unfilteredCA'])
    for n, d in R['vote'].items():
        print('--- vote', n)
        for k, v in d.items(): line(k, v)
    for n, d in R['vote_nonoverlap'].items():
        print('--- vote nonoverlap', n)
        for k, v in d.items(): line(k, v)
    print('--- hybrid')
    for n, v in R['hybrid'].items(): line(n, v)


if __name__ == '__main__':
    main()
