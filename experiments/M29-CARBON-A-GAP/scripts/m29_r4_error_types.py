"""R4: error-type decomposition of unmatched multi-segment primary references (CPU).

For every primary reference (frozen ruler, 22,688) that the tool's output does not reproduce exactly, assign ONE category,
in this priority order (rules are structural only, no score is used):
  confidence_filtered : Carbon-A only - the unfiltered output contains the exact chain, the 0.1-or-NA filter removed it
  other_isoform       : an overlapping same-strand prediction equals another assessable isoform of the same gene (ruler = miss)
  absent_no_call      : no prediction overlaps the reference CDS (either strand)
  wrong_strand_only   : only opposite-strand predictions overlap the CDS
  merge               : best-overlap prediction also overlaps >=60 CDS bp of a different primary reference (same strand)
  split               : >=2 same-strand predictions each overlap >=60 CDS bp of this reference
  endpoint_only       : same exon count, identical internal junctions; only first-exon start / last-exon end differ
  exon_missed         : prediction lacks >=1 reference exon (rest identical; outer terminal boundary may differ)
  exon_extra          : prediction has >=1 extra exon (reference exons all present)
  single_junction_shift : same exon count, exactly one junction replaced by one junction (donor and/or acceptor moved)
  complex_chain       : any other overlapping-but-different chain
"""
import json
from bisect import bisect_left
from collections import Counter, defaultdict
import m29_common as M
from m29_load_preds import tool_chains

OUT = M.M29 / 'R4'
OUT.mkdir(parents=True, exist_ok=True)
CATS = ['confidence_filtered', 'other_isoform', 'absent_no_call', 'wrong_strand_only', 'merge', 'split', 'endpoint_only',
        'exon_missed', 'exon_extra', 'single_junction_shift', 'complex_chain']
MIN_BP = 60


def ov_bp(e1, e2):
    i = j = tot = 0
    while i < len(e1) and j < len(e2):
        a = max(e1[i][0], e2[j][0]); b = min(e1[i][1], e2[j][1])
        if b > a:
            tot += b - a
        if e1[i][1] < e2[j][1]:
            i += 1
        else:
            j += 1
    return tot


class Index:
    def __init__(self, chains):
        self.by = defaultdict(list)
        for c in chains:
            self.by[c['strand']].append(c)
        self.starts, self.pmax = {}, {}
        for s, lst in self.by.items():
            lst.sort(key=lambda c: c['exons'][0][0])
            self.starts[s] = [c['exons'][0][0] for c in lst]
            m, pm = -1, []
            for c in lst:
                m = max(m, c['exons'][-1][1]); pm.append(m)
            self.pmax[s] = pm

    def query(self, strand, lo, hi):
        lst = self.by.get(strand, [])
        if not lst:
            return []
        i = bisect_left(self.starts[strand], hi) - 1
        out = []
        pm = self.pmax[strand]
        while i >= 0 and pm[i] > lo:
            c = lst[i]
            if c['exons'][-1][1] > lo and c['exons'][0][0] < hi:
                out.append(c)
            i -= 1
        return out


def classify_local(Er, Eq, strand):
    nR, nQ = len(Er), len(Eq)
    JR = {(Er[i][1], Er[i + 1][0]) for i in range(nR - 1)}
    JQ = {(Eq[i][1], Eq[i + 1][0]) for i in range(nQ - 1)}
    A, B = JR - JQ, JQ - JR
    if nR == nQ and not A and not B:
        s_diff, e_diff = Er[0][0] != Eq[0][0], Er[-1][1] != Eq[-1][1]
        if strand == '+':
            sub = 'both' if s_diff and e_diff else 'start_codon_side' if s_diff else 'stop_side'
        else:
            sub = 'both' if s_diff and e_diff else 'stop_side' if s_diff else 'start_codon_side'
        return 'endpoint_only', sub
    SR, ER_ = set(Er), {e[1] for e in Er}
    ST = {e[0] for e in Er}

    def represented(Ea, Eb):
        """exons of Ea represented in Eb (exact, or Eb terminal with matching inner boundary)."""
        Sb = set(Eb)
        rep = set()
        for i, e in enumerate(Ea):
            if e in Sb:
                rep.add(i)
        # relaxed terminals
        if Eb[0] not in Sb or True:
            for i, e in enumerate(Ea):
                if Eb[0][1] == e[1] and (i < len(Ea) - 1 or len(Ea) == 1):
                    rep.add(i) if len(Eb) > 1 or Eb[0] in Sb else None
                if Eb[-1][0] == e[0] and (i > 0 or len(Ea) == 1):
                    rep.add(i) if len(Eb) > 1 or Eb[-1] in Sb else None
        return rep

    def contained(Esmall, Ebig):
        """every exon of Esmall present in Ebig exactly, except terminals whose inner boundary matches."""
        Sb = set(Ebig); ends = {e[1] for e in Ebig}; starts = {e[0] for e in Ebig}
        n = len(Esmall)
        for i, e in enumerate(Esmall):
            if e in Sb:
                continue
            if n == 1:
                return False
            if i == 0 and e[1] in ends:
                continue
            if i == n - 1 and e[0] in starts:
                continue
            return False
        return True

    if nQ < nR and contained(Eq, Er):
        # which reference exons are not represented
        rep = set()
        Sq = set(Eq)
        for i, e in enumerate(Er):
            if e in Sq or (Eq[0][1] == e[1] and nQ > 1 and i < nR - 1) or (Eq[-1][0] == e[0] and nQ > 1 and i > 0):
                rep.add(i)
        missing = [i for i in range(nR) if i not in rep]
        if missing:
            sub = 'terminal_exon' if (0 in missing or nR - 1 in missing) else 'internal_exon'
            if len(missing) > 1:
                sub += '_x%d' % (len(missing) if len(missing) < 4 else 4)
            return 'exon_missed', sub
    if nQ > nR and contained(Er, Eq):
        Sr = set(Er)
        rep = set()
        for i, e in enumerate(Eq):
            if e in Sr or (Er[0][1] == e[1] and nR > 1 and i < nQ - 1) or (Er[-1][0] == e[0] and nR > 1 and i > 0):
                rep.add(i)
        extra = [i for i in range(nQ) if i not in rep]
        if extra:
            sub = 'terminal_exon' if (0 in extra or nQ - 1 in extra) else 'internal_exon'
            if len(extra) > 1:
                sub += '_x%d' % (len(extra) if len(extra) < 4 else 4)
            return 'exon_extra', sub
    if nR == nQ and len(A) == 1 and len(B) == 1:
        (d1, a1), (d2, a2) = list(A)[0], list(B)[0]
        sub = 'both_sides' if d1 != d2 and a1 != a2 else 'donor_only' if d1 != d2 else 'acceptor_only'
        if Er[0][0] != Eq[0][0] or Er[-1][1] != Eq[-1][1]:
            sub += '+endpoint'
        return 'single_junction_shift', sub
    return 'complex_chain', f'nQ-nR={min(max(nQ - nR, -3), 3):+d},junction_diff={min(len(A) + len(B), 6)}'


def nseg_bin(n):
    return '2' if n == 2 else '3-5' if n <= 5 else '6-9' if n <= 9 else '10+'


def main():
    tools = {'E_lora': 'E_lora', 'CarbonA_filtered': 'CarbonA_unfiltered', 'CarbonA_unfiltered': 'CarbonA_unfiltered',
             'EM_pool': 'EM_pool', 'ANNEVO': 'ANNEVO', 'Tiberius': 'Tiberius'}
    res = {t: dict(rows=[]) for t in tools}
    hits = {t: set() for t in tools}
    allrefs = []
    for sp, seq, L in M.SCOPE:
        ref = M.build_reference(seq)
        refs = ref['refs']
        iso_by_gene = defaultdict(set)
        for t in ref['isoforms']:
            iso_by_gene[t['gene_id']].add(M.key(t))
        pk_ref = {M.key(t) for t in refs}
        ref_index = Index([dict(strand=t['strand'], exons=[(a, b) for a, b, _ in t['CDS']], id=t['id']) for t in refs])
        preds = {}
        for t, base in tools.items():
            ch = tool_chains(base, seq)
            if t == 'CarbonA_filtered':
                ch = {k: v for k, v in ch.items() if v['na'] or v['conf'] >= 0.1}
            preds[t] = ch
        unf = tool_chains('CarbonA_unfiltered', seq)
        idxs = {t: Index(list(preds[t].values())) for t in tools}
        opp = {t: idxs[t] for t in tools}
        for r in refs:
            k = M.key(r)
            Er = [(a, b) for a, b, _ in r['CDS']]
            nR = len(Er)
            rid = (seq, r['id'])
            allrefs.append((sp, seq, r['id'], nR, k))
            for t in tools:
                if k in preds[t]:
                    hits[t].add(rid)
                    continue
                if nR < 2:
                    res[t]['rows'].append((sp, seq, r['id'], nR, 'single_exon_miss', ''))
                    continue
                strand = r['strand']
                lo, hi = Er[0][0], Er[-1][1]
                cands = idxs[t].query(strand, lo, hi)
                cands = [(ov_bp(c['exons'], Er), c) for c in cands]
                cands = [(o, c) for o, c in cands if o > 0]
                if t == 'CarbonA_filtered' and k in unf:
                    cat, sub = 'confidence_filtered', ''
                elif any(c['key'] in iso_by_gene[r['gene_id']] and c['key'] != k for o, c in cands):
                    cat, sub = 'other_isoform', ''
                elif not cands:
                    other = '-' if strand == '+' else '+'
                    oc = [c for c in idxs[t].query(other, lo, hi) if ov_bp(c['exons'], Er) > 0]
                    cat, sub = ('wrong_strand_only', '') if oc else ('absent_no_call', '')
                else:
                    cands.sort(key=lambda x: -x[0])
                    best = cands[0][1]
                    others = [x for x in ref_index.query(strand, best['exons'][0][0], best['exons'][-1][1])
                              if x['id'] != r['id'] and ov_bp(x['exons'], best['exons']) >= MIN_BP]
                    big = [c for o, c in cands if o >= MIN_BP]
                    if others:
                        cat, sub = 'merge', f'neighbors={min(len(others), 3)}'
                    elif len(big) >= 2:
                        cat, sub = 'split', f'pieces={min(len(big), 4)}'
                    else:
                        cat, sub = classify_local(Er, best['exons'], strand)
                res[t]['rows'].append((sp, seq, r['id'], nR, cat, sub))
        print(seq, 'done', flush=True)
    # ---- aggregate
    nsegs = {(a, b): n for (sp, a, b, n, k) in allrefs}
    out = {'n_ref': len(allrefs), 'n_multi': sum(1 for x in allrefs if x[3] >= 2), 'n_single': sum(1 for x in allrefs if x[3] == 1)}
    ms = Counter(nseg_bin(x[3]) for x in allrefs if x[3] >= 2)
    out['multi_by_nseg'] = dict(ms)
    out['tools'] = {}
    for t in tools:
        rows = res[t]['rows']
        d = {'hits_total': len(hits[t]), 'misses_total': len(rows), 'single_exon_miss': sum(1 for r in rows if r[4] == 'single_exon_miss')}
        d['hit_recall_by_nseg'] = {}
        tot = Counter(nseg_bin(x[3]) if x[3] >= 2 else '1' for x in allrefs)
        hh = Counter()
        for sp, seq, rid, n, k in allrefs:
            if (seq, rid) in hits[t]:
                hh[nseg_bin(n) if n >= 2 else '1'] += 1
        d['hit_recall_by_nseg'] = {b: dict(ref=tot[b], hit=hh[b], recall=hh[b] / tot[b]) for b in ['1', '2', '3-5', '6-9', '10+']}
        multi = [r for r in rows if r[4] != 'single_exon_miss']
        d['multi_miss_total'] = len(multi)
        d['categories'] = dict(Counter(r[4] for r in multi))
        d['categories_sub'] = {c: dict(Counter(r[5] for r in multi if r[4] == c)) for c in CATS if any(r[4] == c for r in multi)}
        d['by_species'] = {s: dict(Counter(r[4] for r in multi if r[0] == s)) for s in ('arabidopsis_thaliana', 'oryza_sativa')}
        d['by_nseg'] = {b: dict(Counter(r[4] for r in multi if nseg_bin(r[3]) == b)) for b in ['2', '3-5', '6-9', '10+']}
        d['by_species_nseg'] = {f'{s}|{b}': dict(Counter(r[4] for r in multi if r[0] == s and nseg_bin(r[3]) == b))
                                for s in ('arabidopsis_thaliana', 'oryza_sativa') for b in ['2', '3-5', '6-9', '10+']}
        out['tools'][t] = d
        with open(OUT / f'misses_{t}.tsv', 'w') as f:
            f.write('species\tseqid\tref_id\tnseg\tcategory\tsub\n')
            for r in rows:
                f.write('\t'.join(map(str, r)) + '\n')
    # ---- paired: where E misses but Carbon-A (filtered) hits and vice versa
    def tab(a, b):
        sa = {(r[1], r[2]): r for r in res[a]['rows'] if r[4] != 'single_exon_miss'}
        return sa
    pair = {}
    for a, b in (('E_lora', 'CarbonA_filtered'), ('CarbonA_filtered', 'E_lora'), ('E_lora', 'CarbonA_unfiltered'),
                 ('EM_pool', 'CarbonA_filtered'), ('E_lora', 'ANNEVO')):
        # categories of a's miss among refs b hits (multi-seg only)
        rows = [r for r in res[a]['rows'] if r[4] != 'single_exon_miss' and (r[1], r[2]) in hits[b]]
        both = [r for r in res[a]['rows'] if r[4] != 'single_exon_miss' and (r[1], r[2]) not in hits[b]]
        pair[f'{a}_miss_given_{b}_hit'] = dict(n=len(rows), categories=dict(Counter(r[4] for r in rows)),
                                               by_species=dict(Counter(r[0] for r in rows)),
                                               by_nseg=dict(Counter(nseg_bin(r[3]) for r in rows)))
        pair[f'{a}_miss_given_{b}_miss'] = dict(n=len(both), categories=dict(Counter(r[4] for r in both)))
    out['paired'] = pair
    # single-exon: hits for both
    sx = {t: sum(1 for x in allrefs if x[3] == 1 and (x[1], x[2]) in hits[t]) for t in tools}
    out['single_exon_hits'] = sx
    with open(OUT / 'summary.json', 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: out[k] for k in ('n_ref', 'n_multi', 'n_single', 'multi_by_nseg', 'single_exon_hits')}))
    for t in tools:
        d = out['tools'][t]
        print(t, 'hits', d['hits_total'], 'multi_miss', d['multi_miss_total'], d['categories'])
    print(json.dumps(pair, indent=1))


if __name__ == '__main__':
    main()
