"""R4 error-type decomposition (same rules/code as m29_r4_error_types) for A1 and comparators on the six TRAIN contigs (CPU).
Usage: a1_r4.py <a1_run_dir> <gate_run_dir> <out.json>    A1 variants: raw(cds96), +Carbon-A conf>=0.1, +own posterior (LOO thresholds from <run>/eval_a1.json)
"""
import json, sys
from collections import Counter, defaultdict
import numpy as np
import torch
import a1_common as K
sys.path.insert(0, str(K.CA_SRC))
import m29_common as M                       # noqa: E402
from m29_load_preds import tool_chains       # noqa: E402
import m29_r4_error_types as R4              # noqa: E402
import a1_eval as E                          # noqa: E402


def run_chains(run, seq, score=None, thr=None, kind='raw'):
    ch = E.load_chains(run, seq)
    if kind == 'conf':
        from constants import CONFIDENCE_NETWORK_PATH
        from confidence_network import load_confidence_network
        cn = E.__dict__.setdefault('_CN', load_confidence_network(CONFIDENCE_NETWORK_PATH, torch.device('cpu')))
        sc = E.conf_scores(seq, ch, cn, E.__dict__.setdefault('_TC', {}))
        ch = {k: v for k, v in ch.items() if sc[k] >= 0.1}
    elif kind == 'post':
        sc = E.posterior_scores(run, seq, ch)
        ch = {k: v for k, v in ch.items() if sc[k] >= thr}
    return {k: dict(key=k, strand=v['strand'], exons=v['exons'], id=str(i), conf=1., na=False) for i, (k, v) in enumerate(ch.items())}


def main():
    a1, gate, outp = sys.argv[1], sys.argv[2], sys.argv[3]
    ev = json.load(open(K.Path(a1) / 'eval_a1.json'))
    thr_loo = ev['TRAIN'].get('posterior_LOO_chromosome', {}).get('thresholds')
    tools = ['A1_raw', 'A1_conf', 'A1_post', 'Gate', 'CarbonA_unfiltered', 'CarbonA_filtered', 'B_freshjoint', 'E_lora']
    if thr_loo is None:
        tools.remove('A1_post')
    res = {t: dict(rows=[]) for t in tools}
    hits = {t: set() for t in tools}
    allrefs = []
    for sp, seq, L in M.SCOPE:
        ref = M.build_reference(seq)
        refs = ref['refs']
        iso_by_gene = defaultdict(set)
        for t in ref['isoforms']:
            iso_by_gene[t['gene_id']].add(M.key(t))
        ref_index = R4.Index([dict(strand=t['strand'], exons=[(a, b) for a, b, _ in t['CDS']], id=t['id']) for t in refs])
        preds = {}
        preds['A1_raw'] = run_chains(a1, seq)
        preds['A1_conf'] = run_chains(a1, seq, kind='conf')
        if 'A1_post' in tools:
            preds['A1_post'] = run_chains(a1, seq, thr=thr_loo[seq], kind='post')
        preds['Gate'] = run_chains(gate, seq)
        ca = tool_chains('CarbonA_unfiltered', seq)
        preds['CarbonA_unfiltered'] = ca
        preds['CarbonA_filtered'] = {k: v for k, v in ca.items() if v['na'] or v['conf'] >= 0.1}
        preds['B_freshjoint'] = tool_chains('B_freshjoint', seq)
        preds['E_lora'] = tool_chains('E_lora', seq)
        idxs = {t: R4.Index(list(preds[t].values())) for t in tools}
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
                cands = [(R4.ov_bp(c['exons'], Er), c) for c in idxs[t].query(strand, lo, hi)]
                cands = [(o, c) for o, c in cands if o > 0]
                if t == 'CarbonA_filtered' and k in ca:
                    cat, sub = 'confidence_filtered', ''
                elif t in ('A1_conf', 'A1_post') and k in preds['A1_raw']:
                    cat, sub = 'confidence_filtered', ''
                elif any(c['key'] in iso_by_gene[r['gene_id']] and c['key'] != k for o, c in cands):
                    cat, sub = 'other_isoform', ''
                elif not cands:
                    other = '-' if strand == '+' else '+'
                    oc = [c for c in idxs[t].query(other, lo, hi) if R4.ov_bp(c['exons'], Er) > 0]
                    cat, sub = ('wrong_strand_only', '') if oc else ('absent_no_call', '')
                else:
                    cands.sort(key=lambda x: -x[0])
                    best = cands[0][1]
                    others = [x for x in ref_index.query(strand, best['exons'][0][0], best['exons'][-1][1])
                              if x['id'] != r['id'] and R4.ov_bp(x['exons'], best['exons']) >= R4.MIN_BP]
                    big = [c for o, c in cands if o >= R4.MIN_BP]
                    if others:
                        cat, sub = 'merge', ''
                    elif len(big) >= 2:
                        cat, sub = 'split', ''
                    else:
                        cat, sub = R4.classify_local(Er, best['exons'], strand)
                res[t]['rows'].append((sp, seq, r['id'], nR, cat, sub))
        print(seq, 'done', flush=True)
    out = {'n_ref': len(allrefs), 'tools': {}}
    cats = R4.CATS
    for t in tools:
        rows = res[t]['rows']
        multi = [r for r in rows if r[4] != 'single_exon_miss']
        tot = Counter(R4.nseg_bin(x[3]) if x[3] >= 2 else '1' for x in allrefs)
        hh = Counter()
        for sp, seq, rid, n, k in allrefs:
            if (seq, rid) in hits[t]:
                hh[R4.nseg_bin(n) if n >= 2 else '1'] += 1
        out['tools'][t] = dict(hits_total=len(hits[t]), multi_miss_total=len(multi), categories=dict(Counter(r[4] for r in multi)),
                               by_species={s: dict(Counter(r[4] for r in multi if r[0] == s)) for s in ('arabidopsis_thaliana', 'oryza_sativa')},
                               recall_by_nseg={b: dict(ref=tot[b], hit=hh[b], recall=hh[b] / tot[b]) for b in ['1', '2', '3-5', '6-9', '10+']})
    json.dump(out, open(outp, 'w'), indent=1)
    for t in tools:
        d = out['tools'][t]
        print(t, 'hits', d['hits_total'], 'multi_miss', d['multi_miss_total'], d['categories'])


if __name__ == '__main__':
    main()
