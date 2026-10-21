"""Full comparison tables: exact-chain F1 (pooled / species / exon-count strata), CDS-base F1 and background burden for A1 and comparators.
Base level: union of all annotated protein-coding CDS bases (primary refs + other assessable isoforms, strand-aware) vs union of unique predicted chain CDS bases.
Background burden: predicted CDS bases (either strand) outside every annotated gene span (M26 'all annotated gene span complement'), in kb and % of the background length.
Usage: a1_compare.py <a1_run> <gate_run> <out.json>
"""
import json, sys
import numpy as np
import torch
import a1_common as K
sys.path.insert(0, str(K.CA_SRC))
import m29_common as M                      # noqa: E402
from m29_load_preds import tool_chains      # noqa: E402
from m29_r3_posterior_filter import carbon_conf   # noqa: E402
import a1_eval as E                         # noqa: E402

TRAIN, DEV = E.TRAIN, E.DEV
refs = E.Refs()
DEVGFF = {'CarbonA': {'NC_003074.8': M.M29 / 'G0/ara_pred/prediction.gff', 'NC_089041.1': M.M29 / 'G0/rice_pred/prediction.gff'},
          'B_freshjoint': {s: M.OUTS / f'M28-FULL-GRAMMAR-JOINT-CONTROLS-DEV-RAW-R1/FRESH_JOINT24576/{s}/prediction/full_grammar.gff3' for s in DEV}}


def cdict(lst):
    return {k: dict(nex=len(v['exons']), exons=v['exons'], strand=v['strand']) for k, v in lst.items()}


def carbon_dev(seq, filtered):
    gff = DEVGFF['CarbonA'][seq]
    conf = carbon_conf(gff)
    out = {}
    for tx in M.load_txs(gff, seq):
        c = conf[tx['gene_id']]
        if filtered and not (c is None or c >= 0.1):
            continue
        ex = [(a, b) for a, b, _ in tx['CDS']]
        out.setdefault(M.key(tx), dict(nex=len(ex), exons=ex, strand=tx['strand']))
    return out


def gff_chains(path, seq):
    out = {}
    for tx in M.load_txs(path, seq):
        ex = [(a, b) for a, b, _ in tx['CDS']]
        out.setdefault(M.key(tx), dict(nex=len(ex), exons=ex, strand=tx['strand']))
    return out


_MASK = {}


def truth_masks(seq):
    if seq in _MASK:
        return _MASK[seq]
    L = M.LENGTH[seq]
    ref = M.build_reference(seq)
    cds = {'+': np.zeros(L, bool), '-': np.zeros(L, bool)}
    for t in list(ref['refs']) + list(ref['isoforms']):
        for a, b, _ in t['CDS']:
            cds[t['strand']][a:b] = True
    gene = np.zeros(L, bool)
    for a, b in ref['gene_spans']:
        gene[max(0, a - 1):b] = True
    _MASK[seq] = (cds, gene)
    return _MASK[seq]


def base_metrics(chains_by_seq, seqs):
    tp = pp = rp = bg = bgsize = 0
    for s in seqs:
        cds, gene = truth_masks(s)
        L = M.LENGTH[s]
        pred = {'+': np.zeros(L, bool), '-': np.zeros(L, bool)}
        for v in chains_by_seq[s].values():
            for a, b in v['exons']:
                pred[v['strand']][a:b] = True
        for st in '+-':
            tp += int((pred[st] & cds[st]).sum()); pp += int(pred[st].sum()); rp += int(cds[st].sum())
        anyp = pred['+'] | pred['-']
        bg += int((anyp & ~gene).sum()); bgsize += int((~gene).sum())
    return dict(base_P=100 * tp / pp, base_R=100 * tp / rp, base_F1=100 * 2 * tp / (pp + rp), bg_pred_kb=bg / 1000, bg_burden_pct=100 * bg / bgsize)


def main():
    a1, gate, outp = sys.argv[1], sys.argv[2], sys.argv[3]
    ev = json.load(open(K.Path(a1) / 'eval_a1.json'))
    thr_loo = ev['TRAIN']['posterior_LOO_chromosome']['thresholds']
    thr_all = ev['TRAIN']['posterior_threshold_all_TRAIN']
    from constants import CONFIDENCE_NETWORK_PATH
    from confidence_network import load_confidence_network
    cn = load_confidence_network(CONFIDENCE_NETWORK_PATH, torch.device('cpu'))
    tc = {}
    methods = {}
    for scope, seqs in (('TRAIN', TRAIN), ('DEV', DEV)):
        m = {}
        a1c = {s: E.load_chains(a1, s) for s in seqs}
        gc = {s: E.load_chains(gate, s) for s in seqs}
        m['A1_unfiltered'] = a1c
        m['Gate_unfiltered'] = gc
        ca1 = {s: E.conf_scores(s, a1c[s], cn, tc) for s in seqs}
        cg = {s: E.conf_scores(s, gc[s], cn, tc) for s in seqs}
        m['A1_carbonA_conf0.1'] = {s: {k: v for k, v in a1c[s].items() if ca1[s][k] >= 0.1} for s in seqs}
        m['Gate_carbonA_conf0.1'] = {s: {k: v for k, v in gc[s].items() if cg[s][k] >= 0.1} for s in seqs}
        post = {s: E.posterior_scores(a1, s, a1c[s]) for s in seqs}
        m['A1_own_posterior'] = {s: {k: v for k, v in a1c[s].items() if post[s][k] >= (thr_loo[s] if scope == 'TRAIN' else thr_all)} for s in seqs}
        if scope == 'TRAIN':
            ca = {s: tool_chains('CarbonA_unfiltered', s) for s in seqs}
            m['CarbonA_unfiltered'] = {s: cdict(ca[s]) for s in seqs}
            m['CarbonA_filtered'] = {s: cdict({k: v for k, v in ca[s].items() if v['na'] or v['conf'] >= .1}) for s in seqs}
            for t in ('B_freshjoint', 'E_lora', 'ANNEVO', 'Tiberius'):
                m[t] = {s: cdict(tool_chains(t, s)) for s in seqs}
        else:
            m['CarbonA_unfiltered'] = {s: carbon_dev(s, False) for s in seqs}
            m['CarbonA_filtered'] = {s: carbon_dev(s, True) for s in seqs}
            m['B_freshjoint'] = {s: gff_chains(DEVGFF['B_freshjoint'][s], s) for s in seqs}
        methods[scope] = {}
        for name, ch in m.items():
            r = E.summarize(ch, refs, seqs)
            r['base'] = base_metrics(ch, seqs)
            methods[scope][name] = r
            print(scope, name, round(r['pooled']['F1'], 3), {k: round(v, 2) for k, v in r['base'].items()}, flush=True)
    json.dump(methods, open(outp, 'w'), indent=1)


if __name__ == '__main__':
    main()
