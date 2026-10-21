"""Scoring of A1 (and gate) predictions with the frozen ruler (M29 `common_ruler` wrapper), exactly as in the gate / G0 scoring.
For each run: unique exact CDS chains (CDS>=96 bp rule on, as in Carbon-A and the gate score 75.76), pooled/per-species/per-contig F1,
stratified by CDS-exon count (1 / 2-5 / >=6), plus the two filters:
  (i)  own gene posterior log P(g): threshold chosen on TRAIN (leave-one-chromosome-out for TRAIN; ALL six TRAIN chromosomes -> DEV, never tuned on DEV)
  (ii) Carbon-A confidence network on our chains, fixed threshold 0.1 (no tuning)
Usage: a1_eval.py <label> <run_dir> [--no-posterior]   -> <run_dir>/eval_<label>.json
"""
import json, sys
from collections import defaultdict
import numpy as np
import torch
import a1_common as K
sys.path.insert(0, str(K.CA_SRC))
import m29_common as M            # noqa: E402
from m29_r3_posterior_filter import curve, best_threshold_pick, auroc_ap   # noqa: E402

TRAIN = [s for _, s, _ in M.SCOPE]
DEV = [s for _, s, _ in M.DEV_SCOPE]
SPECIES = M.SEQ2SPECIES
STRATA = (('1', lambda n: n == 1), ('2-5', lambda n: 2 <= n <= 5), ('>=6', lambda n: n >= 6))


def strat(n):
    return '1' if n == 1 else '2-5' if n <= 5 else '>=6'


def load_chains(run, seq, cds96=True):
    """unique exact chains from the run GFF -> {key: dict(nex, gain, exons, strand)}"""
    gff = K.K_path if False else (K.Path(run) / seq / 'prediction' / 'full_grammar.gff3')
    out = {}
    for tx in M.load_txs(gff, seq):
        ex = [(a, b) for a, b, _ in tx['CDS']]
        if cds96 and sum(b - a for a, b in ex) < 96:
            continue
        k = M.key(tx)
        out.setdefault(k, dict(nex=len(ex), gain=float(tx['score']) if tx.get('score') not in (None, '.') else 0., exons=ex, strand=tx['strand']))
    return out


def counts(pred, refkeys, nex_ref):
    tp = {s: 0 for s, _ in STRATA}; pr = {s: 0 for s, _ in STRATA}; rf = {s: 0 for s, _ in STRATA}
    for k, v in pred.items():
        pr[strat(v['nex'])] += 1
        if k in refkeys:
            tp[strat(v['nex'])] += 1
    for k in refkeys:
        rf[strat(nex_ref[k])] += 1
    return tp, pr, rf


def f1(tp, p, r):
    return 100 * 2 * tp / (p + r) if p + r else 0.


class Refs:
    def __init__(self):
        self.keys, self.nex = {}, {}
        for _, seq, _ in M.SCOPE + M.DEV_SCOPE:
            ref = M.build_reference(seq)
            self.keys[seq] = {M.key(t) for t in ref['refs']}
            for t in ref['refs']:
                self.nex[M.key(t)] = len(t['CDS'])


def summarize(chains_by_seq, refs, seqs):
    """chains_by_seq: seq -> {key: info};  returns pooled / per-species / per-contig / strata dicts."""
    res = {'per_contig': {}, 'pooled': None, 'by_species': {}, 'by_strata': {}}
    agg = defaultdict(lambda: [0, 0, 0])
    sp_agg = defaultdict(lambda: defaultdict(lambda: [0, 0, 0]))
    st_agg = defaultdict(lambda: [0, 0, 0])
    for s in seqs:
        pred = chains_by_seq[s]
        tp = sum(1 for k in pred if k in refs.keys[s])
        res['per_contig'][s] = dict(tp=tp, pred=len(pred), ref=len(refs.keys[s]), F1=f1(tp, len(pred), len(refs.keys[s])))
        for a, v in zip(range(3), (tp, len(pred), len(refs.keys[s]))):
            agg['all'][a] += v; sp_agg[SPECIES[s]]['x'][a] += v
        ctp, cpr, crf = counts(pred, refs.keys[s], refs.nex)
        for st, _ in STRATA:
            st_agg[st][0] += ctp[st]; st_agg[st][1] += cpr[st]; st_agg[st][2] += crf[st]
    t = agg['all']
    res['pooled'] = dict(tp=t[0], pred=t[1], ref=t[2], F1=f1(*t))
    res['by_species'] = {sp: dict(tp=v['x'][0], pred=v['x'][1], ref=v['x'][2], F1=f1(*v['x'])) for sp, v in sp_agg.items()}
    res['by_strata'] = {st: dict(tp=v[0], pred=v[1], ref=v[2], F1=f1(*v), recall=100 * v[0] / v[2] if v[2] else None,
                                 precision=100 * v[0] / v[1] if v[1] else None) for st, v in st_agg.items()}
    return res


def posterior_scores(run, seq, chains):
    """attach log_p to every kept chain from <run>/posterior/<seq>_{+,-}.json (matched by oriented exons)"""
    lp = {}
    for strand in '+-':
        f = K.Path(run) / 'posterior' / f'{seq}_{strand}.json'
        for g in json.load(open(f))['genes']:
            ex = tuple((a, b) for a, b in g['exons'])
            lp[(seq, strand, ex)] = g['log_p']
    return {k: lp[k] for k in chains}


def conf_scores(seq, chains, cn, tracks_cache):
    if seq not in tracks_cache:
        import a1_build_offsets as BO
        import m29_gate_parity as GP
        rid, sequence, pos, neg = BO.load_dev(seq) if seq in BO.DEV_RAW else GP.load_tracks(seq)
        tracks_cache[seq] = (len(sequence), {'+': np.ascontiguousarray(pos), '-': np.ascontiguousarray(neg[::-1])})
    n, tr = tracks_cache[seq]
    out = {}
    for strand in '+-':
        ks = [k for k in chains if chains[k]['strand'] == strand]
        if not ks:
            continue
        ex = [[(a, b - 1) for a, b in chains[k]['exons']] if strand == '+' else [(n - b, n - 1 - a) for a, b in reversed(chains[k]['exons'])] for k in ks]
        conf = cn.inference(tr[strand], ex, list(range(len(ks))))
        for k, v in zip(ks, conf):
            out[k] = float(v)
    return out


def filtered(chains, score, thr):
    return {k: v for k, v in chains.items() if score[k] >= thr}


def main():
    label, run = sys.argv[1], sys.argv[2]
    with_post = '--no-posterior' not in sys.argv
    torch.set_num_threads(4)
    refs = Refs()
    run_p = K.Path(run)
    allseq = TRAIN + DEV
    chains = {s: load_chains(run, s) for s in allseq if (run_p / s / 'prediction' / 'full_grammar.gff3').exists()}
    out = dict(label=label, run=str(run), cds96_rule=True)
    raw = {s: load_chains(run, s, cds96=False) for s in chains}
    for scope, seqs in (('TRAIN', TRAIN), ('DEV', DEV)):
        if not all(s in chains for s in seqs):
            continue
        out[scope] = dict(unfiltered_no_cds96=summarize(raw, refs, seqs)['pooled'], unfiltered=summarize(chains, refs, seqs))
    # (ii) Carbon-A confidence network
    from constants import CONFIDENCE_NETWORK_PATH
    from confidence_network import load_confidence_network
    cn = load_confidence_network(CONFIDENCE_NETWORK_PATH, torch.device('cpu'))
    tracks_cache = {}
    conf = {s: conf_scores(s, chains[s], cn, tracks_cache) for s in chains}
    for scope, seqs in (('TRAIN', TRAIN), ('DEV', DEV)):
        if scope in out:
            out[scope]['carbonA_conf>=0.1'] = summarize({s: filtered(chains[s], conf[s], 0.1) for s in seqs}, refs, seqs)
            lab = [(conf[s][k], int(k in refs.keys[s])) for s in seqs for k in chains[s]]
            out[scope]['carbonA_conf_ranking'] = dict(zip(('AUROC', 'AP'), auroc_ap([a for a, _ in lab], [b for _, b in lab])))
    # (i) own posterior
    if with_post and all((run_p / 'posterior' / f'{s}_{st}.json').exists() for s in chains for st in '+-'):
        post = {s: posterior_scores(run, s, chains[s]) for s in chains}
        rows = {s: [(post[s][k], int(k in refs.keys[s])) for k in chains[s]] for s in chains}
        nref = {s: len(refs.keys[s]) for s in chains}
        if all(s in chains for s in TRAIN):
            # leave-one-chromosome-out on TRAIN
            kept = {}
            thr_used = {}
            for s in TRAIN:
                tr_rows = [r for t in TRAIN if t != s for r in rows[t]]
                thr, _ = best_threshold_pick(tr_rows, sum(nref[t] for t in TRAIN if t != s))
                thr_used[s] = thr
                kept[s] = filtered(chains[s], post[s], thr)
            out['TRAIN']['posterior_LOO_chromosome'] = dict(summarize(kept, refs, TRAIN), thresholds=thr_used)
            thr_all, _ = best_threshold_pick([r for t in TRAIN for r in rows[t]], sum(nref[t] for t in TRAIN))
            out['TRAIN']['posterior_ranking'] = dict(zip(('AUROC', 'AP'), auroc_ap([a for s in TRAIN for a, _ in rows[s]], [b for s in TRAIN for _, b in rows[s]])))
            out['TRAIN']['posterior_threshold_all_TRAIN'] = thr_all
            out['TRAIN']['posterior_posthoc_best'] = None
            allr = [r for s in TRAIN for r in rows[s]]
            t_, k_, tp_, f_ = curve([a for a, _ in allr], [b for _, b in allr], sum(nref[t] for t in TRAIN))
            i = int(np.argmax(f_)); out['TRAIN']['posterior_posthoc_best'] = dict(thr=float(t_[i]), kept=int(k_[i]), tp=int(tp_[i]), F1=100 * float(f_[i]), note='oracle, not achievable')
            if all(s in chains for s in DEV):
                out['DEV']['posterior_thr_from_TRAIN'] = dict(summarize({s: filtered(chains[s], post[s], thr_all) for s in DEV}, refs, DEV), threshold=thr_all)
                dr = [r for s in DEV for r in rows[s]]
                out['DEV']['posterior_ranking'] = dict(zip(('AUROC', 'AP'), auroc_ap([a for a, _ in dr], [b for _, b in dr])))
    (run_p / f'eval_{label}.json').write_text(json.dumps(out, indent=1))
    for scope in ('TRAIN', 'DEV'):
        for k, v in out.get(scope, {}).items():
            if isinstance(v, dict) and 'pooled' in v:
                print(scope, k, {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in v['pooled'].items()},
                      {s: round(d['F1'], 2) for s, d in v['by_species'].items()}, flush=True)


if __name__ == '__main__':
    main()
