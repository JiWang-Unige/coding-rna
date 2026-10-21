"""Evaluation-set (B) scoring of one species: zero-shot arms + adapted settings, with gene-unit paired bootstrap (protocol §6, §8).
Usage: m30_evalB.py <species> <arm,arm,...>   arms: CarbonA_unf CarbonA_filt Gate Gate_post <setting dirs under outputs/M30/evalB/>  (each adapted dir also adds <name>+conf0.1)
-> outputs/M30/results_<species>.json"""
import json, sys
from pathlib import Path
import numpy as np
import torch
import m30_common as M
import a1_common as K
sys.path.insert(0, str(K.CA_SRC))
import m30_score as SC
import m30_offsets as OF
from constants import CONFIDENCE_NETWORK_PATH     # noqa: E402
from confidence_network import load_confidence_network   # noqa: E402

torch.set_num_threads(4)
_TR = {}


def conf_scores(contig, chains, cn):
    if contig not in _TR:
        rid, seq, pos, neg = OF.load(contig)
        _TR[contig] = (len(seq), {'+': np.ascontiguousarray(pos), '-': np.ascontiguousarray(neg[::-1])})
    n, tr = _TR[contig]
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


def main():
    species, arms = sys.argv[1], sys.argv[2].split(',')
    B = M.SPECIES[species]['B']
    ev = SC.Evaluator(B)
    cn = load_confidence_network(CONFIDENCE_NETWORK_PATH, torch.device('cpu'))
    res = {}
    thr = json.load(open(M.OUT / 'gate_posterior_threshold.json')) if 'Gate_post' in arms else None
    for arm in arms:
        if arm == 'CarbonA_unf':
            ch = {c: SC.load_gff_chains(str(M.OUT / 'carbonA_raw' / species / c / 'prediction.gff'), c, cds96=False) for c in B}
        elif arm == 'CarbonA_filt':
            ch = {c: SC.load_gff_chains(str(M.OUT / 'carbonA_raw' / species / c / 'prediction.gff'), c, cds96=False, conf_filter='carbon') for c in B}
        elif arm == 'Gate':
            ch = {c: SC.load_gff_chains(str(M.OUT / 'evalB' / 'gate' / c / 'prediction' / 'full_grammar.gff3'), c) for c in B}
        elif arm == 'Gate_post':
            ch = {}
            for c in B:
                full = SC.load_gff_chains(str(M.OUT / 'evalB' / 'gate' / c / 'prediction' / 'full_grammar.gff3'), c)
                lp = {}
                for st in '+-':
                    for g in json.load(open(M.OUT / 'posterior_gate_tau0.03_B' / f'{c}_{st}.json'))['genes']:
                        lp[(c, st, tuple(tuple(e) for e in g['exons']))] = g['log_p']
                ch[c] = {k: v for k, v in full.items() if lp[k] >= thr['threshold']}
        else:
            ch = {c: SC.load_gff_chains(str(M.OUT / 'evalB' / arm / c / 'prediction' / 'full_grammar.gff3'), c) for c in B}
        res[arm] = ev.add(arm, ch)
        print(arm, res[arm], flush=True)
        if arm not in ('CarbonA_unf', 'CarbonA_filt', 'Gate_post'):
            cf = {}
            for c in B:
                sc = conf_scores(c, ch[c], cn)
                cf[c] = {k: v for k, v in ch[c].items() if sc[k] >= 0.1}
            res[arm + '+conf0.1'] = ev.add(arm + '+conf0.1', cf)
            print(arm + '+conf0.1', res[arm + '+conf0.1'], flush=True)
    names = list(ev.methods)
    ref_arm = 'CarbonA_filt' if 'CarbonA_filt' in names else names[0]
    pairs = [(n, ref_arm) for n in names if n != ref_arm]
    for a in names:                      # U - F pairs by naming convention  <prefix>_U<suffix> / <prefix>_F<suffix>
        for suf in ('', '+conf0.1'):
            if a.endswith('_U' + suf) and '_U' in a:
                b2 = a[:len(a) - len('_U' + suf)] + '_F' + suf
                if b2 in names:
                    pairs.append((a, b2))
    have = set(pairs)
    extra = [(n, 'Gate') for n in names if n not in ('Gate', ref_arm) and (n, 'Gate') not in have]
    extra += [('Gate_post', 'Gate'), ('Gate_post', 'Gate+conf0.1'), ('Gate+conf0.1', 'Gate')]
    ks = [n for n in names if '_k100_FULL' in n]
    for a in ks:                                    # k=100% vs k=25% U / F, and vs +conf0.1 counterparts
        suf = '+conf0.1' if a.endswith('+conf0.1') else ''
        for tr in ('U', 'F'):
            b2 = a.replace('_k100_FULL', '_k25_' + tr).replace('+conf0.1', '') + suf
            if b2 in names:
                extra.append((a, b2))
        if suf == '' and a.replace('_FULL', '_FULL+conf0.1') in names:
            pass
    for tr in ('U', 'F'):                           # k25 +conf0.1 vs Gate+conf0.1
        for n in names:
            if n.endswith(f'_k25_{tr}+conf0.1'):
                extra.append((n, 'Gate+conf0.1'))
    pairs = pairs + [p for p in dict.fromkeys(extra) if p not in have and p[0] in names and p[1] in names]
    cis, diffs = ev.bootstrap(pairs)
    out = dict(species=species, B=B, tables={n: ev.table(n) for n in names}, bootstrap_F1=cis, bootstrap_diffs=diffs, reference_arm=ref_arm)
    if thr:
        out['gate_posterior_threshold'] = thr
    json.dump(out, open(M.OUT / f'results_{species}_final.json', 'w'), indent=1)
    for n in names:
        t = out['tables'][n]
        sts = ', '.join(k + ':' + format(v['F1'], '.1f') for k, v in t['strata'].items())
        print(n + ': F1 ' + format(t['pooled']['F1'], '.2f') + ' CI [' + format(cis[n]['lo'], '.2f') + ',' + format(cis[n]['hi'], '.2f') + '] strata {' + sts + '} base ' + format(t['base']['base_F1'], '.2f') + ' bg ' + format(t['base']['bg_burden_pct'], '.2f') + '%')
    for k, v in diffs.items():
        print(f"  {k}: {v['point']:+.2f} [{v['lo']:+.2f},{v['hi']:+.2f}]")


if __name__ == '__main__':
    main()
