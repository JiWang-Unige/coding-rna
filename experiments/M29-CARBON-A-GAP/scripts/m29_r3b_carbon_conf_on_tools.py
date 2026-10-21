"""R3b: use Carbon-A's own confidence network (14,433 params) as a gene-level filter on OTHER tools' chains.

For every tool, every predicted chain is converted to oriented inclusive exons on its strand track, the strand's decoded roles are
built from that tool's whole chain set, and  confidence_network.inference(p_CDS_track, exons, scored=all)  returns P(structure exactly right).
Input track = Carbon-A raw probabilities (the same 1.2B-encoder signal all our heads use), CPU float32.
Compared on the SAME prediction pool with our posterior (EM pool only) and with Carbon-A's own confidence.
Output: outputs/M29-CARBON-A-GAP/R3b/{conf_<tool>_<seq>.json, summary.json}
"""
import json
import os
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'carbon_a_pipeline_copy' / 'src'))
import torch  # noqa: E402
torch.set_num_threads(int(os.environ.get('M29_THREADS', '3')))
import m29_common as M  # noqa: E402
from m29_load_preds import tool_chains  # noqa: E402
from m29_r3_posterior_filter import auroc_ap, curve  # noqa: E402
from m29_gate_parity import load_tracks   # noqa: E402  (imports the Carbon-A copy)
from constants import CONFIDENCE_NETWORK_PATH  # noqa: E402
from confidence_network import load_confidence_network  # noqa: E402

OUT = M.M29 / 'R3b'
OUT.mkdir(parents=True, exist_ok=True)
TOOLS = ['E_lora', 'EM_pool', 'ANNEVO', 'Tiberius', 'B_freshjoint']


def oriented_exons(chain, n):
    ex = chain['exons']            # genomic half-open, ascending
    if chain['strand'] == '+':
        return [(a, b - 1) for a, b in ex]
    return [(n - b, n - 1 - a) for a, b in reversed(ex)]


def main():
    tools = sys.argv[1].split(',') if len(sys.argv) > 1 else TOOLS
    seqs = sys.argv[2].split(',') if len(sys.argv) > 2 else [s for _, s, _ in M.SCOPE]
    cn = load_confidence_network(CONFIDENCE_NETWORK_PATH, torch.device('cpu'))
    for seq in seqs:
        rid, sequence, pos, neg = load_tracks(seq)
        n = len(sequence)
        tracks = {'+': np.ascontiguousarray(pos), '-': np.ascontiguousarray(neg[::-1])}
        for tool in tools:
            f = OUT / f'conf_{tool}_{seq}.json'
            if f.exists():
                continue
            chains = list(tool_chains(tool, seq).values())
            res = {}
            for strand in '+-':
                sel = [c for c in chains if c['strand'] == strand]
                if not sel:
                    continue
                ex = [oriented_exons(c, n) for c in sel]
                conf = cn.inference(tracks[strand], ex, list(range(len(sel))))
                for c, v in zip(sel, conf):
                    res[json.dumps(c['key'][2:])] = (c['strand'], float(v))
            f.write_text(json.dumps(res))
            print(seq, tool, len(res), flush=True)


def summarize():
    refs = {}
    for _, seq, _ in M.SCOPE:
        refs[seq] = {M.key(t) for t in M.build_reference(seq)['refs']}
    nref = sum(len(v) for v in refs.values())
    out = {}
    for tool in TOOLS:
        score, label, ca_score = [], [], []
        for _, seq, _ in M.SCOPE:
            f = OUT / f'conf_{tool}_{seq}.json'
            if not f.exists():
                break
            d = json.loads(f.read_text())
            for k, (strand, v) in d.items():
                key = (seq, strand, tuple(tuple(x) for x in json.loads(k)[0]))
                score.append(v); label.append(int(key in refs[seq]))
        else:
            auc, ap = auroc_ap(score, label)
            t, kept, tp, f1 = curve(score, label, nref)
            thr = {}
            for th in (0.05, 0.1, 0.2, 0.3, 0.5):
                m = np.asarray(score) >= th
                tpv = int(np.asarray(label)[m].sum()); kv = int(m.sum())
                thr[str(th)] = dict(kept=kv, tp=tpv, F1=100 * 2 * tpv / (kv + nref))
            k = int(np.argmax(f1))
            out[tool] = dict(n=len(score), hits=int(sum(label)), F1_unfiltered=100 * 2 * sum(label) / (len(score) + nref), auroc=auc, ap=ap,
                             best_posthoc=dict(thr=float(t[k]), kept=int(kept[k]), tp=int(tp[k]), F1=100 * float(f1[k])),
                             at_thresholds=thr)
    (OUT / 'summary.json').write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'summarize':
        summarize()
    else:
        main()
