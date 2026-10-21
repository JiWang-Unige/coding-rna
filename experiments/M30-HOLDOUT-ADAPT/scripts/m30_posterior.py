"""Gene posterior log P(g in gene set) of decoded chains under E_tau = tau*(offset+head) (head optional), whole contig/strand.
Usage: m30_posterior.py <run_dir> <contig> <strand> <tau> <out_dir>      (run_dir/chains/<contig>_<strand>.json; run_dir/head/ used if present)
Works for the M30 contigs and for the legacy Arabidopsis/rice contigs (A1 offsets)."""
import json, re, sys, time
from unittest.mock import patch
import numpy as np
import torch
import m30_common as M
import a1_common as K
import a1_offset as O
M29 = M.register()
K.repoint()
import src.m28.full_grammar_gene_posterior as G
from src.m28.full_grammar_scaled_mass import _scaled_prefix_queries, _scaled_suffix_queries


def main():
    run, seq, strand, tau, outd = K.Path(sys.argv[1]), sys.argv[2], sys.argv[3], float(sys.argv[4]), K.Path(sys.argv[5])
    torch.set_num_threads(1)
    outd.mkdir(parents=True, exist_ok=True)
    out = outd / f'{seq}_{strand}.json'
    if out.exists():
        return
    p = M.offset_path(seq, strand)
    if not p.exists():
        p = K.A1 / 'offset' / f'{seq}_{strand}.npy'
    c = np.load(p, mmap_mode='r')
    L = c.shape[1]
    hp = run / 'head' / f'{seq}_{strand}.npy'
    lg = np.load(hp) if hp.exists() else np.zeros((L, 15), np.float32)
    for a in range(0, L, 2_000_000):
        b = min(L, a + 2_000_000)
        lg[a:b] += O.expand(c, a, b)
    lg *= np.float32(tau)
    dna = re.sub('[^ACGT]', 'N', M29.read_sequence(seq))
    dna = dna if strand == '+' else O.rc(dna)
    cj = json.load(open(run / 'chains' / f'{seq}_{strand}.json'))
    chains = []
    for ch in cj['chains']:
        ex = [tuple(e) for e in ch['exons']]
        chains.append(tuple(ex) if strand == '+' else tuple((L - b, L - a) for a, b in reversed(ex)))
    order = sorted(range(len(chains)), key=lambda i: chains[i][0][0])
    t0 = time.perf_counter()
    with patch.multiple(G, _prefix_queries=_scaled_prefix_queries, _suffix_queries=_scaled_suffix_queries):
        r = G.complete_gene_set_posteriors(dna, torch.from_numpy(lg), [chains[i] for i in order], 128)
    rows = [dict(exons=cj['chains'][i]['exons'], gain=cj['chains'][i]['gain'], log_p=float(lp)) for i, lp in zip(order, r['gene_log_probabilities'])]
    json.dump(dict(seq=seq, strand=strand, tau=tau, genes=rows, log_partition=r['log_partition'], seconds=time.perf_counter() - t0), open(out, 'w'))
    print(seq, strand, len(rows), f'{time.perf_counter()-t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
