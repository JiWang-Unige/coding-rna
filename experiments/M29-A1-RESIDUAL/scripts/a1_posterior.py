"""Gene posterior log P(g in gene set) of every decoded chain under the A1 energy (head + fixed offset), whole contig / strand.
Uses the frozen `complete_gene_set_posteriors` with the qualified scaled-mass messages (as the frozen TRAIN readout does).
Usage: a1_posterior.py <run_dir> <seq> <strand>   -> <run_dir>/posterior/<seq>_<strand>.json   (run_dir/head must exist; gate: head absent)
"""
import json, sys, time
from unittest.mock import patch
import numpy as np
import torch
import a1_common as K
import a1_offset as O
import m29_common as M
K.repoint()
import src.m28.full_grammar_gene_posterior as G
from src.m28.full_grammar_scaled_mass import _scaled_prefix_queries, _scaled_suffix_queries


def main():
    run, seq, strand = K.Path(sys.argv[1]), sys.argv[2], sys.argv[3]
    torch.set_num_threads(1)
    out = run / 'posterior' / f'{seq}_{strand}.json'
    out.parent.mkdir(exist_ok=True)
    if out.exists():
        return
    import re
    chains_json = json.load(open(run / 'chains' / f'{seq}_{strand}.json'))
    c = O.load_tracks(seq)[strand]
    L = c.shape[1]
    hp = run / 'head' / f'{seq}_{strand}.npy'
    lg = np.load(hp) if hp.exists() else np.zeros((L, 15), np.float32)
    for a in range(0, L, 2_000_000):
        b = min(L, a + 2_000_000)
        lg[a:b] += O.expand(c, a, b)
    dna = re.sub('[^ACGT]', 'N', M.read_sequence(seq))
    dna = dna if strand == '+' else O.rc(dna)
    chains = []
    for ch in chains_json['chains']:
        ex = [tuple(e) for e in ch['exons']]
        chains.append(tuple(ex) if strand == '+' else tuple((L - b, L - a) for a, b in reversed(ex)))
    order = sorted(range(len(chains)), key=lambda i: chains[i][0][0])
    t0 = time.perf_counter()
    with patch.multiple(G, _prefix_queries=_scaled_prefix_queries, _suffix_queries=_scaled_suffix_queries):
        r = G.complete_gene_set_posteriors(dna, torch.from_numpy(lg), [chains[i] for i in order], 128)
    rows = [dict(exons=chains_json['chains'][i]['exons'], gain=chains_json['chains'][i]['gain'], log_p=float(lp), energy=float(e))
            for i, lp, e in zip(order, r['gene_log_probabilities'], r['chain_gains'])]
    json.dump(dict(seq=seq, strand=strand, genes=rows, log_partition=r['log_partition'], reverse_log_partition=r['reverse_log_partition'],
                   seconds=time.perf_counter() - t0), open(out, 'w'))
    print(seq, strand, len(rows), f'{time.perf_counter()-t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
