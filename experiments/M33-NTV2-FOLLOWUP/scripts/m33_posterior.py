"""Gene-set posterior log P(g in gene set) of model B's decoded chains under E_tau = tau * (C1 head potentials), whole contig / strand
(frozen `complete_gene_set_posteriors` + scaled-mass queries, exactly as M30 `m30_posterior.py` / M29-A2; only the energy source differs:
there energy = tau*(offset+head), here energy = tau*pot because model B has no external offset).
Usage: m33_posterior.py <run_dir> <seq> <strand> <taus comma-separated> <out_dir>
  reads  <run_dir>/pot/<seq>_<strand>.npy  (oriented [L,15] float32 kept by m33_infer.py)  and  <run_dir>/chains/<seq>_<strand>.json
  writes <out_dir>/tau<tau>/<seq>_<strand>.json   {seq,strand,tau,genes:[{exons,gain,log_p}],log_partition,seconds}   (resumable per tau)
"""
import json, re, sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
EXP = HERE.parents[0].parent
for p in (str(EXP / 'M31-C1-NTV2' / 'scripts'), str(EXP / 'M30-HOLDOUT-ADAPT' / 'scripts')):
    if p not in sys.path:
        sys.path.insert(0, p)
from unittest.mock import patch
import numpy as np
import torch
import m30_common as M30
import a1_common as K
M29 = M30.register()
K.repoint()
import manifest as MF
import src.m28.full_grammar_gene_posterior as G
from src.m28.full_grammar_scaled_mass import _scaled_prefix_queries, _scaled_suffix_queries


def oriented_chains(chain_rows, strand, L):
    """genomic exon lists (as stored by m33_infer) -> oriented (start,end) tuples, left to right in the oriented frame"""
    out = []
    for ch in chain_rows:
        ex = [tuple(e) for e in ch['exons']]
        out.append(tuple(ex) if strand == '+' else tuple((L - b, L - a) for a, b in reversed(ex)))
    return out


def posterior(dna_oriented, lg_scaled, chains):
    order = sorted(range(len(chains)), key=lambda i: chains[i][0][0])
    with patch.multiple(G, _prefix_queries=_scaled_prefix_queries, _suffix_queries=_scaled_suffix_queries):
        r = G.complete_gene_set_posteriors(dna_oriented, lg_scaled, [chains[i] for i in order], 128)
    return order, r


def main():
    run, seq, strand = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
    taus = [float(x) for x in sys.argv[4].split(',')]
    outd = Path(sys.argv[5])
    torch.set_num_threads(1)
    cj = json.load(open(run / 'chains' / f'{seq}_{strand}.json'))
    lg0 = np.load(run / 'pot' / f'{seq}_{strand}.npy', mmap_mode='r')
    L = lg0.shape[0]
    dna = re.sub('[^ACGT]', 'N', M29.read_sequence(seq))
    assert len(dna) == L
    dna = dna if strand == '+' else MF.rc(dna)
    chains = oriented_chains(cj['chains'], strand, L)
    for tau in taus:
        out = outd / f'tau{tau:g}' / f'{seq}_{strand}.json'
        if out.exists():
            continue
        out.parent.mkdir(parents=True, exist_ok=True)
        t0 = time.perf_counter()
        lg = torch.from_numpy(np.ascontiguousarray(lg0 * np.float32(tau)))
        order, r = posterior(dna, lg, chains)
        del lg
        rows = [dict(exons=cj['chains'][i]['exons'], gain=cj['chains'][i]['gain'], log_p=float(lp)) for i, lp in zip(order, r['gene_log_probabilities'])]
        json.dump(dict(seq=seq, strand=strand, tau=tau, genes=rows, log_partition=r['log_partition'], seconds=time.perf_counter() - t0), open(str(out) + '.tmp', 'w'))
        Path(str(out) + '.tmp').replace(out)
        print(seq, strand, f'tau={tau:g}', len(rows), f'{time.perf_counter() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
