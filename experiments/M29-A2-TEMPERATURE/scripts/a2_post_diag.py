"""A2 step 1b (CPU): gene-level posterior of the GATE model's decoded chains at inverse temperature tau, region-limited.
log P(g in gene set) under E_tau = tau*offset (head = 0), DNA/energies restricted to the region [a,b) of one contig, one strand.
Chains = gate decode (A1 gate_ownpath) fully inside the region with 50 kb margins; label = exact match to a primary reference chain.
Usage: a2_post_diag.py seq strand a_kb b_kb taus   -> A2/posterior_diag/<seq>_<strand>_<a>_<b>.json
"""
import json, re, sys, time
from unittest.mock import patch
import numpy as np
import torch
import a2_common as C2
import a1_common as K
import a1_offset as O
import m29_common as M
K.repoint()
import src.m28.full_grammar_gene_posterior as G
from src.m28.full_grammar_scaled_mass import _scaled_prefix_queries, _scaled_suffix_queries

def main():
    seq, strand, a, b = sys.argv[1], sys.argv[2], int(sys.argv[3]) * 1000, int(sys.argv[4]) * 1000
    taus = [float(x) for x in sys.argv[5].split(',')]
    torch.set_num_threads(1)
    out = C2.A2 / 'posterior_diag' / f'{seq}_{strand}_{a // 1000}_{b // 1000}.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    c = O.load_tracks(seq)[strand]
    L = c.shape[1]
    # oriented region
    ra, rb = (a, b) if strand == '+' else (L - b, L - a)
    lg = O.expand(c, ra, rb)
    dna = re.sub('[^ACGT]', 'N', M.read_sequence(seq))
    dna = dna if strand == '+' else O.rc(dna)
    dna = dna[ra:rb]
    refk = {M.key(t) for t in M.build_reference(seq)['refs']}
    chains, labels, gen = [], [], []
    for ch in json.load(open(K.A1 / 'gate_ownpath' / 'chains' / f'{seq}_{strand}.json'))['chains']:
        ex = [tuple(e) for e in ch['exons']]
        if ex[0][0] < a + 50_000 or ex[-1][1] > b - 50_000: continue
        if sum(y - x for x, y in ex) < 96: continue
        oc = tuple((x - a, y - a) for x, y in ex) if strand == '+' else tuple((L - y - ra, L - x - ra) for x, y in reversed(ex))
        chains.append(oc); gen.append(ex)
        labels.append(int((seq, strand, tuple(ex)) in refk))
    order = sorted(range(len(chains)), key=lambda i: chains[i][0][0])
    chains = [chains[i] for i in order]; labels = [labels[i] for i in order]
    res = dict(seq=seq, strand=strand, region=[a, b], n_chains=len(chains), n_hits=sum(labels), taus={})
    for tau in taus:
        t0 = time.time()
        with patch.multiple(G, _prefix_queries=_scaled_prefix_queries, _suffix_queries=_scaled_suffix_queries):
            r = G.complete_gene_set_posteriors(dna, torch.from_numpy((tau * lg).astype(np.float32)), chains, 128)
        res['taus'][str(tau)] = dict(log_p=[float(x) for x in r['gene_log_probabilities']], seconds=time.time() - t0)
        res['labels'] = labels
        json.dump(res, open(out, 'w'))
        print(seq, strand, 'tau', tau, f'{time.time()-t0:.0f}s', flush=True)

if __name__ == '__main__':
    main()
