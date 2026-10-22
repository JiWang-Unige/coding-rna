"""Per-base no-gene marginal vs brute force.  CPU only, synthetic DNA; no real data.
 (1) exhaustive enumeration of every complete gene-set path of the frozen automaton on a short sequence (cap 3): log Z and every per-base P(role 0)
     agree with the forward/backward query to 1e-9;
 (2) the same for the production cap 128 (identical on the short sequence);
 (3) log Z equals the frozen `complete_gene_set_posteriors` log partition, and for every decoded chain g, P(no gene) <= 1 - P(g) inside g;
 (4) segment helpers round-trip between oriented and genomic frames.
Run: pytest -q test_nogene_marginal.py   (or python test_nogene_marginal.py)"""
import math
import sys

import numpy as np
import torch

import m36_common as W
import m36_nogene_marginal as NG
from src.m28.full_grammar_boundary_messages import BACKGROUND, _transitions


def enumerate_paths(seq, cap):
    """every role sequence of a B ... B path (= every complete, nonoverlapping gene set with explicit exon structure)"""
    out = []
    sys.setrecursionlimit(10000)

    def rec(i, state, roles):
        if i == len(seq):
            if state == BACKGROUND:
                out.append(list(roles))
            return
        for target, role in _transitions(state, seq[i], cap):
            roles.append(role)
            rec(i + 1, target, roles)
            roles.pop()
    rec(0, BACKGROUND, [])
    return out


def brute_force(seq, logits, cap):
    gain = (logits.double() - logits[:, :1].double()).numpy()
    paths = enumerate_paths(seq, cap)
    E = np.array([sum(gain[i, r] for i, r in enumerate(p)) for p in paths])
    m = E.max()
    logz = m + math.log(np.exp(E - m).sum())
    R = np.array(paths)                       # [paths, n]
    w = np.exp(E - logz)
    p_nogene = np.array([w[R[:, i] == 0].sum() for i in range(len(seq))])
    return logz, p_nogene, len(paths), int((R.sum(1) > 0).sum())


def planted(seed=0, n=72, two_genes=False):
    rng = np.random.default_rng(seed)
    dna = list(rng.choice(list('ACGT'), n))
    gene = 'ATGGCTGAAGTAAGTTTCAGGCTAAATAA'            # ATG GCT GA | GT AAGTTTC AG | A GCT AAA TAA  (exon1 8 bp, intron 11 bp, exon2 10 bp)
    i = 20
    dna[i:i + len(gene)] = list(gene)
    if two_genes:                                      # a second, single-exon gene downstream -> gene sets with 0, 1 or 2 genes
        g2 = 'ATGAAAGCTGAATAA'
        dna[i + len(gene) + 6:i + len(gene) + 6 + len(g2)] = list(g2)
    dna = ''.join(dna)
    lg = rng.normal(0, 1, (n, 15)).astype(np.float32)
    lg[:, 0] += 0.3
    return dna, torch.from_numpy(lg)


def test_brute_force_cap3_and_cap128():
    for seed, two in ((0, False), (1, False), (2, False), (3, True), (4, True)):
        dna, lg = planted(seed, n=72 if not two else 96, two_genes=two)
        logz, p_bf, n_paths, n_gene_paths = brute_force(dna, lg, 3)
        assert n_gene_paths > 0, 'planted sequence produced no gene path'
        for cap in (3, 128):
            r = NG.no_gene_marginal(dna, lg, max_exons=cap)
            assert abs(r['log_partition'] - logz) < 1e-9 * max(1, abs(logz)), (cap, r['log_partition'], logz)
            assert np.allclose(r['p_nogene'], p_bf, atol=1e-9), (cap, np.abs(r['p_nogene'] - p_bf).max())
            assert (r['p_nogene'] <= 1 + 1e-12).all() and (r['p_nogene'] >= 0).all()
        R = np.array(enumerate_paths(dna, 3))
        n_two = int(((R == 13).sum(1) >= 6).sum())    # paths with two start codons = two-gene sets
        print(f'seed {seed}: paths {n_paths} (with genes {n_gene_paths}, two-gene sets {n_two}) logZ {logz:.6f}  min P(no gene) {p_bf.min():.4f}')
        if two:
            assert n_two > 0


def test_against_frozen_posterior_and_chains():
    import src.m28.full_grammar_gene_posterior as G
    from m29_skyline_reuse import decode_full_grammar_skyline
    dna, lg = NG.synthetic(6000, seed=3, bg=0.2)
    r = NG.no_gene_marginal(dna, lg)
    z = G.complete_gene_set_posteriors(dna, lg, [], W.MAX_EXONS)
    assert abs(z['log_partition'] - r['log_partition']) < 1e-6
    d = decode_full_grammar_skyline(dna, lg, W.MAX_EXONS)
    chains = [tuple(map(tuple, c)) for c in d['chains']]
    assert chains, 'no chain decoded on the synthetic problem'
    post = G.complete_gene_set_posteriors(dna, lg, sorted(chains, key=lambda c: c[0][0]), W.MAX_EXONS)
    for c, pg in zip(post['chains'], post['gene_probabilities']):
        a, b = c[0][0], c[-1][1]
        assert (r['p_nogene'][a:b] <= 1 - pg + 1e-9).all(), (c, pg, r['p_nogene'][a:b].max())
    print('chains', len(chains), 'max P(gene)', max(post['gene_probabilities']), 'mean P(no gene)', float(r['p_nogene'].mean()))


def test_segments_round_trip():
    p = np.array([0.995, 0.999, 0.2, 0.99, 0.99, 0.5, 1.0])
    seg = NG.segments(p, 0.99)
    assert seg == [(0, 2), (3, 5), (6, 7)]
    assert NG.segments(p, 0.99, min_len=2) == [(0, 2), (3, 5)]
    g = NG.to_genomic(seg, '-', 7)
    assert g == [(0, 1), (2, 4), (5, 7)]
    assert NG.to_genomic(g, '-', 7) == seg
    import m36_seed_labels as SL
    assert SL.pseudo_background(p, '-', 7, 0.99) == g


if __name__ == '__main__':
    W.K.repoint()
    test_brute_force_cap3_and_cap128()
    test_against_frozen_posterior_and_chains()
    test_segments_round_trip()
    print('ALL OK')
