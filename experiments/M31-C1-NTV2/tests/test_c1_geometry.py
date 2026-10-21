"""CPU tests of the NT-v2 adapter geometry: tokenisation alignment, nearest-centre ownership, and training-path == inference-path potentials."""
import numpy as np
import pytest
import torch
import ntv2_adapter as A
import m31_common as M
from src.m28.contig_emission_lattice import stitch_contig_emissions
import manifest as MF


class FakeVocab(A.NTVocab):
    def __init__(self):
        import itertools
        v = {''.join(p): 6 + i for i, p in enumerate(itertools.product('ACGT', repeat=6))}
        v.update({'<cls>': 3, 'N': 4106})
        super().__init__(v)


def test_tokenisation_alignment():
    v = FakeVocab()
    s = 'ACGTAC' * 3 + 'ACGNAC' + 'NNNNNN' + 'ACG'          # 3 normal chunks, 2 N-containing chunks, short tail
    ids = v.ids(s, width=60)
    assert len(ids) == 1 + 10 and ids[0] == 3
    assert ids[1] == ids[2] == ids[3] == v.v['ACGTAC']
    assert ids[4] == ids[5] == ids[6] == ids[7] == v.n or ids[4] == v.n         # N-containing chunk -> single N id
    assert ids[5] == v.n and ids[6] == v.n and ids[-1] == v.n                    # tail 'ACG'+N padding contains N -> N id
    full = v.ids('ACGTAC' * 2048)
    assert len(full) == 2049
    with pytest.raises(ValueError):
        v.ids('A' * (A.NT_WIDTH + 1))


@pytest.mark.parametrize('length', [3000, 12288, 12289, 20000, 49152, 30001])
def test_nt_geometry_covers_and_owns_nearest_centre(length):
    g = A.nt_geometry(0, length)
    wins = g['local_windows']
    assert wins[0][0] == 0 and wins[-1][1] == length
    assert all(b - a <= A.NT_WIDTH for a, b in wins)
    segs = g['base_owned_segments']
    assert segs[0][1] == 0 and segs[-1][2] == length
    assert all(s[2] == t[1] for s, t in zip(segs, segs[1:]))
    # brute force nearest centre (ties -> lower start)
    own = np.empty(length, int)
    for i, (_, l, r) in enumerate(segs):
        own[l:r] = i
    centres = np.array([(a + b) / 2 for a, b in wins])
    for p in np.random.default_rng(0).integers(0, length, 300):
        d = np.abs(centres - (p + .5))
        cand = [i for i, (a, b) in enumerate(wins) if a <= p < b]
        best = min(cand, key=lambda i: (d[i], i))
        assert own[p] == best


@pytest.mark.parametrize('strand', ['+', '-'])
def test_training_ownership_equals_inference_stitch(strand):
    """Whole-contig context: frozen `owned_context_logits` (training) == frozen `stitch_contig_emissions` (inference)."""
    L = 40000
    rng = np.random.default_rng(1)
    seq = ''.join(rng.choice(list('ACGT'), L))
    g = A.nt_geometry(0, L)
    oriented = []
    for a, b in g['genomic_windows']:
        s = seq[a:b]
        oriented.append(s if strand == '+' else MF.rc(s))
    pot = [torch.from_numpy(rng.standard_normal((len(s), 15)).astype(np.float32)) for s in oriented]
    train = M.G.owned_context_logits(pot, g, strand)
    recs = [dict(window_index=i, sequence=s, local_potentials=p) for i, (s, p) in enumerate(zip(oriented, pot))]
    infer = stitch_contig_emissions(seq, list(g['genomic_windows']), strand, recs)['local_potentials']
    assert train.shape == infer.shape == (L, 15)
    assert torch.equal(train, infer)
