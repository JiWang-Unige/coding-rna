"""CPU tests of the PlantCAD2 adapter: window geometry / nearest-centre ownership, six-base pooling, RC orientation, and
training-path == inference-path potentials (with an RC-equivariant fake encoder)."""
import numpy as np
import pytest
import torch
import pc2_adapter as A
import m35_common as M
import c1_model as MD
import c1_infer as I
from src.m28.contig_emission_lattice import stitch_contig_emissions
import manifest as MF


@pytest.mark.parametrize('length', [3000, 8192, 8193, 20000, 49152, 30001])
def test_pc_geometry_covers_and_owns_nearest_centre(length):
    g = A.pc_geometry(0, length)
    wins = g['local_windows']
    assert wins[0][0] == 0 and wins[-1][1] == length
    assert all(b - a <= A.PC_WIDTH for a, b in wins)
    segs = g['base_owned_segments']
    assert segs[0][1] == 0 and segs[-1][2] == length
    assert all(s[2] == t[1] for s, t in zip(segs, segs[1:]))
    own = np.empty(length, int)
    for i, (_, l, r) in enumerate(segs):
        own[l:r] = i
    centres = np.array([(a + b) / 2 for a, b in wins])
    for p in np.random.default_rng(0).integers(0, length, 300):
        d = np.abs(centres - (p + .5))
        cand = [i for i, (a, b) in enumerate(wins) if a <= p < b]
        best = min(cand, key=lambda i: (d[i], i))
        assert own[p] == best


def test_pool_tokens_mean_over_valid_bases():
    rng = np.random.default_rng(0)
    for n in (8192, 8190, 8191, 4097, 6, 1):
        f = torch.from_numpy(rng.standard_normal((n, 5)).astype(np.float32))
        t = A.pool_tokens(f).float()
        assert t.shape == (A.HEAD_TOKENS, 5)
        k = (n + 5) // 6
        ref = torch.stack([f[6 * i:6 * i + 6].mean(0) for i in range(k)])
        assert torch.allclose(t[:k], ref, atol=3e-2, rtol=1e-2)       # BF16 rounding of the output
        assert float(t[k:].abs().sum()) == 0.


def test_orient_is_involution_and_minus_flips_position_and_channel():
    f = torch.arange(12.).view(4, 3)
    assert torch.equal(A.orient(f, '+'), f)
    m = A.orient(f, '-')
    assert torch.equal(m, f.flip(0).flip(1)) and torch.equal(A.orient(m, '-'), f)


class FakeEnc:
    """Per-base RC-equivariant fake encoder: f(RC(x)) = flip(f(x), pos & channel) holds exactly (A/T and C/G embeddings are channel-flips of each other,
    N is channel-palindromic) plus a position-independent context term (window mean) that is itself RC-symmetric."""
    def __init__(self, seed=0):
        g = torch.Generator().manual_seed(seed)
        a, c = torch.randn(A.HIDDEN, generator=g), torch.randn(A.HIDDEN, generator=g)
        n = torch.randn(A.HIDDEN, generator=g)
        n = n + n.flip(0)
        self.e = {'A': a, 'T': a.flip(0), 'C': c, 'G': c.flip(0), 'N': n}


def fake_encode(enc):
    def f(model, sequences, layer=-1, batch=4, device='cpu'):
        out = []
        for s in sequences:
            x = torch.stack([enc.e[ch] for ch in s])
            # local smoothing keeps equivariance (symmetric 3-tap kernel)
            xp = torch.nn.functional.pad(x.T[None], (1, 1), mode='replicate')
            y = (xp[..., :-2] + 2 * xp[..., 1:-1] + xp[..., 2:])[0].T / 4
            out.append(y.to(torch.bfloat16))
        return out
    return f


def test_fake_encoder_is_rc_equivariant():
    enc = FakeEnc()
    s = ''.join(np.random.default_rng(3).choice(list('ACGTN'), 500))
    f = fake_encode(enc)(None, [s])[0]
    fr = fake_encode(enc)(None, [MF.rc(s)])[0]
    assert torch.equal(fr, A.orient(f, '-'))


@pytest.mark.parametrize('strand', ['+', '-'])
def test_training_path_equals_inference_stitch(strand, monkeypatch):
    L = 40001
    rng = np.random.default_rng(1)
    seq = ''.join(rng.choice(list('ACGTN'), L, p=[.245, .245, .245, .245, .02]))
    monkeypatch.setattr(A, 'encode_windows', fake_encode(FakeEnc()))
    torch.manual_seed(0)
    model = MD.new_head().eval()
    # training path (what c1_train does for one draw)
    g = A.pc_geometry(0, L)
    gen = MD.genomic_window_strings(seq, g)
    oriented = MD.orient_strings(gen, strand)
    toks = A.oriented_tokens(A.encode_windows(None, gen), strand)
    with torch.no_grad():
        train = MD.head_context_logits(model, toks, oriented, g, strand, device='cpu')
    # inference path (what c1_infer does for a region): both strands from one forward
    pots, n = I.region_potentials(model, None, -1, seq, batch=3, device='cpu')
    assert n == len(g['genomic_windows'])
    assert pots[strand].shape == train.shape == (L, 15)
    assert torch.allclose(pots[strand], train, atol=2e-5, rtol=1e-4), float((pots[strand] - train).abs().max())


def test_both_strands_features_consistent_with_separate_forward(monkeypatch):
    """Deriving the minus-strand features from the genomic forward == encoding the reverse-complement window directly."""
    monkeypatch.setattr(A, 'encode_windows', fake_encode(FakeEnc(1)))
    s = ''.join(np.random.default_rng(2).choice(list('ACGT'), 8192))
    direct = A.encode_windows(None, [MF.rc(s)])[0]
    derived = A.orient(A.encode_windows(None, [s])[0], '-')
    assert torch.equal(direct, derived)
