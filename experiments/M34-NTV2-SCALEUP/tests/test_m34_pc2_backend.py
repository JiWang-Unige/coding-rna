"""PlantCAD2 backend (M34_BACKEND=pc2, env ~/envs/plantcad2) CPU tests with an RC-equivariant pseudo-encoder: window geometry, genomic/oriented strings,
training path (window_pots + frozen owned_context_logits) == inference path (region_potentials), identity of both strands from one forward, gradients reach the head."""
import types

import numpy as np
import pytest
import torch

import m34_train_common as T

pytestmark = pytest.mark.skipif(T.BACKEND != 'pc2', reason='run with M34_BACKEND=pc2 under ~/envs/plantcad2')


CTX = types.SimpleNamespace(enc=None, layer=-1)


def fake_features(strings, device='cpu'):
    rng = np.random.default_rng(7)
    X = {b: torch.from_numpy(rng.standard_normal(768).astype(np.float32)) for b in 'ACGT'}
    comp = dict(zip('ACGT', 'TGCA'))
    Y = {b: X[comp[b]].flip(0) for b in 'ACGT'}                      # v(comp b) = reverse(v(b)): exactly RC-equivariant under flip(position).flip(channel)
    out = []
    for s in strings:
        rows = [torch.cat([X[c], Y[c]]) if c in X else torch.zeros(1536) for c in s]
        out.append(torch.stack(rows).to(torch.bfloat16))
    return out


@pytest.fixture(autouse=True)
def patch(monkeypatch):
    import m34_backend as BK
    monkeypatch.setattr(BK, '_pc2_encode', lambda enc, strings, layer=-1, batch=8, grad=False, device='cpu': fake_features(strings))
    monkeypatch.setattr(BK.A, 'encode_windows', lambda enc, strings, layer=-1, batch=4, device='cpu': fake_features(strings))
    monkeypatch.setattr(BK.MD.G, 'one_hot', BK.MD.G.one_hot)


def test_geometry_and_strings():
    import m34_backend as BK
    rng = np.random.default_rng(1)
    chrom = ''.join(rng.choice(list('ACGT'), 60000))
    geo = BK.geometry(5000, 49000)
    assert geo['genomic_windows'][0][0] == 5000 and geo['genomic_windows'][-1][1] == 49000
    assert all(b - a <= 8192 for a, b in geo['genomic_windows'])
    plus, minus = BK.item_windows(chrom, geo, '+'), BK.item_windows(chrom, geo, '-')
    assert plus['seqs'] == plus['gen'] and minus['gen'] == plus['gen']
    assert minus['seqs'] == [T.Q.M.rc(s) for s in plus['gen']]


@pytest.mark.parametrize('strand', ['+', '-'])
def test_train_path_equals_inference_path(strand):
    import m34_backend as BK
    rng = np.random.default_rng(2)
    L = 30000
    seq = ''.join(rng.choice(list('ACGTN'), L, p=[.245, .245, .245, .245, .02]))
    head = BK.MD.new_head(0).eval()
    geo = BK.geometry(0, L)
    it = dict(BK.item_windows(seq, geo, strand), strand=strand)
    with torch.no_grad():
        pots = BK.window_pots(head, CTX, it, list(range(len(it['seqs']))), False, device='cpu')
        ctx = BK.T.G.owned_context_logits(pots, geo, strand)
        out, n = BK.region_potentials(head, CTX, seq, [strand], device='cpu')
    ref = out[strand]
    assert n == len(geo['genomic_windows'])
    assert ref.shape == ctx.shape == (L, 15) and torch.allclose(ref, ctx.cpu(), atol=2e-5, rtol=1e-4)


def test_grad_reaches_head_and_pool_shapes():
    import m34_backend as BK
    rng = np.random.default_rng(3)
    seq = ''.join(rng.choice(list('ACGT'), 9000))
    head = BK.MD.new_head(0)
    geo = BK.geometry(0, 9000)
    it = dict(BK.item_windows(seq, geo, '-'), strand='-')
    pots = BK.window_pots(head, CTX, it, [0, 1], True, device='cpu')
    assert [p.shape for p in pots] == [(len(it['seqs'][i]), 15) for i in (0, 1)]
    sum(p.pow(2).sum() for p in pots).backward()
    assert all(p.grad is not None for p in head.parameters())
