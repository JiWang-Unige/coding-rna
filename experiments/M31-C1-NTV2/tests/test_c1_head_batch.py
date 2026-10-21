"""Batched head call == per-window head calls (full and short/edge windows, with N)."""
import numpy as np
import torch
import ntv2_adapter as A
import c1_model as MD


def test_batched_equals_single():
    torch.manual_seed(0)
    layers = (21, 27)
    model = MD.new_head(layers).eval()
    rng = np.random.default_rng(0)
    lens = [A.NT_WIDTH, A.NT_WIDTH, 9001, 4090, 6]
    seqs = [''.join(rng.choice(list('ACGTN'), n, p=[.24, .24, .24, .24, .04])) for n in lens]
    feats = [torch.randn(A.NT_TOKENS, A.HIDDEN * 2).to(torch.bfloat16) for _ in lens]
    with torch.no_grad():
        batched = MD.window_potentials(model, feats, seqs, device='cpu')
        for f, s, b in zip(feats, seqs, batched):
            dna = MD.G.one_hot(s, A.NT_WIDTH)[None]
            single = model(f.float()[None], dna, [len(s)])['segmentation_logits'][0, :len(s)]
            assert b.shape == (len(s), 15)
            assert torch.allclose(single, b, atol=2e-5, rtol=1e-4), float((single - b).abs().max())
