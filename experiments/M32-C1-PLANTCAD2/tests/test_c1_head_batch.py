"""Batched head call == per-window head calls (full and short/edge windows, with N)."""
import numpy as np
import torch
import pc2_adapter as A
import c1_model as MD


def test_batched_equals_single():
    torch.manual_seed(0)
    model = MD.new_head().eval()
    rng = np.random.default_rng(0)
    lens = [A.PC_WIDTH, A.PC_WIDTH, 7001, 4090, 6]
    seqs = [''.join(rng.choice(list('ACGTN'), n, p=[.24, .24, .24, .24, .04])) for n in lens]
    feats = [torch.randn(A.HEAD_TOKENS, A.HIDDEN).to(torch.bfloat16) for _ in lens]
    with torch.no_grad():
        batched = MD.window_potentials(model, feats, seqs, device='cpu')
        for f, s, b in zip(feats, seqs, batched):
            dna = MD.G.one_hot(s, A.HEAD_W)[None]
            single = model(f.float()[None], dna, [len(s)])['segmentation_logits'][0, :len(s)]
            assert b.shape == (len(s), 15)
            assert torch.allclose(single, b, atol=2e-5, rtol=1e-4), float((single - b).abs().max())
