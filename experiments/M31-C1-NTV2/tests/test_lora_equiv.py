"""GPU test: custom block-wise forward == HF hidden_states path; zero-initialised LoRA is the identity; LoRA state round-trips; gradients reach only LoRA."""
import torch
import ntv2_adapter as A
import numpy as np


def test_custom_forward_and_lora():
    esm, vocab = A.load_ntv2()
    rng = np.random.default_rng(0)
    seqs = [''.join(rng.choice(list('ACGT'), A.NT_WIDTH)) for _ in range(2)]
    ref = A.encode_windows(esm, vocab, seqs, (21, 27), norm=True)
    cus = A.encode_windows(esm, vocab, seqs, (21, 27), norm=True, custom=True)
    for a, b in zip(ref, cus):
        assert torch.equal(a, b) or float((a.float() - b.float()).abs().max()) < 1e-2
    params = A.add_lora(esm, list(range(12, 27)))
    assert len(params) == 15 * 4
    same = A.encode_windows(esm, vocab, seqs, (21, 27), norm=True, custom=True)
    for a, b in zip(cus, same):
        assert torch.equal(a, b)
    f = A.encode_windows(esm, vocab, seqs[:1], (21, 27), norm=True, grad=True)
    f[0].float().pow(2).sum().backward()
    assert all(p.grad is not None for p in params)
    assert esm.embeddings.word_embeddings.weight.grad is None
    with torch.no_grad():
        for p in params[1::2]:
            p.add_(torch.randn_like(p) * .01)
    st = A.lora_state(esm)
    out1 = A.encode_windows(esm, vocab, seqs, (21, 27), norm=True, custom=True)
    esm2, _ = A.load_ntv2()
    A.load_lora_state(esm2, st)
    out2 = A.encode_windows(esm2, vocab, seqs, (21, 27), norm=True, custom=True)
    for a, b in zip(out1, out2):
        assert torch.equal(a, b)
    assert not torch.equal(out1[0], same[0])
