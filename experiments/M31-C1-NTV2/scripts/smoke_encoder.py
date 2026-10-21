"""Encoder smoke: load time, speed/memory of 12,288-bp windows (BF16), hidden_states layout, BF16-vs-FP32 agreement, N handling."""
import json, time
import numpy as np, torch
import ntv2_adapter as A

t0 = time.time()
esm, vocab = A.load_ntv2(fast=False)
ref_ids = torch.tensor([A.NTVocab.ids(vocab, 'ACGT' * 3072)], device='cuda')
with torch.inference_mode(), torch.autocast('cuda', dtype=torch.bfloat16):
    ref = esm(input_ids=ref_ids, attention_mask=torch.ones_like(ref_ids), output_hidden_states=True).hidden_states
A.patch_sdpa(esm)
with torch.inference_mode(), torch.autocast('cuda', dtype=torch.bfloat16):
    new = esm(input_ids=ref_ids, attention_mask=torch.ones_like(ref_ids), output_hidden_states=True).hidden_states
for k in (10, 20, 29):
    a, b = ref[k][:, 1:].float().flatten(0, 1), new[k][:, 1:].float().flatten(0, 1)
    print('eager-vs-sdpa layer', k, 'cos', round(float(torch.nn.functional.cosine_similarity(a, b, dim=-1).mean()), 5), 'relerr', round(float((a - b).norm() / a.norm()), 4), flush=True)
print('load_s', round(time.time() - t0, 1), 'gpu_alloc_GB', torch.cuda.memory_allocated() / 2**30, flush=True)
rng = np.random.default_rng(0)
seqs = [''.join(rng.choice(list('ACGT'), A.NT_WIDTH)) for _ in range(8)]
seqs[1] = seqs[1][:5000] + 'N' * 300 + seqs[1][5300:]
seqs[2] = seqs[2][:9001]                      # short window (right-padded with N)
ids = torch.tensor([vocab.ids(seqs[1])], device='cuda')
print('ids shape', tuple(ids.shape), 'N-chunks', int((ids == vocab.n).sum()), flush=True)
with torch.inference_mode():
    with torch.autocast('cuda', dtype=torch.bfloat16):
        o = esm(input_ids=ids, attention_mask=torch.ones_like(ids), output_hidden_states=True)
    hs = o.hidden_states
    print('n_hidden', len(hs), 'last==hs[-1]', bool(torch.equal(o.last_hidden_state, hs[-1])), tuple(hs[-1].shape))
    print('hs std per layer', [round(float(h.float().std()), 2) for h in (hs[0], hs[10], hs[20], hs[28], hs[29])])
for layers in ((29,),):
    for batch in (1, 4):
        torch.cuda.synchronize(); torch.cuda.reset_peak_memory_stats(); t = time.perf_counter()
        for _ in range(3):
            f = A.encode_windows(esm, vocab, seqs, layers, batch=batch)
        torch.cuda.synchronize(); dt = (time.perf_counter() - t) / 3 / len(seqs)
        print(json.dumps(dict(layers=layers, batch=batch, s_per_window=round(dt, 4), kbp_per_s=round(A.NT_WIDTH / 1000 / dt, 1),
                              peak_GB=round(torch.cuda.max_memory_allocated() / 2**30, 2), out_shape=tuple(f[0].shape))), flush=True)
# BF16 vs FP32
from transformers import AutoModelForMaskedLM
m32 = esm  # weights are FP32 already; the reference forward below runs without autocast
with torch.inference_mode():
    ids = torch.tensor([vocab.ids(seqs[0]), vocab.ids(seqs[3])], device='cuda')
    h32 = m32(input_ids=ids, attention_mask=torch.ones_like(ids), output_hidden_states=True).hidden_states
    with torch.autocast('cuda', dtype=torch.bfloat16):
        h16 = esm(input_ids=ids, attention_mask=torch.ones_like(ids), output_hidden_states=True).hidden_states
    for k in (10, 20, 29):
        a, b = h32[k][:, 1:].float().flatten(0, 1), h16[k][:, 1:].float().flatten(0, 1)
        print('layer', k, 'cos', round(float(torch.nn.functional.cosine_similarity(a, b, dim=-1).mean()), 5), 'relerr', round(float((a - b).norm() / a.norm()), 4))
