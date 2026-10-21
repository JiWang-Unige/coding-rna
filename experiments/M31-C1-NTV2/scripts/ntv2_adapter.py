"""M31-C1-NTV2 adapter: frozen NT-v2-500M-multi-species features -> M28 C0 full-grammar role head.

NT-v2 (self-supervised only, no annotation) is a non-overlapping 6-mer model with rotary positions and at most 2,048 content
tokens (+1 CLS) = 12,288 bp.  The M28 head (`SharedContext`) already expects exactly one feature vector per six bases, so the
alignment is the same as for GENERanno; only the window changes:

    window  NT_WIDTH  = 12,288 bp  (2,048 tokens)        [M28/GENERanno: 24,576 bp, 4,096 tokens]
    stride  NT_STRIDE =  6,144 bp  (50 % overlap, nearest-center ownership via the frozen `context_segments`)

Tokens are built by hand (vocabulary lookup of every non-overlapping 6-mer chunk) because the stock tokenizer splits an
N-containing chunk into single characters, which would break the 1-token-per-6-bases alignment.  Any chunk containing a
non-ACGT character gets the vocabulary's single-character `N` id.  Short windows are right-padded with N (after orientation)
exactly like M28.  The CLS token is prepended (as in pre-training) and its output is dropped.

No frozen M28/Carbon-A file is edited; this module only imports.
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]                         # ~/coding-rna-reorg
OUTS = ROOT / 'outputs'                        # -> ~/coding-rna/outputs
M31 = OUTS / 'M31-C1-NTV2'
for p in (str(ROOT / 'pipeline/m28'), str(ROOT / 'pipeline/m28/scripts'), str(ROOT / 'experiments/M29-A1-RESIDUAL/scripts'),
          str(ROOT / 'experiments/M29-CARBON-A-GAP/scripts'), str(HERE), os.path.expanduser('~/envs/m29_extra')):
    if p not in sys.path:
        sys.path.insert(0, p)

import numpy as np
import torch

MODEL = 'InstaDeepAI/nucleotide-transformer-v2-500m-multi-species'
KMER = 6
NT_WIDTH = 12288
NT_STRIDE = 6144
NT_TOKENS = NT_WIDTH // KMER                   # 2048
N_LAYERS = 29
HIDDEN = 1024
CPUS = [2, 4, 6, 8, 10, 12, 14] + list(range(16, 25))


def load_ntv2(device='cuda', dtype=torch.float32, fast=True):
    """Frozen NT-v2 encoder (EsmModel) + the 6-mer vocabulary.  Weights stay FP32 and the forward runs under BF16 autocast: the
    remote modelling code builds its rotary tables in FP32, so a pure-BF16 module fails on dtype mixing (and autocast keeps the
    rotary / softmax / LayerNorm arithmetic in FP32, as in the historical M9-M23 NT-v2 runs)."""
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    tok = AutoTokenizer.from_pretrained(MODEL, trust_remote_code=True)
    esm = AutoModelForMaskedLM.from_pretrained(MODEL, trust_remote_code=True, torch_dtype=dtype).esm
    esm = esm.eval().requires_grad_(False).to(device)
    if fast:
        patch_sdpa(esm)
    return esm, NTVocab(tok.get_vocab())


def _sdpa_forward(self, hidden_states, attention_mask=None, head_mask=None, encoder_hidden_states=None,
                  encoder_attention_mask=None, past_key_value=None, output_attentions=False):
    """Same arithmetic as the remote `EsmSelfAttention.forward` for our use (self-attention, no padding, no dropout, rotary), but
    with fused scaled-dot-product attention.  The query scale commutes with the (linear) rotary rotation, so SDPA's default
    1/sqrt(d) scale reproduces it; rotary itself stays FP32 (module buffers) and q/k are cast to the value dtype afterwards."""
    import torch.nn.functional as F
    if attention_mask is not None and bool((attention_mask != 0).any()):
        raise ValueError('Padding is not supported by the fused NT-v2 attention path')
    q = self.transpose_for_scores(self.query(hidden_states))
    k = self.transpose_for_scores(self.key(hidden_states))
    v = self.transpose_for_scores(self.value(hidden_states))
    q, k = self.rotary_embeddings(q, k)
    out = F.scaled_dot_product_attention(q.to(v.dtype), k.to(v.dtype), v)
    out = out.permute(0, 2, 1, 3).reshape(*hidden_states.shape[:2], self.all_head_size)
    return (out,)


def patch_sdpa(esm):
    import types
    n = 0
    for m in esm.modules():
        if type(m).__name__ == 'EsmSelfAttention':
            m.forward = types.MethodType(_sdpa_forward, m)
            n += 1
    assert n == N_LAYERS, n
    return n


class NTVocab:
    def __init__(self, vocab):
        self.v = {k: i for k, i in vocab.items() if len(k) == KMER and set(k) <= set('ACGT')}
        assert len(self.v) == 4096
        self.cls, self.n = vocab['<cls>'], vocab['N']

    def ids(self, sequence, width=NT_WIDTH):
        """sequence (<= width bp, ACGTN) -> [1+width/6] ids (CLS + one id per 6 bp chunk; right padded with N)."""
        if not 0 < len(sequence) <= width:
            raise ValueError('Nonempty window of at most NT_WIDTH bp required')
        s = sequence + 'N' * (width - len(sequence))
        get, n = self.v.get, self.n
        return [self.cls] + [get(s[i:i + KMER], n) for i in range(0, width, KMER)]


class LoRALinear(torch.nn.Module):
    """y = W x + (alpha/r) B A x  (A kaiming, B zero => identity at init); W frozen, A/B FP32 trainable."""
    def __init__(self, base, r=8, alpha=16):
        super().__init__()
        import math
        self.base, self.scale = base, alpha / r
        self.A = torch.nn.Parameter(torch.empty(r, base.in_features))
        torch.nn.init.kaiming_uniform_(self.A, a=math.sqrt(5))
        self.B = torch.nn.Parameter(torch.zeros(base.out_features, r))

    def forward(self, x):
        F = torch.nn.functional
        return self.base(x) + F.linear(F.linear(x, self.A), self.B) * self.scale


def add_lora(esm, layer_ids, r=8, alpha=16):
    """LoRA on query and value of encoder layers `layer_ids` (0-based indices of esm.encoder.layer).  Returns the trainable params."""
    params = []
    for i in layer_ids:
        att = esm.encoder.layer[i].attention.self
        att.query, att.value = LoRALinear(att.query, r, alpha).to(att.query.weight.device), LoRALinear(att.value, r, alpha).to(att.value.weight.device)
        params += [att.query.A, att.query.B, att.value.A, att.value.B]
    return params


def lora_state(esm):
    return {f'{i}.{n}.{w}': getattr(getattr(esm.encoder.layer[i].attention.self, n), w).detach().cpu()
            for i in range(N_LAYERS) if isinstance(esm.encoder.layer[i].attention.self.query, LoRALinear) for n in ('query', 'value') for w in ('A', 'B')}


def load_lora_state(esm, state, r=8, alpha=16):
    add_lora(esm, sorted({int(k.split('.')[0]) for k in state}), r, alpha)
    for k, v in state.items():
        i, n, w = k.split('.')
        getattr(getattr(esm.encoder.layer[int(i)].attention.self, n), w).data.copy_(v)


def forward_layers(esm, ids, layers):
    """Run embeddings + the first max(layers) blocks only (hidden_states[k] = output of block k); layers 1..28 (29 would carry the final LayerNorm)."""
    assert max(layers) < N_LAYERS
    h = esm.embeddings(input_ids=ids, attention_mask=torch.ones_like(ids))
    out = {}
    for i in range(max(layers)):
        h = esm.encoder.layer[i](h)[0]
        if i + 1 in layers:
            out[i + 1] = h
    return [out[k] for k in layers]


def encode_windows(esm, vocab, sequences, layers=(N_LAYERS,), device='cuda', batch=4, norm=False, custom=False, grad=False):
    """Frozen (or LoRA-adapted when custom/grad) BF16-autocast NT-v2 features for oriented windows.  Returns a list of
    [2048, HIDDEN*len(layers)] BF16 tensors on `device`.  Layer ids are hidden_states indices (29 = last_hidden_state after the final LayerNorm);
    layers are concatenated in the given order, optionally each through a parameter-free LayerNorm."""
    out = []
    with (torch.enable_grad() if grad else torch.no_grad()):      # no_grad, not inference_mode: the rotary tables are cached across calls and must stay usable under autograd
        for i in range(0, len(sequences), batch):
            ids = torch.tensor([vocab.ids(s) for s in sequences[i:i + batch]], dtype=torch.long, device=device)
            with torch.autocast('cuda', dtype=torch.bfloat16):
                if custom or grad:
                    hs = forward_layers(esm, ids, layers)
                else:
                    allhs = esm(input_ids=ids, attention_mask=torch.ones_like(ids), output_hidden_states=True).hidden_states
                    hs = [allhs[k] for k in layers]
            if norm:      # parameter-free LayerNorm per layer so concatenated layers share one scale
                f = torch.cat([torch.nn.functional.layer_norm(h.float(), (HIDDEN,)) for h in hs], dim=-1)[:, 1:]
            else:
                f = torch.cat(hs, dim=-1)[:, 1:]
            if not bool(torch.isfinite(f).all()):
                raise ValueError('Finite NT-v2 features required')
            out.extend(x.to(torch.bfloat16).contiguous() for x in f)
    return out


def nt_geometry(a, b):
    """NT window grid inside a genomic context [a, b) (all windows lie in the context; coordinates local to `a`).
    Same dict shape that the frozen `owned_context_logits` expects."""
    from src.m28.contig_emission_lattice import context_segments
    import manifest as M
    L = b - a
    local = tuple((s, min(s + NT_WIDTH, L)) for s in M.starts(L, NT_WIDTH, NT_STRIDE))
    return dict(genomic_context=[a, b], genomic_windows=tuple((a + x, a + y) for x, y in local), local_windows=local,
                base_owned_segments=context_segments(L, local))


def new_head(feature_dim=HIDDEN, seed=0):
    """C0 head (frozen-pipeline class) with the default PyTorch initialisation (the M28 first fit used a fresh seed-0 C0)."""
    from src.m28.core import M28Core
    torch.manual_seed(seed)
    return M28Core('C0', feature_dim, 128)
