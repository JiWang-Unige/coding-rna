"""M35-PLANTCAD2-LORA: rank-r LoRA on the Mamba2 `in_proj` / `out_proj` of PlantCAD2-Small (Caduceus-style bidirectional Mamba2 with RCPS).

Design (see docs/review_20261016/21_M35_plantcad2_lora.md):
* PlantCAD2 has no query/value projections.  Each block's mixer is RCPSWrapper(BiMambaWrapper(mamba_fwd, mamba_rev)); both Mamba2 modules own the two big
  linear maps `in_proj` (768 -> 3,224) and `out_proj` (1,536 -> 768).  They are the analogue of M31's q/v projections.
* Weight tying: the released checkpoint stores `mamba_rev.{in,out}_proj` as separate copies of `mamba_fwd.{in,out}_proj` (tied in the original model; the HF
  loader breaks the tie).  Here the frozen base weights are re-tied (after asserting that the two copies are identical) AND `mamba_fwd` and `mamba_rev` share ONE
  LoRA pair per projection, so the effective weights stay tied as in pre-training.
* Reverse-complement equivariance: RCPSWrapper applies the SAME BiMambaWrapper to the forward-strand half and to the flipped RC half of the channels, so any
  parameter inside the mixer (including these shared LoRA matrices) acts identically on both views and f(RC(x)) = flip(f(x), position and channel) stays exact.
* `out_proj` is fused inside `mamba_split_conv1d_scan_combined`, which bypasses the module; the patched forward calls the fused kernel WITHOUT its out-projection
  and applies `out_proj` (base + LoRA) explicitly, so out_proj LoRA is a real additive adapter (no bf16 rounding of a merged W + BA).
* Training memory: blocks of the adapted range are run through torch.utils.checkpoint (non-reentrant) when autograd is on; lower, frozen blocks build no graph.
"""
import math
import types

import torch
from torch import nn
from torch.nn import functional as F

N_LAYERS = 24
LORA_LAYERS = tuple(range(12, 24))            # upper half of the 24 blocks (M31 adapted blocks 13..27 of 29, i.e. the upper ~half)
RANK = 8
ALPHA = 16
NAMES = ('in_proj', 'out_proj')

scan_fn = None                                 # injectable fused scan (CPU tests inject the pure-torch reference)


def _scan():
    if scan_fn is not None:
        return scan_fn
    from mamba_ssm.ops.triton.ssd_combined import mamba_split_conv1d_scan_combined
    return mamba_split_conv1d_scan_combined


class LoRALinear(nn.Module):
    """y = W x + (alpha/r) B (A x).  W frozen (the wrapped Linear), A (r x in) and B (out x r) trainable FP32 Parameters that may be SHARED between wrappers."""

    def __init__(self, base, A, B, scale):
        super().__init__()
        self.base, self.A, self.B, self.scale = base, A, B, scale

    @property
    def weight(self):              # Mamba2 reads in_proj.weight / out_proj.weight for dtype/device bookkeeping
        return self.base.weight

    @property
    def bias(self):
        return self.base.bias

    def forward(self, x):
        return self.base(x) + F.linear(F.linear(x, self.A), self.B) * self.scale


def _mamba2_forward(self, u, seqlen=None, seq_idx=None, cu_seqlens=None, inference_params=None):
    """Mamba2.forward (mamba-ssm 2.3.0, training/prefill path) with the output projection applied OUTSIDE the fused kernel."""
    from einops import rearrange
    assert seqlen is None and seq_idx is None and cu_seqlens is None and inference_params is None
    zxbcdt = self.in_proj(u)
    A = -torch.exp(self.A_log.float())
    y = _scan()(zxbcdt, rearrange(self.conv1d.weight, 'd 1 w -> d w'), self.conv1d.bias, self.dt_bias, A,
                D=rearrange(self.D, '(h p) -> h p', p=self.headdim) if self.D_has_hdim else self.D, chunk_size=self.chunk_size,
                activation=self.activation, rmsnorm_weight=self.norm.weight if self.rmsnorm else None,
                rmsnorm_eps=self.norm.eps if self.rmsnorm else 1e-6, outproj_weight=None, outproj_bias=None,
                headdim=None if self.D_has_hdim else self.headdim, ngroups=self.ngroups, norm_before_gate=self.norm_before_gate)
    return self.out_proj(y)


def mixers(model):
    """[(block_index, BiMambaWrapper)] of a Caduceus backbone (`model` = AutoModel Caduceus or a model with .backbone)."""
    return [(i, layer.mixer.submodule) for i, layer in enumerate(model.backbone.layers)]


def tie_check_and_retie(model):
    """Assert that the HF-loaded `mamba_rev` in/out projections equal `mamba_fwd`'s (the release stores the tied weights twice), then re-tie them (one
    Parameter object).  Returns the number of tied tensors."""
    n = 0
    for i, bi in mixers(model):
        for name in NAMES:
            f, r = getattr(bi.mamba_fwd, name), getattr(bi.mamba_rev, name)
            if f.weight is r.weight:
                n += 1
                continue
            assert torch.equal(f.weight, r.weight), f'block {i} {name}: mamba_rev differs from mamba_fwd (not a tied checkpoint)'
            assert f.bias is None and r.bias is None
            r.weight = f.weight
            n += 1
    return n


def install_lora(model, layers=LORA_LAYERS, r=RANK, alpha=ALPHA, seed=0):
    """Re-tie the frozen projections, then wrap in_proj/out_proj of `mamba_fwd` and `mamba_rev` of the given blocks with ONE shared LoRA pair per projection and
    patch their forward.  A ~ U(-1/sqrt(in), 1/sqrt(in)) (kaiming_uniform a=sqrt5, as M31), B = 0 (identity at init).  Returns the unique trainable
    parameters in the fixed order [block][in.A, in.B, out.A, out.B]."""
    tie_check_and_retie(model)
    params = []
    g = torch.Generator().manual_seed(seed)
    for i, bi in mixers(model):
        if i not in layers:
            continue
        for name in NAMES:
            base_f = getattr(bi.mamba_fwd, name)
            assert isinstance(base_f, nn.Linear) and base_f.bias is None
            w = base_f.weight
            A = nn.Parameter(torch.empty(r, w.shape[1], device=w.device, dtype=torch.float32))
            bound = 1 / math.sqrt(w.shape[1])
            with torch.no_grad():
                A.copy_((torch.rand(A.shape, generator=g) * 2 - 1).to(A.device) * bound)
            B = nn.Parameter(torch.zeros(w.shape[0], r, device=w.device, dtype=torch.float32))
            for m in (bi.mamba_fwd, bi.mamba_rev):
                setattr(m, name, LoRALinear(getattr(m, name), A, B, alpha / r))
            params += [A, B]
        for m in (bi.mamba_fwd, bi.mamba_rev):
            m.forward = types.MethodType(_mamba2_forward, m)
    return params


def lora_layers(model):
    return [i for i, bi in mixers(model) if isinstance(bi.mamba_fwd.in_proj, LoRALinear)]


def lora_params(model):
    out = []
    for i, bi in mixers(model):
        if isinstance(bi.mamba_fwd.in_proj, LoRALinear):
            for name in NAMES:
                m = getattr(bi.mamba_fwd, name)
                out += [m.A, m.B]
    return out


def lora_state(model):
    """{ '<block>.<proj>.<A|B>': cpu tensor } from the mamba_fwd copy (mamba_rev shares the same Parameters)."""
    st = {}
    for i, bi in mixers(model):
        if isinstance(bi.mamba_fwd.in_proj, LoRALinear):
            for name in NAMES:
                m = getattr(bi.mamba_fwd, name)
                st[f'{i}.{name}.A'], st[f'{i}.{name}.B'] = m.A.detach().cpu().clone(), m.B.detach().cpu().clone()
    return st


def load_lora_values(model, state):
    """Copy values into an already installed LoRA (same layers)."""
    with torch.no_grad():
        for i, bi in mixers(model):
            if isinstance(bi.mamba_fwd.in_proj, LoRALinear):
                for name in NAMES:
                    m = getattr(bi.mamba_fwd, name)
                    m.A.copy_(state[f'{i}.{name}.A']); m.B.copy_(state[f'{i}.{name}.B'])


def load_lora_state(model, state, r=RANK, alpha=ALPHA):
    """Install LoRA on the blocks present in `state` and load its values."""
    layers = sorted({int(k.split('.')[0]) for k in state})
    params = install_lora(model, tuple(layers), r, alpha)
    load_lora_values(model, state)
    return params


def enable_layer_checkpoint(model, layers=LORA_LAYERS):
    """Run the given blocks through torch.utils.checkpoint (non-reentrant) whenever autograd is enabled; no effect under no_grad / inference."""
    from torch.utils.checkpoint import checkpoint
    for i, layer in enumerate(model.backbone.layers):
        if i in layers and not getattr(layer, '_m35_ckpt', False):
            orig = layer.forward

            def fwd(hidden_states, residual=None, inference_params=None, _o=orig):
                if torch.is_grad_enabled():
                    return checkpoint(lambda h, r: _o(h, r, inference_params=inference_params), hidden_states, residual,
                                      use_reentrant=False, preserve_rng_state=False)
                return _o(hidden_states, residual, inference_params=inference_params)
            layer.forward = fwd
            layer._m35_ckpt = True


def encode_grad(model, sequences, device='cuda'):
    """Per-base genomic-orientation features [n, 1536] (BF16, autograd to the LoRA parameters) of ACGTN windows of EQUAL length (one batch)."""
    import pc2_adapter as A
    n = len(sequences[0])
    assert all(len(s) == n for s in sequences) and 0 < n <= A.PC_WIDTH
    ids = torch.stack([A.token_ids(s, device) for s in sequences])
    with torch.enable_grad(), torch.autocast('cuda', dtype=torch.bfloat16):
        h = model(ids).last_hidden_state
    if h.shape != (len(sequences), n, A.HIDDEN):
        raise ValueError(f'unexpected feature shape {tuple(h.shape)}')
    if not bool(torch.isfinite(h).all()):
        raise ValueError('Finite PlantCAD2 features required')
    return [x.to(torch.bfloat16) for x in h]


# ---- M34 backend contract (added 2026-10-22 after the M35 decision): add_lora(model, cfg) -> params, lora_state(model) -> {k: tensor},
# load_lora_state(model, state, cfg) and a JSON-able CFG.  Pure wrappers over install_lora / load_lora_values; the M35 (r, alpha) call form
# of load_lora_state keeps working (dispatch on the type of the third argument).
CFG = dict(layers=list(LORA_LAYERS), r=RANK, alpha=ALPHA, seed=0)


def add_lora(model, cfg=None):
    c = dict(CFG, **(cfg or {}))
    return install_lora(model, tuple(c['layers']), c['r'], c['alpha'], c.get('seed', 0))


_load_lora_state_r_alpha = load_lora_state


def load_lora_state(model, state, cfg_or_r=RANK, alpha=ALPHA):
    if isinstance(cfg_or_r, dict):
        if lora_layers(model):                       # already adapted (resume / after add_lora): copy values only
            load_lora_values(model, state)
            return lora_params(model)
        c = dict(CFG, **cfg_or_r)
        params = install_lora(model, tuple(c['layers']), c['r'], c['alpha'], c.get('seed', 0))
        load_lora_values(model, state)
        return params
    return _load_lora_state_r_alpha(model, state, cfg_or_r, alpha)
