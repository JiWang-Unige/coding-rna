"""M32-C1-PLANTCAD2 adapter: frozen PlantCAD2-Small (self-supervised MLM, no label supervision) -> M28 C0 full-grammar role head.

PlantCAD2 (Caduceus-style bidirectional Mamba2 with RCPS reverse-complement parameter sharing) is single-nucleotide and was pre-trained on
8,192-bp windows only; Small degrades above ~8 kb (report 16 section 3.4).  Window design here:

    window  PC_WIDTH  = 8,192 bp      stride PC_STRIDE = 6,144 bp (25 % overlap; nearest-centre ownership via the frozen `context_segments`)

One forward over a GENOMIC window gives the features of both strands: the model is exactly RC-equivariant, f(RC(x)) = flip(f(x), dims=(position, channel))
(report 16 section 3.2, re-checked in `smoke_encoder.py`).  So the minus-strand (oriented) features are `feat.flip(0).flip(-1)`.

The frozen M28 C0 head (`SharedContext`) takes ONE feature vector per six bases (GRU over tokens + DNA conv for the position inside the token).
To keep the head, loss, targets, sampler and decoder identical to the NT-v2 pilot (M31) we aggregate the per-base 1,536-d vectors into six-base
tokens by the mean over the (<=6) valid bases of each aligned group, counted from the start of the ORIENTED window.  The head width is 8,196
(= 1,366 tokens x 6) so that an 8,192-bp window has a ragged 2-base last token (valid length is passed to the head).

No frozen M28 / Carbon-A file is edited; this module only imports.
"""
import os
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]                         # ~/coding-rna-reorg
OUTS = ROOT / 'outputs'                        # -> ~/coding-rna/outputs
M32 = OUTS / 'M32-C1-PLANTCAD2'
EXTRA = os.path.expanduser('~/envs/m32_extra')   # numba/scipy/pandas/pytest built for the py3.11 plantcad2 env (the env itself is untouched)
for p in (str(ROOT / 'pipeline/m28'), str(ROOT / 'pipeline/m28/scripts'), str(ROOT / 'experiments/M29-A1-RESIDUAL/scripts'),
          str(ROOT / 'experiments/M29-CARBON-A-GAP/scripts'), EXTRA, str(HERE)):    # HERE last => first on the path (shim `a1_common`)
    if p in sys.path:
        sys.path.remove(p)
    sys.path.insert(0, p)                       # (re-)insert so that HERE really ends up first even when it was already on the path

import numpy as np
import torch

MODEL_DIR = os.path.expanduser('~/plantcad2_work/models/PlantCAD2-Small-l24-d0768')
MODEL = 'kuleshov-group/PlantCAD2-Small-l24-d0768'
KMER = 6
PC_WIDTH = 8192
PC_STRIDE = 6144
HEAD_W = 8196                                   # 1,366 tokens x 6 bases
HEAD_TOKENS = HEAD_W // KMER
D_MODEL = 768
HIDDEN = 2 * D_MODEL                            # RCPS: forward view | reverse-complement view
N_LAYERS = 24
LUT = np.full(256, 2, dtype=np.int64)           # N / anything else -> [UNK] = 2 (tokenizer ids a,c,g,t = 3..6)
for _c, _i in zip('acgt', (3, 4, 5, 6)):
    LUT[ord(_c)] = _i
    LUT[ord(_c.upper())] = _i


def load_pc2(device='cuda'):
    """Frozen PlantCAD2-Small.  Weights stay FP32 (`torch_dtype=bfloat16` is silently ignored by this remote code); forward under BF16 autocast.
    Note: the checkpoint stores the tied reverse projections twice; irrelevant for frozen inference."""
    from transformers import AutoModel
    m = AutoModel.from_pretrained(MODEL_DIR, trust_remote_code=True)
    return m.eval().requires_grad_(False).to(device)


def token_ids(sequence, device='cuda'):
    return torch.from_numpy(LUT[np.frombuffer(sequence.encode('ascii'), dtype=np.uint8)]).to(device)


@torch.inference_mode()
def encode_windows(model, sequences, layer=-1, batch=4, device='cuda'):
    """Genomic-orientation per-base features of ACGTN windows (any length <= PC_WIDTH; no padding is added): list of [len, 1536] BF16 tensors
    on `device`.  layer = -1 is last_hidden_state (final norm); otherwise a hidden_states index (0 = input, n_layer = last)."""
    out = [None] * len(sequences)
    by_len = {}
    for i, s in enumerate(sequences):
        if not 0 < len(s) <= PC_WIDTH:
            raise ValueError('Nonempty window of at most PC_WIDTH bp required')
        by_len.setdefault(len(s), []).append(i)
    for n, idx in by_len.items():
        for j in range(0, len(idx), batch):
            part = idx[j:j + batch]
            ids = torch.stack([token_ids(sequences[i], device) for i in part])
            with torch.autocast('cuda', dtype=torch.bfloat16):
                if layer == -1:
                    h = model(ids).last_hidden_state
                else:
                    h = model(ids, output_hidden_states=True).hidden_states[layer]
            if h.shape != (len(part), n, HIDDEN):
                raise ValueError(f'unexpected feature shape {tuple(h.shape)}')
            if not bool(torch.isfinite(h).all()):
                raise ValueError('Finite PlantCAD2 features required')
            for k, i in enumerate(part):
                out[i] = h[k].to(torch.bfloat16).contiguous()
    return out


def orient(feat, strand):
    """Genomic per-base features [n,1536] -> features of the ORIENTED window (exact RC equivariance)."""
    return feat if strand == '+' else feat.flip(0).flip(-1)


def pool_tokens(feat):
    """Oriented per-base features [n<=HEAD_W, C] -> [HEAD_TOKENS, C] BF16: mean over the valid bases of every aligned six-base group, zero padded."""
    n, c = feat.shape
    if n > HEAD_W:
        raise ValueError('window longer than the head width')
    g = torch.zeros(HEAD_W, c, dtype=torch.float32, device=feat.device)
    g[:n] = feat.float()
    cnt = (n - KMER * torch.arange(HEAD_TOKENS, device=feat.device)).clamp(0, KMER).clamp(min=1).to(torch.float32)
    return (g.view(HEAD_TOKENS, KMER, c).sum(1) / cnt[:, None]).to(torch.bfloat16)


def oriented_tokens(feats_genomic, strand):
    return [pool_tokens(orient(f, strand)) for f in feats_genomic]


def pc_geometry(a, b):
    """Window grid inside a genomic context [a, b) (coordinates local to `a`), same dict shape the frozen `owned_context_logits` expects."""
    from src.m28.contig_emission_lattice import context_segments
    import manifest as M
    L = b - a
    local = tuple((s, min(s + PC_WIDTH, L)) for s in M.starts(L, PC_WIDTH, PC_STRIDE))
    return dict(genomic_context=[a, b], genomic_windows=tuple((a + x, a + y) for x, y in local), local_windows=local,
                base_owned_segments=context_segments(L, local))


def new_head(feature_dim=HIDDEN, seed=0):
    """C0 head (frozen-pipeline class) with the default PyTorch initialisation (same as M31 / the M28 first fit: fresh seed-0 C0)."""
    from src.m28.core import M28Core
    torch.manual_seed(seed)
    return M28Core('C0', feature_dim, 128)
