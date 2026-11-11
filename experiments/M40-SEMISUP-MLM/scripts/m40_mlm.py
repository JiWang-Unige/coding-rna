"""M40 masked-language-model pieces (pure numpy / torch; no data paths, unit-tested on CPU).

PlantCAD2-Small ships its MLM head (checked 2026-11-11 on ont, protocol 32 section 2): `lm_head.lm_head.weight` (8 x 768) is element-wise identical to the
input embedding `backbone.embeddings.word_embeddings.embedding.weight` (tied), and the head is the Caduceus `RCPSLMHead`:
    logits(h) = h[:, :768] @ E^T  +  flip(h[:, 768:], channels) @ E[complement]^T          (h = last_hidden_state, 1,536 = forward | RC view)
`AutoModel` (what M32/M34/M35 load) returns exactly that last_hidden_state, so the pre-trained MLM head is reproduced from the frozen embedding of the
loaded backbone; no extra weights and no new head.  Vocabulary: [PAD]=0 [MASK]=1 [UNK]=2 a=3 c=4 g=5 t=6 (7 unused, pad to 8).

Masking of one draw (genomic coordinates of the draw's encoder context [ca, cb)):
  * selection is LABEL-BLIND: every ACGT base is selected with probability `rate` (default 5 %), from a per-draw RNG; the CRF targets are never consulted,
    so the position of [MASK] tokens carries no information about which loci are unknown;
  * corruption 80 / 10 / 10 ([MASK] / random base / unchanged), applied once in context coordinates, hence identically in every overlapping encoder window;
  * the loss is counted at selected bases that are in SCOPE (no CRF label: unknown core positions + context flanks outside the core), once per base, in the
    window that owns the base (nearest-centre ownership of the frozen geometry); weight 1, or `repeat_weight` for soft-masked (repeat) bases.
"""
import numpy as np
import torch
import torch.nn.functional as F

PAD, MASK, UNK = 0, 1, 2
NUC = (3, 4, 5, 6)
COMPLEMENT = (0, 1, 2, 6, 5, 4, 3, 7)            # config.json complement_map of PlantCAD2-Small
LUT = np.full(256, UNK, dtype=np.int64)          # identical to pc2_adapter.LUT (asserted in the tests)
for _c, _i in zip('acgt', NUC):
    LUT[ord(_c)] = _i
    LUT[ord(_c.upper())] = _i


def ids_of(s):
    return LUT[np.frombuffer(s.encode('ascii'), dtype=np.uint8)]


# ---- schedule of the auxiliary weight ---------------------------------------------------------------------------------------------------------
def lam_at(step, total, lam_max, start_frac, ramp_frac):
    """lambda for 0-based optimiser step `step`: 0 for steps < s0 = round(start_frac*total); then lam_max * min(1, (step - s0 + 1) / r) with
    r = round(ramp_frac*total) (r = 0: lam_max at once).  Protocol: s0 = 307, r = 307 -> first non-zero step 307 (lam_max/307), lam_max from step 613 on."""
    s0, r = int(round(start_frac * total)), int(round(ramp_frac * total))
    if step < s0:
        return 0.
    if r <= 0:
        return float(lam_max)
    return float(lam_max) * min(1., (step - s0 + 1) / r)


# ---- per-draw mask construction (CPU, prefetch thread) -----------------------------------------------------------------------------------------
def scope_genomic(n_ctx, core_lo, core_hi, unknown_oriented, strand):
    """Boolean [n_ctx] over the genomic context: True = no CRF label.  core_lo/core_hi: core in context-local GENOMIC coordinates; unknown_oriented: [core_hi-core_lo]
    booleans in ORIENTED order (index 0 = oriented start of the core) -- the trainer's targets are oriented."""
    sc = np.ones(n_ctx, dtype=bool)
    u = np.asarray(unknown_oriented, dtype=bool)
    assert len(u) == core_hi - core_lo
    sc[core_lo:core_hi] = u if strand == '+' else u[::-1]
    return sc


def draw_masks(ctx_ids, scope, repeat, local_windows, owned, rng, rate, mask_frac=0.8, rand_frac=0.1, repeat_weight=1e-4):
    """ctx_ids [n] clean token ids of the genomic context; scope / repeat [n] bool; local_windows [(a, b)] window bounds in context coordinates (genomic order);
    owned [(left, right)] the base range each window owns (context coordinates).  -> dict(
        ids=[int64 [b-a] corrupted window ids], pos=[int64 local loss positions], tgt=[int64 true ids], w=[float32 weights],
        n_selected, n_scope_selected, wsum)."""
    n = len(ctx_ids)
    nuc = ctx_ids >= NUC[0]
    sel = np.flatnonzero((rng.random(n) < rate) & nuc)
    kind = rng.random(len(sel))
    rand_base = rng.integers(NUC[0], NUC[-1] + 1, size=len(sel))
    cor = ctx_ids.copy()
    m = kind < mask_frac
    cor[sel[m]] = MASK
    r = (kind >= mask_frac) & (kind < mask_frac + rand_frac)
    cor[sel[r]] = rand_base[r]
    loss_pos = sel[scope[sel]]
    weight = np.where(repeat[loss_pos], repeat_weight, 1.).astype(np.float32)
    out = dict(ids=[], pos=[], tgt=[], w=[], n_selected=int(len(sel)), n_scope_selected=int(len(loss_pos)))
    for (a, b), (lo, hi) in zip(local_windows, owned):
        out['ids'].append(cor[a:b])
        k = (loss_pos >= lo) & (loss_pos < hi)
        p = loss_pos[k]
        out['pos'].append((p - a).astype(np.int64))
        out['tgt'].append(ctx_ids[p].astype(np.int64))
        out['w'].append(weight[k])
    out['wsum'] = float(sum(x.sum() for x in out['w']))
    return out


def owned_ranges(geometry):
    """[(left, right)] per window from the frozen geometry's base_owned_segments (context-local coordinates, genomic order)."""
    return [(left, right) for _, left, right in geometry['base_owned_segments']]


# ---- loss (GPU / CPU) -----------------------------------------------------------------------------------------------------------------------------
def lm_params(enc):
    """(E [V, D] frozen tied embedding, complement index [V]) of a loaded Caduceus backbone (AutoModel)."""
    we = enc.backbone.embeddings.word_embeddings
    return we.embedding.weight, we.complement_map


def lm_logits(h, E, comp):
    """Caduceus RCPSLMHead on rows h [..., 2D] -> [..., V] (float32)."""
    D = E.shape[1]
    E = E.float()
    h = h.float()
    return F.linear(h[..., :D], E) + F.linear(torch.flip(h[..., D:], dims=[-1]), E[comp])


def window_mlm(h, pos, tgt, w, E, comp):
    """Weighted CE sum, weight sum, correct count of one window.  h [len, 2D] per-base features (genomic orientation); pos/tgt/w torch tensors on h's device."""
    if len(pos) == 0:
        z = h.new_zeros((), dtype=torch.float32)
        return z, 0., 0
    lg = lm_logits(h[pos], E, comp)
    ce = F.cross_entropy(lg, tgt, reduction='none')
    return (ce * w).sum(), float(w.sum()), int((lg.argmax(-1) == tgt).sum())
