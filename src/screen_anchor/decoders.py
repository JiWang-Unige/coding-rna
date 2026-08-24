"""Structured decoders for TA-DECODER-M3 (Track A focused batch), on top of the per-base
emission backbone (TiberiusLike -> (B,W,C) emissions). Three mechanisms:

  - LinearChainCRF   : learned KxK transition matrix + forward-NLL + Viterbi (label-sequence
                       coherence; replaces per-base softmax argmax).
  - SemiCRF          : bounded semi-Markov CRF (segment-level scoring up to max_seg_len) +
                       segment DP; models segments not per-base.
  - constrained_decode: deterministic post-processing of per-base argmax (drop sub-min CDS
                       fragments, fill small gaps) to reduce gene-count fragmentation.

All operate on 3 classes {intergenic(0), CDS(1), gene-body-noncoding(2)} with a padding mask.
"""
import numpy as np
import torch
import torch.nn as nn

from .data import NUM_CLASSES, CLASS_CDS, CLASS_INTERGENIC, CLASS_GENEBODY_NC

NEG = -1e4  # finite "−inf" for masked transitions (avoids NaN in logsumexp)


class LinearChainCRF(nn.Module):
    """Standard linear-chain CRF. emissions (B,W,C), tags (B,W) long with pad=-100,
    mask (B,W) bool. trans[i,j] = score of transitioning tag i -> j."""

    def __init__(self, num_tags=NUM_CLASSES):
        super().__init__()
        self.num_tags = num_tags
        self.trans = nn.Parameter(torch.randn(num_tags, num_tags) * 0.1)
        self.start = nn.Parameter(torch.randn(num_tags) * 0.1)
        self.end = nn.Parameter(torch.randn(num_tags) * 0.1)

    def _forward_alg(self, emissions, mask):
        B, W, C = emissions.shape
        alpha = self.start.unsqueeze(0) + emissions[:, 0]          # (B,C)
        for t in range(1, W):
            # (B,C_prev,1) + (1,C_prev,C_next) + (B,1,C_next)
            score = alpha.unsqueeze(2) + self.trans.unsqueeze(0) + emissions[:, t].unsqueeze(1)
            new_alpha = torch.logsumexp(score, dim=1)             # (B,C)
            m = mask[:, t].unsqueeze(1)
            alpha = torch.where(m, new_alpha, alpha)
        alpha = alpha + self.end.unsqueeze(0)
        return torch.logsumexp(alpha, dim=1)                       # (B,)

    def _gold_score(self, emissions, tags, mask):
        B, W, C = emissions.shape
        tags_safe = tags.clone()
        tags_safe[~mask] = 0
        idx0 = tags_safe[:, 0]
        score = self.start[idx0] + emissions[:, 0].gather(1, idx0.unsqueeze(1)).squeeze(1)
        for t in range(1, W):
            emit = emissions[:, t].gather(1, tags_safe[:, t:t + 1]).squeeze(1)
            tr = self.trans[tags_safe[:, t - 1], tags_safe[:, t]]
            score = score + torch.where(mask[:, t], emit + tr, torch.zeros_like(emit))
        last = mask.sum(1).clamp(min=1) - 1                        # (B,) last valid idx
        last_tags = tags_safe.gather(1, last.unsqueeze(1)).squeeze(1)
        score = score + self.end[last_tags]
        return score

    def nll(self, emissions, tags, mask):
        # PER-TOKEN normalized so the CRF NLL matches a per-base CE scale (the summed-over-W NLL
        # ~1000s otherwise swamps the aux class-weighted CE and emissions never learn CDS).
        ntok = mask.sum(1).clamp(min=1).float()
        return ((self._forward_alg(emissions, mask) - self._gold_score(emissions, tags, mask)) / ntok).mean()

    @torch.no_grad()
    def viterbi(self, emissions, mask):
        B, W, C = emissions.shape
        score = self.start.unsqueeze(0) + emissions[:, 0]         # (B,C)
        backptr = []
        for t in range(1, W):
            s = score.unsqueeze(2) + self.trans.unsqueeze(0)     # (B,C_prev,C_next)
            best, bp = s.max(dim=1)                                # (B,C_next)
            cand = best + emissions[:, t]
            m = mask[:, t].unsqueeze(1)
            score = torch.where(m, cand, score)
            backptr.append(torch.where(m, bp, torch.arange(C, device=emissions.device).unsqueeze(0)))
        score = score + self.end.unsqueeze(0)
        best_last = score.argmax(dim=1)                           # (B,)
        # backtrack
        paths = torch.zeros(B, W, dtype=torch.long, device=emissions.device)
        cur = best_last
        paths[:, W - 1] = cur
        for t in range(W - 2, -1, -1):
            cur = backptr[t].gather(1, cur.unsqueeze(1)).squeeze(1)
            paths[:, t] = cur
        return paths                                              # (B,W)


def _logbmm(a, b):
    """Log-space (logsumexp,+) matrix product over the last two dims, batched on leading dims.
    a (...,C,K), b (...,K,C2) -> out (...,C,C2): out[...,i,j]=logsumexp_k(a[...,i,k]+b[...,k,j])."""
    return torch.logsumexp(a.unsqueeze(-1) + b.unsqueeze(-3), dim=-2)


class LinearChainCRFVec(nn.Module):
    """Vectorized linear-chain CRF — drop-in for LinearChainCRF (same nll/viterbi interface).
    PARTITION via a log-space ASSOCIATIVE SCAN (O(log W) tree of log-matmuls) instead of an
    O(W) python loop -> tractable backward at W=2048. GOLD score fully vectorized (gather +
    masked sum, no loop). VITERBI kept sequential (predict-only, runs once). Padded steps use
    an identity transition (diag 0, off -inf) so masked positions are no-ops in the scan."""

    def __init__(self, num_tags=NUM_CLASSES):
        super().__init__()
        self.num_tags = num_tags
        self.trans = nn.Parameter(torch.randn(num_tags, num_tags) * 0.1)
        self.start = nn.Parameter(torch.randn(num_tags) * 0.1)
        self.end = nn.Parameter(torch.randn(num_tags) * 0.1)

    def _forward_alg(self, emissions, mask):
        B, W, C = emissions.shape
        dev = emissions.device
        alpha0 = self.start.unsqueeze(0) + emissions[:, 0]                 # (B,C)
        if W == 1:
            return torch.logsumexp(alpha0 + self.end.unsqueeze(0), dim=1)
        # per-step transition-with-emission matrices T_t[b,i,j] = trans[i,j] + emit_t[b,j], t=1..W-1
        T = self.trans.view(1, 1, C, C) + emissions[:, 1:].unsqueeze(2)    # (B,W-1,C,C)
        # identity (no-op) matrix for padded steps
        ident = torch.full((C, C), NEG, device=dev)
        ident[range(C), range(C)] = 0.0
        m = mask[:, 1:].view(B, W - 1, 1, 1)                               # valid step?
        T = torch.where(m, T, ident.view(1, 1, C, C))
        # associative reduction over the step dim (tree of log-matmuls)
        cur = T
        while cur.shape[1] > 1:
            n = cur.shape[1]
            if n % 2 == 1:
                last = cur[:, -1:]
                cur = cur[:, :-1]
            else:
                last = None
            a = cur[:, 0::2]
            b = cur[:, 1::2]
            cur = _logbmm(a, b)                                            # (B, n//2, C, C)
            if last is not None:
                cur = torch.cat([cur, last], dim=1)
        M = cur[:, 0]                                                      # (B,C,C) full product
        alpha_last = torch.logsumexp(alpha0.unsqueeze(2) + M, dim=1)       # (B,C)
        return torch.logsumexp(alpha_last + self.end.unsqueeze(0), dim=1)  # (B,)

    def _gold_score(self, emissions, tags, mask):
        B, W, C = emissions.shape
        tags_safe = tags.clone()
        tags_safe[~mask] = 0
        emit = emissions.gather(2, tags_safe.unsqueeze(2)).squeeze(2)      # (B,W)
        emit = (emit * mask).sum(1)                                        # masked sum
        # transitions for consecutive valid pairs t-1->t (both valid)
        prev, nxt = tags_safe[:, :-1], tags_safe[:, 1:]
        tr = self.trans[prev, nxt]                                         # (B,W-1)
        pair_mask = mask[:, :-1] & mask[:, 1:]
        tr = (tr * pair_mask).sum(1)
        start_s = self.start[tags_safe[:, 0]]
        last = mask.sum(1).clamp(min=1) - 1
        end_s = self.end[tags_safe.gather(1, last.unsqueeze(1)).squeeze(1)]
        return start_s + emit + tr + end_s

    def nll(self, emissions, tags, mask):
        # PER-TOKEN normalized so the CRF NLL matches a per-base CE scale (the summed-over-W NLL
        # ~1000s otherwise swamps the aux class-weighted CE and emissions never learn CDS).
        ntok = mask.sum(1).clamp(min=1).float()
        return ((self._forward_alg(emissions, mask) - self._gold_score(emissions, tags, mask)) / ntok).mean()

    @torch.no_grad()
    def viterbi(self, emissions, mask):
        # sequential (predict-only); identical to the reference implementation
        B, W, C = emissions.shape
        score = self.start.unsqueeze(0) + emissions[:, 0]
        backptr = []
        for t in range(1, W):
            s = score.unsqueeze(2) + self.trans.unsqueeze(0)
            best, bp = s.max(dim=1)
            cand = best + emissions[:, t]
            mt = mask[:, t].unsqueeze(1)
            score = torch.where(mt, cand, score)
            backptr.append(torch.where(mt, bp, torch.arange(C, device=emissions.device).unsqueeze(0)))
        score = score + self.end.unsqueeze(0)
        cur = score.argmax(dim=1)
        paths = torch.zeros(B, W, dtype=torch.long, device=emissions.device)
        paths[:, W - 1] = cur
        for t in range(W - 2, -1, -1):
            cur = backptr[t].gather(1, cur.unsqueeze(1)).squeeze(1)
            paths[:, t] = cur
        return paths


class SemiCRF(nn.Module):
    """Bounded semi-Markov CRF: scores SEGMENTS (contiguous same-label runs) up to
    max_seg_len, with a label-transition matrix. Segment emission = sum of per-base emissions
    over the segment for that label. O(W * L * C^2). Bounded by max_seg_len to stay tractable.
    For W>max_seg_len, long genes are represented as consecutive same-label segments."""

    def __init__(self, num_tags=NUM_CLASSES, max_seg_len=64):
        super().__init__()
        self.num_tags = num_tags
        self.max_seg_len = max_seg_len
        self.trans = nn.Parameter(torch.randn(num_tags, num_tags) * 0.1)  # label i -> j (j!=i)
        self.start = nn.Parameter(torch.randn(num_tags) * 0.1)

    def _cumemit(self, emissions):
        # prefix sums for O(1) segment-emission queries: csum[b,t,c]=sum_{i<t} emis[b,i,c]
        B, W, C = emissions.shape
        z = torch.zeros(B, 1, C, device=emissions.device)
        return torch.cat([z, torch.cumsum(emissions, dim=1)], dim=1)      # (B,W+1,C)

    def _seg_emit(self, csum, s, e):
        return csum[:, e] - csum[:, s]                                    # (B,C) sum over [s,e)

    def nll(self, emissions, tags, mask):
        # Partition (forward over segments) and gold-path score; both masked by valid length.
        B, W, C = emissions.shape
        csum = self._cumemit(emissions)
        lengths = mask.sum(1)                                             # (B,)
        # ---- forward (log partition) ----
        # alpha[b,t,c] = logsumexp over segmentations of [0,t) ending in a segment of label c
        nll_terms = []
        for b in range(B):
            L = int(lengths[b].item())
            if L <= 0:
                nll_terms.append(emissions.new_zeros(()))
                continue
            alpha = emissions.new_full((L + 1, C), NEG)
            # base: first segment [0,e)
            for e in range(1, min(self.max_seg_len, L) + 1):
                seg = csum[b, e] - csum[b, 0]                             # (C,)
                alpha[e] = torch.logsumexp(torch.stack([alpha[e], self.start + seg]), dim=0)
            for t in range(1, L + 1):
                for e in range(t + 1, min(t + self.max_seg_len, L) + 1):
                    seg = csum[b, e] - csum[b, t]                         # (C,)
                    # prev label c' -> new label c : alpha[t,c'] + trans[c',c] + seg[c]
                    prev = alpha[t].unsqueeze(1) + self.trans + seg.unsqueeze(0)  # (C',C)
                    contrib = torch.logsumexp(prev, dim=0)               # (C,)
                    alpha[e] = torch.logsumexp(torch.stack([alpha[e], contrib]), dim=0)
            logZ = torch.logsumexp(alpha[L], dim=0)
            # ---- gold score ----
            tg = tags[b, :L]
            gold = self._gold_path_score(csum[b], tg)
            nll_terms.append(logZ - gold)
        return torch.stack(nll_terms).mean()

    def _gold_path_score(self, csum_b, tags_b):
        # segment the gold tag sequence into maximal same-label runs (capped at max_seg_len)
        L = tags_b.shape[0]
        score = csum_b.new_zeros(())
        prev_label = None
        i = 0
        while i < L:
            j = i
            while j < L and tags_b[j] == tags_b[i] and (j - i) < self.max_seg_len:
                j += 1
            lab = int(tags_b[i].item())
            seg = csum_b[j] - csum_b[i]                                   # (C,)
            score = score + seg[lab]
            score = score + (self.start[lab] if prev_label is None else self.trans[prev_label, lab])
            prev_label = lab
            i = j
        return score

    @torch.no_grad()
    def decode(self, emissions, mask):
        B, W, C = emissions.shape
        csum = self._cumemit(emissions)
        lengths = mask.sum(1)
        out = torch.zeros(B, W, dtype=torch.long, device=emissions.device)
        for b in range(B):
            L = int(lengths[b].item())
            if L <= 0:
                continue
            dp = emissions.new_full((L + 1, C), NEG)
            bp = [[None] * C for _ in range(L + 1)]                       # (start_idx, prev_label)
            for e in range(1, min(self.max_seg_len, L) + 1):
                seg = csum[b, e] - csum[b, 0]
                for c in range(C):
                    val = self.start[c] + seg[c]
                    if val > dp[e, c]:
                        dp[e, c] = val; bp[e][c] = (0, -1)
            for t in range(1, L + 1):
                for e in range(t + 1, min(t + self.max_seg_len, L) + 1):
                    seg = csum[b, e] - csum[b, t]
                    for c in range(C):
                        prev = dp[t] + self.trans[:, c]                  # (C',)
                        best_prev = int(torch.argmax(prev).item())
                        val = prev[best_prev] + seg[c]
                        if val > dp[e, c]:
                            dp[e, c] = val; bp[e][c] = (t, best_prev)
            # backtrack from (L, best c)
            c = int(torch.argmax(dp[L]).item())
            e = L
            while e > 0:
                s, prevc = bp[e][c]
                out[b, s:e] = c
                e, c = s, prevc
        return out


def constrained_decode(pred_by_seqid, min_cds_len=30, max_fill_gap=20):
    """Deterministic fragmentation-reducing post-processing on per-base argmax labels.
    For each seqid array: (1) drop CDS runs shorter than min_cds_len (-> gene-body-nc);
    (2) within gene-body (label>0) regions, fill short intergenic gaps (<=max_fill_gap)
    that sit between CDS so adjacent CDS merge into one gene. Returns new {seqid: array}."""
    out = {}
    for sid, arr in pred_by_seqid.items():
        a = arr.copy()
        # (1) drop tiny CDS fragments
        is_cds = a == CLASS_CDS
        i = 0
        n = len(a)
        while i < n:
            if is_cds[i]:
                j = i
                while j < n and is_cds[j]:
                    j += 1
                if (j - i) < min_cds_len:
                    a[i:j] = CLASS_GENEBODY_NC
                i = j
            else:
                i += 1
        # (2) fill short intergenic gaps flanked by gene-body (merge fragmented genes)
        gene = a > 0
        i = 0
        while i < n:
            if not gene[i]:
                j = i
                while j < n and not gene[j]:
                    j += 1
                if 0 < i and j < n and (j - i) <= max_fill_gap:
                    a[i:j] = CLASS_GENEBODY_NC
                i = j
            else:
                i += 1
        out[sid] = a
    return out
