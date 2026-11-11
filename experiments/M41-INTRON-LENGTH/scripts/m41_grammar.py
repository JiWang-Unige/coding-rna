"""M41 complete-CDS grammar with Tiberius-style short / long intron states (decision `intron-length-states`, protocol docs/review_20261016/33).

What changes vs. the frozen automaton (`pipeline/m28/src/m28/full_grammar_boundary_messages._transitions`):
  the single intron-interior stage I3 (self-loop, one per (exon count, actual coding carry) = one per phase and history) is replaced by
    * a SHORT branch  I3S1 .. I3Sr  (r = `short_stages`, default 1) and
    * a LONG  branch  I3L,
  entered from I2 and both leaving through the unchanged AG stage.  Every other stage, role, motif, the 5-bp minimum intron and the exon cap are the
  frozen ones (the frozen `_transitions` is called for every non-intron-interior stage).  Each edge carries a TYPE; the edge score is
  `role gain (logit - background logit)` + `bonus[type]` (type 0 = no bonus, the frozen edges).  With `IntronLength(enabled=False)` the state table,
  roles and every score are the frozen ones (bit-identical logZ, gradient and decode -- tested).

Edge types (bonus vector of length N_TYPES, built from three learnable scalars by `m41_prior.bonus`):
  1 I2 -> I3S1 (enter short)   2 I2 -> I3L (enter long)
  3 I3Si -> I3Si (short self)  4 I3Si -> I3S(i+1) (short advance, r > 1 only)  5 I3Sr -> AG (short exit)
  6 I3L -> I3L (long self)     7 I3L -> AG (long exit)
Intron length n = 5 + k, k = number of bases consumed by I3* after the mandatory I2 base; the latent branch path is summed in logZ and maximised in decode.

Partition function / gradient: dense FP64 checkpointed sum-product (same algorithm as `full_grammar_boundary_compiled` with the M29 alpha arena), plus the
gradient wrt the bonus vector.  Decode: the M29 skyline (exact max with same-(stage, carry) lower-exon-count dominance), edge scores incl. bonus.
"""
from dataclasses import dataclass
from functools import lru_cache
import math
import sys
from pathlib import Path

import numpy as np
from numba import njit
import torch
from torch.autograd.function import once_differentiable

_ROOT = Path(__file__).resolve().parents[3]
for _p in (str(_ROOT / 'pipeline/m28'),):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from src.m28.full_grammar_boundary_messages import BACKGROUND, _transitions as _frozen_transitions   # noqa: E402
from src.m28.full_grammar_observed_partition import InfeasibleObservation, allowed_roles_from_local_targets  # noqa: E402

ALPHABET = 'ACGTN'
N_TYPES = 8
T_ENTER_S, T_ENTER_L, T_SELF_S, T_ADV_S, T_EXIT_S, T_SELF_L, T_EXIT_L = 1, 2, 3, 4, 5, 6, 7


@dataclass(frozen=True)
class IntronLength:
    enabled: bool = False
    short_stages: int = 1

    def key(self):
        return (bool(self.enabled), int(self.short_stages) if self.enabled else 0)


OFF = IntronLength(False)


def transitions(state, base, cap, il):
    """(target, role, edge type) for one consumed base.  Frozen transitions everywhere except I2 / the intron-interior stages when enabled."""
    stage, count, carry = state
    if not il.enabled:
        for target, role in _frozen_transitions(state, base, cap):
            yield target, role, 0
        return
    r = il.short_stages
    role = 4 + (len(carry) - 1) % 3 if stage[0] in 'DIA' else None
    if stage == 'I2':
        yield ('I3S1', count, carry), role, T_ENTER_S
        yield ('I3L', count, carry), role, T_ENTER_L
    elif stage.startswith('I3S'):
        i = int(stage[3:])
        yield state, role, T_SELF_S
        if i < r:
            yield (f'I3S{i + 1}', count, carry), role, T_ADV_S
        elif base == 'A':
            yield ('AG', count, carry), role, T_EXIT_S
    elif stage == 'I3L':
        yield state, role, T_SELF_L
        if base == 'A':
            yield ('AG', count, carry), role, T_EXIT_L
    else:
        for target, role0 in _frozen_transitions(state, base, cap):
            assert target[0] != 'I3', 'frozen I2 -> I3 edge must be replaced'
            yield target, role0, 0


@lru_cache(maxsize=8)
def state_table(cap, il_key=(False, 0)):
    """-> states, destinations[5,S,3], roles[5,S,3], types[5,S,3], background index.  Structural closure from background (no path pruning)."""
    if cap < 1:
        raise ValueError('A positive explicit exon cap is required')
    il = IntronLength(*il_key) if il_key[0] else OFF
    if il.enabled and il.short_stages < 1:
        raise ValueError('short_stages >= 1 required')
    seen, pending = {BACKGROUND}, [BACKGROUND]
    while pending:
        state = pending.pop()
        for base in ALPHABET:
            for target, _, _ in transitions(state, base, cap, il):
                if target not in seen:
                    seen.add(target); pending.append(target)
    states = tuple(sorted(seen))
    index = {s: i for i, s in enumerate(states)}
    dest = np.full((5, len(states), 3), -1, dtype=np.int32)
    roles = np.zeros_like(dest, dtype=np.int8)
    types = np.zeros_like(dest, dtype=np.int8)
    for b, base in enumerate(ALPHABET):
        for i, s in enumerate(states):
            for e, (t, role, ty) in enumerate(transitions(s, base, cap, il)):
                dest[b, i, e], roles[b, i, e], types[b, i, e] = index[t], role, ty
    return states, dest, roles, types, index[BACKGROUND]


def encode(sequence):
    lookup = np.full(256, 4, dtype=np.int8)
    for i, base in enumerate(ALPHABET):
        lookup[ord(base)] = i
    return lookup[np.frombuffer(sequence.upper().encode('ascii'), dtype=np.uint8)]


# ---- sum-product -----------------------------------------------------------------------------------------------------------------
@njit(cache=False, fastmath=False, parallel=False)
def _add_mass(a, b):
    if a == -np.inf:
        return b
    if b == -np.inf:
        return a
    return max(a, b) + math.log1p(math.exp(-abs(a - b)))


@njit(cache=False, fastmath=False, parallel=False)
def _advance(alpha, base, gain, allowed, bonus, dest, roles, types):
    following = np.full(len(alpha), -np.inf, dtype=np.float64)
    for state in range(len(alpha)):
        value = alpha[state]
        if value == -np.inf:
            continue
        for edge in range(3):
            target = dest[base, state, edge]
            if target < 0:
                break
            role = roles[base, state, edge]
            if allowed[role]:
                following[target] = _add_mass(following[target], value + gain[role] + bonus[types[base, state, edge]])
    return following


@njit(cache=False, fastmath=False, parallel=False)
def _forward(codes, gain, allowed, bonus, initial, positions, dest, roles, types):
    checkpoints = np.empty((len(positions), len(initial)), dtype=np.float64)
    checkpoints[0] = initial
    alpha = initial.copy()
    c = 1
    for i in range(len(codes)):
        alpha = _advance(alpha, codes[i], gain[i], allowed[i], bonus, dest, roles, types)
        if i + 1 == positions[c]:
            checkpoints[c] = alpha
            c += 1
    return alpha, checkpoints


@njit(cache=False, fastmath=False, parallel=False)
def _reverse(codes, gain, allowed, bonus, positions, checkpoints, upstream, dest, roles, types):
    gradient = np.zeros_like(gain)
    gbonus = np.zeros(len(bonus), dtype=np.float64)
    adjoint = upstream.copy()
    longest = 0
    for block in range(len(positions) - 1):
        longest = max(longest, positions[block + 1] - positions[block])
    arena = np.empty((longest + 1, len(adjoint)), dtype=np.float64)
    peak = 0
    for block in range(len(positions) - 2, -1, -1):
        a, b = positions[block], positions[block + 1]
        alpha = arena[:b - a + 1]
        peak = max(peak, alpha.nbytes)
        alpha[0] = checkpoints[block]
        for i in range(a, b):
            alpha[i - a + 1] = _advance(alpha[i - a], codes[i], gain[i], allowed[i], bonus, dest, roles, types)
        for i in range(b - 1, a - 1, -1):
            previous = np.zeros_like(adjoint)
            for state in range(len(adjoint)):
                value = alpha[i - a, state]
                if value == -np.inf:
                    continue
                for edge in range(3):
                    target = dest[codes[i], state, edge]
                    if target < 0:
                        break
                    role = roles[codes[i], state, edge]
                    g = adjoint[target]
                    if g == 0. or not allowed[i, role]:
                        continue
                    ty = types[codes[i], state, edge]
                    contribution = g * math.exp(value + gain[i, role] + bonus[ty] - alpha[i - a + 1, target])
                    previous[state] += contribution
                    gradient[i, role] += contribution
                    gbonus[ty] += contribution
            adjoint = previous
    return gradient, gbonus, adjoint, peak


class _Apply(torch.autograd.Function):
    @staticmethod
    def forward(ctx, logits, bonus, initial, codes, allowed, cap, il_key, checkpoint_bp, stats):
        _, dest, roles, types, _ = state_table(cap, il_key)
        gain = (logits.detach().double() - logits.detach()[:, :1].double()).numpy()
        bon = bonus.detach().double().numpy().copy()
        bon[0] = 0.
        init = initial.detach().numpy().copy()
        mask = allowed.numpy().copy()
        positions = np.array([0] + list(range(checkpoint_bp, len(codes), checkpoint_bp)) + [len(codes)], dtype=np.int64)
        final, checkpoints = _forward(codes, gain, mask, bon, init, positions, dest, roles, types)
        ctx.saved = (codes, gain, mask, bon, positions, checkpoints, final, dest, roles, types)
        ctx.dtype, ctx.bdtype, ctx.stats = logits.dtype, bonus.dtype, stats
        stats.update(state_count=dest.shape[1], checkpoint_rows=len(positions), checkpoint_bytes=checkpoints.nbytes)
        return torch.from_numpy(final)

    @staticmethod
    @once_differentiable
    def backward(ctx, grad_output):
        codes, gain, mask, bon, positions, checkpoints, final, dest, roles, types = ctx.saved
        upstream = grad_output.detach().numpy().copy()
        if not np.isfinite(upstream).all() or np.any((final == -np.inf) & (upstream != 0)):
            raise ValueError('Only finite, supported boundary-state readouts have defined derivatives')
        gradient, gbonus, ginit, peak = _reverse(codes, gain, mask, bon, positions, checkpoints, upstream, dest, roles, types)
        gradient[:, 0] -= gradient.sum(axis=1)          # every role score is its logit minus the background logit
        gbonus[0] = 0.
        ctx.stats['maximum_recomputed_alpha_bytes'] = peak
        return (torch.from_numpy(gradient).to(ctx.dtype), torch.from_numpy(gbonus).to(ctx.bdtype), torch.from_numpy(ginit),
                None, None, None, None, None, None)


def zero_bonus():
    return torch.zeros(N_TYPES, dtype=torch.float64)


def log_partition(sequence, logits, max_exons, bonus=None, il=OFF, allowed=None, cuts=(), checkpoint_bp=256):
    """Global complete gene-set logZ (FP64) over the M41 automaton; differentiable wrt logits and bonus.  Cuts carry exact messages (as the frozen
    `compiled_boundary_full_grammar_log_partition`)."""
    n = len(sequence)
    cuts = tuple(cuts)
    if not n or tuple(logits.shape) != (n, 15) or cuts != tuple(sorted(set(cuts))) or any(not 0 < x < n for x in cuts):
        raise ValueError('Aligned nonempty global DNA and ordered distinct interior cuts required')
    if logits.device.type != 'cpu' or logits.dtype not in (torch.float32, torch.float64) or not bool(torch.isfinite(logits).all()):
        raise ValueError('Finite CPU FP32/FP64 potentials required')
    if allowed is None:
        allowed = torch.ones((n, 15), dtype=torch.bool)
    if tuple(allowed.shape) != (n, 15) or allowed.dtype != torch.bool:
        raise ValueError('Boolean allowed roles aligned to every DNA row required')
    if bonus is None:
        bonus = zero_bonus()
    if tuple(bonus.shape) != (N_TYPES,) or not bool(torch.isfinite(bonus).all()):
        raise ValueError(f'finite bonus vector of length {N_TYPES} required')
    if not il.enabled and bool((bonus.detach()[1:] != 0).any()):
        raise ValueError('nonzero intron-length bonus with the states disabled')
    states, _, _, _, background = state_table(max_exons, il.key())
    message = torch.full((len(states),), -float('inf'), dtype=torch.float64)
    message[background] = 0.
    codes_all = encode(sequence)
    stats = []
    bounds = (0,) + cuts + (n,)
    for a, b in zip(bounds, bounds[1:]):
        st = {}
        message = _Apply.apply(logits[a:b], bonus, message, codes_all[a:b], allowed[a:b], max_exons, il.key(), checkpoint_bp, st)
        stats.append(st)
    value = message[background]
    if not bool(torch.isfinite(value)):
        raise InfeasibleObservation('No complete global gene set satisfies the supplied observations')
    return dict(log_partition=value, statistics=stats)


def coarse_channel_nll(sequence, logits, max_exons, core_start, core_end, targets, epsilon, bonus=None, il=OFF, *, cuts=(), checkpoint_bp=256):
    """Same loss as the frozen `full_grammar_coarse_channel.coarse_channel_nll` (exact fine rows, epsilon coarse channel, unknown rows free) on the M41
    automaton.  The bonus enters numerator and denominator alike (the observation mask only restricts roles)."""
    n = core_end - core_start
    length = len(sequence)
    if not 0 <= core_start < core_end <= length:
        raise ValueError('Original core must lie inside the supplied fixed DNA view')
    original = allowed_roles_from_local_targets(targets, n)
    fine = torch.as_tensor(targets['seg'], device='cpu')[:n]
    region = torch.as_tensor(targets['region'], device='cpu')[:n]
    epsilon = torch.as_tensor(epsilon, dtype=torch.float64, device='cpu')
    if epsilon.ndim or not bool(torch.isfinite(epsilon)) or not 0 < float(epsilon.detach()) < 1:
        raise ValueError('One finite scalar observation error probability strictly between zero and one required')
    exact = fine >= 0
    coarse = (fine < 0) & (region >= 0)
    allowed = torch.ones((length, 15), dtype=torch.bool)
    allowed[core_start:core_end][exact] = original[exact]
    logq = torch.zeros((length, 15), dtype=torch.float64)
    positions = torch.where(coarse)[0] + core_start
    logq[positions] = torch.where(original[coarse], torch.log1p(-epsilon), torch.log(epsilon / 2))
    den = log_partition(sequence, logits, max_exons, bonus, il, cuts=cuts, checkpoint_bp=checkpoint_bp)
    wei = log_partition(sequence, logits.double() + logq, max_exons, bonus, il, allowed=allowed, cuts=cuts, checkpoint_bp=checkpoint_bp)
    observed = wei['log_partition'] + logq[:, 0].sum()
    return dict(loss=den['log_partition'] - observed, log_partition_all=den['log_partition'], log_partition_observed=observed,
                statistics_all=den['statistics'], statistics_observed=wei['statistics'])


# ---- skyline max decode (M29 `m29_skyline_reuse` + edge bonus) ------------------------------------------------------------------
@lru_cache(maxsize=8)
def _skyline_groups(cap, il_key):
    states = state_table(cap, il_key)[0]
    residuals = sorted({(stage, carry) for stage, _, carry in states})
    index = {r: i for i, r in enumerate(residuals)}
    return np.array([index[stage, carry] for stage, _, carry in states], dtype=np.int32), len(residuals)


@njit(cache=False, fastmath=False, parallel=False)
def _prune(following, candidates, count, groups, best, active):
    candidates[:count].sort()
    best.fill(-np.inf)
    kept = 0
    for i in range(count):
        state = candidates[i]
        group = groups[state]
        value = following[state]
        if value > best[group]:
            best[group] = value
            active[kept] = state
            kept += 1
        else:
            following[state] = -np.inf
    return kept, count - kept


@njit(cache=False, fastmath=False, parallel=False)
def _step(alpha, active, active_count, base, row, bonus, dest, roles, types, background,
          following, predecessor, gain, candidates, groups, best, next_active):
    for role in range(15):
        gain[role] = np.float64(row[role]) - np.float64(row[0])
    following[background] = alpha[background]
    predecessor[background] = background * 3 + 1
    candidates[0] = background
    candidate_count, proposals = 1, 0
    for i in range(active_count):
        state = active[i]
        for edge in range(3):
            target = dest[base, state, edge]
            if target < 0:
                break
            proposals += 1
            value = alpha[state] + gain[roles[base, state, edge]] + bonus[types[base, state, edge]]
            if value > following[target]:
                if following[target] == -np.inf:
                    candidates[candidate_count] = target
                    candidate_count += 1
                following[target] = value
                predecessor[target] = state * 3 + edge + 1
    kept, discarded = _prune(following, candidates, candidate_count, groups, best, next_active)
    for i in range(active_count):
        alpha[active[i]] = -np.inf
    return kept, proposals, discarded


@njit(cache=False, fastmath=False, parallel=False)
def _sky_forward(codes, logits, bonus, positions, dest, roles, types, background, groups, group_count):
    size = dest.shape[1]
    alpha = np.full(size, -np.inf, dtype=np.float64)
    alpha[background] = 0.
    following = np.full(size, -np.inf, dtype=np.float64)
    active = np.empty(size, dtype=np.int32)
    next_active = np.empty(size, dtype=np.int32)
    active[0] = background
    active_count = 1
    candidates = np.empty(size, dtype=np.int32)
    predecessor = np.empty(size, dtype=np.uint32)
    gain = np.empty(15, dtype=np.float64)
    best = np.empty(group_count, dtype=np.float64)
    checkpoints = np.empty((len(positions), size), dtype=np.float64)
    checkpoints[0] = alpha
    cursor = 1
    visited, proposals, discarded, maximum = 0, 0, 0, 1
    for i in range(len(codes)):
        visited += active_count
        kept, edges, removed = _step(alpha, active, active_count, codes[i], logits[i], bonus, dest, roles, types, background,
                                     following, predecessor, gain, candidates, groups, best, next_active)
        proposals += edges
        discarded += removed
        maximum = max(maximum, kept)
        alpha, following = following, alpha
        active, next_active = next_active, active
        active_count = kept
        if i + 1 == positions[cursor]:
            checkpoints[cursor] = alpha
            cursor += 1
    return alpha[background], checkpoints, (visited, proposals, discarded, maximum)


@njit(cache=False, fastmath=False, parallel=False)
def _sky_trace(codes, logits, bonus, positions, checkpoints, dest, roles, types, background, groups, group_count):
    size = dest.shape[1]
    path = np.empty(len(codes), dtype=np.uint8)
    ptype = np.empty(len(codes), dtype=np.uint8)
    tail = background
    longest = 0
    for block in range(len(positions) - 1):
        longest = max(longest, positions[block + 1] - positions[block])
    arena = np.empty((longest, size), dtype=np.uint32)
    gain = np.empty(15, dtype=np.float64)
    candidates = np.empty(size, dtype=np.int32)
    next_active = np.empty(size, dtype=np.int32)
    best = np.empty(group_count, dtype=np.float64)
    visited, proposals, discarded, maximum = 0, 0, 0, 1
    for block in range(len(positions) - 2, -1, -1):
        a, b = positions[block], positions[block + 1]
        predecessors = arena[:b - a]
        alpha = checkpoints[block].copy()
        following = np.full(size, -np.inf, dtype=np.float64)
        active = np.empty(size, dtype=np.int32)
        active_count = 0
        for state in range(size):
            if alpha[state] != -np.inf:
                active[active_count] = state
                active_count += 1
        for i in range(a, b):
            visited += active_count
            kept, edges, removed = _step(alpha, active, active_count, codes[i], logits[i], bonus, dest, roles, types, background,
                                         following, predecessors[i - a], gain, candidates, groups, best, next_active)
            proposals += edges
            discarded += removed
            maximum = max(maximum, kept)
            alpha, following = following, alpha
            active, next_active = next_active, active
            active_count = kept
        if alpha[tail] != checkpoints[block + 1, tail]:
            raise ValueError('Skyline selected boundary score changed under recomputation')
        for i in range(b - 1, a - 1, -1):
            token = predecessors[i - a, tail]
            if token == 0:
                raise ValueError('Missing skyline predecessor on selected path')
            state, edge = (token - 1) // 3, (token - 1) % 3
            if dest[codes[i], state, edge] != tail:
                raise ValueError('Invalid skyline selected predecessor edge')
            path[i] = roles[codes[i], state, edge]
            ptype[i] = types[codes[i], state, edge]
            tail = state
    if tail != background:
        raise ValueError('Skyline selected path does not begin in background')
    return path, ptype, (visited, proposals, discarded, maximum)


@njit(cache=False, fastmath=False, parallel=False)
def _selected_scores(logits, bonus, path, ptype, starts, ends):
    """Same accumulation order as the forward pass: ((total + role gain) + bonus).  Background continuation adds nothing."""
    total = 0.
    baseline = 0.
    cursor = 0
    gains = np.empty(len(starts), dtype=np.float64)
    inside = False
    for i in range(len(path)):
        if cursor < len(starts) and i == starts[cursor]:
            baseline = total
            inside = True
        if inside:
            total = total + (np.float64(logits[i, path[i]]) - np.float64(logits[i, 0])) + bonus[ptype[i]]
        if cursor < len(ends) and i + 1 == ends[cursor]:
            gains[cursor] = total - baseline
            cursor += 1
            inside = False
    if cursor != len(starts):
        raise ValueError('Incomplete selected gene score intervals')
    return total, gains


@torch.no_grad()
def decode_skyline(sequence, logits, max_exons, bonus=None, il=OFF, checkpoint_bp=8192):
    """Exact max over (complete gene set, latent intron branch path).  Returns the frozen decoder's dict (chains, gain, chain_gains, counts) plus
    `introns` = per selected intron (branch S|L).  With il=OFF and bonus=0 the result equals the frozen `m29_skyline_reuse.decode_full_grammar_skyline`."""
    sequence = sequence.upper()
    n = len(sequence)
    if (not n or set(sequence) - set('ACGTN') or checkpoint_bp < 1 or logits.device.type != 'cpu'
            or logits.dtype not in (torch.float32, torch.float64) or tuple(logits.shape) != (n, 15) or not bool(torch.isfinite(logits).all())):
        raise ValueError('Normalized DNA, aligned finite CPU roles and positive checkpoint span required')
    bon = np.zeros(N_TYPES) if bonus is None else np.asarray(bonus.detach().double().numpy() if torch.is_tensor(bonus) else bonus, dtype=np.float64).copy()
    bon[0] = 0.
    if not il.enabled and np.any(bon != 0):
        raise ValueError('nonzero intron-length bonus with the states disabled')
    states, dest, roles, types, background = state_table(max_exons, il.key())
    if len(states) * 3 >= 2**32:
        raise ValueError('State table exceeds explicit uint32 predecessor representation')
    groups, group_count = _skyline_groups(max_exons, il.key())
    codes = encode(sequence)
    positions = np.array([0] + list(range(checkpoint_bp, n, checkpoint_bp)) + [n], dtype=np.int64)
    values = logits.detach().numpy()
    score, checkpoints, fstats = _sky_forward(codes, values, bon, positions, dest, roles, types, background, groups, group_count)
    path, ptype, tstats = _sky_trace(codes, values, bon, positions, checkpoints, dest, roles, types, background, groups, group_count)
    starts = np.flatnonzero((path == 13) & np.r_[True, path[:-1] != 13])
    ends = np.flatnonzero((path == 14) & np.r_[path[1:] != 14, True]) + 1
    if len(starts) != len(ends) or np.any(starts >= ends):
        raise ValueError('Incomplete skyline selected gene boundary roles')
    coding = (path != 0) & ((path < 4) | (path > 6))
    opening = np.r_[True, ~coding[:-1]]
    opening[starts] = True
    closing = np.r_[~coding[1:], True]
    closing[ends - 1] = True
    exon_starts = np.flatnonzero(coding & opening)
    exon_ends = np.flatnonzero(coding & closing) + 1
    chains, cursor = [], 0
    for a, b in zip(starts, ends):
        chain = []
        while cursor < len(exon_starts) and exon_starts[cursor] < b:
            x, y = int(exon_starts[cursor]), int(exon_ends[cursor])
            if not a <= x < y <= b:
                raise ValueError('Skyline selected exon outside its complete gene')
            chain.append((x, y))
            cursor += 1
        if not 1 <= len(chain) <= max_exons:
            raise ValueError('Skyline selected chain exceeds unchanged exon cap')
        chains.append(tuple(chain))
    if cursor != len(exon_starts):
        raise ValueError('Skyline selected coding rows outside complete genes')
    reconstructed, gains = _selected_scores(values, bon, path, ptype, starts, ends)
    if reconstructed != score:
        raise ValueError('Skyline selected per-base score differs from optimum')
    if fstats != tstats:
        raise ValueError('Skyline recomputation changed traversal statistics')
    branch = []
    for x in np.flatnonzero((ptype == T_ENTER_S) | (ptype == T_ENTER_L)):
        branch.append('S' if ptype[x] == T_ENTER_S else 'L')
    return dict(chains=tuple(chains), gain=float(score), chain_gains=tuple(float(v) for v in gains), intron_branches=tuple(branch),
                counts=dict(state_count=len(states), residual_group_count=group_count, checkpoint_rows=len(positions),
                            checkpoint_bytes=checkpoints.nbytes, checkpoint_bp=checkpoint_bp, forward_active_state_visits=int(fstats[0]),
                            forward_edge_proposals=int(fstats[1]), forward_dominated_states=int(fstats[2]), maximum_active_states=int(fstats[3])),
                intron_length=dict(enabled=il.enabled, short_stages=il.short_stages, bonus=[float(v) for v in bon]))
