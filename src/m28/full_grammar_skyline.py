"""Exact max-only exon-budget skyline on the complete boundary automaton.

At a fixed DNA position, a no-worse score with the same stage and actual
coding carry but fewer used exons dominates: every future of the higher
count remains available with the lower count. This is NOT a logZ reduction.
Only reachable, nondominated states are traversed. The dense checkpoint and
single-block predecessor layout is retained; neither closes interior genes.
"""
from functools import lru_cache

import numpy as np
from numba import njit
import torch

from .full_grammar_boundary_compiled import _state_table
from .full_grammar_viterbi import _selected_scores


@lru_cache(maxsize=4)
def _skyline_groups(cap):
    states, _, _, _ = _state_table(cap)
    residuals = sorted({(stage, carry) for stage, _, carry in states})
    index = {residual: i for i, residual in enumerate(residuals)}
    return np.array([index[stage, carry] for stage, _, carry in states], dtype=np.int32), len(residuals)


@njit(cache=False, fastmath=False, parallel=False)
def _prune(following, candidates, count, groups, best, active):
    # State-table order is stage, exon count, carry: each residual is visited
    # in increasing count even though different carries are interleaved.
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
def _step(alpha, active, active_count, base, row, destinations, roles, background,
          following, predecessor, gain, candidates, groups, best, next_active):
    for role in range(15):
        gain[role] = np.float64(row[role]) - np.float64(row[0])
    # Background-first tie preference preserves the all-zero null path.
    following[background] = alpha[background]
    predecessor[background] = background * 3 + 1
    candidates[0] = background
    candidate_count, proposals = 1, 0
    for i in range(active_count):
        state = active[i]
        for edge in range(3):
            target = destinations[base, state, edge]
            if target < 0:
                break
            proposals += 1
            value = alpha[state] + gain[roles[base, state, edge]]
            if value > following[target]:
                if following[target] == -np.inf:
                    candidates[candidate_count] = target
                    candidate_count += 1
                following[target] = value
                predecessor[target] = state * 3 + edge + 1
    kept, discarded = _prune(following, candidates, candidate_count, groups, best, next_active)
    # Clear just the current frontier before using this array as scratch.
    for i in range(active_count):
        alpha[active[i]] = -np.inf
    return kept, proposals, discarded


@njit(cache=False, fastmath=False, parallel=False)
def _forward(codes, logits, positions, destinations, roles, background, groups, group_count):
    size = destinations.shape[1]
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
        kept, edges, removed = _step(alpha, active, active_count, codes[i], logits[i],
            destinations, roles, background, following, predecessor, gain, candidates, groups, best, next_active)
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
def _trace(codes, logits, positions, checkpoints, destinations, roles, background, groups, group_count):
    size = destinations.shape[1]
    path = np.empty(len(codes), dtype=np.uint8)
    tail, maximum_bytes = background, 0
    gain = np.empty(15, dtype=np.float64)
    candidates = np.empty(size, dtype=np.int32)
    next_active = np.empty(size, dtype=np.int32)
    best = np.empty(group_count, dtype=np.float64)
    visited, proposals, discarded, maximum = 0, 0, 0, 1
    for block in range(len(positions) - 2, -1, -1):
        a, b = positions[block], positions[block + 1]
        # Unreached cells are never read: selected traceback follows only
        # explicitly written tokens of surviving states.
        predecessors = np.empty((b - a, size), dtype=np.uint32)
        maximum_bytes = max(maximum_bytes, predecessors.nbytes)
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
            kept, edges, removed = _step(alpha, active, active_count, codes[i], logits[i],
                destinations, roles, background, following, predecessors[i - a], gain,
                candidates, groups, best, next_active)
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
            if destinations[codes[i], state, edge] != tail:
                raise ValueError('Invalid skyline selected predecessor edge')
            path[i] = roles[codes[i], state, edge]
            tail = state
    if tail != background:
        raise ValueError('Skyline selected path does not begin in background')
    return path, maximum_bytes, (visited, proposals, discarded, maximum)


@torch.no_grad()
def decode_full_grammar_skyline(sequence, logits, max_exons, checkpoint_bp=8192):
    sequence = sequence.upper()
    n = len(sequence)
    if (not n or set(sequence) - set('ACGTN') or checkpoint_bp < 1
            or logits.device.type != 'cpu' or logits.dtype not in (torch.float32, torch.float64)
            or tuple(logits.shape) != (n, 15) or not bool(torch.isfinite(logits).all())):
        raise ValueError('Normalized DNA, aligned finite CPU roles and positive checkpoint span required')
    states, destinations, roles, background = _state_table(max_exons)
    if len(states) * 3 >= 2**32:
        raise ValueError('State table exceeds explicit uint32 predecessor representation')
    groups, group_count = _skyline_groups(max_exons)
    lookup = np.full(256, 4, dtype=np.int8)
    for i, base in enumerate('ACGTN'):
        lookup[ord(base)] = i
    codes = lookup[np.frombuffer(sequence.encode('ascii'), dtype=np.uint8)]
    positions = np.array([0] + list(range(checkpoint_bp, n, checkpoint_bp)) + [n], dtype=np.int64)
    values = logits.detach().numpy()
    score, checkpoints, forward_stats = _forward(codes, values, positions,
        destinations, roles, background, groups, group_count)
    path, trace_bytes, trace_stats = _trace(codes, values, positions, checkpoints,
        destinations, roles, background, groups, group_count)
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
    reconstructed, gains = _selected_scores(values, path, starts, ends)
    if reconstructed != score:
        raise ValueError('Skyline selected per-base score differs from optimum')
    if forward_stats != trace_stats:
        raise ValueError('Skyline recomputation changed traversal statistics')
    return dict(chains=tuple(chains), gain=float(score), chain_gains=tuple(float(v) for v in gains),
        counts=dict(state_count=len(states), residual_group_count=group_count,
            checkpoint_rows=len(positions), checkpoint_bytes=checkpoints.nbytes, checkpoint_bp=checkpoint_bp,
            maximum_block_predecessor_bytes=int(trace_bytes), selected_role_path_bytes=path.nbytes,
            input_logit_bytes=logits.numel() * logits.element_size(), forward_and_recomputed_bases=2*n,
            forward_active_state_visits=int(forward_stats[0]), forward_edge_proposals=int(forward_stats[1]),
            forward_dominated_states=int(forward_stats[2]), maximum_active_states=int(forward_stats[3])),
        reference_used=False, pruning='exact_same_stage_actual_carry_lower_exon_budget_dominance_max_only',
        tie_policy='background_continuation_then_state_table_arrival_lower_count_dominates_equal_score',
        score_accumulation='FP64_per_base_forward_order')
