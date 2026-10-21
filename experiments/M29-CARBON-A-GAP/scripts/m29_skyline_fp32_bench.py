"""G0b: FP32 vs FP64 skyline forward-pass microbenchmark (copy of src/m28/full_grammar_skyline.py _forward with float32 alpha/gain).
Real oriented strand + Carbon-A fixed-offset logits (as in m29_gate_parity), 1 Mb slice, CPU single thread.  Forward pass only
(the trace pass repeats it); reports seconds, maximum active states and whether the optimum score / forward boundary score agree."""
import sys, time
from functools import lru_cache
import numpy as np
from numba import njit
import torch
from src.m28.full_grammar_boundary_compiled import _state_table

@lru_cache(maxsize=4)
def _skyline_groups(cap):
    states, _, _, _ = _state_table(cap)
    residuals = sorted({(stage, carry) for stage, _, carry in states})
    index = {residual: i for i, residual in enumerate(residuals)}
    return np.array([index[stage, carry] for stage, _, carry in states], dtype=np.int32), len(residuals)


@njit(cache=False, fastmath=False, parallel=False)
def _prune64(following, candidates, count, groups, best, active):
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
def _step64(alpha, active, active_count, base, row, destinations, roles, background,
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
    kept, discarded = _prune64(following, candidates, candidate_count, groups, best, next_active)
    # Clear just the current frontier before using this array as scratch.
    for i in range(active_count):
        alpha[active[i]] = -np.inf
    return kept, proposals, discarded


@njit(cache=False, fastmath=False, parallel=False)
def _forward64(codes, logits, positions, destinations, roles, background, groups, group_count):
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
        kept, edges, removed = _step64(alpha, active, active_count, codes[i], logits[i],
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
def _prune32(following, candidates, count, groups, best, active):
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
def _step32(alpha, active, active_count, base, row, destinations, roles, background,
          following, predecessor, gain, candidates, groups, best, next_active):
    for role in range(15):
        gain[role] = np.float32(row[role]) - np.float32(row[0])
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
    kept, discarded = _prune32(following, candidates, candidate_count, groups, best, next_active)
    # Clear just the current frontier before using this array as scratch.
    for i in range(active_count):
        alpha[active[i]] = -np.inf
    return kept, proposals, discarded


@njit(cache=False, fastmath=False, parallel=False)
def _forward32(codes, logits, positions, destinations, roles, background, groups, group_count):
    size = destinations.shape[1]
    alpha = np.full(size, -np.inf, dtype=np.float32)
    alpha[background] = 0.
    following = np.full(size, -np.inf, dtype=np.float32)
    active = np.empty(size, dtype=np.int32)
    next_active = np.empty(size, dtype=np.int32)
    active[0] = background
    active_count = 1
    candidates = np.empty(size, dtype=np.int32)
    predecessor = np.empty(size, dtype=np.uint32)
    gain = np.empty(15, dtype=np.float32)
    best = np.empty(group_count, dtype=np.float32)
    checkpoints = np.empty((len(positions), size), dtype=np.float32)
    checkpoints[0] = alpha
    cursor = 1
    visited, proposals, discarded, maximum = 0, 0, 0, 1
    for i in range(len(codes)):
        visited += active_count
        kept, edges, removed = _step32(alpha, active, active_count, codes[i], logits[i],
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




def bench(seq, logits, cap=128):
    states, destinations, roles, background = _state_table(cap)
    groups, gc = _skyline_groups(cap)
    lookup = np.full(256, 4, dtype=np.int8)
    for i, base in enumerate('ACGTN'):
        lookup[ord(base)] = i
    codes = lookup[np.frombuffer(seq.encode('ascii'), dtype=np.uint8)]
    positions = np.array([0] + list(range(8192, len(seq), 8192)) + [len(seq)], dtype=np.int64)
    out = {}
    for name, fn, lg in (('fp64', _forward64, logits.numpy()), ('fp32', _forward32, logits.numpy())):
        fn(codes[:2000], lg[:2000], np.array([0, 2000], dtype=np.int64), destinations, roles, background, groups, gc)  # compile
        t = time.perf_counter()
        score, ck, stats = fn(codes, lg, positions, destinations, roles, background, groups, gc)
        out[name] = dict(seconds=time.perf_counter() - t, score=float(score), visited=int(stats[0]), proposals=int(stats[1]), max_active=int(stats[3]))
        print(name, out[name], flush=True)
    return out


if __name__ == '__main__':
    import m29_common as M
    from m29_gate_parity import load_tracks, build_logits
    seq = sys.argv[1] if len(sys.argv) > 1 else 'NC_089044.1'
    a, n = 6_000_000, 1_000_000
    rid, sequence, pos, neg = load_tracks(seq)
    tt = np.load(M.CACHE / f'tt_{seq}_+.npy')
    s = sequence[a:a + n]
    lg = torch.from_numpy(build_logits(s, pos[a:a + n], tt[a:a + n], 'A_full'))
    res = bench(s, lg)
    import json
    (M.M29 / 'G0b' / f'skyline_fp32_bench_{seq}.json').write_text(json.dumps(res, indent=1))
