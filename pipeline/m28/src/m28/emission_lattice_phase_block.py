"""Isolated artificial A42 preparation: array blocks across exon-count groups.

Only array dispatch granularity changes from emission_lattice_batched.
Frozen/research callers and the previous batch implementation are unchanged.
"""
from bisect import bisect_right
from collections import defaultdict

import numpy as np
import torch

from .candidates import Budget, ORFIndex, STOP
from .emission_lattice import _Live, _Margins, _Solution, _retain
from .exon_end_bound import exon_event_range
from .long_exon_phase import long_exon_phase


def _phase_block_winners(groups, deltas):
    """Keep each exon-count row independent; add FP64 before strict compare.

Actual carry columns are in the original lexical arrival order. Missing
slots are -inf placeholders, never paths. Iterate carry columns with a
strict comparison to retain the first tied exit. Working arrays are 2-D;
there is no exon-count by carry by endpoint intermediate.
"""
    carries = sorted({row.carry for rows in groups for row in rows})
    columns = {carry: i for i, carry in enumerate(carries)}
    entries = np.full((len(groups), len(carries)), -np.inf, dtype=np.float64)
    lookup = [[None] * len(carries) for _ in groups]
    for i, rows in enumerate(groups):
        for row in rows:
            col = columns[row.carry]
            entries[i, col] = row.score
            lookup[i][col] = row
    delta = np.asarray(deltas, dtype=np.float64)[None, :]
    values = entries[:, :1] + delta
    winners = np.zeros(values.shape, dtype=np.int64)
    for col in range(1, len(carries)):
        candidate = entries[:, col:col + 1] + delta
        better = candidate > values
        values[better] = candidate[better]
        winners[better] = col
    return lookup, winners.tolist(), values.tolist(), len(carries)


@torch.no_grad()
def decode_phase_block_emission_lattice(sequence, logits, max_exons):
    """Original grammar/FP64/ties with independent count rows in each block."""
    sequence = sequence.upper()
    length = len(sequence)
    if max_exons < 1:
        raise ValueError('A positive explicit exon cap is required')
    if logits.ndim != 2 or tuple(logits.shape) != (length, 15) or not length:
        raise ValueError('One valid-length 15-role row per oriented DNA base required')
    if logits.dtype != torch.float32 or logits.device.type != 'cpu' or not bool(torch.isfinite(logits).all()):
        raise ValueError('Existing finite CPU FP32 emission input required')
    scores, orf = _Margins(logits), ORFIndex(sequence)
    starts = {i for i in range(length - 2) if sequence[i:i + 3] == 'ATG'}
    acceptors = {i + 2 for i in range(length - 2) if sequence[i:i + 2] == 'AG'}
    donors = {i for i in range(1, length - 1) if sequence[i:i + 2] in ('GT', 'GC')}
    stops = {i + 3 for i in range(length - 2) if sequence[i:i + 3] in STOP}
    events = sorted([(d, False) for d in donors] + [(b, True) for b in stops])
    gap = Budget().min_intron_bp
    maturities = {d + gap: d for d in donors if d + gap <= length}
    positions = sorted(starts | acceptors | stops | set(maturities))
    pending, active, completed = defaultdict(dict), {}, {}
    first_cache, short_cache, long_cache = {}, {}, {}
    best = _Solution()
    counts = dict(first_exon_edges=0, short_actual_carry_edges=0, shared_long_exon_edges=0,
        logical_exon_extensions=0, logical_terminal_extensions=0,
        logical_cap_nonterminal_extensions=0, state_update_proposals=0,
        donor_state_proposals=0, donor_occupied_key_comparisons=0,
        mature_intron_key_comparisons=0, maximum_active_states=0,
        short_actual_carry_extensions=0, long_history_groups=0,
        long_group_history_extensions=0, long_exit_proposals=0,
        long_array_blocks=0, long_carry_column_batches=0, long_nonempty_history_groups=0)

    def choices(a, carry, kind):
        cache = first_cache if kind == 'first' else short_cache if kind == 'short' else long_cache
        key = a if kind == 'first' else (a, carry) if kind == 'short' else (a, len(carry))
        if key not in cache:
            lo, hi = exon_event_range(orf, events, a, carry)
            if kind == 'short':
                hi = min(hi, bisect_right(events, (a + 2, True)))
            else:
                lo = max(lo, bisect_right(events, (a + 2, True)))
            result = []
            for b, terminal in events[lo:hi]:
                if kind == 'short' and terminal:
                    continue
                remainder = orf.extend(a, b, carry, terminal)
                if remainder is not None:
                    delta = scores.exon_gain(a, b, len(carry), kind == 'first', terminal)
                    result.append((b, terminal, remainder, delta))
            count_key = 'first_exon_edges' if kind == 'first' else 'short_actual_carry_edges' if kind == 'short' else 'shared_long_exon_edges'
            counts[count_key] += len(result)
            cache[key] = (result, np.asarray([row[3] for row in result], dtype=np.float64))
        return cache[key]

    def update(a, row, b, terminal, carry, score):
        counts['state_update_proposals'] += 1
        exon_count = len(row.chain) + 1
        if terminal:
            old = completed.get(b)
            if old is None or score > old.score:
                chain = row.chain + ((a, b),)
                completed[b] = _Solution(score, row.preceding.chains + (chain,),
                    row.preceding.gains + (score - row.preceding.score,))
        elif exon_count < max_exons:
            counts['donor_state_proposals'] += 1
            table, key = pending[b], (exon_count, carry)
            old = table.get(key)
            counts['donor_occupied_key_comparisons'] += old is not None
            if old is None or score > old.score:
                table[key] = _Live(score, row.preceding, row.chain + ((a, b),), carry)

    def logical_count(row, edges, histories=1):
        counts['logical_exon_extensions'] += histories * len(edges)
        counts['logical_terminal_extensions'] += histories * sum(edge[1] for edge in edges)
        if len(row.chain) + 1 >= max_exons:
            counts['logical_cap_nonterminal_extensions'] += histories * sum(not edge[1] for edge in edges)

    for a in positions:
        solution = completed.pop(a, None)
        if solution is not None and solution.score > best.score:
            best = solution
        d = maturities.get(a)
        if d is not None:
            for key, row in sorted(pending.pop(d, {}).items()):
                intron_role = 4 + (len(row.carry) - 1) % 3
                normalized = _Live(row.score - scores.prefix[intron_role][d],
                                   row.preceding, row.chain, row.carry)
                counts['mature_intron_key_comparisons'] += _retain(active, key, normalized)
            counts['maximum_active_states'] = max(counts['maximum_active_states'], len(active))
        if a in starts:
            row = _Live(best.score, best, (), '')
            edges, _ = choices(a, '', 'first')
            logical_count(row, edges)
            for b, terminal, carry, delta in edges:
                update(a, row, b, terminal, carry, row.score + delta)
        if a not in acceptors:
            continue
        groups = {}
        for _, original in sorted(active.items()):
            intron_role = 4 + (len(original.carry) - 1) % 3
            row = _Live(original.score + scores.prefix[intron_role][a],
                        original.preceding, original.chain, original.carry)
            short, _ = choices(a, row.carry, 'short')
            logical_count(row, short)
            counts['short_actual_carry_extensions'] += len(short)
            for b, terminal, carry, delta in short:
                update(a, row, b, terminal, carry, row.score + delta)
            phase = long_exon_phase(sequence, a, row.carry)
            if phase is not None:
                groups.setdefault((len(row.chain), phase), []).append(row)
        phases = {}
        for (_, phase), rows in groups.items():
            phases.setdefault(phase, []).append(rows)
        # Within each phase count rows keep ascending original count order.
        # Different phases cannot compete at a given terminal endpoint or
        # donor DP key. The original strict updates and first exons remain.
        for block in phases.values():
            edges, deltas = choices(a, block[0][0].carry, 'long')
            for rows in block:
                logical_count(rows[0], edges, len(rows))
                counts['long_history_groups'] += 1
                counts['long_group_history_extensions'] += len(rows) * len(edges)
                counts['long_exit_proposals'] += len(edges)
            if edges:
                lookup, winners, values, width = _phase_block_winners(block, deltas)
                counts['long_array_blocks'] += 1
                counts['long_carry_column_batches'] += width
                counts['long_nonempty_history_groups'] += len(block)
                for i in range(len(block)):
                    for edge, winner, value in zip(edges, winners[i], values[i]):
                        b, terminal, carry, _ = edge
                        update(a, lookup[i][winner], b, terminal, carry, value)
    return dict(chains=best.chains, gain=best.score, chain_gains=best.gains,
                counts=counts, reference_used=False)
