"""Isolated artificial A42 preparation: batch long-exon history comparisons.

Frozen/research callers are unchanged. The supported input here is the
existing CPU FP32 15-column emission format, with original FP64 arithmetic.
"""
from bisect import bisect_right
from collections import defaultdict

import numpy as np
import torch

from .candidates import Budget, ORFIndex, STOP
from .emission_lattice import _Live, _Margins, _Solution, _retain
from .exon_end_bound import exon_event_range
from .long_exon_phase import long_exon_phase


def _long_exit_winners(rows, deltas):
    """Add each ORIGINAL entry score before first-arrival argmax per exit.

Do not select the largest entry score before the common delta: FP64 may
round unequal entries to a tied exit score. Row order is original arrival.
"""
    values = np.asarray([row.score for row in rows], dtype=np.float64)[:, None] + np.asarray(deltas, dtype=np.float64)[None, :]
    winners = values.argmax(axis=0)
    return winners.tolist(), values[winners, np.arange(len(deltas))].tolist()


@torch.no_grad()
def decode_batched_emission_lattice(sequence, logits, max_exons):
    """Same complete-CDS target and first-arrival chain identity, not a beam.

Short middle exons retain actual carry and individual updates. Long exits
share legality/role deltas only within eligible phase and completed-exon
count. Every history still contributes its independently rounded exit sum.
Counts distinguish logical extensions from actual state-update proposals.
"""
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
        long_group_history_extensions=0, long_exit_proposals=0)

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
        # Group insertion retains ascending exon-count order. For a fixed
        # terminal endpoint exactly one phase can finish; at donor endpoints
        # different phases have distinct outgoing carry lengths/DP keys.
        # Short and long endpoints are disjoint. These updates cannot reorder
        # a competing equal-score arrival across different phase groups.
        for rows in groups.values():
            edges, deltas = choices(a, rows[0].carry, 'long')
            logical_count(rows[0], edges, len(rows))
            counts['long_history_groups'] += 1
            counts['long_group_history_extensions'] += len(rows) * len(edges)
            counts['long_exit_proposals'] += len(edges)
            if edges:
                winners, values = _long_exit_winners(rows, deltas)
                for edge, winner, value in zip(edges, winners, values):
                    b, terminal, carry, _ = edge
                    update(a, rows[winner], b, terminal, carry, value)
    return dict(chains=best.chains, gain=best.score, chain_gains=best.gains,
                counts=counts, reference_used=False)
