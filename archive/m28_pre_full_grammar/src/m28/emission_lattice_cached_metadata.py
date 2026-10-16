"""Isolated artificial A42 preparation: DNA metadata and integer bookkeeping.

Scores, state updates, traceback and phase-block exit comparisons are unchanged.
No frozen/research caller uses this preparation implementation.
"""
from bisect import bisect_right
from collections import defaultdict

import numpy as np
import torch

from .candidates import Budget, ORFIndex, STOP
from .emission_lattice import _Live, _Margins, _Solution, _retain
from .emission_lattice_phase_block import _phase_block_winners
from .exon_end_bound import exon_event_range
from .long_exon_phase import long_exon_phase


@torch.no_grad()
def decode_cached_metadata_emission_lattice(sequence, logits, max_exons):
    """Same grammar/FP64/ties; aggregate only integer diagnostic counters."""
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
        long_array_blocks=0, long_carry_column_batches=0, long_nonempty_history_groups=0,
        long_phase_metadata_evaluations=0, logical_counter_batches=0)

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
            cache[key] = (result, np.asarray([row[3] for row in result], dtype=np.float64),
                          len(result), sum(row[1] for row in result))
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

    def account(histories, capped_histories, edge_count, terminals):
        # Integer sums commute; no score, history winner or path is merged.
        counts['logical_exon_extensions'] += histories * edge_count
        counts['logical_terminal_extensions'] += histories * terminals
        counts['logical_cap_nonterminal_extensions'] += capped_histories * (edge_count - terminals)
        counts['logical_counter_batches'] += 1

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
            edges, _, edge_count, terminals = choices(a, '', 'first')
            account(1, int(max_exons <= 1), edge_count, terminals)
            for b, terminal, carry, delta in edges:
                update(a, row, b, terminal, carry, row.score + delta)
        if a not in acceptors:
            continue
        groups, carry_metadata = {}, {}
        for _, original in sorted(active.items()):
            intron_role = 4 + (len(original.carry) - 1) % 3
            row = _Live(original.score + scores.prefix[intron_role][a],
                        original.preceding, original.chain, original.carry)
            if row.carry not in carry_metadata:
                # A list exists even when phase is None or short edges empty.
                # It is local to this acceptor, so DNA identity is implicit.
                carry_metadata[row.carry] = [choices(a, row.carry, 'short'),
                    long_exon_phase(sequence, a, row.carry), 0, 0]
                counts['long_phase_metadata_evaluations'] += 1
            metadata = carry_metadata[row.carry]
            metadata[2] += 1
            metadata[3] += len(row.chain) + 1 >= max_exons
            short, _, _, _ = metadata[0]
            for b, terminal, carry, delta in short:
                update(a, row, b, terminal, carry, row.score + delta)
            phase = metadata[1]
            if phase is not None:
                groups.setdefault((len(row.chain), phase), []).append(row)
        for metadata, _, histories, capped in carry_metadata.values():
            _, _, edge_count, terminals = metadata
            account(histories, capped, edge_count, terminals)
            counts['short_actual_carry_extensions'] += histories * edge_count
        phases = {}
        for (_, phase), rows in groups.items():
            phases.setdefault(phase, []).append(rows)
        for block in phases.values():
            edges, deltas, edge_count, terminals = choices(a, block[0][0].carry, 'long')
            histories = sum(len(rows) for rows in block)
            capped = sum(len(rows) for rows in block if len(rows[0].chain) + 1 >= max_exons)
            account(histories, capped, edge_count, terminals)
            counts['long_history_groups'] += len(block)
            counts['long_group_history_extensions'] += histories * edge_count
            counts['long_exit_proposals'] += len(block) * edge_count
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
