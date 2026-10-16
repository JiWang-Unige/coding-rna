"""Isolated same-partition/gradient recurrence with necessary endpoint indexing.

All retained events still pass the original ORF test. Original sum-product,
score operation order, state keys and research callers remain unchanged.
"""
from collections import defaultdict

import torch

from .candidates import Budget, ORFIndex, STOP
from .exon_end_bound import exon_event_range
from .full_grammar_partition import _TensorMargins


def indexed_full_grammar_log_partition(sequence, logits, max_exons):
    """FP64 logZ and gradient of unit-role gain over all legal gene sets.

    Motifs, actual unfinished codons, ambiguity policy, intron minimum and
    explicit exon cap match decode_emission_lattice. Donor histories with the
    same (exon count, actual carry) are summed after intron-prefix subtraction,
    never winner-pruned. The null prefix enters once; earlier completed sets
    also enter once. Adjacent genes are allowed, overlapping spans are not.
    No observed-mask or gold-event likelihood is defined by this function.
    """
    sequence = sequence.upper()
    length = len(sequence)
    if max_exons < 1:
        raise ValueError('A positive explicit exon cap is required')
    if not length or logits.ndim != 2 or tuple(logits.shape) != (length, 15):
        raise ValueError('One valid-length 15-role row per oriented DNA base required')
    if not logits.is_floating_point() or not bool(torch.isfinite(logits).all()):
        raise ValueError('Finite floating local potentials required')
    scores, orf = _TensorMargins(logits), ORFIndex(sequence)
    starts = {i for i in range(length - 2) if sequence[i:i + 3] == 'ATG'}
    acceptors = {i + 2 for i in range(length - 2) if sequence[i:i + 2] == 'AG'}
    donors = {i for i in range(1, length - 1) if sequence[i:i + 2] in ('GT', 'GC')}
    stops = {i + 3 for i in range(length - 2) if sequence[i:i + 3] in STOP}
    events = sorted([(d, False) for d in donors] + [(b, True) for b in stops])
    gap = Budget().min_intron_bp
    maturities = {d + gap: d for d in donors if d + gap <= length}
    positions = sorted(starts | acceptors | stops | set(maturities))
    pending, active, completed, edge_cache = defaultdict(dict), {}, {}, {}
    prefix_mass = logits[0, 0].double() * 0
    counts = dict(legal_exon_edges=0, state_exon_extensions=0,
                  donor_state_proposals=0, mature_state_proposals=0,
                  terminal_state_proposals=0, maximum_active_states=0,
                  indexed_cache_keys=0, indexed_events_inspected=0,
                  original_all_event_loop_iterations=0)

    def add(table, key, value):
        old = table.get(key)
        table[key] = value if old is None else torch.logaddexp(old, value)

    def legal_ends(a, carry, first):
        key = (a, carry, first)
        if key not in edge_cache:
            choices = []
            lo, hi = exon_event_range(orf, events, a, carry)
            counts['indexed_cache_keys'] += 1
            counts['indexed_events_inspected'] += hi - lo
            counts['original_all_event_loop_iterations'] += len(events)
            for b, terminal in events[lo:hi]:
                if b <= a or ((first or terminal) and b - a < 3):
                    continue
                remainder = orf.extend(a, b, carry, terminal)
                if remainder is not None:
                    choices.append((b, terminal, remainder,
                                    scores.exon_gain(a, b, len(carry), first, terminal)))
            edge_cache[key] = choices
            counts['legal_exon_edges'] += len(choices)
        return edge_cache[key]

    for a in positions:
        mass = completed.pop(a, None)
        if mass is not None:
            prefix_mass = torch.logaddexp(prefix_mass, mass)
        d = maturities.get(a)
        if d is not None:
            for key, value in sorted(pending.pop(d, {}).items()):
                _, carry = key
                role = 4 + (len(carry) - 1) % 3
                add(active, key, value - scores.prefix[role, d])
                counts['mature_state_proposals'] += 1
            counts['maximum_active_states'] = max(counts['maximum_active_states'], len(active))
        entries = [(0, '', prefix_mass)] if a in starts else []
        if a in acceptors:
            for (count, carry), value in sorted(active.items()):
                role = 4 + (len(carry) - 1) % 3
                entries.append((count, carry, value + scores.prefix[role, a]))
        for count, carry, value in entries:
            for b, terminal, remainder, delta in legal_ends(a, carry, count == 0):
                counts['state_exon_extensions'] += 1
                proposal = value + delta
                if terminal:
                    add(completed, b, proposal)
                    counts['terminal_state_proposals'] += 1
                elif count + 1 < max_exons:
                    add(pending[b], (count + 1, remainder), proposal)
                    counts['donor_state_proposals'] += 1
    return dict(log_partition=prefix_mass, counts=counts, reference_used=False,
                enumerated_candidate_pool_used=False, Viterbi_winner_pruning=False)
