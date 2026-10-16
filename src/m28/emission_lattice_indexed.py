"""Isolated A42 grammar variant that skips necessarily illegal end events.

The frozen emission_lattice module and its callers are unchanged. This
variant shares the original role arithmetic and DP states, not model weights.
"""
from collections import defaultdict

import torch

from .candidates import Budget, ORFIndex, STOP
from .emission_lattice import _Live, _Margins, _Solution, _retain
from .exon_end_bound import exon_event_range


@torch.no_grad()
def decode_indexed_emission_lattice(sequence, logits, max_exons):
    """Same grammar, scores, ordered LEGAL edges, and ties as the original DP.

Only the necessary start/end event range is indexed. Every retained event
still passes ORFIndex.extend, including terminal STOP and short carry exits.
No gene length, motif top-k, beam, zero shortcut, or scoring term is added.
"""
    sequence = sequence.upper()
    length = len(sequence)
    if max_exons < 1:
        raise ValueError('A positive explicit exon cap is required')
    if logits.ndim != 2 or tuple(logits.shape) != (length, 15) or not length:
        raise ValueError('One valid-length 15-role row per oriented DNA base required')
    if not bool(torch.isfinite(logits).all()):
        raise ValueError('Finite local potentials required')
    scores, orf = _Margins(logits), ORFIndex(sequence)
    starts = {i for i in range(length - 2) if sequence[i:i + 3] == 'ATG'}
    acceptors = {i + 2 for i in range(length - 2) if sequence[i:i + 2] == 'AG'}
    donors = {i for i in range(1, length - 1) if sequence[i:i + 2] in ('GT', 'GC')}
    stops = {i + 3 for i in range(length - 2) if sequence[i:i + 3] in STOP}
    events = sorted([(d, False) for d in donors] + [(b, True) for b in stops])
    gap = Budget().min_intron_bp
    maturities = {d + gap: d for d in donors if d + gap <= length}
    positions = sorted(starts | acceptors | stops | set(maturities))
    pending, active, completed, edge_cache = defaultdict(dict), {}, {}, {}
    best = _Solution()
    counts = dict(legal_exon_edges=0, exon_extensions=0, donor_states_created=0,
                  donor_states_dominated=0, intron_states_dominated=0,
                  completed_extensions=0, exon_cap_nonterminal_edges=0,
                  maximum_active_states=0)

    def legal_ends(a, carry, first):
        key = (a, carry, first)
        if key not in edge_cache:
            choices = []
            lo, hi = exon_event_range(orf, events, a, carry)
            for b, terminal in events[lo:hi]:
                if (first or terminal) and b - a < 3:
                    continue
                remainder = orf.extend(a, b, carry, terminal)
                if remainder is not None:
                    delta = scores.exon_gain(a, b, len(carry), first, terminal)
                    choices.append((b, terminal, remainder, delta))
            edge_cache[key] = choices
            counts['legal_exon_edges'] += len(choices)
        return edge_cache[key]

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
                counts['intron_states_dominated'] += _retain(active, key, normalized)
            counts['maximum_active_states'] = max(counts['maximum_active_states'], len(active))
        entries = []
        if a in starts:
            entries.append(_Live(best.score, best, (), ''))
        if a in acceptors:
            for _, row in sorted(active.items()):
                intron_role = 4 + (len(row.carry) - 1) % 3
                entries.append(_Live(row.score + scores.prefix[intron_role][a],
                                     row.preceding, row.chain, row.carry))
        for row in entries:
            first = not row.chain
            for b, terminal, carry, delta in legal_ends(a, row.carry, first):
                counts['exon_extensions'] += 1
                score = row.score + delta
                exon_count = len(row.chain) + 1
                if terminal:
                    counts['completed_extensions'] += 1
                    old = completed.get(b)
                    if old is None or score > old.score:
                        chain = row.chain + ((a, b),)
                        completed[b] = _Solution(score, row.preceding.chains + (chain,),
                            row.preceding.gains + (score - row.preceding.score,))
                elif exon_count < max_exons:
                    counts['donor_states_created'] += 1
                    table, key = pending[b], (exon_count, carry)
                    old = table.get(key)
                    counts['donor_states_dominated'] += old is not None
                    if old is None or score > old.score:
                        table[key] = _Live(score, row.preceding, row.chain + ((a, b),), carry)
                else:
                    counts['exon_cap_nonterminal_edges'] += 1
    return dict(chains=best.chains, gain=best.score, chain_gains=best.gains,
                counts=counts, reference_used=False)
