"""Isolated differentiable sum-product over the existing complete-CDS grammar.

Preparation only: no loss, reference constraint, training or inference caller.
This sums all nonoverlapping complete-gene sets including null; it does not
reuse Viterbi winners, an enumerated candidate pool or beam-pruned support.
"""
from collections import defaultdict

import torch

from .candidates import Budget, ORFIndex, STOP


class _TensorMargins:
    def __init__(self, logits):
        self.gain = logits.double() - logits[:, :1].double()
        self.prefix = torch.cat((self.gain.new_zeros((1, 15)), self.gain.cumsum(0))).T
        positions = torch.arange(len(logits), device=logits.device)
        coding = torch.stack([self.gain[positions, 1 + (positions + offset) % 3]
                              for offset in range(3)], dim=1)
        self.coding = torch.cat((coding.new_zeros((1, 3)), coding.cumsum(0))).T

    def interval(self, a, b, role):
        return self.prefix[role, b] - self.prefix[role, a]

    def coding_interval(self, a, b, phase):
        offset = (phase - a) % 3
        return self.coding[offset, b] - self.coding[offset, a]

    def exon_gain(self, a, b, phase, first, terminal):
        size = b - a
        value = self.coding_interval(a, b, phase)
        if first:
            value = value + (self.interval(a, a + 3, 13) - self.coding_interval(a, a + 3, 0))
        elif not (terminal and size == 3):
            value = value + (self.gain[a, 10 + phase] - self.gain[a, 1 + phase])
        if terminal:
            value = value + (self.interval(b - 3, b, 14) - self.coding_interval(
                b - 3, b, (phase + size - 3) % 3))
            if first and size == 3:
                value = value - (self.interval(a, b, 13) - self.coding_interval(a, b, 0))
        elif not (first and size == 3):
            end_phase = (phase + size - 1) % 3
            prior = 10 + phase if not first and size == 1 else 1 + end_phase
            value = value + (self.gain[b - 1, 7 + end_phase] - self.gain[b - 1, prior])
        return value


def full_grammar_log_partition(sequence, logits, max_exons):
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
                  terminal_state_proposals=0, maximum_active_states=0)

    def add(table, key, value):
        old = table.get(key)
        table[key] = value if old is None else torch.logaddexp(old, value)

    def legal_ends(a, carry, first):
        key = (a, carry, first)
        if key not in edge_cache:
            choices = []
            for b, terminal in events:
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
