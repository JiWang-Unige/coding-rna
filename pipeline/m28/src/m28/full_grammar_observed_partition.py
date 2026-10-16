"""Partial role observations over the same complete, nonoverlapping gene sets.

Artificial preparation only. Unknown rows admit all roles. No boundary repair,
R4 input access, loss normalization, optimizer, or training caller is supplied.
"""
from collections import defaultdict

import numpy as np
import torch
from torch.autograd.function import once_differentiable

from .candidates import Budget, ORFIndex, STOP
from .emission_lattice_compact import _CompactMargins
from .exon_end_bound import exon_event_range
from .full_grammar_partition_compiled_reverse import (
    _compiled_reverse, compiled_reverse_full_grammar_log_partition)
from .full_grammar_partition_reverse import _PartitionTape
from .labels import REGION_COLUMNS


class InfeasibleObservation(ValueError):
    """No complete gene set satisfies the supplied role observations."""


def allowed_roles_from_local_targets(targets, valid_bases):
    """Match local_loss's fine-first, coarse-second, otherwise-unknown policy.

    Padding must already be ignored. This does not manufacture background from
    positive chains, mask boundary genes, or select a policy for infeasible rows.
    """
    fine = torch.as_tensor(targets['seg'], device='cpu')
    region = torch.as_tensor(targets['region'], device='cpu')
    n = int(valid_bases)
    if (fine.ndim != 1 or region.shape != fine.shape or not 0 < n <= len(fine)
            or fine.is_floating_point() or region.is_floating_point()):
        raise ValueError('Original integer fine/coarse rows and valid length required')
    if bool((fine[n:] >= 0).any()) or bool((region[n:] >= 0).any()):
        raise ValueError('Padded target positions must already be ignored')
    fine, region = fine[:n], region[:n]
    if bool((fine >= 15).any()) or bool(((fine < 0) & (region >= 3)).any()):
        raise ValueError('Expected original 15-role and three-region labels')
    allowed = torch.ones((n, 15), dtype=torch.bool)
    known = fine >= 0
    allowed[known] = False
    rows = torch.arange(n)[known]
    allowed[rows, fine[known].long()] = True
    for label, columns in enumerate(REGION_COLUMNS):
        rows = torch.where((fine < 0) & (region == label))[0]
        allowed[rows] = False
        for column in columns:
            allowed[rows, column] = True
    return allowed


class _Observation:
    def __init__(self, allowed):
        bad = (~allowed).numpy().astype(np.int64)
        self.bad = bad
        self.prefix = np.concatenate((np.zeros((1, 15), dtype=np.int64),
                                      bad.cumsum(axis=0)), axis=0).T
        positions = np.arange(len(bad))
        coding = np.stack([bad[positions, 1 + (positions + offset) % 3]
                           for offset in range(3)], axis=1)
        self.coding = np.concatenate((np.zeros((1, 3), dtype=np.int64),
                                      coding.cumsum(axis=0)), axis=0).T

    def interval(self, a, b, role):
        return int(self.prefix[role, b] - self.prefix[role, a])

    def coding_interval(self, a, b, phase):
        offset = (phase - a) % 3
        return int(self.coding[offset, b] - self.coding[offset, a])

    def exon_allowed(self, a, b, phase, first, terminal):
        # Integer analogue of the original exon role precedence, not a score.
        size = b - a
        bad = self.coding_interval(a, b, phase)
        if first:
            bad += self.interval(a, a + 3, 13) - self.coding_interval(a, a + 3, 0)
        elif not (terminal and size == 3):
            bad += int(self.bad[a, 10 + phase]) - int(self.bad[a, 1 + phase])
        if terminal:
            bad += (self.interval(b - 3, b, 14)
                    - self.coding_interval(b - 3, b, (phase + size - 3) % 3))
            if first and size == 3:
                bad -= self.interval(a, b, 13) - self.coding_interval(a, b, 0)
        elif not (first and size == 3):
            end_phase = (phase + size - 1) % 3
            prior = 10 + phase if not first and size == 1 else 1 + end_phase
            bad += int(self.bad[b - 1, 7 + end_phase]) - int(self.bad[b - 1, prior])
        return bad == 0


def _observed_forward(sequence, logits, max_exons, allowed, counts):
    scores, orf, tape = _CompactMargins(logits), ORFIndex(sequence), _PartitionTape()
    observation = _Observation(allowed)
    length = len(sequence)
    starts = {i for i in range(length - 2) if sequence[i:i + 3] == 'ATG'}
    acceptors = {i + 2 for i in range(length - 2) if sequence[i:i + 2] == 'AG'}
    donors = {i for i in range(1, length - 1) if sequence[i:i + 2] in ('GT', 'GC')}
    stops = {i + 3 for i in range(length - 2) if sequence[i:i + 3] in STOP}
    events = sorted([(d, False) for d in donors] + [(b, True) for b in stops])
    gap = Budget().min_intron_bp
    maturities = {d + gap: d for d in donors if d + gap <= length}
    positions = sorted(starts | acceptors | stops | set(maturities))
    pending, active, completed, edge_cache = defaultdict(dict), {}, {}, {}
    prefix_mass, previous = 0, 0

    def add(table, key, value):
        old = table.get(key)
        table[key] = value if old is None else tape.combine(old, value)

    def shift(node, role, position, sign):
        delta = float(scores.prefix[role, position])
        return tape.node(2 if sign > 0 else -2, node, role * (length + 1) + position,
                         tape.values[node] + sign * delta)

    def legal_ends(a, carry, first):
        key = (a, carry, first)
        if key not in edge_cache:
            choices = []
            lo, hi = exon_event_range(orf, events, a, carry)
            for b, terminal in events[lo:hi]:
                if b <= a or ((first or terminal) and b - a < 3):
                    continue
                remainder = orf.extend(a, b, carry, terminal)
                if remainder is not None:
                    if not observation.exon_allowed(a, b, len(carry), first, terminal):
                        counts['rejected_observed_exon_edges'] += 1
                        continue
                    edge = len(tape.edges)
                    tape.edges.append((a, b, len(carry), first, terminal))
                    choices.append((b, terminal, remainder, edge,
                                    scores.exon_gain(a, b, len(carry), first, terminal)))
            edge_cache[key] = choices
            counts['legal_exon_edges'] += len(choices)
        return edge_cache[key]

    for a in positions:
        # Previously completed sets extend through background, not through a
        # required coding/intron observation. New completions at a are separate.
        if observation.interval(previous, a, 0):
            prefix_mass = None
        for key in list(active):
            role = 4 + (len(key[1]) - 1) % 3
            if observation.interval(previous, a, role):
                del active[key]
                counts['invalidated_observed_intron_states'] += 1
        mass = completed.pop(a, None)
        if mass is not None:
            prefix_mass = mass if prefix_mass is None else tape.combine(prefix_mass, mass)
        d = maturities.get(a)
        if d is not None:
            for key, value in sorted(pending.pop(d, {}).items()):
                _, carry = key
                role = 4 + (len(carry) - 1) % 3
                if observation.interval(d, a, role):
                    counts['rejected_observed_mature_states'] += 1
                    continue
                add(active, key, shift(value, role, d, -1))
                counts['mature_state_proposals'] += 1
            counts['maximum_active_states'] = max(counts['maximum_active_states'], len(active))
        entries = [(0, '', prefix_mass)] if a in starts and prefix_mass is not None else []
        if a in acceptors:
            for (count, carry), value in sorted(active.items()):
                role = 4 + (len(carry) - 1) % 3
                entries.append((count, carry, shift(value, role, a, 1)))
        for count, carry, value in entries:
            for b, terminal, remainder, edge, delta in legal_ends(a, carry, count == 0):
                counts['state_exon_extensions'] += 1
                if terminal or count + 1 < max_exons:
                    proposal = tape.node(1, value, edge, tape.values[value] + delta)
                    if terminal:
                        add(completed, b, proposal)
                        counts['terminal_state_proposals'] += 1
                    else:
                        add(pending[b], (count + 1, remainder), proposal)
                        counts['donor_state_proposals'] += 1
        previous = a
    if observation.interval(previous, length, 0):
        prefix_mass = None
    if prefix_mass is None:
        raise InfeasibleObservation('No complete gene set matches the observations; no boundary repair applied')
    counts['reverse_tape_nodes'] = len(tape.values)
    counts['reverse_tape_numeric_bytes'] = sum(len(x) * x.itemsize for x in
        (tape.values, tape.left, tape.right, tape.kind))
    return tape, prefix_mass


class _ObservedPartition(torch.autograd.Function):
    @staticmethod
    def forward(ctx, logits, sequence, max_exons, allowed, counts):
        tape, result = _observed_forward(sequence, logits, max_exons, allowed, counts)
        ctx.tape, ctx.result, ctx.length, ctx.dtype = tape, result, len(sequence), logits.dtype
        return torch.tensor(tape.values[result], dtype=torch.float64)

    @staticmethod
    @once_differentiable
    def backward(ctx, grad_output):
        gradient = torch.from_numpy(_compiled_reverse(ctx.tape, ctx.result, ctx.length)) * grad_output
        return gradient.to(dtype=ctx.dtype), None, None, None, None


def observed_full_grammar_log_partition(sequence, logits, max_exons, allowed_roles):
    """Log mass of exactly the complete gene sets allowed by explicit row masks."""
    sequence = sequence.upper()
    length = len(sequence)
    if max_exons < 1:
        raise ValueError('A positive explicit exon cap is required')
    if not length or logits.ndim != 2 or tuple(logits.shape) != (length, 15):
        raise ValueError('One valid-length 15-role row per oriented DNA base required')
    if (logits.device.type != 'cpu' or logits.dtype not in (torch.float32, torch.float64)
            or not bool(torch.isfinite(logits).all())):
        raise ValueError('Finite CPU FP32/FP64 local potentials required')
    if (allowed_roles.device.type != 'cpu' or allowed_roles.dtype != torch.bool
            or tuple(allowed_roles.shape) != (length, 15)):
        raise ValueError('Explicit valid-length CPU boolean 15-role observation rows required')
    counts = dict(legal_exon_edges=0, state_exon_extensions=0, donor_state_proposals=0,
        mature_state_proposals=0, terminal_state_proposals=0, maximum_active_states=0,
        rejected_observed_exon_edges=0, rejected_observed_mature_states=0,
        invalidated_observed_intron_states=0)
    value = _ObservedPartition.apply(logits, sequence, max_exons, allowed_roles, counts)
    return dict(log_partition=value, counts=counts, observation_constraints_used=True,
                enumerated_candidate_pool_used=False, Viterbi_winner_pruning=False)


def partial_observation_nll(sequence, logits, max_exons, allowed_roles):
    """Unnormalized -log P(observations); no coefficient or denominator chosen."""
    total = compiled_reverse_full_grammar_log_partition(sequence, logits, max_exons)
    observed = observed_full_grammar_log_partition(sequence, logits, max_exons, allowed_roles)
    return dict(loss=total['log_partition'] - observed['log_partition'],
                all_log_partition=total['log_partition'],
                observed_log_partition=observed['log_partition'],
                all_counts=total['counts'], observed_counts=observed['counts'])
