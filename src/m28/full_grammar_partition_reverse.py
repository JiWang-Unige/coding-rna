"""Isolated first-order CPU logZ with a compact grammar-specific reverse tape.

Same complete path support and scalar recurrence; no per-proposal Torch graph.
Only artificial preparation callers use this branch. No observed likelihood.
"""
from array import array
from collections import defaultdict
import math

import numpy as np
import torch
from torch.autograd.function import once_differentiable

from .candidates import Budget, ORFIndex, STOP
from .emission_lattice_compact import _CompactMargins
from .exon_end_bound import exon_event_range


class _PartitionTape:
    """Only null, exon gain, signed intron prefix and logaddexp nodes."""
    def __init__(self):
        self.values = array('d', [0.])
        self.left = array('q', [0])
        self.right = array('q', [0])
        self.kind = array('b', [0])
        self.edges = []

    def node(self, kind, left, right, value):
        index = len(self.values)
        self.values.append(value)
        self.left.append(left)
        self.right.append(right)
        self.kind.append(kind)
        return index

    def combine(self, left, right):
        a, b = self.values[left], self.values[right]
        value = max(a, b) + math.log1p(math.exp(-abs(a - b)))
        return self.node(3, left, right, value)

    def reverse(self, result, length):
        adjoint = np.zeros(len(self.values), dtype=np.float64)
        adjoint[result] = 1.
        edge_gradient = np.zeros(len(self.edges), dtype=np.float64)
        prefix = np.zeros((15, length + 1), dtype=np.float64)
        coding = np.zeros((3, length + 1), dtype=np.float64)
        gain = np.zeros((length, 15), dtype=np.float64)
        flat_prefix = prefix.reshape(-1)
        for i in range(len(self.values) - 1, 0, -1):
            g = float(adjoint[i])
            if g == 0.:
                continue
            left, right, kind = self.left[i], self.right[i], self.kind[i]
            if kind == 3:
                adjoint[left] += g * math.exp(self.values[left] - self.values[i])
                adjoint[right] += g * math.exp(self.values[right] - self.values[i])
            else:
                adjoint[left] += g
                if kind == 1:
                    edge_gradient[right] += g
                else:
                    flat_prefix[right] += g if kind == 2 else -g

        def interval(a, b, role, g):
            prefix[role, b] += g
            prefix[role, a] -= g

        def coding_interval(a, b, phase, g):
            offset = (phase - a) % 3
            coding[offset, b] += g
            coding[offset, a] -= g

        for (a, b, phase, first, terminal), g in zip(self.edges, edge_gradient):
            if g == 0.:
                continue
            size = b - a
            coding_interval(a, b, phase, g)
            if first:
                interval(a, a + 3, 13, g)
                coding_interval(a, a + 3, 0, -g)
            elif not (terminal and size == 3):
                gain[a, 10 + phase] += g
                gain[a, 1 + phase] -= g
            if terminal:
                interval(b - 3, b, 14, g)
                coding_interval(b - 3, b, (phase + size - 3) % 3, -g)
                if first and size == 3:
                    interval(a, b, 13, -g)
                    coding_interval(a, b, 0, g)
            elif not (first and size == 3):
                end_phase = (phase + size - 1) % 3
                prior = 10 + phase if not first and size == 1 else 1 + end_phase
                gain[b - 1, 7 + end_phase] += g
                gain[b - 1, prior] -= g
        gain += np.cumsum(prefix[:, 1:][:, ::-1], axis=1)[:, ::-1].T
        coding_gradient = np.cumsum(coding[:, 1:][:, ::-1], axis=1)[:, ::-1]
        positions = np.arange(length)
        for offset in range(3):
            gain[positions, 1 + (positions + offset) % 3] += coding_gradient[offset]
        gain[:, 0] -= gain.sum(axis=1)
        return gain


def _forward(sequence, logits, max_exons, counts):
    scores, orf, tape = _CompactMargins(logits), ORFIndex(sequence), _PartitionTape()
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
    prefix_mass = 0

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
                    edge = len(tape.edges)
                    tape.edges.append((a, b, len(carry), first, terminal))
                    choices.append((b, terminal, remainder, edge,
                                    scores.exon_gain(a, b, len(carry), first, terminal)))
            edge_cache[key] = choices
            counts['legal_exon_edges'] += len(choices)
        return edge_cache[key]

    for a in positions:
        mass = completed.pop(a, None)
        if mass is not None:
            prefix_mass = tape.combine(prefix_mass, mass)
        d = maturities.get(a)
        if d is not None:
            for key, value in sorted(pending.pop(d, {}).items()):
                _, carry = key
                role = 4 + (len(carry) - 1) % 3
                add(active, key, shift(value, role, d, -1))
                counts['mature_state_proposals'] += 1
            counts['maximum_active_states'] = max(counts['maximum_active_states'], len(active))
        entries = [(0, '', prefix_mass)] if a in starts else []
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
    counts['reverse_tape_nodes'] = len(tape.values)
    counts['reverse_tape_numeric_bytes'] = sum(len(x) * x.itemsize for x in
        (tape.values, tape.left, tape.right, tape.kind))
    return tape, prefix_mass


class _FullGrammarReverse(torch.autograd.Function):
    @staticmethod
    def forward(ctx, logits, sequence, max_exons, counts):
        tape, result = _forward(sequence, logits, max_exons, counts)
        ctx.tape, ctx.result, ctx.length, ctx.dtype = tape, result, len(sequence), logits.dtype
        return torch.tensor(tape.values[result], dtype=torch.float64)

    @staticmethod
    @once_differentiable
    def backward(ctx, grad_output):
        gradient = torch.from_numpy(ctx.tape.reverse(ctx.result, ctx.length)) * grad_output
        return gradient.to(dtype=ctx.dtype), None, None, None


def reverse_full_grammar_log_partition(sequence, logits, max_exons):
    """Same all-path logZ; CPU FP32/FP64 and first-order gradients only."""
    sequence = sequence.upper()
    length = len(sequence)
    if max_exons < 1:
        raise ValueError('A positive explicit exon cap is required')
    if not length or logits.ndim != 2 or tuple(logits.shape) != (length, 15):
        raise ValueError('One valid-length 15-role row per oriented DNA base required')
    if (logits.device.type != 'cpu' or logits.dtype not in (torch.float32, torch.float64)
            or not bool(torch.isfinite(logits).all())):
        raise ValueError('Finite CPU FP32/FP64 local potentials required')
    counts = dict(legal_exon_edges=0, state_exon_extensions=0, donor_state_proposals=0,
                  mature_state_proposals=0, terminal_state_proposals=0, maximum_active_states=0)
    value = _FullGrammarReverse.apply(logits, sequence, max_exons, counts)
    return dict(log_partition=value, counts=counts, reference_used=False,
                enumerated_candidate_pool_used=False, Viterbi_winner_pruning=False)
