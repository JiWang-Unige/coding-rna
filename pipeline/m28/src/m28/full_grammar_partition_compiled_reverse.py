"""Isolated sequential nopython reverse loop; original forward untouched."""
import math

import numpy as np
import torch
from numba import njit
from torch.autograd.function import once_differentiable

from .full_grammar_partition_reverse import _FullGrammarReverse


@njit(cache=False, fastmath=False, parallel=False)
def _accumulate_reverse_nodes(values, left_nodes, right_nodes, kinds, adjoint, edge_gradient, flat_prefix):
    for i in range(len(values) - 1, 0, -1):
        g = float(adjoint[i])
        if g == 0.:
            continue
        left, right, kind = left_nodes[i], right_nodes[i], kinds[i]
        if kind == 3:
            adjoint[left] += g * math.exp(values[left] - values[i])
            adjoint[right] += g * math.exp(values[right] - values[i])
        else:
            adjoint[left] += g
            if kind == 1:
                edge_gradient[right] += g
            else:
                flat_prefix[right] += g if kind == 2 else -g

def _compiled_reverse(tape, result, length):
    adjoint = np.zeros(len(tape.values), dtype=np.float64)
    adjoint[result] = 1.
    edge_gradient = np.zeros(len(tape.edges), dtype=np.float64)
    prefix = np.zeros((15, length + 1), dtype=np.float64)
    coding = np.zeros((3, length + 1), dtype=np.float64)
    gain = np.zeros((length, 15), dtype=np.float64)
    flat_prefix = prefix.reshape(-1)
    _accumulate_reverse_nodes(np.frombuffer(tape.values, dtype=np.float64),
        np.frombuffer(tape.left, dtype=np.int64), np.frombuffer(tape.right, dtype=np.int64),
        np.frombuffer(tape.kind, dtype=np.int8), adjoint, edge_gradient, flat_prefix)

    def interval(a, b, role, g):
        prefix[role, b] += g
        prefix[role, a] -= g

    def coding_interval(a, b, phase, g):
        offset = (phase - a) % 3
        coding[offset, b] += g
        coding[offset, a] -= g

    for (a, b, phase, first, terminal), g in zip(tape.edges, edge_gradient):
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


class _CompiledFullGrammarReverse(_FullGrammarReverse):
    @staticmethod
    @once_differentiable
    def backward(ctx, grad_output):
        gradient = torch.from_numpy(_compiled_reverse(ctx.tape, ctx.result, ctx.length)) * grad_output
        return gradient.to(dtype=ctx.dtype), None, None, None


def compiled_reverse_full_grammar_log_partition(sequence, logits, max_exons):
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
    value = _CompiledFullGrammarReverse.apply(logits, sequence, max_exons, counts)
    return dict(log_partition=value, counts=counts, reference_used=False,
                enumerated_candidate_pool_used=False, Viterbi_winner_pruning=False)
