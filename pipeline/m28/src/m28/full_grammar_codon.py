"""Complete spliced-CDS codon energies on the unchanged gene automaton.

Every completed codon, including ATG and the terminal stop, is scored once.
Introns preserve actual unfinished coding bases. Two group means (sense and
stop) are removed from 64 raw coefficients; zero coefficients give the parent
model. This is a scoring experiment, not a new structural gene grammar.
"""
from functools import lru_cache
from itertools import product

import numpy as np
from numba import njit
import torch
from torch.autograd.function import once_differentiable

from .candidates import STOP
from .full_grammar_boundary_compiled import (
    ALPHABET, CompiledBoundaryMessage, _state_table, _add_mass)
from .full_grammar_observed_partition import InfeasibleObservation, allowed_roles_from_local_targets

CODONS = tuple(''.join(bases) for bases in product('ACGT', repeat=3))
CODON_INDEX = {codon: i for i, codon in enumerate(CODONS)}
STOP_MASK = torch.tensor([codon in STOP for codon in CODONS], dtype=torch.bool)


def centered_codon_scores(raw_scores):
    if (tuple(raw_scores.shape) != (64,) or raw_scores.device.type != 'cpu'
            or raw_scores.dtype not in (torch.float32, torch.float64)
            or not bool(torch.isfinite(raw_scores).all())):
        raise ValueError('64 finite CPU FP32/FP64 raw codon coefficients required')
    scores = raw_scores.double()
    return torch.where(STOP_MASK, scores - scores[STOP_MASK].mean(),
                       scores - scores[~STOP_MASK].mean())


@lru_cache(maxsize=4)
def _codon_edge_indices(cap):
    states, destinations, roles, _ = _state_table(cap)
    indices = np.full(destinations.shape, -1, dtype=np.int16)
    for b, base in enumerate(ALPHABET):
        for i, (stage, _, carry) in enumerate(states):
            for edge, target in enumerate(destinations[b, i]):
                if target < 0:
                    break
                role = roles[b, i, edge]
                codon = None
                if stage == 'S2' and role == 13:
                    codon = 'ATG'
                elif stage in ('TA', 'TG') and role == 14:
                    codon = stage + base
                elif stage in ('F', 'E', 'X') and len(carry) == 2 and role in (1, 2, 3, 7, 8, 9, 10, 11, 12):
                    codon = carry + base
                if codon is not None:
                    indices[b, i, edge] = CODON_INDEX[codon]
    return indices


@njit(cache=False, fastmath=False, parallel=False)
def _advance(alpha, base, gain, allowed, destinations, roles, indices, scores):
    following = np.full(len(alpha), -np.inf, dtype=np.float64)
    for state in range(len(alpha)):
        value = alpha[state]
        if value == -np.inf:
            continue
        for edge in range(3):
            target = destinations[base, state, edge]
            if target < 0:
                break
            role = roles[base, state, edge]
            if allowed[role]:
                index = indices[base, state, edge]
                bonus = scores[index] if index >= 0 else 0.
                following[target] = _add_mass(following[target], value + (gain[role] + bonus))
    return following


@njit(cache=False, fastmath=False, parallel=False)
def _forward(codes, gain, allowed, initial, positions, destinations, roles, indices, scores):
    checkpoints = np.empty((len(positions), len(initial)), dtype=np.float64)
    checkpoints[0] = initial
    alpha = initial.copy()
    cursor = 1
    for i in range(len(codes)):
        alpha = _advance(alpha, codes[i], gain[i], allowed[i], destinations, roles, indices, scores)
        if i + 1 == positions[cursor]:
            checkpoints[cursor] = alpha
            cursor += 1
    return alpha, checkpoints


@njit(cache=False, fastmath=False, parallel=False)
def _reverse(codes, gain, allowed, positions, checkpoints, upstream, destinations, roles, indices, scores):
    gradient = np.zeros_like(gain)
    codon_gradient = np.zeros(64, dtype=np.float64)
    adjoint = upstream.copy()
    maximum_bytes = 0
    for block in range(len(positions) - 2, -1, -1):
        a, b = positions[block], positions[block + 1]
        alpha = np.empty((b - a + 1, len(adjoint)), dtype=np.float64)
        maximum_bytes = max(maximum_bytes, alpha.nbytes)
        alpha[0] = checkpoints[block]
        for i in range(a, b):
            alpha[i - a + 1] = _advance(alpha[i - a], codes[i], gain[i], allowed[i],
                                       destinations, roles, indices, scores)
        for i in range(b - 1, a - 1, -1):
            previous = np.zeros_like(adjoint)
            for state in range(len(adjoint)):
                value = alpha[i - a, state]
                if value == -np.inf:
                    continue
                for edge in range(3):
                    target = destinations[codes[i], state, edge]
                    if target < 0:
                        break
                    role = roles[codes[i], state, edge]
                    g = adjoint[target]
                    if g == 0. or not allowed[i, role]:
                        continue
                    index = indices[codes[i], state, edge]
                    bonus = scores[index] if index >= 0 else 0.
                    contribution = g * np.exp(value + (gain[i, role] + bonus) - alpha[i - a + 1, target])
                    previous[state] += contribution
                    gradient[i, role] += contribution
                    if index >= 0:
                        codon_gradient[index] += contribution
            adjoint = previous
    return gradient, codon_gradient, adjoint, maximum_bytes


class _CodonBoundaryApply(torch.autograd.Function):
    @staticmethod
    def forward(ctx, logits, scores, initial, codes, allowed, cap, checkpoint_bp, statistics):
        states, destinations, roles, _ = _state_table(cap)
        indices = _codon_edge_indices(cap)
        gain = (logits.detach().double() - logits.detach()[:, :1].double()).numpy()
        coefficients = scores.detach().numpy().copy()
        mask = allowed.numpy().copy()
        positions = np.array([0] + list(range(checkpoint_bp, len(codes), checkpoint_bp)) + [len(codes)], dtype=np.int64)
        final, checkpoints = _forward(codes, gain, mask, initial.detach().numpy().copy(), positions,
                                      destinations, roles, indices, coefficients)
        ctx.codes, ctx.gain, ctx.allowed, ctx.positions = codes, gain, mask, positions
        ctx.checkpoints, ctx.final = checkpoints, final
        ctx.destinations, ctx.roles, ctx.indices, ctx.scores = destinations, roles, indices, coefficients
        ctx.dtype, ctx.statistics = logits.dtype, statistics
        statistics.update(state_count=len(states), checkpoint_rows=len(positions),
                          checkpoint_bytes=checkpoints.nbytes, checkpoint_bp=checkpoint_bp)
        return torch.from_numpy(final)

    @staticmethod
    @once_differentiable
    def backward(ctx, grad_output):
        upstream = grad_output.detach().numpy().copy()
        if (not np.isfinite(upstream).all() or np.any((ctx.final == -np.inf) & (upstream != 0))):
            raise ValueError('Only finite, supported boundary readouts have defined derivatives')
        gradient, codon_gradient, initial_gradient, peak = _reverse(ctx.codes, ctx.gain, ctx.allowed,
            ctx.positions, ctx.checkpoints, upstream, ctx.destinations, ctx.roles, ctx.indices, ctx.scores)
        gradient[:, 0] -= gradient.sum(axis=1)
        ctx.statistics['maximum_recomputed_alpha_bytes'] = peak
        return (torch.from_numpy(gradient).to(ctx.dtype), torch.from_numpy(codon_gradient),
                torch.from_numpy(initial_gradient), None, None, None, None, None)


def codon_boundary_forward(sequence, logits, raw_scores, message, allowed=None, checkpoint_bp=256):
    sequence = sequence.upper()
    n = len(sequence)
    states = _state_table(message.max_exons)[0]
    scores = centered_codon_scores(raw_scores)
    if (tuple(logits.shape) != (n, 15) or logits.device.type != 'cpu'
            or logits.dtype not in (torch.float32, torch.float64)
            or not bool(torch.isfinite(logits).all()) or checkpoint_bp < 1):
        raise ValueError('Aligned finite CPU potentials and positive checkpoint span required')
    initial = message.log_masses
    if (tuple(initial.shape) != (len(states),) or initial.dtype != torch.float64
            or initial.device.type != 'cpu' or bool(torch.isnan(initial).any())
            or bool((initial == float('inf')).any())):
        raise ValueError('Incoming FP64 masses must follow the fixed-cap state table')
    if allowed is None:
        allowed = torch.ones((n, 15), dtype=torch.bool)
    if tuple(allowed.shape) != (n, 15) or allowed.dtype != torch.bool or allowed.device.type != 'cpu':
        raise ValueError('Aligned CPU Boolean role observations required')
    if not n:
        return CompiledBoundaryMessage(message.max_exons, initial, {})
    codes = np.array([ALPHABET.find(base) if base in ALPHABET else 4 for base in sequence], dtype=np.int8)
    statistics = {}
    masses = _CodonBoundaryApply.apply(logits, scores, initial, codes, allowed,
                                      message.max_exons, checkpoint_bp, statistics)
    return CompiledBoundaryMessage(message.max_exons, masses, statistics)


def codon_full_grammar_log_partition(sequence, logits, max_exons, raw_scores,
                                     allowed=None, cuts=(), checkpoint_bp=256):
    n, cuts = len(sequence), tuple(cuts)
    if (not n or tuple(logits.shape) != (n, 15) or cuts != tuple(sorted(set(cuts)))
            or any(not 0 < x < n for x in cuts)):
        raise ValueError('Nonempty aligned DNA and ordered distinct interior cuts required')
    if allowed is not None and tuple(allowed.shape) != (n, 15):
        raise ValueError('Observation rows must retain their original geometry')
    states, _, _, background = _state_table(max_exons)
    initial = torch.full((len(states),), -float('inf'), dtype=torch.float64)
    initial[background] = 0.
    message, statistics = CompiledBoundaryMessage(max_exons, initial), []
    bounds = (0,) + cuts + (n,)
    for a, b in zip(bounds, bounds[1:]):
        message = codon_boundary_forward(sequence[a:b], logits[a:b], raw_scores, message,
            None if allowed is None else allowed[a:b], checkpoint_bp)
        statistics.append(message.statistics)
    value = message.log_masses[background]
    if not bool(torch.isfinite(value)):
        raise InfeasibleObservation('No complete gene set satisfies the supplied observations')
    return dict(log_partition=value, final_message=message, statistics=statistics,
                reference_used=False, enumerated_candidate_pool_used=False,
                Viterbi_winner_pruning=False, interior_edges_closed=False,
                global_complete_gene_grammar=True, first_order_CPU_only=True,
                structural_grammar_changed=False, codon_scores_centered_by_sense_and_stop=True)


def codon_coarse_channel_nll(sequence, logits, max_exons, raw_scores, core_start, core_end,
                             targets, epsilon, *, cuts=(), checkpoint_bp=256):
    n, length = core_end - core_start, len(sequence)
    if not 0 <= core_start < core_end <= length:
        raise ValueError('Original core must lie inside supplied fixed DNA')
    original = allowed_roles_from_local_targets(targets, n)
    fine = torch.as_tensor(targets['seg'], device='cpu')[:n]
    region = torch.as_tensor(targets['region'], device='cpu')[:n]
    epsilon = torch.as_tensor(epsilon, dtype=torch.float64, device='cpu')
    if epsilon.ndim or not bool(torch.isfinite(epsilon)) or not 0 < float(epsilon.detach()) < 1:
        raise ValueError('One observation error probability strictly between zero and one required')
    exact, coarse = fine >= 0, (fine < 0) & (region >= 0)
    allowed = torch.ones((length, 15), dtype=torch.bool)
    allowed[core_start:core_end][exact] = original[exact]
    logq = torch.zeros((length, 15), dtype=torch.float64)
    positions = torch.where(coarse)[0] + core_start
    logq[positions] = torch.where(original[coarse], torch.log1p(-epsilon), torch.log(epsilon / 2))
    denominator = codon_full_grammar_log_partition(sequence, logits, max_exons, raw_scores,
                                                  cuts=cuts, checkpoint_bp=checkpoint_bp)
    weighted = codon_full_grammar_log_partition(sequence, logits.double() + logq, max_exons,
        raw_scores, allowed=allowed, cuts=cuts, checkpoint_bp=checkpoint_bp)
    correction = logq[:, 0].sum()
    observed = weighted['log_partition'] + correction
    return dict(loss=denominator['log_partition'] - observed,
        log_partition_all=denominator['log_partition'], log_partition_observed=observed,
        weighted_relative_log_partition=weighted['log_partition'], background_correction=correction,
        fine_positions=int(exact.sum()), coarse_positions=int(coarse.sum()),
        fully_unknown_positions=n - int(exact.sum()) - int(coarse.sum()), epsilon=epsilon,
        statistics_all=denominator['statistics'], statistics_observed=weighted['statistics'],
        original_labels_changed=False, coarse_observation_rows_dropped=False,
        inference_grammar_changed=False, fit_error_probability_selected=False)
