"""Read-only complete-CDS event marginals under the unchanged gene-set grammar.

The supplied genes are queries, never the partition's candidate universe.
Forward/backward messages are conditioned only on oriented DNA and role scores.
No observation labels, temperature, fitting or posterior threshold enters here.
"""
import math

import numpy as np
from numba import njit
import torch

from .full_grammar_boundary_compiled import ALPHABET, _add_mass, _advance, _state_table
from .full_grammar_boundary_messages import BACKGROUND, _transitions
from .phase_evidence_chain import candidate_roles


@njit(cache=False, fastmath=False, parallel=False)
def _prefix_queries(codes, gain, positions, destinations, roles, background):
    alpha = np.full(destinations.shape[1], -np.inf, dtype=np.float64)
    alpha[background] = 0.
    output = np.empty(len(positions), dtype=np.float64)
    cursor = 0
    if positions[0] == 0:
        output[0] = 0.; cursor = 1
    allowed = np.ones(15, dtype=np.bool_)
    for i in range(len(codes)):
        alpha = _advance(alpha, codes[i], gain[i], allowed, destinations, roles)
        if cursor < len(positions) and i+1 == positions[cursor]:
            output[cursor] = alpha[background]; cursor += 1
    return alpha[background], output


@njit(cache=False, fastmath=False, parallel=False)
def _suffix_queries(codes, gain, positions, destinations, roles, background):
    beta = np.full(destinations.shape[1], -np.inf, dtype=np.float64)
    beta[background] = 0.
    output = np.empty(len(positions), dtype=np.float64)
    cursor = len(positions)-1
    if positions[cursor] == len(codes):
        output[cursor] = 0.; cursor -= 1
    for i in range(len(codes)-1, -1, -1):
        previous = np.full(len(beta), -np.inf, dtype=np.float64)
        for state in range(len(beta)):
            for edge in range(3):
                target = destinations[codes[i], state, edge]
                if target < 0:
                    break
                if beta[target] != -np.inf:
                    role = roles[codes[i], state, edge]
                    previous[state] = _add_mass(previous[state], gain[i, role]+beta[target])
        beta = previous
        if cursor >= 0 and i == positions[cursor]:
            output[cursor] = beta[background]; cursor -= 1
    return beta[background], output


def _gene_roles(sequence, chain, max_exons):
    if not chain or len(chain) > max_exons:
        raise ValueError('Each query must be a nonempty complete CDS within the fixed exon cap')
    start, end = chain[0][0], chain[-1][1]
    if not 0 <= start < end <= len(sequence):
        raise ValueError('Query coordinates must remain inside the supplied oriented DNA')
    local = tuple((a-start, b-start) for a, b in chain)
    result = candidate_roles(local, end-start).numpy()
    states = {BACKGROUND}
    for base, role in zip(sequence[start:end], result):
        states = {target for state in states for target, emitted in _transitions(state, base, max_exons)
                  if emitted == role}
        if not states:
            raise ValueError('Query is not a legal literal-DNA complete CDS in the unchanged grammar')
    if BACKGROUND not in states:
        raise ValueError('Query must close a complete gene, not an edge-truncated fragment')
    return result


def _log_probability(value):
    # Only cancellation-sized overshoot is rounded, never a substantive P>1.
    if not math.isfinite(value) or value > 1e-8:
        raise ValueError('Finite supported event mass cannot exceed the complete partition')
    return min(0., value)


@torch.no_grad()
def complete_gene_set_posteriors(sequence, logits, chains, max_exons):
    """Gene-event and exact-set probabilities for a supplied legal nonoverlapping set.

For g=[a,b), log P(g present)=alpha_B(a)+E(g)+beta_B(b)-logZ.
Alpha_B and beta_B close only the events' exterior sets, not the whole
partition at these query points. The empty set has energy zero. No
gradient or biological calibration is claimed; retain logs for ranking
when exp(log probability) underflows in floating point.
"""
    sequence = sequence.upper()
    n = len(sequence)
    if (not n or tuple(logits.shape) != (n, 15) or logits.device.type != 'cpu'
            or logits.dtype not in (torch.float32, torch.float64)
            or not bool(torch.isfinite(logits).all())):
        raise ValueError('Finite CPU FP32/FP64 roles aligned to nonempty oriented DNA required')
    chains = tuple(tuple((int(a), int(b)) for a, b in chain) for chain in chains)
    queries = []
    previous = 0
    for chain in chains:
        gene_roles = _gene_roles(sequence, chain, max_exons)
        if chain[0][0] < previous:
            raise ValueError('Supplied set must be ordered and nonoverlapping on this orientation')
        queries.append(gene_roles); previous = chain[-1][1]
    states, destinations, roles, background = _state_table(max_exons)
    codes = np.array([ALPHABET.find(base) if base in ALPHABET else 4 for base in sequence], dtype=np.int8)
    gain = (logits.detach().double()-logits.detach()[:, :1].double()).numpy()
    starts = np.array(sorted({0, *(chain[0][0] for chain in chains)}), dtype=np.int64)
    ends = np.array(sorted({n, *(chain[-1][1] for chain in chains)}), dtype=np.int64)
    log_z, prefix = _prefix_queries(codes, gain, starts, destinations, roles, background)
    reverse_z, suffix = _suffix_queries(codes, gain, ends, destinations, roles, background)
    if not math.isclose(log_z, reverse_z, rel_tol=1e-12, abs_tol=1e-8) or log_z < -1e-8:
        raise ValueError('Forward/backward total mass differs or omits the zero-energy empty set')
    prefix = dict(zip(starts.tolist(), prefix.tolist()))
    suffix = dict(zip(ends.tolist(), suffix.tolist()))
    energies, log_probabilities = [], []
    for chain, gene_roles in zip(chains, queries):
        a, b = chain[0][0], chain[-1][1]
        energy = float(gain[np.arange(a, b), gene_roles].sum())
        energies.append(energy)
        log_probabilities.append(_log_probability(prefix[a]+energy+suffix[b]-log_z))
    set_log_probability = _log_probability(sum(energies)-log_z)
    return dict(chains=chains, chain_gains=tuple(energies), log_partition=float(log_z),
        reverse_log_partition=float(reverse_z), gene_log_probabilities=tuple(log_probabilities),
        gene_probabilities=tuple(math.exp(v) for v in log_probabilities),
        queried_set_log_probability=set_log_probability, queried_set_probability=math.exp(set_log_probability),
        empty_set_log_probability=-float(log_z), empty_set_probability=math.exp(-float(log_z)),
        state_count=len(states), stored_boundary_scalars=len(starts)+len(ends),
        reference_used=False, observation_model_used=False, enumerated_candidate_pool_used=False,
        Viterbi_winner_pruning=False, optimizer_steps=0, posterior_threshold_applied=False,
        interpretation='model_structural_probability_not_empirically_calibrated_biological_correctness')
