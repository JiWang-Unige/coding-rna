"""Artificial finite-state preparation for exact complete-gene chunk composition.

Interior cuts carry scored messages, not uniform boundary priors or clipped
genes. Only global sequence edges are closed. No research caller uses this.
"""
from dataclasses import dataclass

import torch

from .candidates import STOP
from .full_grammar_observed_partition import InfeasibleObservation

BACKGROUND = ('B', 0, '')


@dataclass(frozen=True)
class BoundaryMessage:
    """FP64 log masses, including autograd history, under one fixed exon cap."""
    max_exons: int
    log_masses: dict


def _coding_carry(carry, base):
    if base not in 'ACGT':
        return None
    codon = carry + base
    if len(codon) == 3:
        return None if codon in STOP else ''
    return codon


def _transitions(state, base, cap):
    """Consume exactly one DNA base and emit its final, precedence-correct role."""
    stage, count, carry = state
    if stage == 'B':
        yield BACKGROUND, 0
        if base == 'A':
            yield ('S1', 1, ''), 13
    elif stage == 'S1':
        if base == 'T':
            yield ('S2', 1, ''), 13
    elif stage == 'S2':
        if base == 'G':
            yield ('F', 1, ''), 13
    elif stage in ('F', 'E', 'X'):
        phase = len(carry)
        # F follows the three start bases. Its length-three first exon can
        # end without replacing the final start role with a donor role.
        if stage == 'F' and count < cap and base == 'G':
            yield ('D1', count, carry), 4 + (phase - 1) % 3
        updated = _coding_carry(carry, base)
        if updated is not None:
            yield ('E', count, updated), (10 if stage == 'X' else 1) + phase
            if count < cap:
                # A length-one middle exon's donor overrides its acceptor.
                yield ('D0', count, updated), 7 + phase
        if not carry and base == 'T':
            # The three stop roles override a length-three terminal exon's
            # acceptor. STOP is consumed within one exon, not across splicing.
            yield ('T', count, ''), 14
    elif stage in ('D0', 'D1', 'I2', 'I3', 'AG'):
        role = 4 + (len(carry) - 1) % 3
        if stage == 'D0' and base == 'G':
            yield ('D1', count, carry), role
        elif stage == 'D1' and base in ('T', 'C'):
            yield ('I2', count, carry), role
        elif stage == 'I2':
            yield ('I3', count, carry), role
        elif stage == 'I3':
            yield state, role
            if base == 'A':
                yield ('AG', count, carry), role
        elif stage == 'AG' and base == 'G' and count < cap:
            # At least donor2 + interior1 + acceptor2 bases have been read.
            yield ('X', count + 1, carry), role
    elif stage == 'T':
        if base in ('A', 'G'):
            yield ('T' + base, count, ''), 14
    elif stage == 'TA':
        if base in ('A', 'G'):
            yield BACKGROUND, 14
    elif stage == 'TG':
        if base == 'A':
            yield BACKGROUND, 14
    else:
        raise ValueError('Unknown complete-gene boundary state')


def boundary_forward(sequence, logits, message, allowed=None):
    """Apply a chunk to an incoming scored message; never close its right edge.

States retain actual unfinished codon strings, gene exon counts, partial ATG
and STOP literals, and donor/maturation/acceptor stages. A missing state has
zero mass. Empty chunks are the identity. No annotation chooses a state.
"""
    sequence = sequence.upper()
    length = len(sequence)
    cap = message.max_exons
    if cap < 1:
        raise ValueError('A positive fixed exon cap is required')
    if (tuple(logits.shape) != (length, 15) or logits.device.type != 'cpu'
            or logits.dtype not in (torch.float32, torch.float64)
            or not bool(torch.isfinite(logits).all())):
        raise ValueError('Finite CPU FP32/FP64 potentials aligned to chunk DNA required')
    if allowed is not None and (tuple(allowed.shape) != (length, 15)
            or allowed.dtype != torch.bool or allowed.device.type != 'cpu'):
        raise ValueError('CPU Boolean allowed roles aligned to every chunk base required')
    for (stage, count, carry), value in message.log_masses.items():
        if (not isinstance(value, torch.Tensor) or value.ndim != 0
                or value.dtype != torch.float64 or value.device.type != 'cpu'
                or not bool(torch.isfinite(value))):
            raise ValueError('Incoming scored boundary masses must be finite CPU FP64 scalars')
        if ((stage == 'B' and (count != 0 or carry))
                or (stage != 'B' and not 1 <= count <= cap)
                or len(carry) > 2 or set(carry) - set('ACGT')):
            raise ValueError('Boundary gene count or actual coding carry is outside its fixed cap')
    gain = logits.double() - logits[:, :1].double()
    masses = message.log_masses
    for position, base in enumerate(sequence):
        following = {}
        for state, value in sorted(masses.items()):
            for target, role in _transitions(state, base, cap):
                if allowed is not None and not bool(allowed[position, role]):
                    continue
                proposal = value + gain[position, role]
                previous = following.get(target)
                following[target] = proposal if previous is None else torch.logaddexp(previous, proposal)
        masses = following
    return BoundaryMessage(cap, masses)


def boundary_full_grammar_log_partition(sequence, logits, max_exons, allowed=None, cuts=()):
    """The original complete, nonoverlapping gene-set universe, split anywhere.

Global background-only initialization and closure are explicit. Interior
messages come from preceding DNA/potentials, not a chosen open-boundary prior.
This native Torch implementation is currently artificial preparation only.
"""
    length = len(sequence)
    cuts = tuple(cuts)
    if not length or cuts != tuple(sorted(set(cuts))) or any(not 0 < x < length for x in cuts):
        raise ValueError('Nonempty global DNA and ordered distinct interior cuts required')
    if tuple(logits.shape) != (length, 15):
        raise ValueError('One 15-role row per global DNA base required')
    if allowed is not None and tuple(allowed.shape) != (length, 15):
        raise ValueError('Global allowed roles must retain every original DNA row')
    message = BoundaryMessage(max_exons, {BACKGROUND: logits.double().sum() * 0.})
    bounds = (0,) + cuts + (length,)
    for a, b in zip(bounds, bounds[1:]):
        local_allowed = None if allowed is None else allowed[a:b]
        message = boundary_forward(sequence[a:b], logits[a:b], message, local_allowed)
    value = message.log_masses.get(BACKGROUND)
    if value is None:
        raise InfeasibleObservation('No complete global gene set satisfies the supplied observations')
    return dict(log_partition=value, final_message=message, chunks=len(bounds)-1,
        reference_used=False, enumerated_candidate_pool_used=False, Viterbi_winner_pruning=False,
        interior_edges_closed=False, global_complete_gene_grammar=True)


def boundary_core_observation_nll(sequence, logits, max_exons, core_start, core_end, allowed_core, cuts=()):
    """Global logZ minus core-observed global logZ; outside roles stay unknown.

DNA/potentials outside the core remain real inputs and retain gradients.
This does not fetch context, select its extent, truncate genes, normalize the
loss, change target rows, or define a likelihood for missing outside DNA.
"""
    length = len(sequence)
    if not 0 <= core_start < core_end <= length:
        raise ValueError('A nonempty core inside the supplied global DNA required')
    if (tuple(allowed_core.shape) != (core_end-core_start, 15)
            or allowed_core.dtype != torch.bool or allowed_core.device.type != 'cpu'):
        raise ValueError('All original Boolean core role observations required')
    allowed = torch.ones((length, 15), dtype=torch.bool)
    allowed[core_start:core_end] = allowed_core
    all_mass = boundary_full_grammar_log_partition(sequence, logits, max_exons, cuts=cuts)
    observed = boundary_full_grammar_log_partition(sequence, logits, max_exons, allowed, cuts)
    return dict(loss=all_mass['log_partition'] - observed['log_partition'],
        log_partition_all=all_mass['log_partition'], log_partition_observed=observed['log_partition'],
        outside_roles_unknown=True, supplied_DNA_and_potentials_used=True,
        loss_normalized=False, context_extent_selected=False)
