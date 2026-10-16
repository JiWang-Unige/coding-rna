"""Checkpointed nopython application of the unchanged complete-gene automaton.

Dense boundary messages preserve all reachable states and first-order CPU
gradients through incoming messages. No research caller uses this preparation.
"""
from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from numba import njit
import torch
from torch.autograd.function import once_differentiable

from .full_grammar_boundary_messages import BACKGROUND, _transitions
from .full_grammar_observed_partition import InfeasibleObservation

ALPHABET = 'ACGTN'


@lru_cache(maxsize=4)
def _state_table(cap):
    """Structural closure over the original transitions, not path pruning."""
    if cap < 1:
        raise ValueError('A positive explicit exon cap is required')
    seen, pending = {BACKGROUND}, [BACKGROUND]
    while pending:
        state = pending.pop()
        for base in ALPHABET:
            for target, _ in _transitions(state, base, cap):
                if target not in seen:
                    seen.add(target); pending.append(target)
    states = tuple(sorted(seen))
    index = {state:i for i,state in enumerate(states)}
    destinations = np.full((5,len(states),3),-1,dtype=np.int32)
    roles = np.zeros_like(destinations,dtype=np.int8)
    for b,base in enumerate(ALPHABET):
        for i,state in enumerate(states):
            for edge,(target,role) in enumerate(_transitions(state,base,cap)):
                destinations[b,i,edge] = index[target]
                roles[b,i,edge] = role
    return states,destinations,roles,index[BACKGROUND]


@njit(cache=False,fastmath=False,parallel=False)
def _add_mass(a,b):
    if a == -np.inf:
        return b
    if b == -np.inf:
        return a
    return max(a,b) + math.log1p(math.exp(-abs(a-b)))


@njit(cache=False,fastmath=False,parallel=False)
def _advance(alpha,base,gain,allowed,destinations,roles):
    following = np.full(len(alpha),-np.inf,dtype=np.float64)
    for state in range(len(alpha)):
        value = alpha[state]
        if value == -np.inf:
            continue
        for edge in range(3):
            target = destinations[base,state,edge]
            if target < 0:
                break
            role = roles[base,state,edge]
            if allowed[role]:
                following[target] = _add_mass(following[target],value+gain[role])
    return following


@njit(cache=False,fastmath=False,parallel=False)
def _forward_messages(codes,gain,allowed,initial,positions,destinations,roles):
    checkpoints = np.empty((len(positions),len(initial)),dtype=np.float64)
    checkpoints[0] = initial
    alpha = initial.copy()
    checkpoint = 1
    for i in range(len(codes)):
        alpha = _advance(alpha,codes[i],gain[i],allowed[i],destinations,roles)
        if i+1 == positions[checkpoint]:
            checkpoints[checkpoint] = alpha
            checkpoint += 1
    return alpha,checkpoints


@njit(cache=False,fastmath=False,parallel=False)
def _reverse_messages(codes,gain,allowed,positions,checkpoints,upstream,destinations,roles):
    gradient = np.zeros_like(gain)
    adjoint = upstream.copy()
    maximum_block_bytes = 0
    for block in range(len(positions)-2,-1,-1):
        a,b = positions[block],positions[block+1]
        alpha = np.empty((b-a+1,len(adjoint)),dtype=np.float64)
        maximum_block_bytes = max(maximum_block_bytes,alpha.nbytes)
        alpha[0] = checkpoints[block]
        for i in range(a,b):
            alpha[i-a+1] = _advance(alpha[i-a],codes[i],gain[i],allowed[i],destinations,roles)
        for i in range(b-1,a-1,-1):
            previous = np.zeros_like(adjoint)
            for state in range(len(adjoint)):
                value = alpha[i-a,state]
                if value == -np.inf:
                    continue
                for edge in range(3):
                    target = destinations[codes[i],state,edge]
                    if target < 0:
                        break
                    role = roles[codes[i],state,edge]
                    g = adjoint[target]
                    if g == 0. or not allowed[i,role]:
                        continue
                    contribution = g*math.exp(value+gain[i,role]-alpha[i-a+1,target])
                    previous[state] += contribution
                    gradient[i,role] += contribution
            adjoint = previous
    return gradient,adjoint,maximum_block_bytes


@dataclass(frozen=True)
class CompiledBoundaryMessage:
    """FP64 tensor columns follow the original state table under a fixed cap."""
    max_exons: int
    log_masses: torch.Tensor
    statistics: dict = None


class _CompiledBoundaryApply(torch.autograd.Function):
    @staticmethod
    def forward(ctx,logits,initial,codes,allowed,cap,checkpoint_bp,statistics):
        states,destinations,roles,_ = _state_table(cap)
        gain = (logits.detach().double()-logits.detach()[:,:1].double()).numpy()
        # Retain the forward inputs used by recomputation, even if a caller
        # reuses an incoming message or target buffer before backward.
        initial_values = initial.detach().numpy().copy()
        mask = allowed.numpy().copy()
        positions = np.array([0]+list(range(checkpoint_bp,len(codes),checkpoint_bp))+[len(codes)],dtype=np.int64)
        final,checkpoints = _forward_messages(codes,gain,mask,initial_values,positions,destinations,roles)
        ctx.codes,ctx.gain,ctx.allowed,ctx.positions = codes,gain,mask,positions
        ctx.checkpoints,ctx.final,ctx.destinations,ctx.roles = checkpoints,final,destinations,roles
        ctx.dtype,ctx.statistics = logits.dtype,statistics
        statistics.update(state_count=len(states),checkpoint_rows=len(positions),
            checkpoint_bytes=checkpoints.nbytes,checkpoint_bp=checkpoint_bp)
        return torch.from_numpy(final)

    @staticmethod
    @once_differentiable
    def backward(ctx,grad_output):
        upstream = grad_output.detach().numpy().copy()
        if (not np.isfinite(upstream).all()
                or np.any((ctx.final == -np.inf)&(upstream != 0))):
            raise ValueError('Only finite, supported boundary-state readouts have defined derivatives')
        gradient,initial_gradient,peak = _reverse_messages(ctx.codes,ctx.gain,ctx.allowed,
            ctx.positions,ctx.checkpoints,upstream,ctx.destinations,ctx.roles)
        # Every role score is its logit minus the background logit.
        gradient[:,0] -= gradient.sum(axis=1)
        ctx.statistics['maximum_recomputed_alpha_bytes'] = peak
        return torch.from_numpy(gradient).to(ctx.dtype),torch.from_numpy(initial_gradient),None,None,None,None,None


def compiled_boundary_forward(sequence,logits,message,allowed=None,checkpoint_bp=256):
    """Apply oriented DNA without closing the chunk or detaching its message."""
    sequence = sequence.upper()
    n = len(sequence)
    states,_,_,_ = _state_table(message.max_exons)
    if (tuple(logits.shape) != (n,15) or logits.device.type != 'cpu'
            or logits.dtype not in (torch.float32,torch.float64)
            or not bool(torch.isfinite(logits).all()) or checkpoint_bp < 1):
        raise ValueError('Finite aligned CPU FP32/FP64 potentials and positive checkpoint span required')
    initial = message.log_masses
    if (tuple(initial.shape) != (len(states),) or initial.dtype != torch.float64
            or initial.device.type != 'cpu' or bool(torch.isnan(initial).any())
            or bool((initial == float('inf')).any())):
        raise ValueError('Incoming FP64 masses must follow the fixed-cap state table; absence is minus infinity')
    if allowed is None:
        allowed = torch.ones((n,15),dtype=torch.bool)
    if (tuple(allowed.shape) != (n,15) or allowed.dtype != torch.bool or allowed.device.type != 'cpu'):
        raise ValueError('Boolean allowed roles aligned to every DNA row required')
    if not n:
        return CompiledBoundaryMessage(message.max_exons,initial,{})
    codes = np.array([ALPHABET.find(base) if base in ALPHABET else 4 for base in sequence],dtype=np.int8)
    statistics = {}
    masses = _CompiledBoundaryApply.apply(logits,initial,codes,allowed,message.max_exons,checkpoint_bp,statistics)
    return CompiledBoundaryMessage(message.max_exons,masses,statistics)


def compiled_boundary_full_grammar_log_partition(sequence,logits,max_exons,allowed=None,cuts=(),checkpoint_bp=256):
    """Same global complete gene sets; internal messages remain differentiable."""
    n = len(sequence); cuts = tuple(cuts)
    if (not n or tuple(logits.shape) != (n,15) or cuts != tuple(sorted(set(cuts)))
            or any(not 0 < x < n for x in cuts)):
        raise ValueError('Aligned nonempty global DNA and ordered distinct interior cuts required')
    if allowed is not None and tuple(allowed.shape) != (n,15):
        raise ValueError('Global observation rows must retain their original geometry')
    states,_,_,background = _state_table(max_exons)
    initial = torch.full((len(states),),-float('inf'),dtype=torch.float64)
    initial[background] = 0.
    message = CompiledBoundaryMessage(max_exons,initial)
    statistics = []
    bounds = (0,)+cuts+(n,)
    for a,b in zip(bounds,bounds[1:]):
        message = compiled_boundary_forward(sequence[a:b],logits[a:b],message,
            None if allowed is None else allowed[a:b],checkpoint_bp)
        statistics.append(message.statistics)
    value = message.log_masses[background]
    if not bool(torch.isfinite(value)):
        raise InfeasibleObservation('No complete global gene set satisfies the supplied observations')
    return dict(log_partition=value,final_message=message,statistics=statistics,
        reference_used=False,enumerated_candidate_pool_used=False,Viterbi_winner_pruning=False,
        interior_edges_closed=False,global_complete_gene_grammar=True,first_order_CPU_only=True)


def compiled_boundary_core_observation_nll(sequence,logits,max_exons,core_start,core_end,allowed_core,cuts=(),checkpoint_bp=256):
    """Unnormalized global/core-observed probability with unchanged core rows."""
    n = len(sequence)
    if (not 0 <= core_start < core_end <= n or tuple(allowed_core.shape) != (core_end-core_start,15)
            or allowed_core.dtype != torch.bool or allowed_core.device.type != 'cpu'):
        raise ValueError('All original Boolean core observations within supplied DNA required')
    allowed = torch.ones((n,15),dtype=torch.bool); allowed[core_start:core_end] = allowed_core
    all_mass = compiled_boundary_full_grammar_log_partition(sequence,logits,max_exons,cuts=cuts,checkpoint_bp=checkpoint_bp)
    observed = compiled_boundary_full_grammar_log_partition(sequence,logits,max_exons,allowed,cuts,checkpoint_bp)
    return dict(loss=all_mass['log_partition']-observed['log_partition'],
        log_partition_all=all_mass['log_partition'],log_partition_observed=observed['log_partition'],
        statistics_all=all_mass['statistics'],statistics_observed=observed['statistics'],
        loss_normalized=False,context_extent_selected=False,outside_roles_unknown=True)
