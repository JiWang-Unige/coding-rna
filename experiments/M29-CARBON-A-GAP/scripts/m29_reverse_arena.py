"""M29 copy of src/m28/full_grammar_boundary_compiled._reverse_messages with ONE change: a single reusable alpha arena per call
(the original allocates a fresh (block+1) x states FP64 array for every 256-bp checkpoint block).  Use by monkeypatching
src.m28.full_grammar_boundary_compiled._reverse_messages = _reverse_messages_arena (frozen file untouched)."""
import math
import numpy as np
from numba import njit
from src.m28.full_grammar_boundary_compiled import _advance

@njit(cache=False,fastmath=False,parallel=False)
def _reverse_messages_arena(codes,gain,allowed,positions,checkpoints,upstream,destinations,roles):
    gradient = np.zeros_like(gain)
    adjoint = upstream.copy()
    maximum_block_bytes = 0
    longest = 0
    for block in range(len(positions)-1):
        longest = max(longest,positions[block+1]-positions[block])
    arena = np.empty((longest+1,len(adjoint)),dtype=np.float64)   # M29: one reusable alpha arena
    for block in range(len(positions)-2,-1,-1):
        a,b = positions[block],positions[block+1]
        alpha = arena[:b-a+1]
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


