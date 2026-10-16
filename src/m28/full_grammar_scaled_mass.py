"""Blockwise gauge changes of the same all-state forward/backward log masses.

Messages are shifted by their finite background component, with compensated
offset accumulation. Shifts are restored at every query and at the global
boundary. No state, path, role or gene-set mass is selected or pruned.
"""
import numpy as np
from numba import njit

from .full_grammar_boundary_compiled import _add_mass, _advance


@njit(cache=False, fastmath=False, parallel=False)
def _recenter(message, background, offset, correction):
    shift = message[background]
    if shift != 0.:
        total = offset + shift
        if abs(offset) >= abs(shift):
            correction += (offset-total) + shift
        else:
            correction += (shift-total) + offset
        offset = total
        for state in range(len(message)):
            message[state] -= shift
    return offset, correction


@njit(cache=False, fastmath=False, parallel=False)
def _scaled_prefix_queries(codes, gain, positions, destinations, roles, background, block_bp=256):
    alpha = np.full(destinations.shape[1], -np.inf, dtype=np.float64)
    alpha[background] = 0.
    output = np.empty(len(positions), dtype=np.float64)
    cursor, offset, correction = 0, 0., 0.
    if positions[0] == 0:
        output[0] = 0.; cursor = 1
    allowed = np.ones(15, dtype=np.bool_)
    for i in range(len(codes)):
        alpha = _advance(alpha, codes[i], gain[i], allowed, destinations, roles)
        if (i+1) % block_bp == 0:
            offset, correction = _recenter(alpha, background, offset, correction)
        if cursor < len(positions) and i+1 == positions[cursor]:
            output[cursor] = offset + (correction + alpha[background]); cursor += 1
    return offset + (correction + alpha[background]), output


@njit(cache=False, fastmath=False, parallel=False)
def _scaled_suffix_queries(codes, gain, positions, destinations, roles, background, block_bp=256):
    beta = np.full(destinations.shape[1], -np.inf, dtype=np.float64)
    beta[background] = 0.
    output = np.empty(len(positions), dtype=np.float64)
    cursor, offset, correction = len(positions)-1, 0., 0.
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
        if (len(codes)-i) % block_bp == 0:
            offset, correction = _recenter(beta, background, offset, correction)
        if cursor >= 0 and i == positions[cursor]:
            output[cursor] = offset + (correction + beta[background]); cursor -= 1
    return offset + (correction + beta[background]), output
