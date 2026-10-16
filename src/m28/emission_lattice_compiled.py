"""Historical A43 endpoint/trace-arena prototype with a stopped execution branch.

Retains the scalar indexed recurrence, actual carry, FP64 operation order and
first-arrival ties. Linked traces replace copied exon/gene tuples, not states.
No fastmath, phase-only merge, pruning or zero shortcut. This module loads no
model weights, FASTA files, reference metadata or observation labels.
The failed A43 R1 remains stopped; the current full-grammar raw interface uses
the separate checkpointed state-automaton Viterbi implementation instead.
"""
import numpy as np
from numba import njit, types
from numba.typed import Dict, List
import torch

from .candidates import Budget

_I = types.int64
_LIVE = types.Tuple((types.float64, types.float64, _I, _I))
_SLOTS = types.DictType(_I, _LIVE)
_SOLUTION = types.Tuple((types.float64, _I))
_EDGE = types.Tuple((_I, types.boolean, _I, types.float64))
_EDGES = types.ListType(_EDGE)
_KEY = types.UniTuple(_I, 3)
_EXON = types.UniTuple(_I, 3)
_GENE = types.Tuple((_I, _I, types.float64))
_EVENT = types.UniTuple(_I, 2)

# Exactly Python's lexical order: '', A, AA, AC, AG, AT, C, CA, ... TT.
_CARRIES = ('', 'A', 'AA', 'AC', 'AG', 'AT', 'C', 'CA', 'CC', 'CG', 'CT',
            'G', 'GA', 'GC', 'GG', 'GT', 'T', 'TA', 'TC', 'TG', 'TT')
_PHASES = np.asarray([len(s) for s in _CARRIES], dtype=np.int64)
_CHARS = np.full((21, 2), -1, dtype=np.int64)
for _i, _carry in enumerate(_CARRIES):
    for _j, _base in enumerate(_carry):
        _CHARS[_i, _j] = 'ACGT'.index(_base)


def emission_arrays(logits):
    """The original torch _Margins operations, without boxing every base as lists."""
    gain = (logits.double() - logits[:, :1].double()).cpu()
    prefix = torch.cat((torch.zeros(1, 15, dtype=torch.float64), gain.cumsum(0))).T.contiguous()
    positions = torch.arange(len(gain))
    coding = torch.stack([gain[positions, 1+(positions+offset) % 3] for offset in range(3)], dim=1)
    coding = torch.cat((torch.zeros(1, 3, dtype=torch.float64), coding.cumsum(0))).T.contiguous()
    return gain.numpy(), prefix.numpy(), coding.numpy()


@njit(fastmath=False)
def _stop(x, y, z):
    return x == 3 and ((y == 0 and (z == 0 or z == 2)) or (y == 2 and z == 0))


@njit(fastmath=False)
def _remainder(dna, bad, stops, a, b, carry, terminal):
    if bad[b] != bad[a]:
        return -1
    phase, size = _PHASES[carry], b-a
    if phase:
        need = 3-phase
        if size < need:
            return -1 if terminal else 2+5*_CHARS[carry, 0]+dna[a]
        if phase == 1:
            stop = _stop(_CHARS[carry, 0], dna[a], dna[a+1])
        else:
            stop = _stop(_CHARS[carry, 0], _CHARS[carry, 1], dna[a])
        if stop:
            return -1
        first = a+need
    else:
        first = a
    leftover = (phase+size) % 3
    number = stops[first % 3, b-2]-stops[first % 3, first] if b >= 2 else 0
    if terminal:
        if leftover or number != 1 or not _stop(dna[b-3], dna[b-2], dna[b-1]):
            return -1
    elif number > 0:
        return -1
    if not leftover:
        return 0
    return 1+5*dna[b-1] if leftover == 1 else 2+5*dna[b-2]+dna[b-1]


@njit(fastmath=False)
def _exon_gain(gain, prefix, coding, a, b, phase, first, terminal):
    size, offset = b-a, (phase-a) % 3
    value = coding[offset, b]-coding[offset, a]
    if first:
        value += (prefix[13, a+3]-prefix[13, a])-(coding[(-a) % 3, a+3]-coding[(-a) % 3, a])
    elif not (terminal and size == 3):
        value += gain[a, 10+phase]-gain[a, 1+phase]
    if terminal:
        terminal_offset = ((phase+size-3) % 3-(b-3)) % 3
        value += (prefix[14, b]-prefix[14, b-3])-(coding[terminal_offset, b]-coding[terminal_offset, b-3])
        if first and size == 3:
            value -= (prefix[13, b]-prefix[13, a])-(coding[(-a) % 3, b]-coding[(-a) % 3, a])
    elif not (first and size == 3):
        end_phase = (phase+size-1) % 3
        prior = 10+phase if not first and size == 1 else 1+end_phase
        value += gain[b-1, 7+end_phase]-gain[b-1, prior]
    return value


@njit(fastmath=False)
def _metadata(dna, gap):
    n = len(dna)
    bad = np.zeros(n+1, dtype=np.int64)
    stops = np.zeros((3, n+1), dtype=np.int64)
    next_bad = np.full(n+1, n, dtype=np.int64)
    next_stop = np.full((3, n+1), n+3, dtype=np.int64)
    starts = np.zeros(n+1, dtype=np.bool_)
    acceptors = np.zeros(n+1, dtype=np.bool_)
    donors = np.zeros(n+1, dtype=np.bool_)
    positions = np.zeros(n+1, dtype=np.bool_)
    ends = List.empty_list(_EVENT)
    for i in range(n):
        bad[i+1] = bad[i]+(dna[i] == 4)
        for frame in range(3):
            stops[frame, i+1] = stops[frame, i]
        if i < n-2:
            if dna[i] == 0 and dna[i+1] == 3 and dna[i+2] == 2:
                starts[i] = True; positions[i] = True
            if dna[i] == 0 and dna[i+1] == 2:
                acceptors[i+2] = True; positions[i+2] = True
            if _stop(dna[i], dna[i+1], dna[i+2]):
                stops[i % 3, i+1] += 1
                ends.append((i+3, 1)); positions[i+3] = True
        if 1 <= i < n-1 and dna[i] == 2 and dna[i+1] in (1, 3):
            donors[i] = True; ends.append((i, 0))
            if i+gap <= n:
                positions[i+gap] = True
    for i in range(n-1, -1, -1):
        next_bad[i] = i if dna[i] == 4 else next_bad[i+1]
        for frame in range(3):
            next_stop[frame, i] = next_stop[frame, i+1]
        if i < n-2 and _stop(dna[i], dna[i+1], dna[i+2]):
            next_stop[i % 3, i] = i
    return bad, stops, next_bad, next_stop, starts, acceptors, donors, positions, sorted(ends)


@njit(fastmath=False)
def _choices(dna, gain, prefix, coding, bad, stops, next_bad, next_stop, ends, a, carry, first):
    edges = List.empty_list(_EDGE)
    phase, n = _PHASES[carry], len(dna)
    limit = next_bad[a]
    entry = a+(3-phase) if phase else a
    if phase and entry <= n:
        if phase == 1:
            is_stop = _stop(_CHARS[carry, 0], dna[a], dna[a+1])
        else:
            is_stop = _stop(_CHARS[carry, 0], _CHARS[carry, 1], dna[a])
        if is_stop:
            limit = min(limit, entry-1)
    if entry <= n:
        limit = min(limit, next_stop[entry % 3, entry]+3)
    # Bisect the original (end, terminal) order; the bound removes only illegal edges.
    lo, hi = 0, len(ends)
    while lo < hi:
        mid = (lo+hi)//2
        if ends[mid][0] <= a:
            lo = mid+1
        else:
            hi = mid
    for j in range(lo, len(ends)):
        b, kind = ends[j]
        if b > limit:
            break
        terminal = kind == 1
        if (first or terminal) and b-a < 3:
            continue
        remainder = _remainder(dna, bad, stops, a, b, carry, terminal)
        if remainder >= 0:
            delta = _exon_gain(gain, prefix, coding, a, b, phase, first, terminal)
            edges.append((b, terminal, remainder, delta))
    return edges


@njit(fastmath=False)
def _solve(dna, gain, prefix, coding, cap, gap):
    bad, stops, next_bad, next_stop, starts, acceptors, donors, positions, ends = _metadata(dna, gap)
    pending = Dict.empty(_I, _SLOTS)
    active = Dict.empty(_I, _LIVE)
    completed = Dict.empty(_I, _SOLUTION)
    cache = Dict.empty(_KEY, _EDGES)
    exons, genes = List.empty_list(_EXON), List.empty_list(_GENE)
    counts = np.zeros(8, dtype=np.int64)
    best_score, best_gene = 0., -1
    for a in range(len(dna)+1):
        if not positions[a]:
            continue
        if a in completed:
            solution = completed[a]
            del completed[a]
            if solution[0] > best_score:
                best_score, best_gene = solution
        d = a-gap
        if d >= 1 and donors[d] and d in pending:
            table = pending[d]
            del pending[d]
            for key in sorted(table.keys()):
                row = table[key]; role = 4+(_PHASES[key % 21]-1) % 3
                value = (row[0]-prefix[role, d], row[1], row[2], row[3])
                if key in active:
                    counts[4] += 1
                    if value[0] > active[key][0]:
                        active[key] = value
                else:
                    active[key] = value
            counts[7] = max(counts[7], len(active))
        entries, entry_keys = List.empty_list(_LIVE), List.empty_list(_I)
        if starts[a]:
            entries.append((best_score, best_score, best_gene, -1)); entry_keys.append(0)
        if acceptors[a]:
            for key in sorted(active.keys()):
                row = active[key]; role = 4+(_PHASES[key % 21]-1) % 3
                entries.append((row[0]+prefix[role, a], row[1], row[2], row[3])); entry_keys.append(key)
        for i in range(len(entries)):
            row, key = entries[i], entry_keys[i]
            first, carry, exon_count = key == 0, key % 21, key//21+1
            edge_key = (a, carry, int(first))
            if edge_key not in cache:
                cache[edge_key] = _choices(dna, gain, prefix, coding, bad, stops, next_bad, next_stop,
                    ends, a, carry, first)
                counts[0] += len(cache[edge_key])
            for b, terminal, remainder, delta in cache[edge_key]:
                counts[1] += 1
                score = row[0]+delta
                if terminal:
                    counts[5] += 1
                    if b not in completed or score > completed[b][0]:
                        exons.append((a, b, row[3]))
                        genes.append((len(exons)-1, row[2], score-row[1]))
                        completed[b] = (score, len(genes)-1)
                elif exon_count < cap:
                    counts[2] += 1
                    if b not in pending:
                        pending[b] = Dict.empty(_I, _LIVE)
                    slot = exon_count*21+remainder; table = pending[b]
                    occupied = slot in table
                    counts[3] += occupied
                    if not occupied or score > table[slot][0]:
                        exons.append((a, b, row[3]))
                        table[slot] = (score, row[1], row[2], len(exons)-1)
                else:
                    counts[6] += 1
    return best_score, best_gene, exons, genes, counts


@torch.no_grad()
def decode_compiled_emission_lattice(sequence, logits, max_exons):
    sequence = sequence.upper()
    if not sequence or set(sequence)-set('ACGTN') or max_exons < 1:
        raise ValueError('Nonempty normalized DNA and positive explicit exon cap required')
    if (logits.device.type != 'cpu' or logits.dtype != torch.float32 or
            tuple(logits.shape) != (len(sequence), 15) or not bool(torch.isfinite(logits).all())):
        raise ValueError('Finite CPU FP32 valid-length fifteen-role rows required')
    dna = np.asarray(['ACGTN'.index(c) for c in sequence], dtype=np.int64)
    gain, prefix, coding = emission_arrays(logits)
    score, tail, exons, genes, counts = _solve(dna, gain, prefix, coding, max_exons, Budget().min_intron_bp)
    if not _solve.nopython_signatures:
        raise ValueError('The compiled preparation requires actual nopython execution, no fallback')
    chains, gains = [], []
    while tail >= 0:
        exon, previous, delta = genes[tail]; chain = []
        while exon >= 0:
            a, b, exon = exons[exon]; chain.append((a, b))
        chains.append(tuple(reversed(chain))); gains.append(delta); tail = previous
    names = ('legal_exon_edges', 'exon_extensions', 'donor_states_created', 'donor_states_dominated',
        'intron_states_dominated', 'completed_extensions', 'exon_cap_nonterminal_edges', 'maximum_active_states')
    return dict(chains=tuple(reversed(chains)), gain=score, chain_gains=tuple(reversed(gains)),
        counts={k: int(v) for k, v in zip(names, counts)}, reference_used=False)
