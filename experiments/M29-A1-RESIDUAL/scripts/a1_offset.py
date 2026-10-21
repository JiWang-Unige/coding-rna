"""Fixed Carbon-A offset (parity-gate construction, `m29_gate_parity.build_logits` variant A_full) in compact form.

Per oriented strand a float32 array [5, L]:  0: d = logit(p) (roles 1,2,3,7..14)   1: enter_gene add (role 13, at ATG)
2: leave_gene add (role 14, at the last base of TAA/TAG/TGA)   3: enter_intron add (roles 7,8,9: tr[e+1,2] (+GC prior) where seq[e+1:e+3] in GT/GC)
4: leave_intron add (roles 10,11,12: tr[i,3]).   `expand` reproduces build_logits bit for bit (tested).
"""
import numpy as np
from pathlib import Path
import a1_common as K

STOPS = ('TAA', 'TAG', 'TGA')
COMP = str.maketrans('ACGTN', 'TGCAN')
D_ROLES = (1, 2, 3, 7, 8, 9, 10, 11, 12, 13, 14)


def rc(s):
    return s.translate(COMP)[::-1]


def compact(seqstr, p16, tt, min_prob, max_prob, gc_prior):
    """seqstr: oriented uppercase sequence; p16: oriented CDS probability (float16/32); tt: [L,4] transition log-odds tracks."""
    L = len(seqstr)
    p = np.clip(p16.astype(np.float64), min_prob, max_prob)
    out = np.zeros((5, L), dtype=np.float32)
    out[0] = (np.log(p) - np.log1p(-p)).astype(np.float32)
    arr = np.frombuffer(seqstr.encode('ascii'), dtype=np.uint8)
    A, T, C, G = (arr == ord(c) for c in 'ATCG')
    atg = np.zeros(L, bool)
    atg[:-2] = A[:-2] & T[1:-1] & G[2:]
    out[1, atg] = tt[atg, 0]
    stop = np.zeros(L, bool)
    for s in STOPS:
        m = np.ones(L - 2, bool)
        for k, ch in enumerate(s):
            m &= arr[k:L - 2 + k] == ord(ch)
        stop[2:] |= m
    out[2, stop] = tt[stop, 1]
    ent = np.zeros(L, np.float32)
    ent[:-1] = tt[1:, 2]
    ok_gt = np.zeros(L, bool); ok_gc = np.zeros(L, bool)
    ok_gt[:-2] = G[1:-1] & T[2:]
    ok_gc[:-2] = G[1:-1] & C[2:]
    ent[ok_gc] += np.float32(gc_prior)
    ok = ok_gt | ok_gc
    out[3, ok] = ent[ok]
    out[4] = tt[:, 3]
    return out


def expand(c, a=0, b=None):
    """[5, L] compact -> float32 [b-a, 15] role potentials, identical to m29_gate_parity.build_logits rows a:b."""
    b = c.shape[1] if b is None else b
    c = np.asarray(c[:, a:b])
    row = np.zeros((b - a, 15), dtype=np.float32)
    for r in D_ROLES:
        row[:, r] = c[0]
    row[:, 13] += c[1]
    row[:, 14] += c[2]
    for r in (7, 8, 9):
        row[:, r] += c[3]
    for r in (10, 11, 12):
        row[:, r] += c[4]
    return row


def context_offset(tracks, strand, a, b, length):
    """Oriented offset for the genomic context [a,b) on `strand`; tracks = {'+': compact, '-': compact} over the full contig."""
    if strand == '+':
        return expand(tracks['+'], a, b)
    return expand(tracks['-'], length - b, length - a)


def offset_path(seq, strand):
    return K.A1 / 'offset' / f'{seq}_{strand}.npy'


def load_tracks(seq, mmap=True):
    return {s: np.load(offset_path(seq, s), mmap_mode='r' if mmap else None) for s in '+-'}
