"""Fixed raw-coordinate neighbor view and differentiable base-owned emissions.

This preparation view does not inspect reference boundaries. It is not yet an
approved fit context or a chromosome-wide complete-gene inference contract.
"""
import torch
from .contig_emission_lattice import context_segments

WIDTH, STRIDE = 24576, 12288


def context_geometry(start, end, length):
    """Original core plus one original-stride neighbor on either side.

Clip neighbor STARTS to the raw chromosome's valid window range; deduplicate
at chromosome ends. All coordinates are genomic half-open, on either strand.
"""
    if not 0 <= start < end == min(start+WIDTH,length) <= length:
        raise ValueError('Original valid core window and raw chromosome length required')
    limit = max(0,length-WIDTH)
    starts = sorted({max(0,start-STRIDE),start,min(limit,start+STRIDE)})
    windows = tuple((a,min(a+WIDTH,length)) for a in starts)
    a,b = windows[0][0],windows[-1][1]
    local = tuple((x-a,y-a) for x,y in windows)
    return dict(genomic_context=[a,b],genomic_windows=windows,local_windows=local,
        base_owned_segments=context_segments(b-a,local))


def owned_context_logits(logits, geometry, strand):
    """Same old nearest-center ownership, retaining every selected derivative.

Inputs are valid-length locally oriented rows in increasing GENOMIC window
order. On minus, reverse piece order, not role columns or local row order.
"""
    windows,segments = geometry['local_windows'],geometry['base_owned_segments']
    if strand not in ('+','-') or len(logits) != len(windows):
        raise ValueError('One aligned potential tensor per window and strand required')
    pieces = []
    for value,(a,b),(_,left,right) in zip(logits,windows,segments):
        if tuple(value.shape) != (b-a,15):
            raise ValueError('Original valid-length fifteen-role rows required')
        pieces.append(value[left-a:right-a] if strand == '+' else value[b-right:b-left])
    return torch.cat(pieces if strand == '+' else pieces[::-1],dim=0)
