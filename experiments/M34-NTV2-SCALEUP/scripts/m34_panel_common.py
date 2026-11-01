"""P5 independent-panel inference: shared registry of the panel sequences (outline paper_outline_20261027 P5, scope 'core 10 species' of report 25 §9;
data prepared and frozen by PANEL-PREP, report 25 §3).  No model code here and no reference annotation is read.

`register(MC)` makes the frozen M29 helpers (`m29_common`: SEQ2SPECIES / LENGTH / read_sequence) aware of the panel main sequences, reading them from
data/panel/prepared/<key>/genome.fa.  It is called on top of whatever registration the caller already did (M34 validation pseudo-contigs, M30 contigs), is idempotent,
and must run in every process that reads a sequence (the inference main process, its spawned decode workers, each posterior task).
"""
import json
import os
from pathlib import Path

PREP = Path(os.environ.get('PANEL_PREP', Path.home() / 'coding-rna' / 'data' / 'panel' / 'prepared'))
ALL12 = ['dryas_octopetala', 'eschscholzia_californica', 'akebia_trifoliata', 'hakea_prostrata', 'chenopodium_album', 'zizania_latifolia', 'themeda_triandra',
         'juncus_effusus', 'rhynchospora_breviuscula', 'acorus_gramineus', 'dipteris_shenzhenensis', 'amborella_trichopoda']
CORE10 = [k for k in ALL12 if k not in ('chenopodium_album', 'dipteris_shenzhenensis')]       # report 25 §9 (2026-10-27): the fixed scope


def lengths(key):
    return {k: int(v) for k, v in json.load(open(PREP / key / 'lengths.json')).items()}


def contigs(key):
    """main sequences of a panel species, longest first (the order inference and the posterior queue use)"""
    L = lengths(key)
    return sorted(L, key=lambda c: (-L[c], c))


def total_bp(key):
    return sum(lengths(key).values())


_REG = {}


def register(MC):
    """Add the panel main sequences of all 12 species to `MC` (the m29_common module).  Idempotent."""
    if _REG.get('done'):
        return MC
    seq2key = {}
    for key in ALL12:
        p = PREP / key / 'lengths.json'
        if not p.exists():
            continue
        for c, n in lengths(key).items():
            if c in MC.LENGTH and int(MC.LENGTH[c]) != n:
                raise ValueError(f'{c}: length {n} in the panel, {MC.LENGTH[c]} already registered')
            seq2key[c] = key
            MC.SEQ2SPECIES[c] = key
            MC.LENGTH[c] = n
    orig_read = MC.read_sequence
    cache = {}

    def read_sequence(seq):
        if seq not in seq2key:
            return orig_read(seq)
        if seq not in cache:
            cache.clear()
            cache[seq] = MC.C.M.read_allowed(PREP / seq2key[seq] / 'genome.fa', {seq: MC.LENGTH[seq]})[seq].upper()
        return cache[seq]

    MC.read_sequence = read_sequence
    _REG.update(done=True, seq2key=seq2key)
    return MC
