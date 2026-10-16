"""Artificial preparation only: carry classes for exon lengths >= 3.

Used only by the isolated artificial batch decoder, not frozen/research
callers. This is a legality/role equivalence, not a safe FP64 history winner
selector; short middle exons keep their actual carry.
"""
from .candidates import STOP


def long_exon_phase(sequence, start, carry):
    """Return the shared phase, or None if no >=3-base exon can use carry.

Inputs are the grammar's valid actual ACGT carries of length 0, 1 or 2.
For nonempty carry the first codon must finish within the new exon. If it
contains ambiguity or is STOP it cannot be an internal coding codon.
Remaining DNA, motifs and terminal validity still require full ORF checks.
"""
    if carry:
        first = carry + sequence[start:start + 3 - len(carry)]
        if len(first) != 3 or set(first) - set('ACGT') or first in STOP:
            return None
    return len(carry)
