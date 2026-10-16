"""Necessary ORF endpoint bounds, not a beam or a change to the CDS grammar.

Unused by the frozen A42 decoder. Future callers still run ORFIndex.extend
on every event in the retained interval; this bound is not a legality test.
"""
from bisect import bisect_left, bisect_right

from .candidates import STOP


def maximum_exon_end(orf, start, carry):
    """No legal extension can end beyond the first ambiguity/internal STOP.

The first complete in-exon STOP may itself be the terminal codon, so its END
is retained. An unfinished STOP prefix may leave through a donor before that
codon completes. A STOP completed by incoming actual carry excludes only
ends that would include its completion, preserving short middle exons.
"""
    limit = bisect_right(orf.bad, orf.bad[start]) - 1
    if carry:
        need = 3 - len(carry)
        if carry + orf.sequence[start:start + need] in STOP:
            return min(limit, start + need - 1)
        first = start + need
    else:
        first = start
    stops = orf.stops[first % 3]
    i = bisect_left(stops, first)
    return min(limit, stops[i] + 3) if i < len(stops) else limit


def exon_event_range(orf, events, start, carry):
    """Bounds in the ORIGINAL sorted (end, terminal) event order, both ties kept."""
    return (bisect_right(events, (start, True)),
            bisect_right(events, (maximum_exon_end(orf, start, carry), True)))
