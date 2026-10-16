"""Reference-relative geometry of an immutable, nonoverlapping predicted gene set."""
from bisect import bisect_left, bisect_right
from collections import defaultdict

CATEGORIES = (
    'exact_primary',
    'same_outer_boundaries_different_internal_CDS',
    'both_boundaries_in_different_predictions',
    'transcriptional_start_only',
    'transcriptional_stop_only',
    'same_strand_span_overlap_neither_boundary',
    'no_same_strand_span_overlap',
)


def key(tx):
    return tx['seqid'], tx['strand'], tuple((a, b) for a, b, _ in tx['CDS'])


def boundaries(tx):
    a, b = tx['CDS'][0][0], tx['CDS'][-1][1]
    return (a, b) if tx['strand'] == '+' else (b, a)


class PredictionIndex:
    def __init__(self, predictions):
        self.keys = {key(tx) for tx in predictions}
        if len(self.keys) != len(predictions):
            raise ValueError('Duplicate predicted CDS geometry')
        groups = defaultdict(list)
        for tx in predictions:
            groups[tx['seqid'], tx['strand']].append(tx)
        self.groups = {}
        for group, rows in groups.items():
            rows.sort(key=lambda tx: tx['CDS'][0][0])
            starts = [tx['CDS'][0][0] for tx in rows]
            ends = [tx['CDS'][-1][1] for tx in rows]
            if any(b > c for b, c in zip(ends, starts[1:])):
                raise ValueError('Span index requires nonoverlapping same-strand predictions')
            self.groups[group] = rows, starts, ends

    def classify(self, reference):
        a, b = reference['CDS'][0][0], reference['CDS'][-1][1]
        rows, starts, ends = self.groups.get((reference['seqid'], reference['strand']), ([], [], []))
        overlapping = rows[bisect_right(ends, a):bisect_left(starts, b)]
        if key(reference) in self.keys:
            category = CATEGORIES[0]
        elif not overlapping:
            category = CATEGORIES[6]
        else:
            target = boundaries(reference)
            candidates = [boundaries(tx) for tx in overlapping]
            start = any(x == target[0] for x, _ in candidates)
            stop = any(y == target[1] for _, y in candidates)
            if target in candidates:
                category = CATEGORIES[1]
            elif start and stop:
                category = CATEGORIES[2]
            elif start:
                category = CATEGORIES[3]
            elif stop:
                category = CATEGORIES[4]
            else:
                category = CATEGORIES[5]
        return dict(category=category, same_strand_span_overlaps=len(overlapping))
