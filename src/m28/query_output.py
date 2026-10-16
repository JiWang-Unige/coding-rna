"""Reference-free A query outputs in the existing chromosome assembly schema.

The caller explicitly supplies an existence-logit operating point. All completed
queries, including nonpositive gains and duplicate chains, stay in the free
record. Unterminated queries remain diagnostics, never assembled CDS genes.
This is an oriented-window bridge, not a raw-FASTA encoder or trained model.
"""
import math
from .query_decode import decode_query_chains


def query_window_candidates(model, features, sequence, window_index, window_id,
                            max_exons, existence_logit_threshold, valid_bases=None):
    threshold = float(existence_logit_threshold)
    if not math.isfinite(threshold):
        raise ValueError('A finite explicit existence-logit threshold is required')
    queries = decode_query_chains(model, features, sequence, max_exons, valid_bases)
    if any(not math.isfinite(row['query_logit']) for row in queries):
        raise ValueError('Nonfinite query existence logit')
    completed = [row for row in queries if row['terminated']]
    gains = [row['query_logit'] - threshold for row in completed]
    return {'window_index': window_index, 'window_id': window_id,
            'chains': [row['chain'] for row in completed], 'score': gains,
            'query_indices': [row['query_index'] for row in completed],
            'raw_queries': queries,
            'existence_logit_threshold': threshold,
            'score_kind': 'query_existence_logit_minus_threshold',
            'reference_used': False,
            'counts': {'queries': len(queries), 'completed_queries': len(completed),
                       'unterminated_queries': len(queries)-len(completed),
                       'completed_nonpositive_gain': sum(value <= 0 for value in gains)}}
