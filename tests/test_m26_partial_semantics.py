import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M26-SAME-SCOPE-MECHANISM'))
import partial_semantics as P


def test_parent_partial_is_not_CDS_partial_and_all_CDS_rows_count():
    lines = [
        'chr\tRefSeq\tgene\t1\t30\t.\t+\t.\tID=g;partial=true\n',
        'chr\tRefSeq\tmRNA\t1\t30\t.\t+\t.\tID=t;Parent=g;partial=true;start_range=.,1\n',
        'chr\tRefSeq\tCDS\t5\t13\t.\t+\t0\tID=c;Parent=t\n',
        'chr\tRefSeq\tCDS\t20\t28\t.\t+\t0\tID=c;Parent=t\n',
    ]
    g, t, c, counts = P.raw_features(lines, 'chr')
    assert P.S.E.is_partial(g['g']) and P.S.E.is_partial(t['t'])
    assert not any(P.S.E.is_partial(a) for a in c['t'])
    assert counts['CDS']['unflagged'] == 2
    lines[-1] = lines[-1].rstrip() + ';end_range=28,.\n'
    assert any(P.S.E.is_partial(a) for a in P.raw_features(lines, 'chr')[2]['t'])


def test_union_splits_CDS_from_span_without_double_counting():
    txs = [{'CDS': [(2,5,'0'), (8,11,'0')]}, {'CDS': [(4,9,'0')]}]
    labels = np.array([5,5,2,2,2,2,2,2,2,2,2,5])
    split, span, cds = P.split_coverage(txs, labels)
    assert split['span'][2] == 9 and split['CDS'][2] == 9
    assert sum(split['span_only']) == 0
    split, _, _ = P.split_coverage(txs[:1], labels)
    assert split['CDS'][2] == 6 and split['span_only'][2] == 3


def test_sequence_compatibility_is_strand_ordered_not_expression_evidence():
    tx = {'strand': '-', 'CDS': [(0,6,'0'), (10,13,'0')]}
    assert P.compatible(tx, 'ATGAAATAA')
    assert not P.compatible(tx, 'ATGTAATAA')
    assert not P.compatible(tx, 'ATGNNNTAA')
    tx['CDS'][-1] = (10,13,'1')
    assert not P.compatible(tx, 'ATGAAATAA')
