"""Label-builder semantics on the Brachypodium validation crop (small, CPU).  Skipped when the prepared data are absent."""
import numpy as np
import pytest

import m34_common as C

pytestmark = pytest.mark.skipif(not (C.PREP / 'brachypodium_distachyon' / 'splits.json').exists(), reason='prepared data not present')


@pytest.fixture(scope='module')
def lab():
    import m34_labels as LB
    sp = C.species_table()['brachypodium_distachyon']
    return LB.Labeler('brachypodium_distachyon', sp['val'], 0.5)


def first_window_with_unlabeled(lab):
    for t in sorted(lab.catalog.values(), key=lambda t: t['span'][0]):
        if t['eligible'] and t['id'] not in lab.labeled:
            c, s = t['seqid'], t['strand']
            a = max(0, t['span'][0] - 3000)
            return c, s, a, min(lab.lengths[c], a + 24576), t
    raise AssertionError


def test_selection_fraction_and_nesting(lab):
    import m34_labels as LB
    l25 = LB.Labeler('brachypodium_distachyon', lab.contigs, 0.25)
    l100 = LB.Labeler('brachypodium_distachyon', lab.contigs, 1.0)
    assert l25.labeled <= lab.labeled <= l100.labeled and len(l100.labeled) == l100.n_eligible
    assert abs(len(lab.labeled) / lab.n_eligible - 0.5) < 0.002


def test_U_vs_F_vs_FULL(lab):
    c, s, a, b, t = first_window_with_unlabeled(lab)
    U, _, _ = lab.targets(c, s, a, b, 'U')
    F, _, _ = lab.targets(c, s, a, b, 'F')
    FULL, _, _ = lab.targets(c, s, a, b, 'FULL')
    n = b - a
    x, y = (t['span'][0] - a, t['span'][1] - a) if s == '+' else (b - t['span'][1], b - t['span'][0])
    x, y = max(0, x), min(n, y)
    assert (U['seg'][x:y] < 0).any() and (F['seg'][x:y] >= 0).all()
    assert (FULL['seg'][x:y] > 0).any() and not (F['seg'][x:y] > 0).any()                 # F: unselected gene is background; FULL: it is labeled
    both = (U['seg'] >= 0) & (F['seg'] >= 0)
    assert (U['seg'][both] == F['seg'][both]).all()                                       # U differs from F only by turning background into unknown
    labeled_roles = F['seg'] > 0
    assert (U['seg'][labeled_roles] == F['seg'][labeled_roles]).all()


def test_full_equals_k1(lab):
    import m34_labels as LB
    l1 = LB.Labeler('brachypodium_distachyon', lab.contigs, 1.0)
    c, s, a, b, _ = first_window_with_unlabeled(lab)
    A, _, _ = l1.targets(c, s, a, b, 'U')
    B, _, _ = l1.targets(c, s, a, b, 'FULL')
    assert np.array_equal(A['seg'], B['seg']) and np.array_equal(A['region'], B['region'])


def test_window_row_schema(lab):
    row = lab.window_row(lab.contigs[0], '+', 0, 24576, 'val')
    for k in ('id', 'species', 'split', 'seqid', 'start', 'end', 'width', 'valid_bases', 'strand', 'positive_ids', 'overlapping_transcript_ids',
              'chain_unknown_intervals', 'fine_label_unknown_intervals', 'no_annotated_gene_either_strand', 'motifs', 'label_policy'):
        assert k in row
    assert row['id'] == f'brachypodium_distachyon|{lab.contigs[0]}|0|+'
