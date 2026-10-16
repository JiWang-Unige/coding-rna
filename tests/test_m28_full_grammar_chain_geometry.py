import unittest

from src.m28.full_grammar_chain_geometry import PredictionIndex


def tx(parts, strand='+', phase='0', seqid='chr'):
    return dict(seqid=seqid, strand=strand, CDS=[(a, b, phase) for a, b in parts])


class ChainGeometryTests(unittest.TestCase):
    def test_transcriptional_start_and_stop_reverse_on_minus(self):
        for strand, category in (('+', 'transcriptional_start_only'), ('-', 'transcriptional_stop_only')):
            result = PredictionIndex([tx([(10, 20), (30, 45)], strand)]).classify(
                tx([(10, 20), (30, 50)], strand))
            self.assertEqual(result['category'], category)

    def test_exact_geometry_ignores_phase_but_internal_CDS_does_not(self):
        ref = tx([(10, 20), (30, 50)])
        self.assertEqual(PredictionIndex([tx([(10, 20), (30, 50)], phase='2')]).classify(ref)['category'],
                         'exact_primary')
        self.assertEqual(PredictionIndex([tx([(10, 22), (32, 50)])]).classify(ref)['category'],
                         'same_outer_boundaries_different_internal_CDS')

    def test_split_boundary_support_is_not_a_single_complete_chain(self):
        for strand in ('+', '-'):
            result = PredictionIndex([tx([(10, 20)], strand), tx([(40, 50)], strand)]).classify(
                tx([(10, 20), (30, 50)], strand))
            self.assertEqual(result, dict(category='both_boundaries_in_different_predictions',
                                          same_strand_span_overlaps=2))

    def test_half_open_spans_opposite_strand_and_overlapping_references(self):
        index = PredictionIndex([tx([(0, 10)]), tx([(50, 60)]), tx([(10, 50)], '-')])
        self.assertEqual(index.classify(tx([(10, 50)]))['category'], 'no_same_strand_span_overlap')
        index = PredictionIndex([tx([(15, 45)])])
        categories = [index.classify(r)['category'] for r in (tx([(10, 50)]), tx([(12, 48)]))]
        self.assertEqual(categories, ['same_strand_span_overlap_neither_boundary'] * 2)
        with self.assertRaisesRegex(ValueError, 'nonoverlapping'):
            PredictionIndex([tx([(10, 30)]), tx([(20, 40)])])


if __name__ == '__main__':
    unittest.main()
