import math
import unittest

from src.m28.gene_posterior_ranking import oriented_chain, ranking_summary, probability_diagnostics


class RankingTests(unittest.TestCase):
    def test_coordinate_roundtrip_and_short_terminal_exons(self):
        chain = ((3, 7), (15, 16), (22, 25))
        reverse = oriented_chain(chain, 31, '-')
        self.assertEqual(reverse, ((6, 9), (15, 16), (24, 28)))
        self.assertEqual(oriented_chain(reverse, 31, '-'), chain)
        self.assertEqual(oriented_chain(chain, 31, '+'), chain)

    def test_exact_ties_and_full_primary_recall_denominator(self):
        rows = [dict(score=s, match=y) for s, y in ((3., True), (2., True), (2., False), (0., False))]
        result = ranking_summary(rows, 'score', 'match', 10)
        self.assertAlmostEqual(result['prediction_pool_AUROC'], .875)
        self.assertAlmostEqual(result['prediction_pool_average_precision'], (1+2/3)/2)
        half = result['retained_fraction_diagnostics']['0.5']
        self.assertEqual(half['retained_predictions'], 3)
        self.assertEqual(half['recall_of_all_primary_references'], .2)
        self.assertEqual(half['actual_fraction'], .75)
        # A permutation inside a score tie must not affect the result.
        self.assertEqual(result, ranking_summary(rows[::-1], 'score', 'match', 10))

    def test_log_underflow_preserves_ranking_and_single_class_is_undefined(self):
        rows = [dict(log_probability=s, match=y) for s, y in ((-1000., True), (-2000., False))]
        self.assertEqual(math.exp(-1000.), math.exp(-2000.))
        result = ranking_summary(rows, 'log_probability', 'match')
        self.assertEqual(result['prediction_pool_AUROC'], 1.)
        all_negative = ranking_summary([dict(score=1., match=False)], 'score', 'match')
        self.assertIsNone(all_negative['prediction_pool_AUROC'])
        self.assertIsNone(all_negative['prediction_pool_average_precision'])

    def test_fixed_probability_bins_and_Brier_without_probability_refit(self):
        rows = [dict(log_probability=math.log(p), match=y) for p, y in ((.25, True), (.75, False), (1., True))]
        result = probability_diagnostics(rows, 'match')
        self.assertEqual(result['reference_relative_Brier'], .375)
        self.assertEqual(result['fixed_probability_bins'][2]['count'], 1)
        self.assertEqual(result['fixed_probability_bins'][9]['count'], 1)
        self.assertFalse(result['empirical_biological_calibration_claimed'])
        self.assertEqual(sum(b['count'] for b in result['fixed_probability_bins']), 3)
