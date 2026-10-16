import sys
import unittest
from pathlib import Path
import torch
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
from analyze_true_prefix_train import probe_query, coordinate_rank
from src.m28.query_beam import decode_query_beam
from test_m28_query_decode import fixture


class TruePrefixTests(unittest.TestCase):
    def test_score_ties_use_coordinate_not_input_order(self):
        self.assertEqual(coordinate_rank([9, 0, 6], torch.tensor([1., 1., 2.]), 9), (3, 1.))

    def test_fifth_true_start_is_locally_excluded_not_reference_injected(self):
        sequence = 'ATGAAATAA'*5
        chain = ((36, 45),)
        model, features = fixture(sequence, chain)
        for a, score in zip(range(0, 45, 9), [6, 5, 4, 3, 2]):
            features[a, 0] = score
        profile = probe_query(model.for_window(features), sequence, chain, 0)
        self.assertEqual(profile['first_local_exclusion'], dict(exon=1, stage='start', rank=5))
        self.assertFalse(profile['local_width4_necessary_conditions_passed'])
        natural = decode_query_beam(model, features, sequence, 128, 4, 4)
        self.assertTrue(all(chain not in [r['chain'] for r in q['completed']] for q in natural))

    def test_multiexon_carry_and_scores_match_actual_returned_path(self):
        sequence, chain = 'ATGCGTTAGTTTAA', ((0, 4), (9, 14))
        model, features = fixture(sequence, chain)
        profile = probe_query(model.for_window(features), sequence, chain, 0)
        self.assertTrue(profile['grammar_supported'])
        self.assertTrue(profile['local_width4_necessary_conditions_passed'])
        self.assertEqual([r['carry_after'] for r in profile['events']], ['C', ''])
        actual = decode_query_beam(model, features, sequence, 128, 4, 4)[0]['completed']
        score = next(r['path_log_score'] for r in actual if r['chain'] == chain)
        self.assertAlmostEqual(score, profile['path_log_score'], places=5)
        self.assertAlmostEqual(profile['path_log_score'], sum(
            r['entry_logp']+r['end_logp']+r['eos_logp'] for r in profile['events']), places=5)

    def test_unsupported_junction_is_retained_not_counted_as_rank_failure(self):
        sequence, chain = 'ATGTGTTAGAAAAATAA', ((0, 4), (9, 17))
        model, features = fixture(sequence, chain)
        profile = probe_query(model.for_window(features), sequence, chain, 0)
        self.assertFalse(profile['grammar_supported'])
        self.assertEqual(profile['unsupported_reason'], 'end_stop')
        self.assertNotIn('local_width4_necessary_conditions_passed', profile)


if __name__ == '__main__':
    unittest.main()
