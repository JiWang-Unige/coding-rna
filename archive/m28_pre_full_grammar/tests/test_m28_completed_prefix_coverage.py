"""Completed-pool prefix distinctions, including mixed tied frontiers."""
import importlib.util
from pathlib import Path
import unittest

PATH = Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART/audit_completed_prefix_coverage.py'
spec = importlib.util.spec_from_file_location('completed_prefix_audit', PATH)
A = importlib.util.module_from_spec(spec); spec.loader.exec_module(A)


class CompletedPrefixTests(unittest.TestCase):
    def test_start_absence_and_empty_pool_are_not_early_termination(self):
        for pool in ([], [[(3, 6)]]):
            r = A.classify([(2, 4), (8, 12)], pool)
            self.assertEqual(r['category'], 'completed_start_absent')
            self.assertEqual(r['same_start_completed_candidate_indices'], [])
        self.assertTrue(A.classify([(2, 4)], [])['empty_completed_pool'])

    def test_exact_chain_not_outer_span_and_first_end_difference(self):
        target = [(2, 4), (8, 12)]
        self.assertEqual(A.classify(target, [[(2, 5), (8, 12)]])['category'], 'first_end_absent')
        self.assertEqual(A.classify(target, [[(2, 4), (9, 12)]])['category'], 'continuation_start_absent')
        self.assertEqual(A.classify(target, [target, [(2, 4), (9, 12)]])['category'], 'natural_hit')

    def test_underextension_and_overextension_are_separate(self):
        target = [(2, 4), (8, 12)]
        self.assertEqual(A.classify(target, [[(2, 4)]])['category'], 'completed_prefix_only')
        r = A.classify(target, [[(2, 4), (8, 12), (15, 19)]])
        self.assertEqual(r['category'], 'target_proper_prefix_of_completed')
        self.assertEqual(r['max_exact_prefix_exons'], 2)

    def test_mixed_deepest_prefix_keeps_all_witnesses_not_one_tie(self):
        target = [(2, 4), (8, 12), (18, 22)]
        r = A.classify(target, [[(2, 4)], [(2, 4), (9, 13)]])
        self.assertEqual(r['category'], 'continuation_start_absent')
        self.assertEqual(r['deepest_prefix_terminated_candidate_indices'], [0])
        self.assertEqual(r['deepest_prefix_continuing_candidate_indices'], [1])
        r = A.classify(target, [[(2, 4)], [(2, 4), (9, 13)], [(2, 4), (8, 14)]])
        self.assertEqual(r['category'], 'continuation_end_absent')
        self.assertEqual(r['deepest_prefix_matching_next_start_candidate_indices'], [2])

    def test_deepest_exact_prefix_beats_shorter_mismatched_paths(self):
        target = [(2, 4), (8, 12), (18, 22)]
        r = A.classify(target, [[(2, 4), (8, 14)], [(2, 4), (8, 12)]])
        self.assertEqual(r['category'], 'completed_prefix_only')
        self.assertEqual(r['max_exact_prefix_exons'], 2)
        self.assertEqual(r['deepest_prefix_candidate_indices'], [1])


if __name__ == '__main__':
    unittest.main()
