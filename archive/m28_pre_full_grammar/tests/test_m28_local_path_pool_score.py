"""Saved geometry joins and unknown-participating selection, without a model/data."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import diagnose_local_path_pool_score_train as D


def value(chains, scores):
    return dict(chains=chains, scores=scores, selected=D.R.select_nonoverlapping(chains, scores))


class LocalPathPoolScoreTests(unittest.TestCase):
    def test_reordered_geometry_separates_score_and_pool_changes(self):
        positive, wrong, extra = [[0, 9]], [[0, 12]], [[20, 29]]
        natural = dict(variants=dict(A3_greedy=value([positive, wrong], [.1, 1.]),
                                    A31_A3=value([wrong, positive, extra], [.5, 2., 1.])))
        saved = dict(labels=dict(A3_greedy=[1, 0], A31_A3=[0, 1, 0]))
        counts, local = D.three_cells(natural, saved)
        self.assertEqual(local['scores'], [2., .5]); self.assertEqual(local['labels'], [1, 0])
        self.assertEqual([counts[a]['selected_correct'] for a in D.ARMS], [0, 1, 1])
        self.assertEqual([counts[a]['selected_known_nonmatch'] for a in D.ARMS], [1, 0, 1])
        refs = {'synthetic': dict(eligible=1)}
        metrics = D.R.proxies({'synthetic': counts}, refs)['synthetic']
        first = D.contrast(counts['greedy_query'], counts['greedy_local'], metrics['greedy_query'], metrics['greedy_local'])
        second = D.contrast(counts['greedy_local'], counts['union_local'], metrics['greedy_local'], metrics['union_local'])
        self.assertAlmostEqual(first['proxy_change_percentage_points']['F1']+
                               second['proxy_change_percentage_points']['F1'], 100*(2/3))

    def test_unknown_still_wins_selection_and_is_not_known_false_positive(self):
        positive, unknown = [[0, 9]], [[0, 20]]
        natural = dict(variants=dict(A3_greedy=value([positive], [1.]),
                                    A31_A3=value([positive, unknown], [1., 2.])))
        counts, _ = D.three_cells(natural, dict(labels=dict(A3_greedy=[1], A31_A3=[1, -1])))
        self.assertEqual(counts['greedy_local']['selected_correct'], 1)
        self.assertEqual(counts['union_local']['selected_correct'], 0)
        self.assertEqual(counts['union_local']['selected_unknown'], 1)
        self.assertEqual(counts['union_local']['selected_known_nonmatch'], 0)

    def test_missing_geometry_or_inconsistent_saved_label_is_not_silently_joined(self):
        natural = dict(variants=dict(A3_greedy=value([[[0, 9]]], [1.]),
                                    A31_A3=value([[[0, 12]]], [1.])))
        saved = dict(labels=dict(A3_greedy=[1], A31_A3=[0]))
        with self.assertRaisesRegex(ValueError, 'subset'):
            D.three_cells(natural, saved)
        natural['variants']['A31_A3'] = value([[[0, 9]]], [1.])
        with self.assertRaisesRegex(ValueError, 'same saved R4'):
            D.three_cells(natural, saved)


if __name__ == '__main__':
    unittest.main()
