import sys
from pathlib import Path
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import diagnose_quality_train as D


class QualityTrainDiagnosticTests(unittest.TestCase):
    def test_unknown_excluded_from_candidate_mean_mass(self):
        self.assertEqual(D.known_mass([1, 0, 0, -1], 3.), {'0': 2., '1': 1.})
        self.assertEqual(D.known_mass([-1], 3.), {'0': 0., '1': 0.})

    def test_overlap_only_known_nonmatch_and_touching_is_not_overlap(self):
        rows = [{'chain': [(0, 9)], 'gain': 100., 'label': -1},
                {'chain': [(0, 9)], 'gain': 99., 'label': 1},
                {'chain': [(9, 18)], 'gain': 98., 'label': 0},
                {'chain': [(3, 12)], 'gain': -2., 'label': 0}]
        self.assertEqual(D.overlapping_max([(0, 9)], rows), -2.)
        self.assertIsNone(D.overlapping_max([(18, 27)], rows))

    def test_zero_gain_is_not_positive_and_empty_group_is_reported(self):
        row = dict(gain=0., additive=1., nonadditive=2., null=3., label=1,
                   max_overlapping_known_nonmatch_gain=0.)
        result = D.profile([row])
        self.assertEqual(result['positive_gain'], 0)
        self.assertEqual(result['nonpositive_gain'], 1)
        self.assertEqual(result['positive_candidates_not_above_max_known_overlap'], 1)
        self.assertIsNone(D.profile([])['gain'])


if __name__ == '__main__':
    unittest.main()
