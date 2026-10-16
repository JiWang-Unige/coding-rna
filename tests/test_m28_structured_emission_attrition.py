from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import diagnose_structured_emission_attrition as A


class StructuredEmissionAttritionTests(unittest.TestCase):
    def test_zero_is_nonpositive_kept_gained_and_positive_unselected_are_disjoint(self):
        chains = [[[0,4]], [[5,9]], [[10,14]], [[6,8]]]
        values = A.transitions(chains, [1]*4, dict(scores=[1.,4.,-1.,-1.], selected=[0,1]),
                               dict(scores=[0.,4.,3.,2.], selected=[1,2]))['1']
        self.assertEqual(dict(values), dict(pool=4, old_selected=2, new_selected=2,
            lost_nonpositive=1, kept=1, gained=1, neither_positive_not_selected=1))

    def test_positive_lost_and_five_exon_subset_are_not_combined_with_null_rejection(self):
        chains = [[[0,1],[2,3],[4,5],[6,7],[8,9]], [[10,14]]]
        result = A.transitions(chains, [1,1], dict(scores=[3.,1.], selected=[0,1]),
                               dict(scores=[2.,-.1], selected=[]))
        self.assertEqual(result['1']['lost_positive_not_selected'], 1)
        self.assertEqual(result['1']['lost_nonpositive'], 1)
        self.assertEqual(result['correct_at_least_5_exons']['pool'], 1)
        self.assertEqual(result['correct_at_least_5_exons']['lost_positive_not_selected'], 1)

    def test_unknown_stays_separate_and_selected_nonpositive_or_bad_label_rejected(self):
        chains = [[[0,4]], [[5,9]]]
        old = dict(scores=[1.,1.], selected=[0,1]); new = dict(scores=[0.,-.1], selected=[])
        values = A.transitions(chains, [0,-1], old, new)
        self.assertEqual(values['0']['lost_nonpositive'], 1)
        self.assertEqual(values['-1']['lost_nonpositive'], 1)
        self.assertEqual(values['1']['pool'], 0)
        with self.assertRaisesRegex(ValueError, 'positive-gain'):
            A.transitions(chains, [0,-1], old, {**new, 'selected':[0]})
        with self.assertRaisesRegex(ValueError, 'candidate scores/labels'):
            A.transitions(chains, [0,2], old, new)


if __name__ == '__main__':
    unittest.main()
