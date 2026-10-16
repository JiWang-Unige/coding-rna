"""Synthetic checks for A22 categories and exact feasible-set accounting."""
import importlib.util
from pathlib import Path
import unittest

path = Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART/audit_chain_displacement.py'
spec = importlib.util.spec_from_file_location('audit_chain_displacement', path)
A = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A)


def candidate(start, end, label, gain, origin='free'):
    return dict(chain=[[start, end]], label=label, gain=gain, origin=origin)


class TestDisplacement(unittest.TestCase):
    def test_coverage_categories_do_not_inject_gold(self):
        local = [candidate(0, 10, 0, 2), candidate(20, 30, 0, 2)]
        original = local+[candidate(20, 29, 1, 100, 'reference_injected')]
        audit = A.audit_window(local, original)
        self.assertEqual(audit['selected_indices'], [0, 1])
        self.assertEqual([r['category'] for r in audit['selected_known_nonmatch']], list(A.CATEGORIES[:2]))

    def test_available_nonpositive_gold_and_touching_intervals(self):
        local = [candidate(0, 10, 1, 0), candidate(0, 9, 0, 1), candidate(10, 20, 1, 3)]
        audit = A.audit_window(local, local)
        self.assertEqual(audit['selected_indices'], [1, 2])
        self.assertEqual(audit['selected_known_nonmatch'][0]['category'], A.CATEGORIES[2])
        self.assertEqual(audit['missed_nonpositive_natural_positive_indices'], [0])
        self.assertFalse(A.overlap(local[1]['chain'], local[2]['chain']))

    def test_two_weaker_candidates_jointly_displace_gold_with_unknown_retained(self):
        local = [candidate(0, 20, 1, 3), candidate(0, 10, 0, 2), candidate(10, 20, -1, 2)]
        audit = A.audit_window(local, local)
        self.assertEqual(audit['selected_indices'], [1, 2])
        self.assertEqual(audit['selected_known_nonmatch'][0]['category'], A.CATEGORIES[3])
        item = audit['positive_gain_missed_exchanges'][0]
        self.assertEqual(item['selected_conflict_indices'], [1, 2])
        self.assertEqual((item['inserted_gain'], item['removed_gain']), (3, 4))
        self.assertEqual(item['conflict_counts_by_label'], {'1': 0, '0': 1, '-1': 1})

    def test_joint_exchange_detects_a_nonoptimal_set(self):
        local = [candidate(0, 10, 1, 3), candidate(10, 20, 1, 3), candidate(0, 20, 0, 4)]
        self.assertEqual(A.exchange(local, [2], [0])['removed_gain'], 4)
        self.assertEqual(A.exchange(local, [2], [1])['removed_gain'], 4)
        with self.assertRaisesRegex(ValueError, 'feasible gold-set exchange'):
            A.exchange(local, [2], [0, 1])
        self.assertEqual(A.audit_window(local, local)['selected_indices'], [0, 1])


if __name__ == '__main__':
    unittest.main()
