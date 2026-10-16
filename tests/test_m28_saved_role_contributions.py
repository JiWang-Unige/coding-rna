import importlib.util
from pathlib import Path
import sys
import unittest
import torch

BASE = Path(__file__).resolve().parents[1]
SCRIPTS = BASE/'scripts/experiments/M28-METHOD-RESTART'
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location('role_parts', SCRIPTS/'diagnose_saved_role_contributions.py')
D = importlib.util.module_from_spec(spec); spec.loader.exec_module(D)


class RoleContributionTests(unittest.TestCase):
    def test_complete_partition_preserves_original_role_gain(self):
        logits = torch.zeros(24, 15)
        for name, value in dict(coding=2., intron=-3., donor=4., acceptor=5., start=6., stop=7.).items():
            logits[:, list(D.GROUPS[name])] = value
        chain = [[1, 7], [13, 19]]; parts = D.components(logits, chain)
        self.assertEqual(parts['groups']['intron']['role_bases'], 6)
        self.assertEqual(parts['groups']['intron']['gain'], -18.)
        self.assertEqual(parts['dominant_negative_group'], 'intron')
        self.assertEqual(sum(x['role_bases'] for x in parts['groups'].values()), 18)
        self.assertEqual(parts['component_sum'], 38.)
        self.assertEqual(parts['component_sum'], float(D.D.R.local_path_scores(logits, [chain])[0]))

    def test_absent_role_is_not_zero_per_base_evidence(self):
        parts = D.components(torch.ones(10, 15), [[1, 7]])
        for name in ('coding', 'intron', 'donor', 'acceptor'):
            self.assertEqual(parts['groups'][name]['role_bases'], 0)
            self.assertIsNone(parts['groups'][name]['gain_per_role_base'])
        self.assertIsNone(parts['dominant_negative_group'])
        self.assertEqual(parts['component_sum'], 0.)

    def test_pooled_per_base_statistic_is_weighted_by_role_bases(self):
        rows = []
        for bases, gain in ((2, 2.), (10, -10.)):
            groups = {name: dict(role_bases=0, gain=0.) for name in D.GROUPS}
            groups['coding'] = dict(role_bases=bases, gain=gain)
            rows.append(dict(components=dict(groups=groups, dominant_negative_group='coding' if gain < 0 else None)))
        stats = D.statistics(rows)
        self.assertEqual(stats['groups']['coding']['role_bases'], 12)
        self.assertAlmostEqual(stats['groups']['coding']['pooled_gain_per_role_base'], -8/12)
        self.assertEqual(stats['dominant_negative_group_counts']['coding'], 1)


if __name__ == '__main__':
    unittest.main()
