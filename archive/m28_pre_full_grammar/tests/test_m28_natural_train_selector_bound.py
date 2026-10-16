"""Artificial lexicographic interval bounds and their interpretation only."""
from copy import deepcopy
from itertools import combinations
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import diagnose_natural_train_selector_bound as M


class NaturalTrainSelectorBounds(unittest.TestCase):
    def test_two_bounds_equal_enumeration_but_individual_maxima_need_not_jointly_occur(self):
        chains = [[[0,2]], [[2,4]], [[0,1],[1,2],[2,3],[3,4],[4,5]], [[8,10]],
                  [[10,11],[11,12],[12,13],[13,14],[14,15]], [[20,21]], [[22,23]]]
        labels = [1,1,1,1,1,-1,0]
        actual = M.bounds(chains, labels)
        feasible = []
        for n in range(len(chains)+1):
            for indices in combinations(range(len(chains)), n):
                if any(labels[i] != 1 for i in indices): continue
                spans = sorted((chains[i][0][0], chains[i][-1][1]) for i in indices)
                if any(a[1] > b[0] for a,b in zip(spans,spans[1:])): continue
                feasible.append((len(indices),sum(len(chains[i]) >= 5 for i in indices)))
        primary, long_only = [actual[mode]['counts'] for mode in M.MODES]
        self.assertEqual((primary['selected_correct'],primary['selected_five_or_more_correct']),max(feasible))
        self.assertEqual(long_only['selected_five_or_more_correct'],max(v[1] for v in feasible))
        self.assertEqual((primary['selected_correct'],primary['selected_five_or_more_correct']), (4,1))
        self.assertEqual(long_only['selected_five_or_more_correct'],2)
        self.assertTrue(any(correct > 2 and long > 1 for correct,long in feasible))
        totals = {'toy': {'A31_A3':dict(selected_correct=2,selected_known_nonmatch=0,
                   selected_five_or_more_correct=1), **{mode:actual[mode]['counts'] for mode in M.MODES}}}
        verdict = M.interpretation(totals)
        self.assertTrue(all(verdict['individual_upper_bound_checks'].values()))
        self.assertEqual(verdict['status'],'individual_bounds_sufficient_but_joint_feasibility_not_established')
        self.assertEqual(primary['pool_unknown'],1); self.assertEqual(primary['selected_unknown'],0)
        for chains0, labels0 in (([],[]),([[[0,1]]],[-1])):
            for value in M.bounds(chains0,labels0).values(): self.assertEqual(value['counts']['selected'],0)
        with self.assertRaisesRegex(ValueError,'protected labels'): M.bounds(chains,[1])

    def test_missing_actual_CPU_proof_rejects_before_saved_labels_are_read(self):
        prepared = dict(completed=True,configuration=M.R.CONFIG,protocol=M.R.PROTOCOL,windows=1491,
                        selected_counts_per_species={},eligible_reference_pairs={})
        checked = deepcopy(prepared)
        checked.update(source_jobs=dict(readout='13398521'),
            actual_saved_natural_pool_and_greedy_subset_identity_verified=True,
            actual_R4_labels_selection_reference_denominators_and_proxy_metrics_recomputed=True)
        M.validate_sources(prepared,checked)
        for proof in ('actual_saved_natural_pool_and_greedy_subset_identity_verified',
                      'actual_R4_labels_selection_reference_denominators_and_proxy_metrics_recomputed'):
            broken = {**checked,proof:False}
            with self.subTest(proof=proof), patch.object(M,'source_jobs'), \
                 patch.object(M.R,'completed',side_effect=[prepared,broken]), patch.object(M.R.C.R,'rows') as rows:
                with self.assertRaisesRegex(ValueError,'CPU reconciliation'): M.main()
                rows.assert_not_called()


if __name__ == '__main__':
    unittest.main()
