"""New A41 diagnostic joins and blockers on artificial saved-score rows only."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
from test_m28_emission_residual_raw import receipts

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import diagnose_emission_residual_attrition as M


class EmissionResidualAttrition(unittest.TestCase):
    def test_components_transitions_joint_blocker_labels_and_full_denominators(self):
        chains = [[[0,6]], [[10,13],[16,19],[22,25],[28,31],[36,39]],
                  [[10,25]], [[25,39]], [[50,56]], [[60,66]], [[70,76]]]
        labels = [1,1,0,-1,1,0,-1]
        old = dict(chains=chains, scores=[2.,2.,0.,0.,-1.,1.,1.], selected=[0,1,5,6])
        saved = []; previous = []; labelled = []
        for species in M.R.C.R.ALLOW:
            window = dict(species=species, split='train', id=species+'|toy|0|+')
            new = dict(window=window, chains=chains, scores=[-1.,5.,3.,3.,2.,0.,1.], selected=[2,3,4,6],
                       frozen_path_gain=old['scores'], residual_path_gain=[-3.,3.,3.,3.,3.,-1.,0.],
                       null=0., reference_used_in_scoring=False)
            saved.append(new)
            previous.append(dict(window=window, variants={'A31_A3':old}, reference_used_in_scoring=False))
            labelled.append(dict(window=window, labels=labels,
                counts=dict(A31_A3=M.R.selected_counts(old, labels), A41=M.R.selected_counts(new, labels)),
                eligible_references=5, eligible_five_or_more=3))
        result = M.census(saved, previous, labelled)
        self.assertEqual((result['windows'], result['natural_candidates']), (2,14))
        self.assertEqual(result['max_original_A31_frozen_component_error'], 0.)
        for species in M.R.C.R.ALLOW:
            groups = result['per_species'][species]
            self.assertEqual(groups['1']['lost_nonpositive'], 1)
            self.assertEqual(groups['1']['lost_positive_not_selected'], 1)
            self.assertEqual(groups['1']['gained'], 1)
            self.assertEqual(groups['-1']['gained'], 1)
            self.assertEqual(result['eligible_reference_pairs'][species], dict(eligible=5, five_or_more=3))
            self.assertEqual(result['lost_positive_correct_blocker_label_sets'][species]
                ['correct_at_least_5_exons'], {'chains':1, '-1,0':1})
        self.assertEqual(result['pooled_lost_positive_correct_blocker_label_sets']['correct']['chains'], 2)
        broken = deepcopy(saved); broken[0]['residual_path_gain'][0] += 1.
        with self.assertRaisesRegex(ValueError, 'FP64'): M.census(broken, previous, labelled)
        broken = deepcopy(previous); broken[0]['variants']['A31_A3']['chains'][0] = [[0,7]]
        with self.assertRaisesRegex(ValueError, 'identity'): M.census(saved, broken, labelled)
        with self.assertRaisesRegex(ValueError, 'Complete'): M.census(saved[:-1], previous, labelled)

    def test_FAIL_is_analyzable_but_missing_proof_precedes_score_access(self):
        fitted, report, checked = receipts()
        for value in (report, checked):
            for arms in value['selected_counts_per_species'].values():
                arms['A41'] = deepcopy(arms['A31_A3'])
            value['predeclared_TRAIN_gate'] = M.D.gate(value['selected_counts_per_species'])
            value['TRAIN_preparation_supported'] = False
        self.assertFalse(all(M.validate_closed(fitted, report, checked).values()))
        for proof in M.I.CPU_PROOFS:
            broken = {**checked, proof:False}
            with self.subTest(proof=proof), patch.object(M.D.F.A.A, 'require_jobs'), \
                 patch.object(M.R, 'completed', side_effect=[fitted, report, broken]), \
                 patch.object(M.R.C.R, 'rows') as rows:
                with self.assertRaisesRegex(ValueError, 'six CPU reconciliations'): M.main()
                rows.assert_not_called()


if __name__ == '__main__':
    unittest.main()
