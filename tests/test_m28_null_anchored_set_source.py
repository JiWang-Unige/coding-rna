"""Regression for two-arm competition versus three-arm closed A38 receipts."""
import copy
from pathlib import Path
import sys
import unittest
from unittest import mock
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import preflight_null_anchored_set as A


class ClosedSourceTests(unittest.TestCase):
    def fixture(self):
        counts = {s: {arm: dict(selected_correct=2, selected_known_nonmatch=1,
                               selected_five_or_more_correct=1)
                      for arm in ('A31_A3', 'A37', 'A38')} for s in A.R.C.R.ALLOW}
        report = dict(configuration=A.L.CONFIG, selected_counts_per_species=counts,
                      eligible_reference_pairs={'toy': {'eligible': 2, 'five_or_more': 1}},
                      TRAIN_preparation_supported=False, predeclared_TRAIN_gate=A.L.gate(counts))
        checked = copy.deepcopy(report)
        checked['source_jobs'] = dict(fit='13412630', readout='13412701')
        for key in ('actual_original_draw_weights_counts_and_4608_updates_verified',
                    'actual_three_model_checkpoint_tensors_and_branch_optimizer_exposure_verified',
                    'actual_saved_natural_identity_and_FP32_score_components_verified',
                    'actual_saved_labels_DP_baseline_denominators_and_five_conditions_recomputed'):
            checked[key] = True
        attrition = copy.deepcopy(report)
        attrition.update(windows=1491, natural_candidates=16321,
                         scientific_status='A38_TRAIN_FAIL_closed_no_rescue')
        competition = copy.deepcopy(attrition)
        competition['selected_counts_per_species'] = {
            s: {arm: value[arm] for arm in ('A31_A3', 'A38')} for s, value in counts.items()}
        competition['all_individual_and_joint_actual_selected_exchanges_verified'] = True
        competition['source_jobs'] = dict(fit='13412630', readout='13412701', CPU='13412702',
                                          attrition='13414584', baseline_readout='13398521')
        return [{}, report, checked, attrition, competition]

    def run_closed(self, values):
        with mock.patch.object(A, 'require_jobs'), mock.patch.object(A.L, 'validate_fit'), \
             mock.patch.object(Path, 'read_text', side_effect=['COMPLETED', A.SCALAR_JOB]), \
             mock.patch.object(A.B.R, 'completed', side_effect=values):
            A.closed_a38()

    def test_exact_common_two_arms_accept_without_fabricating_third_arm(self):
        self.run_closed(self.fixture())

    def test_changed_common_count_or_missing_actual_proof_still_rejects(self):
        values = self.fixture()
        species = next(iter(values[-1]['selected_counts_per_species']))
        values[-1]['selected_counts_per_species'][species]['A38']['selected_correct'] += 1
        with self.assertRaises(ValueError):
            self.run_closed(values)
        values = self.fixture()
        values[2]['actual_saved_natural_identity_and_FP32_score_components_verified'] = False
        with self.assertRaises(ValueError):
            self.run_closed(values)


if __name__ == '__main__':
    unittest.main()
