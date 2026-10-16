"""A38 post-result routing and two-comparator integration, no actual score access."""
import copy
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import test_m28_logistic_context_readout as readout_tests
from test_m28_logistic_context_raw import receipts

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import diagnose_logistic_context_attrition as A


def paired_fixture():
    saved, previous, labels, context = readout_tests.LogisticReadoutTests().fixture()
    labelled, _, _ = A.D.paired_readout(saved, previous, labels, context)
    return saved, previous, labelled, context


class LogisticAttritionTests(unittest.TestCase):
    def test_two_comparisons_keep_unknown_and_full_reference_grain_separate(self):
        result = A.census(*paired_fixture())
        self.assertEqual(result['windows'], 2)
        self.assertEqual(result['natural_candidates'], 6)
        baseline, context = result['pooled']['A31_A3_to_A38'], result['pooled']['A37_to_A38']
        self.assertEqual(baseline['1']['kept'], 2)
        self.assertEqual(context['1']['gained'], 2)
        self.assertEqual(baseline['-1']['kept'], 2)
        self.assertEqual(context['-1']['gained'], 2)
        self.assertEqual(context['0']['pool'], 2)
        self.assertEqual(context['0']['new_selected'], 0)
        self.assertEqual(context['correct_at_least_5_exons']['gained'], 2)
        self.assertEqual(sum(v['eligible'] for v in result['eligible_reference_pairs'].values()), 4)
        self.assertEqual(sum(v['five_or_more'] for v in result['eligible_reference_pairs'].values()), 2)

    def test_new_arm_can_lose_to_stronger_baseline_while_gaining_over_A37(self):
        saved, previous, labelled, context = paired_fixture()
        for new, lab in zip(saved, labelled):
            new['scores'] = [-1., 2., 1.]
            new['selected'] = A.R.select_nonoverlapping(new['chains'], new['scores'])
            lab['counts']['A38'] = A.R.selected_counts(new, lab['labels'])
        result = A.census(saved, previous, labelled, context)
        old, prior = result['pooled']['A31_A3_to_A38'], result['pooled']['A37_to_A38']
        self.assertEqual(old['1']['lost_nonpositive'], 2)
        self.assertEqual(prior['1']['neither_nonpositive'], 2)
        self.assertEqual(prior['0']['gained'], 2)
        self.assertFalse(all(A.D.gate(result['selected_counts_per_species']).values()))

    def test_PASS_and_FAIL_both_allow_diagnosis_but_stale_gate_and_source_do_not(self):
        fitted, report, checked = receipts()
        self.assertTrue(all(A.validate_completed(fitted, report, checked).values()))
        for value in (report, checked):
            values = value['selected_counts_per_species']['arabidopsis_thaliana']
            values['A38']['selected_correct'] = values['A31_A3']['selected_correct']
            value['predeclared_TRAIN_gate'] = A.D.gate(value['selected_counts_per_species'])
            value['TRAIN_preparation_supported'] = False
        self.assertFalse(all(A.validate_completed(fitted, report, checked).values()))
        bad = copy.deepcopy(report); bad['TRAIN_preparation_supported'] = True
        with self.assertRaisesRegex(ValueError, 'five-condition'):
            A.validate_completed(fitted, bad, checked)
        with self.assertRaisesRegex(ValueError, 'source jobs'):
            A.validate_completed(fitted, {**report, 'source_fit_job': '13410061'}, checked)

    def test_incomplete_CPU_proof_blocks_before_any_saved_rows_are_opened(self):
        fitted, report, checked = receipts()
        checked[A.I.CPU_PROOFS[-1]] = False
        with patch.object(A.D.F.A, 'require_jobs'), \
             patch.object(A.R, 'completed', side_effect=[fitted, report, checked]), \
             patch.object(A.R.C.R, 'rows') as rows, patch.object(A.R.C.R, 'write_summary') as write:
            with self.assertRaisesRegex(ValueError, 'CPU reconciliations'):
                A.main()
            rows.assert_not_called(); write.assert_not_called()

    def test_new_four_way_census_rejects_missing_or_reordered_context(self):
        rows = paired_fixture(); rows[3].pop()
        with self.assertRaisesRegex(ValueError, 'four-way'):
            A.census(*rows)
        rows = paired_fixture(); rows[3][0]['chains'] = rows[3][0]['chains'][::-1]
        with self.assertRaisesRegex(ValueError, 'identity'):
            A.census(*rows)


if __name__ == '__main__':
    unittest.main()
