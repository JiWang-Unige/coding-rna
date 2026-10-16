"""New A38 saved-record integration; unchanged A22 exchange tests are not repeated."""
import copy
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
from test_m28_logistic_context_raw import receipts

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import diagnose_logistic_context_competition as C


def fixture(chains, labels, old_scores, new_scores):
    window = dict(id='toy', species='arabidopsis_thaliana', split='train')
    old = dict(chains=chains, scores=old_scores, selected=C.R.select_nonoverlapping(chains, old_scores))
    new = dict(window=window, reference_used_in_scoring=False, chains=chains, scores=new_scores,
               selected=C.R.select_nonoverlapping(chains, new_scores))
    prior = dict(window=window, reference_used_in_scoring=False, variants={'A31_A3': old})
    lab = dict(window=window, labels=labels,
        counts={a: C.R.selected_counts(v, labels) for a, v in (('A31_A3', old), ('A38', new))},
        eligible_references=7, eligible_five_or_more=4)
    return new, prior, lab


class LogisticCompetitionTests(unittest.TestCase):
    def test_adapter_keeps_unknown_and_accounts_for_two_weaker_blockers(self):
        gold = [[0, 2], [3, 4], [5, 6], [7, 8], [9, 20]]
        rows = fixture([gold, [[0, 10]], [[10, 20]]], [1, 0, -1], [3., 1., 1.], [3., 2., 2.])
        record = C.diagnose_window(*rows); item = record['cases'][0]
        self.assertTrue(item['baseline_selected']); self.assertTrue(item['five_or_more'])
        self.assertEqual(item['selected_conflict_indices'], [1, 2])
        self.assertEqual(item['conflict_counts_by_label'], {'1': 0, '0': 1, '-1': 1})
        self.assertTrue(item['all_individual_blockers_strictly_weaker'])
        self.assertEqual((item['inserted_gain'], item['removed_gain']), (3., 4.))
        p = C.profiles(record['cases'])
        self.assertEqual(p['baseline_lost_positive']['blocker_label_combinations'], {'0+-1': 1})
        self.assertEqual(p['baseline_lost_positive_five_or_more']['pairs'], 1)

    def test_baseline_lost_and_never_selected_positive_strata_do_not_merge(self):
        rows = fixture([[[0, 10]], [[0, 9]], [[20, 30]], [[20, 29]]], [1, 0, 1, -1],
                       [3., -1., 1., 2.], [3., 4., 3., 4.])
        p = C.profiles(C.diagnose_window(*rows)['cases'])
        self.assertEqual(p['all_positive_unselected']['pairs'], 2)
        self.assertEqual(p['baseline_lost_positive']['pairs'], 1)
        self.assertEqual(p['all_positive_unselected']['blocker_label_combinations'], {'0': 1, '-1': 1})
        self.assertEqual(p['baseline_lost_positive']['blocker_label_combinations'], {'0': 1})
        new, old, lab = fixture([[[0, 10]], [[0, 9]]], [1, -1], [3., 1.], [0., 2.])
        self.assertEqual(C.diagnose_window(new, old, lab)['cases'], [])

    def test_missing_pair_identity_and_stale_saved_selection_are_rejected(self):
        rows = fixture([[[0, 10]], [[0, 9]]], [1, 0], [3., 1.], [1., 2.])
        with self.assertRaisesRegex(ValueError, 'three-way'):
            C.census([rows[0]], [], [rows[2]])
        broken = copy.deepcopy(rows); broken[1]['variants']['A31_A3']['chains'] = [[[1, 10]], [[0, 9]]]
        with self.assertRaisesRegex(ValueError, 'identity'):
            C.diagnose_window(*broken)
        broken = copy.deepcopy(rows); broken[0]['selected'] = [0]
        with self.assertRaisesRegex(ValueError, 'selection'):
            C.diagnose_window(*broken)

    def test_incomplete_CPU_proof_blocks_before_saved_score_rows(self):
        fitted, report, checked = receipts(); checked[C.I.CPU_PROOFS[-1]] = False
        with patch.object(C.D.F.A, 'require_jobs'), \
             patch.object(C.R, 'completed', side_effect=[fitted, report, checked]), \
             patch.object(C.R.C.R, 'rows') as rows, patch.object(C.R.C.R, 'write_summary') as write:
            with self.assertRaisesRegex(ValueError, 'CPU reconciliations'):
                C.main()
            rows.assert_not_called(); write.assert_not_called()


if __name__ == '__main__':
    unittest.main()
