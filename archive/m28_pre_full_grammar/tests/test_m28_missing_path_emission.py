import importlib.util
from pathlib import Path
import sys
import unittest

BASE = Path(__file__).resolve().parents[1]
SCRIPTS = BASE/'scripts/experiments/M28-METHOD-RESTART'
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location('missing_emission', SCRIPTS/'diagnose_missing_path_emission_train.py')
D = importlib.util.module_from_spec(spec); spec.loader.exec_module(D)


def pool(chains, scores):
    return dict(chains=chains, scores=scores, selected=D.R.select_nonoverlapping(chains, scores))


def ref(chain, cohort='missing'):
    return dict(window_id='old_train', original_candidate_index=1, species='arabidopsis_thaliana',
        chain=chain, exon_class='one', cohort=cohort, mechanism='all_queries_locally_excluded')


class ReferenceAssistedEmissionTests(unittest.TestCase):
    def test_original_null_and_global_interval_competition(self):
        p = pool([[[0, 9]], [[10, 19]]], [2., 2.])
        for gain, state in ((-1., D.STATES[0]), (0., D.STATES[0]), (3., D.STATES[1])):
            rows, joint = D.evaluate_window([ref([[0, 19]])], [gain], p, [0, 0])
            self.assertEqual(rows[0]['reference_assisted_single_state'], state)
            self.assertEqual(joint['injected_references_selected'], 0)
            self.assertEqual(joint['selected_chains'], p['chains'])

    def test_existing_correct_displacement_is_not_hidden(self):
        p = pool([[[0, 20]]], [3.])
        for label, state in ((0, D.STATES[2]), (-1, D.STATES[2]), (1, D.STATES[3])):
            rows, joint = D.evaluate_window([ref([[1, 19]])], [5.], p, [label])
            self.assertEqual(rows[0]['reference_assisted_single_state'], state)
            self.assertEqual(joint['original_selected_dropped_label_counts'][str(label)], 1)
            self.assertEqual(joint['selected_label_count_changes']['1'], 0 if label == 1 else 1)

    def test_individual_recoveries_are_not_summed_as_joint(self):
        p = pool([[[0, 25]]], [3.])
        refs = [ref([[1, 19]]), ref([[2, 20]])]
        rows, joint = D.evaluate_window(refs, [4., 5.], p, [0])
        self.assertTrue(all(r['reference_assisted_single_state'] == D.STATES[2] for r in rows))
        self.assertEqual(joint['injected_missing_reference_geometries'], 2)
        self.assertEqual(joint['injected_references_selected'], 1)
        self.assertEqual(joint['selected_chains'], [[[2, 20]]])

    def test_controls_and_missing_identity_keep_original_scores(self):
        p = pool([[[0, 20]]], [3.])
        rows, joint = D.evaluate_window([ref([[0, 20]], 'natural_hit_control')], [3.], p, [1])
        self.assertTrue(rows[0]['original_control_selected'])
        self.assertEqual(joint['injected_missing_reference_geometries'], 0)
        self.assertEqual(joint['selected_chains'], p['chains'])
        with self.assertRaises(ValueError):
            D.evaluate_window([ref([[0, 20]], 'natural_hit_control')], [4.], p, [1])
        with self.assertRaises(ValueError):
            D.assisted_selection(p, [1], [([[0, 20]], 3.)])


if __name__ == '__main__':
    unittest.main()
