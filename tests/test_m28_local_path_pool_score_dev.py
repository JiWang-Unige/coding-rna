"""Saved DEV geometry mapping and unchanged owner/duplicate/selection behavior."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import diagnose_local_path_pool_score_dev as D


def row(chains, scores, sources=None):
    value = dict(window_index=0, window_id='synthetic|x|0|+', chains=chains, score=scores)
    if sources is not None:
        value['decoder_sources'] = sources
    return value


class SavedDEVPoolScoreTests(unittest.TestCase):
    def test_reordered_geometry_preserves_greedy_duplicates_and_max_rule(self):
        right, wrong, extra = [[0, 9]], [[0, 12]], [[20, 29]]
        old = row([right, wrong, right], [.1, 1., .2])
        union = row([wrong, right, extra], [.5, 2., 1.], [['greedy'], ['greedy', 'beam4'], ['beam4']])
        local = D.paired_window(old, union)
        self.assertEqual(local['score'], [2., .5, 2.])
        selected = [D.assemble_predictions([r], [[0, 100]], '+') for r in (old, local, union)]
        self.assertEqual([x['counts']['selected_unique'] for x in selected], [1, 1, 2])
        self.assertEqual(selected[0]['stages']['selected_unique'], [((0, 12),)])
        self.assertEqual(selected[1]['stages']['selected_unique'], [((0, 9),)])
        self.assertEqual(selected[1]['counts']['duplicate_owner_records'], 1)

    def test_nonreference_candidate_still_competes_without_label_filter(self):
        right, false = [[10, 20]], [[9, 41]]
        old = row([right], [1.])
        union = row([right, false], [1., 2.], [['greedy'], ['beam4']])
        local = D.paired_window(old, union)
        reference = [dict(chain=((10, 20),))]
        middle = D.assemble_predictions([local], [[0, 100]], '+')
        expanded = D.assemble_predictions([union], [[0, 100]], '+')
        self.assertEqual(D.scope_result(middle, reference)['exact_reference_hits']['selected_unique'], 1)
        self.assertEqual(D.scope_result(expanded, reference)['exact_reference_hits']['selected_unique'], 0)
        self.assertEqual(expanded['counts']['selected_unique'], 1)

    def test_missing_extra_greedy_geometry_or_identity_cannot_be_silently_matched(self):
        old = row([[[0, 9]]], [1.])
        union = row([[[0, 12]]], [1.], [['greedy']])
        with self.assertRaisesRegex(ValueError, 'geometry must equal'):
            D.paired_window(old, union)
        union = row([[[0, 9]], [[20, 29]]], [1., 1.], [['greedy'], ['greedy']])
        with self.assertRaisesRegex(ValueError, 'geometry must equal'):
            D.paired_window(old, union)
        union = row([[[0, 9]]], [1.], [['greedy']]); union['window_index'] = 1
        with self.assertRaisesRegex(ValueError, 'identity'):
            D.paired_window(old, union)


if __name__ == '__main__':
    unittest.main()
