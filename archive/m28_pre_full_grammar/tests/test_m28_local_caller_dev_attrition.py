"""Synthetic owner/sign/selection classification and source-closure checks."""
from pathlib import Path
from types import SimpleNamespace
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import diagnose_local_path_caller_dev_attrition as D


def scope(windows, strand, total, unique, owners, owner_unique, selected):
    return dict(windows=windows, strand=strand, assembly_counts=dict(free_records=total,
        free_unique=unique, owner_records=owners, owner_unique=owner_unique, selected_unique=selected))


class NativeAttritionTests(unittest.TestCase):
    def test_zero_negative_and_false_competition_have_distinct_states(self):
        chains = [[[10, 20]], [[30, 40]], [[9, 41]], [[50, 60]], [[70, 80]], [[90, 95]]]
        free, owner, selected = D.scope_pool([dict(window_index=0, chains=chains, score=[1, 1, 5, 0, -1, 2])],
            scope([[0, 100]], '+', 6, 6, 6, 6, 2))
        # The spanning competitor is not a reference. One other correct chain is selected.
        references = [dict(id=str(i), chain=tuple(map(tuple, chain)))
                      for i, chain in enumerate(chains) if i != 2]
        audit = dict(reference=5, not_in_actual_free_candidates=0, lost_at_owner=0,
            lost_at_conflict_or_nonpositive_gain=4, exact_reference_hits=dict(selected_unique=1))
        labels = D.scope_states(references, free, owner, selected, audit)
        self.assertEqual(list(labels.values()), ['positive_unselected', 'positive_unselected',
            'nonpositive_gain', 'nonpositive_gain', 'selected'])
        audit['reference'] = 4
        with self.assertRaisesRegex(ValueError, 'Every reference'):
            D.scope_states(references, free, owner, selected, audit)

    def test_nonowner_only_is_not_rescued_and_owner_duplicate_uses_max(self):
        records = [dict(window_index=0, chains=[[[40, 60]], [[40, 60]]], score=[2, 3]),
                   dict(window_index=1, chains=[[[10, 20]]], score=[10])]
        free, owner, selected = D.scope_pool(records, scope([[0, 100], [50, 150]], '+', 3, 2, 2, 1, 1))
        self.assertEqual(owner, {((40, 60),): 3.})
        self.assertEqual(D.D.state(((60, 70),), free, owner, False), 'lost_owner')
        self.assertEqual(selected, {((40, 60),)})

    def test_minus_coordinates_and_owner_tie_match_original_assembly(self):
        free, owner, selected = D.scope_pool([dict(window_index=0, chains=[[[10, 20], [30, 40]]], score=[2])],
            scope([[0, 100], [50, 150]], '-', 1, 1, 1, 1, 1))
        expected = ((60, 70), (80, 90))
        self.assertEqual(free, {expected}); self.assertEqual(selected, {expected})
        self.assertEqual(owner[expected], 2.)

    def test_incomplete_root_stops_before_any_native_or_reference_output_read(self):
        output = '\n'.join(job+'|'+('RUNNING' if name == 'A23_eval' else 'COMPLETED')+'|0:0'
                           for name, job in D.SOURCE_JOBS.items())
        with patch.object(D.subprocess, 'run', return_value=SimpleNamespace(stdout=output)), \
                patch.object(D.V.I.R, 'completed') as read_native, patch.object(D.B, 'source') as read_refs:
            with self.assertRaisesRegex(ValueError, 'source roots'):
                D.load_sources()
            read_native.assert_not_called(); read_refs.assert_not_called()


if __name__ == '__main__':
    unittest.main()
