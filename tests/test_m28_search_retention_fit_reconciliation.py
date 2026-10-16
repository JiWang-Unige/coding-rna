"""Concrete full-artifact/Slurm-root boundaries; no real weights or TRAIN updates."""
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import reconcile_query_search_retention_fit as T


class FinalFitReconciliationTests(unittest.TestCase):
    def test_live_source_stops_before_checkpoint_or_train_metadata(self):
        with patch.object(T.K, 'require_source_timeout', side_effect=ValueError('TIMEOUT required')), \
                patch.object(T.K, 'restore_checkpoint') as restore, patch.object(T.R.Q, 'read_rows') as draws:
            with self.assertRaisesRegex(ValueError, 'TIMEOUT'):
                T.reconcile()
            restore.assert_not_called(); draws.assert_not_called()

    def test_partial_summary_cannot_use_readout_only_recovery(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)
            fitted = dict(arm='A30', configuration=T.R.CONFIG, protocol=T.R.PROTOCOL,
                optimizer_steps=4352, feature_revision=T.R.A.F.REVISION, completed=True,
                DEV_used_for_selection=False, test_setaria_used=False,
                final_checkpoint=str((T.R.FIT/'step_004608.pt').relative_to(T.R.ROOT)))
            (source/'summary.json').write_text(json.dumps(fitted))
            with patch.object(T.K, 'require_source_timeout', return_value=dict(state='TIMEOUT')), \
                    patch.object(T.K, 'SOURCE', source), patch.object(T.K, 'restore_checkpoint') as restore, \
                    patch.object(T.R.Q, 'read_rows') as draws:
                with self.assertRaisesRegex(ValueError, 'last generator'):
                    T.reconcile()
                restore.assert_not_called(); draws.assert_not_called()

    def test_batch_success_does_not_replace_successful_reconciliation_root(self):
        for state, exit_code in [('TIMEOUT', '0:0'), ('COMPLETED', '1:0'), ('RUNNING', '0:0')]:
            result = SimpleNamespace(stdout=f'123|{state}|{exit_code}\n123.batch|COMPLETED|0:0\n')
            with patch.object(T.subprocess, 'run', return_value=result):
                with self.assertRaisesRegex(ValueError, 'actual reconciliation root'):
                    T.reconciliation_receipt('123')
        with patch.object(T.subprocess, 'run', return_value=SimpleNamespace(stdout='123|COMPLETED|0:0\n')):
            T.reconciliation_receipt('123')


if __name__ == '__main__':
    unittest.main()
