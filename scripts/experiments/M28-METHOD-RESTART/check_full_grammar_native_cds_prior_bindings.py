#!/usr/bin/env python3
"""Synthetic final identity and explicit existing TRAIN evaluation bindings."""
import copy
import os
import runpy
import unittest
from unittest.mock import patch
import fit_full_grammar_native_cds_prior as F
import annotate_full_grammar_native_cds_prior as A
import infer_full_grammar_native_cds_prior_train as I
import evaluate_full_grammar_original_train as E

OUT = F.ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CDS-PRIOR-BINDING-CHECKS'


class BindingChecks(unittest.TestCase):
    def test_fixed_final_identity_matches_fit_and_rejects_old_checkpoint(self):
        metadata = dict(configuration=copy.deepcopy(F.CONFIG), protocol=A.PROTOCOL,
            step=9216, epoch=6, encoder_revision=F.Q.F.REVISION)
        self.assertEqual(A.validate_checkpoint(metadata), F.CONFIG)
        for key, value in (('step', 4608), ('epoch', 3), ('protocol', A.R.PROTOCOL)):
            old = copy.deepcopy(metadata); old[key] = value
            with self.assertRaises(ValueError): A.validate_checkpoint(old)
        old = copy.deepcopy(metadata); old['configuration']['model'] = A.R.IDENTITY['model']
        with self.assertRaises(ValueError): A.validate_checkpoint(old)

    def test_shared_ruler_receives_only_new_explicit_paths_and_identity(self):
        calls = []
        with patch.object(E, 'main', side_effect=lambda *args: calls.append(args)):
            runpy.run_module('evaluate_full_grammar_native_cds_prior_train', run_name='__main__')
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0], (I.OUT, I.CHECKPOINT, A.IDENTITY,
            I.ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CDS-PRIOR-TRAIN-DIAGNOSTIC'))
        self.assertNotEqual(I.OUT, E.I.OUT)
        self.assertNotEqual(I.CHECKPOINT, E.I.CHECKPOINT)


if __name__ == '__main__':
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(BindingChecks))
    F.C.write(OUT/'summary.json', dict(completed=tests.wasSuccessful(), job_id=os.environ['SLURM_JOB_ID'],
        protocol=F.PROTOCOL, tests_run=tests.testsRun, failures=len(tests.failures), errors=len(tests.errors),
        synthetic_metadata_and_mocked_ruler_only=True, research_data_read=False,
        models_or_weights_read=False, native_forward=False, optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False))
    if not tests.wasSuccessful(): raise SystemExit(1)
