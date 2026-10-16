#!/usr/bin/env python3
"""Artificial-only combined-cell fit, strict-final and original-TRAIN bindings."""
import copy
import os
import runpy
import unittest
from unittest.mock import patch
import fit_full_grammar_joint24_native_prior as F
import annotate_full_grammar_joint24_native_cds_prior as A
import infer_full_grammar_joint24_native_prior_train as I
import verify_full_grammar_joint24_native_prior_fit as V
import evaluate_full_grammar_original_train as E


class BindingChecks(unittest.TestCase):
    def test_new_final_identity_matches_fit_and_rejects_old_routes(self):
        state = dict(configuration=copy.deepcopy(F.CONFIG),protocol=F.PROTOCOL,step=9216,epoch=6,
            encoder_revision=F.Q.F.REVISION)
        self.assertEqual(A.validate_checkpoint(state),F.CONFIG)
        for configuration in (A.N.IDENTITY,A.J.IDENTITY):
            with self.assertRaises(ValueError): A.validate_checkpoint(dict(state,configuration=configuration))

    def test_fit_requires_completed_exact_combined_preflight(self):
        row = dict(zero_prior_original_joint_window_logits_exact=3,zero_prior_original_joint_loss_matches=True,
            strict_neural_weight_in_memory_roundtrip=True,updated_raw_adapter_potentials_exact=True)
        value = dict(completed=True,job_id='fixture',configuration=F.CONFIG,protocol=F.PROTOCOL,
            optimizer_steps_this_run=2,checkpoints_written=False,formal_fit_started=False,disposable_updates_only=True,
            pretrained_encoder_and_heads_frozen_no_grad=True,DEV_used=False,sealed_data_used=False,records=[row,row])
        F.require_preflight(value,'COMPLETED','fixture')
        for field,replacement in (('configuration',F.N.CONFIG),('checkpoints_written',True),('DEV_used',True)):
            with self.assertRaises(ValueError):
                F.require_preflight(dict(value,**{field:replacement}),'COMPLETED','fixture')
        with self.assertRaises(ValueError): F.require_preflight(value,'RUNNING','fixture')

    def test_raw_requires_only_combined_verified_final_paths(self):
        fitted = dict(completed=True,job_id='fixture',cumulative_optimizer_steps=9216,optimizer_steps_this_run=4608,
            DEV_used=False,sealed_data_used=False,protocol=F.PROTOCOL,configuration=F.CONFIG,
            final_checkpoint=str(I.CHECKPOINT.relative_to(I.ROOT)))
        verified = dict(completed=True,source_fit_job='fixture',final_checkpoint_strict_load=True,
            configuration=F.CONFIG,protocol=F.PROTOCOL,final_checkpoint=fitted['final_checkpoint'])
        I.require_completed_fit(fitted,verified,'fixture','COMPLETED','COMPLETED')
        for field,value in (('cumulative_optimizer_steps',6144),('configuration',F.N.CONFIG),
                ('final_checkpoint',str((F.N.FIT/'step_009216.pt').relative_to(I.ROOT)))):
            with self.assertRaises(ValueError):
                I.require_completed_fit(dict(fitted,**{field:value}),verified,'fixture','COMPLETED','COMPLETED')
        with self.assertRaises(ValueError):
            I.require_completed_fit(fitted,verified,'fixture','COMPLETED','RUNNING')

    def test_original_exposure_arithmetic_and_unknown_rows_are_enforced(self):
        plan = [dict(draw_ordinal=i,window_id='toy'+str(i),species='toy',importance_weight=2.,
            known_observation_denominator=4,fine_positions=4,coarse_positions=0,fully_unknown_positions=0,
            objective_multiplier=.5) for i in range(1536)]
        rows = [dict(**{k:v for k,v in plan[i%1536].items() if k != 'objective_multiplier'},
            step=4609+i,additional_step=1+i,epoch=4+i//1536,raw_NLL=8.,normalized_NLL=2.,weighted_total=4.,
            gradient_norm_before_clip=1.,prior_weight_l2=.1,seconds=.01) for i in range(4608)]
        V.verify_rows(rows,plan)
        for field,value in (('step',4608),('weighted_total',8.),('importance_weight',1.)):
            bad = list(rows); bad[0] = dict(rows[0],**{field:value})
            with self.assertRaises(ValueError): V.verify_rows(bad,plan)
        plan[0] = dict(plan[0],fine_positions=0,fully_unknown_positions=4)
        bad = list(rows); bad[0] = dict(rows[0],fine_positions=0,fully_unknown_positions=4)
        with self.assertRaises(ValueError): V.verify_rows(bad,plan)

    def test_unchanged_ruler_gets_new_paths_not_old_prior_or_original(self):
        calls = []
        with patch.object(E,'main',side_effect=lambda *args:calls.append(args)):
            runpy.run_module('evaluate_full_grammar_joint24_native_prior_train',run_name='__main__')
        self.assertEqual(calls,[(I.OUT,I.CHECKPOINT,A.IDENTITY,
            I.ROOT/'outputs/M28-FULL-GRAMMAR-JOINT24-NATIVE-PRIOR-TRAIN-DIAGNOSTIC')])
        self.assertNotEqual(I.OUT,E.I.OUT)
        self.assertNotEqual(I.CHECKPOINT,F.N.FIT/'step_009216.pt')


if __name__ == '__main__':
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(BindingChecks))
    F.C.write(F.BINDINGS/'summary.json',dict(completed=result.wasSuccessful(),job_id=os.environ['SLURM_JOB_ID'],
        configuration=F.CONFIG,protocol=F.PROTOCOL,tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        artificial_metadata_and_mocked_ruler_only=True,research_data_read=False,models_or_weights_read=False,
        native_forward=False,optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False))
    if not result.wasSuccessful(): raise SystemExit(1)
