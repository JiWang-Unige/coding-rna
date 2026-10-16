#!/usr/bin/env python3
"""Artificial-only preparation: never read research DEV or actual weights."""
import copy
import json
import os
from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import patch
import infer_full_grammar_objective_contrast_dev as I
import evaluate_full_grammar_dev_baselines as E
import check_full_grammar_dev_preparation as Old

OUT = I.ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-DEV-CHECKS'


class DEVPreparationChecks(unittest.TestCase):
    def test_three_candidate_identities_and_no_reference_before_reconciliation(self):
        seqid, species, length = 'artificial', 'synthetic', 40
        scope = dict(seqid=seqid, windows=[[0, length]], reference_used=False,
            joint_whole_oriented_contig_decode=True, context_policy=I.A.R.CONTEXT,
            pruning='exact_same_stage_actual_carry_lower_exon_budget_dominance_max_only',
            score_accumulation='FP64_per_base_forward_order',
            tie_policy='background_continuation_then_state_table_arrival_lower_count_dominates_equal_score')
        for arm in I.A.ARMS:
            checkpoint = I.final_checkpoint(arm)
            identity = dict(I.A.IDENTITY, arm=arm)
            outer = dict(inference_complete=True, task=0, species=species, seqid=seqid, length=length,
                checkpoint=str(checkpoint), raw_scope=I.I.RAW_SCOPE, DEV_used=True,
                sealed_data_used=False, reference_used=False, labels_or_feature_cache_used=False, optimizer_steps_this_run=0)
            raw = dict(inference_complete=True, checkpoint=str(checkpoint), protocol=I.A.PROTOCOL,
                lengths={seqid:length}, configuration=identity, reference_used=False,
                labels_or_feature_cache_used=False, observation_model_used_at_inference=False,
                optimizer_steps_this_run=0, backend='exact_exon_budget_skyline',
                scopes=[dict(scope, strand='+'), dict(scope, strand='-')])
            args = dict(checkpoint=checkpoint, identity=identity, protocol=I.A.PROTOCOL)
            E.check_candidate(outer, raw, 0, species, seqid, length, **args)
            with self.assertRaises(ValueError):
                E.check_candidate(outer, raw, 0, species, seqid, length)
            changed = copy.deepcopy(raw)
            changed['configuration']['arm'] = I.A.ARMS[(I.A.ARMS.index(arm)+1)%3]
            with self.assertRaises(ValueError):
                E.check_candidate(outer, changed, 0, species, seqid, length, **args)
            with patch.object(E, 'reconcile_inputs', side_effect=ValueError('synthetic incomplete inputs')) as inputs, \
                 patch.object(E.C, 'select_reference') as reference:
                with self.assertRaises(ValueError):
                    E.main(Path('/synthetic'), checkpoint, identity, I.A.PROTOCOL, 'M28_'+arm, Path('/synthetic/out'))
                inputs.assert_called_once_with(Path('/synthetic'), checkpoint, identity, I.A.PROTOCOL)
                reference.assert_not_called()

    def test_all_three_TRAIN_closure_and_missing_branch_prevents_DEV_read(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for arm in I.A.ARMS:
                verify, trained = root/'verify'/arm, root/'train'/arm
                verify.mkdir(parents=True)
                trained.mkdir(parents=True)
                for folder in (verify, trained):
                    (folder/'STATUS').write_text('COMPLETED\n')
                (verify/'summary.json').write_text(json.dumps(dict(completed=True, arm=arm,
                    verified_cumulative_optimizer_steps=9216, final_checkpoint_strict_load=True,
                    final_checkpoint=str(I.final_checkpoint(arm).relative_to(I.ROOT)), DEV_used=False, sealed_data_used=False)))
                (trained/'summary.json').write_text(json.dumps(dict(completed=True, job_id=I.TRAIN_JOBS[arm],
                    checkpoint=str(I.final_checkpoint(arm)), DEV_used=False, sealed_data_used=False,
                    status='original_TRAIN_end_to_end_diagnostic_not_DEV_or_generalization')))
            with patch.object(I, 'VERIFICATION', root/'verify'), patch.object(I, 'TRAIN_EVAL', root/'train'):
                I.require_train_closure()
                (root/'train'/I.A.ARMS[2]/'STATUS').write_text('FAILED exit=1\n')
                with patch.object(I.I.R.M, 'read_allowed') as DNA, patch.object(I.A, 'main') as encoder:
                    with self.assertRaises(ValueError):
                        I.main(I.A.ARMS[0], 0)
                    DNA.assert_not_called()
                    encoder.assert_not_called()

    def test_driver_preserves_IUPAC_without_reference_or_training(self):
        arm, species, seqid = I.A.ARMS[1], 'synthetic', 'artificial'
        DNA = 'acgtRynnacgt'
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            folder = root/arm/seqid
            folder.mkdir(parents=True)
            def source(path, allowed):
                self.assertEqual(path, I.ROOT/'data/m1_screen'/species/'genome.fa')
                self.assertEqual(allowed, {seqid:len(DNA)})
                return {seqid:DNA}
            def caller(fasta, checkpoint, out, selected_arm):
                self.assertEqual(fasta.read_text(), '>'+seqid+'\n'+DNA.upper()+'\n')
                self.assertEqual((checkpoint, out, selected_arm), (I.final_checkpoint(arm), folder/'prediction', arm))
            with patch.object(I, 'OUT', root), patch.object(I.I, 'SCOPE', ((species,seqid,len(DNA)),)), \
                 patch.object(I, 'require_train_closure') as closure, \
                 patch.object(I.I.R.M, 'read_allowed', side_effect=source) as DNA_source, \
                 patch.object(I.A, 'main', side_effect=caller) as raw:
                I.main(arm, 0)
                closure.assert_called_once()
                DNA_source.assert_called_once()
                raw.assert_called_once()
            receipt = json.loads((folder/'summary.json').read_text())
            self.assertTrue(receipt['DEV_used'])
            self.assertFalse(receipt['reference_used'])
            self.assertFalse(receipt['sealed_data_used'])
            self.assertEqual(receipt['optimizer_steps_this_run'], 0)

    def test_affected_original_final4608_DEV_path_regression(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(Old, 'OUT', Path(directory)):
            Old.main()  # All five original adapter checks use synthetic mocks.


if __name__ == '__main__':
    began = time.perf_counter()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(DEVPreparationChecks))
    with (OUT/'summary.json').open('x') as handle:
        json.dump(dict(completed=result.wasSuccessful(), job_id=os.environ['SLURM_JOB_ID'],
            scope='artificial_three_final9216_DEV_preparation_only', tests_run=result.testsRun,
            failures=len(result.failures), errors=len(result.errors), DEV_used=False, sealed_data_used=False,
            research_DNA_or_labels_opened=False, historical_DEV_receipts_or_predictions_opened=False,
            checkpoint_loaded=False, actual_native_GLM_forward=False, optimizer_steps_this_run=0,
            seconds=time.perf_counter()-began), handle, indent=2, allow_nan=False)
        handle.write('\n')
    if not result.wasSuccessful():
        raise SystemExit(1)
