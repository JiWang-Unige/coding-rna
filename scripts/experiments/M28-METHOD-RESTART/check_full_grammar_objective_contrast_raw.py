#!/usr/bin/env python3
"""Artificial checkpoint identity and shared raw caller/evaluator routing."""
import copy
import json
import os
from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import patch
import torch
import fit_full_grammar_objective_contrast as C
import annotate_full_grammar_objective_contrast as A
import evaluate_full_grammar_original_train as E

OUT = C.ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-RAW-CHECKS'


def state(arm):
    return dict(configuration=C.arm_config(arm), protocol=C.PROTOCOL, step=9216, epoch=6,
                encoder_revision=C.Q.F.REVISION, model={})


class RawChecks(unittest.TestCase):
    def test_final_identity_accepts_declared_arms_rejects_drift(self):
        self.assertEqual(C.ARMS, A.ARMS)
        for arm in C.ARMS:
            saved = state(arm)
            self.assertEqual(A.validate_checkpoint(saved, arm), C.arm_config(arm))
            with self.assertRaises(ValueError):
                A.R.validate_checkpoint(saved)  # Old entry point stays final4608-only.
            for key, value in (('step', 7680), ('epoch', 5), ('protocol', A.R.PROTOCOL)):
                changed = copy.deepcopy(saved)
                changed[key] = value
                with self.assertRaises(ValueError):
                    A.validate_checkpoint(changed, arm)
            for key, value in (('arm', C.ARMS[(C.ARMS.index(arm)+1)%3]), ('epsilon', .2), ('null', 1)):
                changed = copy.deepcopy(saved)
                changed['configuration'][key] = value
                with self.assertRaises(ValueError):
                    A.validate_checkpoint(changed, arm)

    def test_raw_main_strict_C0_load_routes_to_same_contig_caller(self):
        class CPUCore(A.R.M28Core):
            def to(self, device):
                if device != 'cuda':
                    raise ValueError('Unexpected execution device')
                return self
        saved = state(C.ARMS[1])
        saved['model'] = CPUCore('C0', 2048, 128).state_dict()
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            fasta = folder/'input.fa'
            fasta.write_text('>artificial\nCCCATGAAATAACCCC\n')
            with (patch.object(A.R.torch, 'load', return_value=saved),
                  patch.object(A.R, 'M28Core', CPUCore),
                  patch.object(A.R, 'load_native', return_value=('backbone', 'tokenizer')),
                  patch.object(A.R, 'annotate_contig', return_value=([], [])) as caller,
                  patch.object(torch.cuda, 'reset_peak_memory_stats'),
                  patch.object(torch.cuda, 'max_memory_allocated', return_value=0),
                  patch.object(torch.cuda, 'max_memory_reserved', return_value=0)):
                A.main(fasta, folder/'final.pt', folder/'prediction', C.ARMS[1])
            self.assertEqual(caller.call_count, 1)
            self.assertEqual(caller.call_args.args[1:],
                             ('backbone', 'tokenizer', 'artificial', 'CCCATGAAATAACCCC', 0, 'cuda'))
            report = json.loads((folder/'prediction/summary.json').read_text())
            self.assertEqual(report['configuration'], C.arm_config(C.ARMS[1]))
            self.assertFalse(report['reference_used'])
            self.assertFalse(report['observation_model_used_at_inference'])
            self.assertEqual(report['backend'], 'exact_exon_budget_skyline')

    def test_shared_reconciliation_uses_explicit_checkpoint_and_identity(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(E.I, 'SCOPE', (('synthetic', 'fixture', 32),)):
            root = Path(directory)
            folder = root/'fixture'
            (folder/'prediction').mkdir(parents=True)
            checkpoint = root/'final.pt'
            identity = dict(A.IDENTITY, arm=C.ARMS[2])
            outer = dict(inference_complete=True, task=0, species='synthetic', seqid='fixture', length=32,
                checkpoint=str(checkpoint), DEV_used=False, sealed_data_used=False, reference_used=False,
                labels_or_feature_cache_used=False, optimizer_steps_this_run=0)
            scopes = [dict(seqid='fixture', strand=strand, windows=[[0, 32]], reference_used=False,
                joint_whole_oriented_contig_decode=True, context_policy=A.R.CONTEXT,
                pruning='exact_same_stage_actual_carry_lower_exon_budget_dominance_max_only',
                score_accumulation='FP64_per_base_forward_order') for strand in ('+', '-')]
            raw = dict(inference_complete=True, checkpoint=str(checkpoint), lengths={'fixture':32},
                reference_used=False, labels_or_feature_cache_used=False, observation_model_used_at_inference=False,
                optimizer_steps_this_run=0, backend='exact_exon_budget_skyline', configuration=identity, scopes=scopes)
            (folder/'STATUS').write_text('COMPLETED\n')
            (folder/'summary.json').write_text(json.dumps(outer))
            (folder/'prediction/summary.json').write_text(json.dumps(raw))
            (folder/'prediction/full_grammar.gff3').write_text('##gff-version 3\n')
            self.assertEqual(set(E.reconcile_inputs(root, checkpoint, identity)), {'fixture'})
            with self.assertRaises(ValueError):
                E.reconcile_inputs(root, checkpoint)  # New weights cannot pass old default identity.
            with self.assertRaises(ValueError):
                E.reconcile_inputs(root, checkpoint, dict(identity, arm=C.ARMS[0]))
            raw['configuration'] = A.R.IDENTITY
            (folder/'prediction/summary.json').write_text(json.dumps(raw))
            self.assertEqual(set(E.reconcile_inputs(root, checkpoint)), {'fixture'})


if __name__ == '__main__':
    torch.set_num_threads(2)
    began = time.perf_counter()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RawChecks))
    C.write(OUT/'summary.json', dict(completed=result.wasSuccessful(), job_id=os.environ['SLURM_JOB_ID'],
        configuration=C.CONFIG, tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        actual_native_GLM_forward=False, research_DNA_or_labels_opened=False, saved_source_weights_read=False,
        optimizer_steps_this_run=0, seconds=time.perf_counter()-began))
    if not result.wasSuccessful():
        raise SystemExit(1)
