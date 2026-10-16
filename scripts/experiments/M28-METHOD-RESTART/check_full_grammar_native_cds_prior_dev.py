#!/usr/bin/env python3
"""Artificial DEV wiring and full-family arithmetic; no research data reads."""
import copy
import json
import os
from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import patch
import infer_full_grammar_native_cds_prior_dev as I
import evaluate_full_grammar_dev_baselines as E
import aggregate_full_grammar_continuation_dev as B

OUT = I.ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CDS-PRIOR-DEV-CHECKS'


def synthetic_reports():
    reports = {}
    starts = dict(zip(B.SPECS, (300, 400, 390, 410, 500)))
    for label in B.SPECS:
        report = dict(completed=True, status='fixed_two_DEV_development_comparison_not_independent_generalization',
            checkpoint=str(B.SPECS[label][1]), DEV_used=True, sealed_data_used=False, optimizer_steps_this_run=0,
            input_records_directly_equal=True, unsupported_reference_strata_retained=True,
            reference_policy='M26_CDS_assessable_longest_eligible_CDS_tie_transcript_ID', species={}, pooled={})
        for index, (species, seqid, length) in enumerate(E.I.SCOPE):
            methods = {}
            for method in (label, *E.N.METHODS):
                tp = starts[label] if method == label else dict(zip(E.N.METHODS, (700, 600, 500)))[method]
                methods[method] = dict(exact_chain=E.C.S.prf(tp, 2000 if index == 0 else 1000, E.EXPECTED[species]),
                                       duplicate_chain_records=0)
            report['species'][species] = dict(seqid=seqid, length=length, reference_primary_chains=E.EXPECTED[species],
                retained_reference_strata={'synthetic_retained_case':1}, methods=methods)
        for method in (label, *E.N.METHODS):
            rows = [report['species'][s]['methods'][method]['exact_chain'] for s in E.EXPECTED]
            report['pooled'][method] = dict(exact_chain=E.C.S.prf(*(sum(r[k] for r in rows) for k in ('tp', 'predicted', 'reference'))),
                macro_species_exact_chain_F1=sum(r['f1'] for r in rows)/2,
                worst_species_exact_chain_F1=min(r['f1'] for r in rows))
        reports[label] = report
    return reports


class PreparedDEVChecks(unittest.TestCase):
    def test_TRAIN_closure_prevents_any_DEV_read_or_encoding(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); verify, train = root/'verify', root/'train'
            verify.mkdir(); train.mkdir()
            for folder in (verify, train): (folder/'STATUS').write_text('COMPLETED\n')
            (verify/'summary.json').write_text(json.dumps(dict(completed=True, job_id='13453180', source_fit_job='13453177',
                verified_old_C0_optimizer_steps=9216, verified_new_prior_optimizer_steps=4608,
                final_checkpoint_strict_load=True, final_checkpoint=str(I.CHECKPOINT.relative_to(I.ROOT)),
                configuration=I.A.IDENTITY, DEV_used=False, sealed_data_used=False)))
            (train/'summary.json').write_text(json.dumps(dict(completed=True, job_id='13453184', checkpoint=str(I.CHECKPOINT),
                status='original_TRAIN_end_to_end_diagnostic_not_DEV_or_generalization', DEV_used=False, sealed_data_used=False)))
            with patch.object(I, 'VERIFY', verify), patch.object(I, 'TRAIN_EVAL', train), \
                 patch.object(I.C, 'require_train_closure') as other:
                I.require_train_closure(); other.assert_called_once()
                (train/'STATUS').write_text('FAILED exit=1\n')
                with patch.object(I.I.R.M, 'read_allowed') as DNA, patch.object(I.A, 'main') as model:
                    with self.assertRaises(ValueError): I.main(0)
                    DNA.assert_not_called(); model.assert_not_called()
            with patch.object(I.C, 'require_train_closure', side_effect=ValueError('synthetic old TRAIN unfinished')), \
                 patch.object(I.I.R.M, 'read_allowed') as DNA, patch.object(I.A, 'main') as model:
                with self.assertRaises(ValueError): I.main(0)
                DNA.assert_not_called(); model.assert_not_called()

    def test_driver_keeps_uppercase_IUPAC_and_fixed_raw_arguments(self):
        species, seqid, DNA = 'synthetic', 'artificial', 'acgtRynnacgt'
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); folder = root/seqid; folder.mkdir()
            def source(path, allowed):
                self.assertEqual(path, I.ROOT/'data/m1_screen'/species/'genome.fa')
                self.assertEqual(allowed, {seqid:len(DNA)})
                return {seqid:DNA}
            def raw(fasta, checkpoint, out):
                self.assertEqual(fasta.read_text(), '>'+seqid+'\n'+DNA.upper()+'\n')
                self.assertEqual((checkpoint, out), (I.CHECKPOINT, folder/'prediction'))
            with patch.object(I, 'OUT', root), patch.object(I.I, 'SCOPE', ((species, seqid, len(DNA)),)), \
                 patch.object(I, 'require_train_closure') as closure, \
                 patch.object(I.I.R.M, 'read_allowed', side_effect=source) as DNA_source, patch.object(I.A, 'main', side_effect=raw) as model:
                I.main(0); closure.assert_called_once(); DNA_source.assert_called_once(); model.assert_called_once()
            receipt = json.loads((folder/'summary.json').read_text())
            self.assertTrue(receipt['DEV_used'] and receipt['all_four_TRAIN_branches_closed'])
            self.assertFalse(receipt['sealed_data_used'] or receipt['reference_used'])

    def test_candidate_identity_and_no_reference_before_input_reconciliation(self):
        species, seqid, length = 'synthetic', 'artificial', 40
        scope = dict(seqid=seqid, windows=[[0, length]], reference_used=False, joint_whole_oriented_contig_decode=True,
            context_policy=I.A.R.CONTEXT, pruning='exact_same_stage_actual_carry_lower_exon_budget_dominance_max_only',
            score_accumulation='FP64_per_base_forward_order',
            tie_policy='background_continuation_then_state_table_arrival_lower_count_dominates_equal_score')
        outer = dict(inference_complete=True, task=0, species=species, seqid=seqid, length=length,
            checkpoint=str(I.CHECKPOINT), raw_scope=I.I.RAW_SCOPE, DEV_used=True, sealed_data_used=False,
            reference_used=False, labels_or_feature_cache_used=False, optimizer_steps_this_run=0)
        raw = dict(inference_complete=True, checkpoint=str(I.CHECKPOINT), protocol=I.A.PROTOCOL,
            lengths={seqid:length}, configuration=I.A.IDENTITY, reference_used=False, labels_or_feature_cache_used=False,
            observation_model_used_at_inference=False, optimizer_steps_this_run=0, backend='exact_exon_budget_skyline',
            scopes=[dict(scope, strand='+'), dict(scope, strand='-')])
        kwargs = dict(checkpoint=I.CHECKPOINT, identity=I.A.IDENTITY, protocol=I.A.PROTOCOL)
        E.check_candidate(outer, raw, 0, species, seqid, length, **kwargs)
        changed = copy.deepcopy(raw); changed['configuration']['native_head_indices'] = [0]
        with self.assertRaises(ValueError): E.check_candidate(outer, changed, 0, species, seqid, length, **kwargs)
        with self.assertRaises(ValueError): E.check_candidate(outer, raw, 0, species, seqid, length)
        with patch.object(E, 'reconcile_inputs', side_effect=ValueError('synthetic incomplete input')) as inputs, \
             patch.object(E.C, 'select_reference') as reference:
            with self.assertRaises(ValueError): E.main(Path('/synthetic'), I.CHECKPOINT, I.A.IDENTITY, I.A.PROTOCOL, I.METHOD, Path('/synthetic/out'))
            inputs.assert_called_once_with(Path('/synthetic'), I.CHECKPOINT, I.A.IDENTITY, I.A.PROTOCOL)
            reference.assert_not_called()

    def test_full_family_denominators_aggregation_contrasts_and_changed_baseline(self):
        reports = synthetic_reports(); value = B.build_comparison(reports)
        self.assertEqual(len(value['exact_chain_table']), 8)
        self.assertEqual(value['pooled_reference_primary_chains'], 7728)
        self.assertEqual(value['matched_primary_control'], B.CONTROL)
        contrast = value['native_prior_differences_percentage_points'][B.CONTROL]
        self.assertAlmostEqual(contrast['pooled_micro'], 100*400/(3000+7728), places=12)
        for index, s in enumerate(E.EXPECTED):
            self.assertAlmostEqual(contrast['per_species'][s], 100*200/((2000 if index == 0 else 1000)+E.EXPECTED[s]), places=12)
        missing = copy.deepcopy(reports); missing.pop('M28_LOCAL_CHANNEL_final9216')
        with self.assertRaises(ValueError): B.build_comparison(missing)
        changed = copy.deepcopy(reports)
        s = next(iter(E.EXPECTED)); old = changed[I.METHOD]['species'][s]['methods']['ANNEVO']['exact_chain']
        changed[I.METHOD]['species'][s]['methods']['ANNEVO']['exact_chain'] = E.C.S.prf(old['tp']+1, old['predicted'], old['reference'])
        with self.assertRaisesRegex(ValueError, 'Fixed baseline'): B.build_comparison(changed)
        changed = copy.deepcopy(reports); changed[I.METHOD]['pooled'][I.METHOD]['exact_chain']['f1'] += .01
        with self.assertRaisesRegex(ValueError, 'aggregation'): B.build_comparison(changed)


if __name__ == '__main__':
    began = time.perf_counter()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(PreparedDEVChecks))
    with (OUT/'summary.json').open('x') as handle:
        json.dump(dict(completed=tests.wasSuccessful(), job_id=os.environ['SLURM_JOB_ID'], plan=I.PLAN,
            scope='artificial_native_prior_DEV_and_complete_frozen_family_only', tests_run=tests.testsRun,
            failures=len(tests.failures), errors=len(tests.errors), research_DNA_or_labels_opened=False,
            historical_DEV_receipts_or_predictions_opened=False, checkpoint_loaded=False, actual_native_GLM_forward=False,
            optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False, seconds=time.perf_counter()-began), handle, indent=2, allow_nan=False)
        handle.write('\n')
    if not tests.wasSuccessful(): raise SystemExit(1)
