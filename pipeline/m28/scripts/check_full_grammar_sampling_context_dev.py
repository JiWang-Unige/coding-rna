#!/usr/bin/env python3
"""Artificial-only preparation for new DEV wiring and the separate extension."""
import copy
import json
import os
from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import patch
import infer_full_grammar_sampling_context_dev as I
import evaluate_full_grammar_sampling_context_dev as V
import aggregate_full_grammar_sampling_context_dev as B
import check_full_grammar_native_cds_prior_dev as Old

OUT = I.ROOT/'outputs/M28-FULL-GRAMMAR-SAMPLING-CONTEXT-DEV-CHECKS'
E = B.E


def synthetic_reports():
    reports = Old.synthetic_reports()
    for branch,tp in zip(I.BRANCHES,(550,600)):
        train,label = I.candidate(branch)
        report = copy.deepcopy(reports['M28_final4608'])
        report['checkpoint'] = str(train.CHECKPOINT)
        report['pooled'][label] = report['pooled'].pop('M28_final4608')
        for index,species in enumerate(E.EXPECTED):
            item = report['species'][species]
            item['methods'][label] = item['methods'].pop('M28_final4608')
            item['methods'][label]['exact_chain'] = E.C.S.prf(tp,2000 if index == 0 else 1000,E.EXPECTED[species])
        rows = [report['species'][s]['methods'][label]['exact_chain'] for s in E.EXPECTED]
        report['pooled'][label] = dict(exact_chain=E.C.S.prf(*(sum(r[k] for r in rows)
            for k in ('tp','predicted','reference'))),macro_species_exact_chain_F1=sum(r['f1'] for r in rows)/2,
            worst_species_exact_chain_F1=min(r['f1'] for r in rows))
        reports[label] = report
    for i,report in enumerate(reports.values()):
        report['job_id'] = str(100+i)
    return reports


class PreparationChecks(unittest.TestCase):
    def test_all_six_TRAIN_closures_before_any_DEV_input(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for branch,job in zip(I.BRANCHES,('13453583','13453780')):
                train,_ = I.candidate(branch)
                folder = root/branch; folder.mkdir()
                (folder/'STATUS').write_text('COMPLETED\n')
                (folder/'summary.json').write_text(json.dumps(dict(completed=True,job_id=job,
                    checkpoint=str(train.CHECKPOINT),status='original_TRAIN_end_to_end_diagnostic_not_DEV_or_generalization',
                    chromosomes={r[1]:{} for r in train.I.SCOPE},DEV_used=False,sealed_data_used=False,optimizer_steps_this_run=0)))
            with patch.object(I.N,'require_train_closure') as original_four, \
                 patch.object(I.F,'source') as fresh_source, patch.object(I.J,'source') as joint_source, \
                 patch.object(I.F,'EVAL',root/I.BRANCHES[0]), patch.object(I.J,'EVAL',root/I.BRANCHES[1]):
                I.require_train_closure()
                original_four.assert_called_once()
                fresh_source.assert_called_once_with('13453573')
                joint_source.assert_called_once_with('13453777')
                for branch in I.BRANCHES:
                    folder = root/branch
                    (folder/'STATUS').write_text('FAILED exit=1\n')
                    with patch.object(I.I.R.M,'read_allowed') as DNA, \
                         patch.object(I.F.A,'main') as fresh, patch.object(I.J.A,'main') as joint:
                        with self.assertRaises(ValueError):
                            I.main(I.BRANCHES[0],0)
                        DNA.assert_not_called(); fresh.assert_not_called(); joint.assert_not_called()
                    (folder/'STATUS').write_text('COMPLETED\n')
            with patch.object(I.N,'require_train_closure',side_effect=ValueError('synthetic original four unfinished')), \
                 patch.object(I.I.R.M,'read_allowed') as DNA, patch.object(I.F.A,'main') as raw:
                with self.assertRaises(ValueError):
                    I.main(I.BRANCHES[0],0)
                DNA.assert_not_called(); raw.assert_not_called()

    def test_both_raw_drivers_keep_original_IUPAC_and_fixed_arguments(self):
        species,seqid,sequence = 'synthetic','artificial','acgtRynnacgt'
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for branch in I.BRANCHES:
                train,_ = I.candidate(branch)
                folder = root/branch/seqid; folder.mkdir(parents=True)
                def source(path,allowed):
                    self.assertEqual(path,I.ROOT/'data/m1_screen'/species/'genome.fa')
                    self.assertEqual(allowed,{seqid:len(sequence)})
                    return {seqid:sequence}
                def raw(fasta,checkpoint,out):
                    self.assertEqual(fasta.read_text(),'>'+seqid+'\n'+sequence.upper()+'\n')
                    self.assertEqual((checkpoint,out),(train.CHECKPOINT,folder/'prediction'))
                with patch.object(I,'OUT',root), patch.object(I.I,'SCOPE',((species,seqid,len(sequence)),)), \
                     patch.object(I,'require_train_closure') as closure, \
                     patch.object(I.I.R.M,'read_allowed',side_effect=source) as DNA, \
                     patch.object(train.A,'main',side_effect=raw) as model:
                    I.main(branch,0)
                    closure.assert_called_once(); DNA.assert_called_once(); model.assert_called_once()
                receipt = json.loads((folder/'summary.json').read_text())
                self.assertTrue(receipt['DEV_used'] and receipt['all_six_TRAIN_branches_closed'])
                self.assertEqual(receipt['branch'],branch)
                self.assertFalse(receipt['sealed_data_used'] or receipt['reference_used'])

    def test_new_identities_and_input_reconciliation_before_reference(self):
        species,seqid,length = 'synthetic','artificial',40
        scope = dict(seqid=seqid,windows=[[0,length]],reference_used=False,joint_whole_oriented_contig_decode=True,
            context_policy=I.I.R.CONTEXT,pruning='exact_same_stage_actual_carry_lower_exon_budget_dominance_max_only',
            score_accumulation='FP64_per_base_forward_order',
            tie_policy='background_continuation_then_state_table_arrival_lower_count_dominates_equal_score')
        for branch in I.BRANCHES:
            train,method = I.candidate(branch)
            outer = dict(inference_complete=True,task=0,species=species,seqid=seqid,length=length,
                checkpoint=str(train.CHECKPOINT),raw_scope=I.I.RAW_SCOPE,DEV_used=True,sealed_data_used=False,
                reference_used=False,labels_or_feature_cache_used=False,optimizer_steps_this_run=0)
            raw = dict(inference_complete=True,checkpoint=str(train.CHECKPOINT),protocol=train.A.PROTOCOL,
                lengths={seqid:length},configuration=train.A.IDENTITY,reference_used=False,labels_or_feature_cache_used=False,
                observation_model_used_at_inference=False,optimizer_steps_this_run=0,backend='exact_exon_budget_skyline',
                scopes=[dict(scope,strand='+'),dict(scope,strand='-')])
            args = dict(checkpoint=train.CHECKPOINT,identity=train.A.IDENTITY,protocol=train.A.PROTOCOL)
            E.check_candidate(outer,raw,0,species,seqid,length,**args)
            changed = copy.deepcopy(raw)
            key,value = ('sampling_seed',2) if branch == 'FRESH_DRAWS' else ('encoder_blocks_bp',6144)
            changed['configuration'][key] = value
            with self.assertRaises(ValueError):
                E.check_candidate(outer,changed,0,species,seqid,length,**args)
            with self.assertRaises(ValueError):
                E.check_candidate(outer,raw,0,species,seqid,length)
            with patch.object(I,'require_train_closure') as closure, \
                 patch.object(E,'reconcile_inputs',side_effect=ValueError('synthetic incomplete inputs')) as inputs, \
                 patch.object(E.C,'select_reference') as reference:
                with self.assertRaises(ValueError):
                    V.main(branch)
                closure.assert_called_once()
                inputs.assert_called_once_with(I.OUT/branch,train.CHECKPOINT,train.A.IDENTITY,train.A.PROTOCOL)
                reference.assert_not_called()

    def test_full_extension_counts_contrasts_and_frozen_eight_preservation(self):
        reports = synthetic_reports()
        original = B.B.build_comparison({k:reports[k] for k in B.B.SPECS})
        result = B.build_comparison(reports)
        self.assertEqual(len(result['exact_chain_table']),10)
        self.assertEqual(result['exact_chain_table'][:8],original['exact_chain_table'])
        self.assertEqual(result['pooled_reference_primary_chains'],7728)
        self.assertEqual(result['matched_primary_control'],B.B.CONTROL)
        self.assertTrue(result['original_eight_method_family_preserved'])
        for label,tp in zip(I.METHODS,(550,600)):
            delta = result['new_candidate_differences_percentage_points'][label][B.B.CONTROL]
            self.assertAlmostEqual(delta['pooled_micro'],100*4*(tp-400)/(3000+7728),places=12)
            desired = []
            for index,s in enumerate(E.EXPECTED):
                wanted = 100*2*(tp-400)/((2000 if index == 0 else 1000)+E.EXPECTED[s])
                desired.append(wanted)
                self.assertAlmostEqual(delta['per_species'][s],wanted,places=12)
            self.assertAlmostEqual(delta['macro_species'],sum(desired)/2,places=12)
            self.assertAlmostEqual(delta['worst_species'],desired[0],places=12)
        missing = copy.deepcopy(reports); missing.pop(I.METHODS[0])
        with self.assertRaises(ValueError):
            B.build_comparison(missing)
        bad = copy.deepcopy(reports)
        species = next(iter(E.EXPECTED))
        row = bad[I.METHODS[0]]['species'][species]['methods']['ANNEVO']['exact_chain']
        bad[I.METHODS[0]]['species'][species]['methods']['ANNEVO']['exact_chain'] = E.C.S.prf(row['tp']+1,row['predicted'],row['reference'])
        with self.assertRaisesRegex(ValueError,'Fixed baseline'):
            B.build_comparison(bad)
        bad = copy.deepcopy(reports); bad[I.METHODS[1]]['pooled'][I.METHODS[1]]['exact_chain']['f1'] += .01
        with self.assertRaisesRegex(ValueError,'aggregation'):
            B.build_comparison(bad)
        bad = copy.deepcopy(reports); bad[I.METHODS[1]]['species'][species]['reference_primary_chains'] -= 1
        with self.assertRaisesRegex(ValueError,'denominator'):
            B.build_comparison(bad)

    def test_real_extension_writer_only_creates_separate_artificial_output(self):
        reports = synthetic_reports()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); out = root/'extension'; out.mkdir()
            unchanged = root/'original_family'; unchanged.mkdir()
            sentinel = unchanged/'comparison.md'; sentinel.write_text('synthetic frozen eight\n')
            specs = {}
            for label,(_,checkpoint) in B.SPECS.items():
                folder = root/label; folder.mkdir()
                (folder/'STATUS').write_text('COMPLETED\n')
                (folder/'summary.json').write_text(json.dumps(reports[label]))
                specs[label] = folder,checkpoint
            with patch.object(B,'SPECS',specs), patch.object(B,'OUT',out), \
                 patch.object(B.B,'OUT',unchanged), patch.object(I,'require_train_closure') as closure:
                B.main(); closure.assert_called_once()
            written = json.loads((out/'summary.json').read_text())
            self.assertEqual(written['methods'],list(B.METHODS))
            self.assertEqual(len(written['measurements']),10)
            self.assertTrue(written['original_eight_method_family_preserved'])
            self.assertEqual(sentinel.read_text(),'synthetic frozen eight\n')
            text = (out/'comparison.md').read_text()
            self.assertTrue(all(method in text for method in B.METHODS))


if __name__ == '__main__':
    began = time.perf_counter()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(PreparationChecks))
    with (OUT/'summary.json').open('x') as handle:
        json.dump(dict(completed=result.wasSuccessful(),job_id=os.environ['SLURM_JOB_ID'],plan=I.PLAN,
            scope='artificial_sampling_context_DEV_and_separate_ten_method_extension_only',tests_run=result.testsRun,
            failures=len(result.failures),errors=len(result.errors),research_DNA_or_labels_opened=False,
            historical_DEV_receipts_or_predictions_opened=False,checkpoint_loaded=False,actual_native_GLM_forward=False,
            optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,seconds=time.perf_counter()-began),
            handle,indent=2,allow_nan=False)
        handle.write('\n')
    if not result.wasSuccessful():
        raise SystemExit(1)
