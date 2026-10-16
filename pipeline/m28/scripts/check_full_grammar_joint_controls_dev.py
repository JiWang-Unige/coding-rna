#!/usr/bin/env python3
"""Artificial-only three-candidate DEV wiring and legacy-default regression."""
import copy
from contextlib import ExitStack
import json
import os
from pathlib import Path
import tempfile
import time
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import torch
import infer_full_grammar_joint_controls_dev as I
import evaluate_full_grammar_joint_controls_dev as V
import aggregate_full_grammar_joint_controls_dev as G
import check_full_grammar_sampling_context_dev as Old

OUT = I.ROOT/'outputs/M28-FULL-GRAMMAR-JOINT-CONTROLS-DEV-CHECKS-R1'
E = G.E


def write_completed(folder, job, **fields):
    folder.mkdir(parents=True)
    (folder/'STATUS').write_text('COMPLETED\n')
    (folder/'JOBID').write_text(job+'\n')
    (folder/'summary.json').write_text(json.dumps(dict(completed=True,job_id=job,
        DEV_used=False,sealed_data_used=False,**fields)))


def synthetic_reports():
    template = Old.Old.synthetic_reports()['M28_final4608']
    reports = {}
    for index, (branch,label) in enumerate(zip(I.BRANCHES,I.METHODS)):
        train, _ = I.candidate(branch)
        report = copy.deepcopy(template)
        report['checkpoint'] = str(train.CHECKPOINT)
        report['job_id'] = str(80000+index)
        for species,item in report['species'].items():
            item['methods'][label] = item['methods'].pop('M28_final4608')
            old = item['methods'][label]['exact_chain']
            item['methods'][label]['exact_chain'] = E.C.S.prf(400+50*index,old['predicted'],old['reference'])
            for row in item['methods'].values():
                counts = row['exact_chain']
                row.update(CDS_base_unstranded_union=E.C.S.prf(*(3*counts[k] for k in ('tp','predicted','reference'))),
                    background=dict(bp=1000000,predicted_span_bp=1000,predicted_CDS_bp=500,
                        wholly_background_unique_chains=2,span_FPR=.001,wholly_background_chains_per_Mb=2),
                    prediction_records=counts['predicted'],exact_other_assessable_isoform_not_primary_TP=7,
                    phase_chain_consistent_prediction_records=counts['predicted'],
                    descriptive_structure_metrics_not_reference_filters=dict(
                        single_CDS_segment=E.C.S.prf(counts['tp']//2,counts['predicted']//2,counts['reference']//2),
                        span_longer_than_native_window=E.C.S.prf(0,0,1)))
        report['pooled'] = {}
        for method in (label,*E.N.METHODS):
            rows = [report['species'][s]['methods'][method] for s in E.EXPECTED]
            f1s = [row['exact_chain']['f1'] for row in rows]
            report['pooled'][method] = dict(E.aggregate(rows),macro_species_exact_chain_F1=sum(f1s)/2,
                worst_species_exact_chain_F1=min(f1s))
        reports[label] = report
    return reports


class JointControlsChecks(unittest.TestCase):
    def test_all_three_TRAIN_closures_block_DEV_and_failed_slurm(self):
        with tempfile.TemporaryDirectory() as directory, ExitStack() as stack:
            root = Path(directory)
            specs = {}
            jobs = set()
            for index,branch in enumerate(I.BRANCHES):
                train,fit,verify,evaluation,paired,_ = I.SPECS[branch]
                for attr in ('FIT','VERIFY','EVAL','OUT'):
                    stack.enter_context(patch.object(train,attr,root/branch/attr))
                stack.enter_context(patch.object(train,'source'))
                pair_root = root/branch/'PAIRED'
                specs[branch] = train,fit,verify,evaluation,paired,pair_root
                write_completed(train.FIT,fit)
                write_completed(train.VERIFY,verify)
                write_completed(train.EVAL,evaluation,checkpoint=str(train.CHECKPOINT),
                    status='original_TRAIN_end_to_end_diagnostic_not_DEV_or_generalization',
                    chromosomes={row[1]:{} for row in I.B.I.SCOPE},optimizer_steps_this_run=0)
                if index == 0:
                    pair = dict(all_four_caller_inputs_and_normalized_original_DNA_directly_equal=True,
                        counts_and_species_pooled_PRF_reconciled=True,
                        source_jobs=dict(F00_original_split='13451625',F10_fresh_split='13453583',
                            F01_original_joint='13453780',F11_fresh_joint=evaluation))
                else:
                    arm,control = ('uniform','informative') if index == 1 else ('codon','baseline')
                    pair = dict(all_paired_caller_inputs_and_normalized_original_DNA_directly_equal=True,
                        counts_and_PRF_reconciled=True,source_jobs={control:'13494898',arm:evaluation},
                        source_fit_jobs={control:'13494895',arm:fit})
                    pair['DNA_comparison_before_'+arm+'_reference_access'] = True
                write_completed(pair_root,paired,optimizer_steps_this_run=0,**pair)
                jobs.update((fit,verify,evaluation,paired))
                for task,(species,seqid,length) in enumerate(I.B.I.SCOPE):
                    job = str(90000+6*index+task)
                    write_completed(train.OUT/seqid,job,inference_complete=True,task=task,species=species,
                        seqid=seqid,length=length,source_fit_job=fit,checkpoint=str(train.CHECKPOINT),
                        reference_used=False,labels_or_feature_cache_used=False,optimizer_steps_this_run=0)
                    jobs.add(job)
            stack.enter_context(patch.object(I,'SPECS',specs))
            rows = '\n'.join(job+suffix+'|COMPLETED|0:0|' for job in sorted(jobs)
                for suffix in ('','.batch','.extern'))
            run = stack.enter_context(patch.object(I.subprocess,'run',return_value=SimpleNamespace(stdout=rows)))
            I.require_train_closure()
            self.assertEqual(len(jobs),30)
            self.assertEqual(set(run.call_args.args[0][2].split(',')),jobs)
            source = specs[I.BRANCHES[0]][-1]/'summary.json'
            saved_source = json.loads(source.read_text())
            for changed in (dict(F00='13451625',F10='13453583',F01='13453780',F11='13494898'),
                            dict(saved_source['source_jobs'],F11_fresh_joint='13451625')):
                source.write_text(json.dumps(dict(saved_source,source_jobs=changed)))
                with patch.object(I.I.R.M,'read_allowed') as DNA,patch.object(I.B.A,'main') as raw:
                    with self.assertRaisesRegex(ValueError,'TRAIN comparison'):
                        I.main(I.BRANCHES[0],0)
                    DNA.assert_not_called(); raw.assert_not_called()
            source.write_text(json.dumps(saved_source))
            for branch in I.BRANCHES:
                train,_ = I.candidate(branch)
                (train.EVAL/'STATUS').write_text('FAILED exit=1\n')
                with patch.object(I.I.R.M,'read_allowed') as DNA,patch.object(train.A,'main') as raw:
                    with self.assertRaises(ValueError):
                        I.main(branch,0)
                    DNA.assert_not_called(); raw.assert_not_called()
                (train.EVAL/'STATUS').write_text('COMPLETED\n')
            run.return_value.stdout = rows.replace('13497947.batch|COMPLETED|0:0|',
                '13497947.batch|FAILED|1:0|')
            with patch.object(I.I.R.M,'read_allowed') as DNA,patch.object(I.D.A,'main') as raw:
                with self.assertRaisesRegex(ValueError,'Slurm terminal'):
                    I.main('CODON_ENERGY',0)
                DNA.assert_not_called(); raw.assert_not_called()

    def test_all_three_raw_drivers_keep_IUPAC_and_fixed_last_weights(self):
        species,seqid,sequence = 'synthetic','artificial','acgtRynNacgt'
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for branch in I.BRANCHES:
                train,_ = I.candidate(branch)
                folder = root/branch/seqid
                folder.mkdir(parents=True)
                def raw(fasta,checkpoint,out):
                    self.assertEqual(fasta.read_text(),'>'+seqid+'\n'+sequence.upper()+'\n')
                    self.assertEqual((checkpoint,out),(train.CHECKPOINT,folder/'prediction'))
                with patch.object(I,'OUT',root),patch.object(I.I,'SCOPE',((species,seqid,len(sequence)),)), \
                     patch.object(I,'require_train_closure') as closure, \
                     patch.object(I.I.R.M,'read_allowed',return_value={seqid:sequence}) as DNA, \
                     patch.object(train.A,'main',side_effect=raw) as model:
                    I.main(branch,0)
                    closure.assert_called_once()
                    DNA.assert_called_once_with(I.ROOT/'data/m1_screen'/species/'genome.fa',{seqid:len(sequence)})
                    model.assert_called_once()
                receipt = json.loads((folder/'summary.json').read_text())
                self.assertEqual(receipt['source_fit_job'],I.SPECS[branch][1])
                self.assertTrue(receipt['all_three_TRAIN_flows_closed'] and receipt['DEV_used'])
                self.assertFalse(receipt['reference_used'] or receipt['sealed_data_used'])

    def test_branch_identity_full_codon_scores_and_coverage(self):
        species,seqid,length = 'synthetic','artificial',40
        scores = torch.arange(64,dtype=torch.float32)
        for branch in I.BRANCHES:
            train,_ = I.candidate(branch)
            outer = dict(inference_complete=True,task=0,species=species,seqid=seqid,length=length,
                branch=branch,checkpoint=str(train.CHECKPOINT),source_fit_job=I.SPECS[branch][1],
                plan=I.PLAN,all_three_TRAIN_flows_closed=True,raw_scope=I.I.RAW_SCOPE,
                DEV_used=True,sealed_data_used=False,reference_used=False,labels_or_feature_cache_used=False,
                optimizer_steps_this_run=0)
            scope = dict(seqid=seqid,windows=[[0,length]],reference_used=False,joint_whole_oriented_contig_decode=True,
                context_policy=I.I.R.CONTEXT,pruning='exact_same_stage_actual_carry_lower_exon_budget_dominance_max_only',
                score_accumulation='FP64_per_base_forward_order',
                tie_policy='background_continuation_then_state_table_arrival_lower_count_dominates_equal_score')
            raw = dict(inference_complete=True,checkpoint=str(train.CHECKPOINT),protocol=train.A.PROTOCOL,
                lengths={seqid:length},configuration=copy.deepcopy(train.A.IDENTITY),reference_used=False,
                labels_or_feature_cache_used=False,observation_model_used_at_inference=False,
                optimizer_steps_this_run=0,backend='exact_exon_budget_skyline')
            if branch == 'CODON_ENERGY':
                scope.update(score_accumulation=I.D.ACCUMULATION,backend=I.D.BACKEND)
                raw.update(backend=I.D.BACKEND,score_kind=I.D.SCORE_KIND,raw_codon_scores=scores.tolist(),
                    centered_codon_scores=I.D.centered_codon_scores(scores).tolist())
            raw['scopes'] = [dict(scope,strand='+'),dict(scope,strand='-')]
            args = dict(codon_scores=scores,checkpoint=train.CHECKPOINT,identity=train.A.IDENTITY,protocol=train.A.PROTOCOL)
            V.check_candidate(branch,outer,raw,0,species,seqid,length,**args)
            bad = dict(outer,source_fit_job='1')
            with self.assertRaises(ValueError):
                V.check_candidate(branch,bad,raw,0,species,seqid,length,**args)
            bad = copy.deepcopy(raw); bad['scopes'][1]['windows'] = [[1,length]]
            with self.assertRaises(ValueError):
                V.check_candidate(branch,outer,bad,0,species,seqid,length,**args)
            bad = copy.deepcopy(raw); bad['configuration']['sampling_seed'] = 99
            with self.assertRaises(ValueError):
                V.check_candidate(branch,outer,bad,0,species,seqid,length,**args)
            if branch == 'CODON_ENERGY':
                for key,value in (('backend','exact_exon_budget_skyline'),('protocol','wrong'),
                                  ('raw_codon_scores',[0.]*64),('centered_codon_scores',[0.]*64)):
                    bad = copy.deepcopy(raw); bad[key] = value
                    with self.assertRaises(ValueError):
                        V.check_candidate(branch,outer,bad,0,species,seqid,length,**args)

    def test_real_evaluator_hook_rejects_before_DNA_or_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            seqid = E.I.SCOPE[0][1]
            folder = root/seqid
            (folder/'prediction').mkdir(parents=True)
            (folder/'STATUS').write_text('COMPLETED\n')
            (folder/'summary.json').write_text('{}')
            (folder/'prediction/summary.json').write_text('{}')
            with patch.object(E,'check_candidate') as legacy, \
                 patch.object(E.C.M,'read_allowed') as DNA,patch.object(E.C,'select_reference') as reference:
                def reject(*args,**kwargs):
                    raise ValueError('synthetic explicit identity rejection')
                with self.assertRaisesRegex(ValueError,'explicit identity'):
                    E.reconcile_inputs(root,candidate_validator=reject)
                legacy.assert_not_called(); DNA.assert_not_called(); reference.assert_not_called()
        for branch in I.BRANCHES:
            with patch.object(I,'require_train_closure',side_effect=ValueError('synthetic not closed')), \
                 patch.object(V.torch,'load') as checkpoint,patch.object(E,'main') as evaluate:
                with self.assertRaises(ValueError):
                    V.main(branch)
                checkpoint.assert_not_called(); evaluate.assert_not_called()

    def test_complete_family_counts_baselines_and_worst_not_min_delta(self):
        reports = synthetic_reports()
        B,U = I.METHODS[:2]
        for label,tps in ((B,(2175,1604)),(U,(4349,1146))):
            report = reports[label]
            for species,tp in zip(E.EXPECTED,tps):
                row = report['species'][species]['methods'][label]
                row['exact_chain'] = E.C.S.prf(tp,E.EXPECTED[species],E.EXPECTED[species])
            rows = [report['species'][s]['methods'][label] for s in E.EXPECTED]
            f1s = [row['exact_chain']['f1'] for row in rows]
            report['pooled'][label] = dict(E.aggregate(rows),macro_species_exact_chain_F1=sum(f1s)/2,
                worst_species_exact_chain_F1=min(f1s))
        result = G.build_comparison(reports)
        self.assertEqual(result['methods'],list(G.METHODS))
        self.assertEqual(len(result['exact_chain_table']),6)
        self.assertEqual(result['pooled_reference_primary_chains'],7728)
        delta = result['candidate_differences_percentage_points'][U][B]
        self.assertGreater(delta['worst_species'],10)
        self.assertLess(min(delta['per_species'].values()),-19)
        self.assertAlmostEqual(delta['macro_species'],sum(delta['per_species'].values())/2)
        changed = copy.deepcopy(reports); changed.pop(I.METHODS[2])
        with self.assertRaises(ValueError):
            G.build_comparison(changed)
        changed = copy.deepcopy(reports)
        changed[U]['species']['arabidopsis_thaliana']['methods'][U]['exact_chain']['f1'] += .01
        with self.assertRaisesRegex(ValueError,'PRF'):
            G.build_comparison(changed)
        changed = copy.deepcopy(reports)
        changed[U]['species']['arabidopsis_thaliana']['methods']['ANNEVO']['prediction_records'] += 1
        with self.assertRaisesRegex(ValueError,'baseline'):
            G.build_comparison(changed)
        changed = copy.deepcopy(reports); changed[U]['pooled'][U]['background']['predicted_span_bp'] += 1
        with self.assertRaisesRegex(ValueError,'aggregation'):
            G.build_comparison(changed)

    def test_actual_writer_preserves_old_families_in_artificial_tree(self):
        reports = synthetic_reports()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            out = root/'new'; out.mkdir()
            sentinel = root/'old_ten_method_summary.json'
            sentinel.write_text('synthetic frozen ten\n')
            for branch,label in zip(I.BRANCHES,I.METHODS):
                folder = root/'comparisons'/branch
                folder.mkdir(parents=True)
                (folder/'STATUS').write_text('COMPLETED\n')
                (folder/'JOBID').write_text(reports[label]['job_id']+'\n')
                (folder/'summary.json').write_text(json.dumps(reports[label]))
            with patch.object(G,'OUT',out),patch.object(I,'COMPARISON',root/'comparisons'), \
                 patch.object(I,'require_train_closure') as closure:
                G.main(); closure.assert_called_once()
            report = json.loads((out/'summary.json').read_text())
            self.assertEqual(report['methods'],list(G.METHODS))
            self.assertEqual(set(report['source_comparison_jobs']),set(I.METHODS))
            self.assertTrue(report['counts_base_background_and_micro_macro_worst_reconciled'])
            self.assertEqual(sentinel.read_text(),'synthetic frozen ten\n')
            self.assertTrue(all(method in (out/'comparison.md').read_text() for method in G.METHODS))


if __name__ == '__main__':
    began = time.perf_counter()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(JointControlsChecks)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        # The failed research launch exposed a real producer/consumer naming mismatch.
        # Check completed TRAIN receipts and actual Slurm only, before any DEV read.
        I.require_train_closure()
    with (OUT/'summary.json').open('x') as handle:
        json.dump(dict(completed=result.wasSuccessful(),job_id=os.environ['SLURM_JOB_ID'],plan=I.PLAN,
            scope='six_changed_interface_fixtures_and_real_completed_TRAIN_closure_no_DEV',tests_run=result.testsRun,
            failures=len(result.failures),errors=len(result.errors),research_DNA_or_labels_opened=False,
            historical_DEV_receipts_or_predictions_opened=False,checkpoint_loaded=False,
            actual_native_GLM_forward=False,DEV_used=False,sealed_data_used=False,
            real_TRAIN_receipts_and_Slurm_read=result.wasSuccessful(),
            all_three_real_TRAIN_flows_closed_before_DEV=result.wasSuccessful(),
            optimizer_steps_this_run=0,seconds=time.perf_counter()-began),handle,indent=2,allow_nan=False)
        handle.write('\n')
    if not result.wasSuccessful():
        raise SystemExit(1)
