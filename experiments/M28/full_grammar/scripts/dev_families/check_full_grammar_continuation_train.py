#!/usr/bin/env python3
"""Artificial-only completeness, denominator and arithmetic checks for the TRAIN family."""
import copy
import json
import os
import time
import unittest
import aggregate_full_grammar_continuation_train as A


def totals(report):
    report['per_species'] = {s:A.E.aggregate([report['chromosomes'][q] for t,q,_ in A.E.I.SCOPE if t == s])
        for s in A.SPECIES}
    report['pooled'] = A.E.aggregate([report['chromosomes'][q] for _,q,_ in A.E.I.SCOPE])
    values = [report['per_species'][s]['exact_chain']['f1'] for s in A.SPECIES]
    report['macro_species_exact_chain_F1'] = sum(values)/len(values)
    report['worst_species_exact_chain_F1'] = min(values)


def fixtures():
    reports = {}
    for index,(label,(_,checkpoint,job,_)) in enumerate(A.SPECS.items()):
        rows = {}
        for species,seqid,length in A.E.I.SCOPE:
            count = A.REFERENCE_COUNTS[seqid]
            rows[seqid] = dict(species=species,seqid=seqid,length=length,
                exact_chain=A.E.C.S.prf(count//2+index,count,count),
                CDS_base_unstranded_union=A.E.C.S.prf(40+index,100,100),
                reference_primary_records=count,reference_assessable_isoform_records=count+2,
                retained_unsupported_reference_strata={'nonstandard_splice':2},
                prediction_records=count,duplicate_chain_records=0,phase_checked_prediction_records=count,
                exact_other_assessable_isoform_not_primary_TP=1,
                background=dict(definition='outside_all_annotated_gene_feature_spans_either_strand',
                    bp=1000,predicted_span_bp=10+index,predicted_CDS_bp=5,wholly_background_unique_chains=index),
                raw_inference_loop_seconds=1,native_encoding_seconds=.5,joint_decoding_seconds=.5)
        report = dict(completed=True,status='original_TRAIN_end_to_end_diagnostic_not_DEV_or_generalization',
            job_id=job,checkpoint=str(checkpoint),
            reference_policy='M26_CDS_assessable_longest_eligible_CDS_tie_transcript_ID',
            reference_used_in_inference=False,training_eligibility_does_not_filter_reference=True,
            unsupported_reference_strata_retained=True,DEV_used=False,sealed_data_used=False,
            optimizer_steps_this_run=0,chromosomes=rows)
        totals(report)
        reports[label] = report
    return reports


class FamilyChecks(unittest.TestCase):
    def test_complete_family_micro_macro_worst_and_percentage_points(self):
        reports = fixtures(); saved = copy.deepcopy(reports)
        result = A.build_comparison(reports)
        self.assertEqual(reports,saved)
        self.assertEqual(result['methods'],list(A.SPECS))
        self.assertEqual(result['pooled_reference_primary_chains'],22688)
        self.assertEqual(result['reference_primary_chains'],{'arabidopsis_thaliana':17739,'oryza_sativa':4949})
        label = 'M28_LOCAL_CHANNEL_final9216'
        self.assertAlmostEqual(result['differences_from_control_percentage_points'][label]['pooled_micro'],
            600/22688,delta=1e-12)
        self.assertFalse(result['differences_from_control_percentage_points'][A.ANCHOR]['cumulative_optimizer_steps_matched'])
        self.assertTrue(result['differences_from_control_percentage_points'][label]['cumulative_optimizer_steps_matched'])

    def test_missing_candidate_or_chromosome_rejected(self):
        reports = fixtures(); del reports['M28_FRESH_DRAWS_final9216']
        with self.assertRaisesRegex(ValueError,'All seven'): A.build_comparison(reports)
        reports = fixtures(); del reports[A.CONTROL]['chromosomes']['NC_089046.1']
        with self.assertRaisesRegex(ValueError,'scope'): A.build_comparison(reports)

    def test_checkpoint_and_job_swap_rejected(self):
        for field in ('checkpoint','job_id'):
            with self.subTest(field=field):
                reports = fixtures(); reports[A.CONTROL][field] = reports[A.ANCHOR][field]
                with self.assertRaisesRegex(ValueError,'identity'): A.build_comparison(reports)

    def test_changed_reference_metadata_rejected_even_if_reaggregated(self):
        for field in ('primary','CDS','background','strata','isoforms'):
            with self.subTest(field=field):
                reports = fixtures(); report = reports[A.CONTROL]; row = report['chromosomes']['NC_003070.9']
                if field == 'primary':
                    row['exact_chain'] = A.E.C.S.prf(row['exact_chain']['tp'],row['exact_chain']['predicted'],7136)
                    row['reference_primary_records'] = 7136
                elif field == 'CDS': row['CDS_base_unstranded_union'] = A.E.C.S.prf(41,100,101)
                elif field == 'background': row['background']['bp'] += 1
                elif field == 'strata': row['retained_unsupported_reference_strata'] = {}
                else: row['reference_assessable_isoform_records'] += 1
                totals(report)
                with self.assertRaisesRegex(ValueError,'reference denominator'): A.build_comparison(reports)

    def test_corrupt_counts_micro_macro_worst_or_phase_rejected(self):
        for field in ('chromosome','species','pooled','macro','worst','phase'):
            with self.subTest(field=field):
                reports = fixtures(); report = reports[A.CONTROL]
                if field == 'chromosome': report['chromosomes']['NC_003070.9']['exact_chain']['f1'] += .01
                elif field == 'species': report['per_species'][A.SPECIES[0]]['exact_chain']['tp'] += 1
                elif field == 'pooled': report['pooled']['exact_chain']['f1'] += .01
                elif field == 'macro': report['macro_species_exact_chain_F1'] += .01
                elif field == 'worst': report['worst_species_exact_chain_F1'] += .01
                else: report['chromosomes']['NC_003070.9']['phase_checked_prediction_records'] -= 1
                with self.assertRaises(ValueError): A.build_comparison(reports)

    def test_data_scope_or_reference_inference_flag_rejected(self):
        for field in ('DEV_used','sealed_data_used','reference_used_in_inference'):
            with self.subTest(field=field):
                reports = fixtures(); reports[A.CONTROL][field] = True
                with self.assertRaisesRegex(ValueError,'identity'): A.build_comparison(reports)


if __name__ == '__main__':
    began = time.perf_counter()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(FamilyChecks))
    with (A.CHECKS/'summary.json').open('x') as handle:
        json.dump(dict(completed=result.wasSuccessful(),job_id=os.environ['SLURM_JOB_ID'],tests_run=result.testsRun,
            failures=len(result.failures),errors=len(result.errors),research_data_read=False,models_or_weights_read=False,
            native_forward_run=False,optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,
            seconds=time.perf_counter()-began),handle,indent=2,allow_nan=False)
        handle.write('\n')
    if not result.wasSuccessful():
        raise SystemExit(1)
