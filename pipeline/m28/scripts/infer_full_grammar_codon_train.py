#!/usr/bin/env python3
"""Strict final9216 codon raw TRAIN, paired DNA check, unchanged CDS ruler."""
import argparse
import json
import os
import re
import time
import torch
import infer_full_grammar_original_train as I
import infer_full_grammar_fresh_joint24_train as B
import annotate_full_grammar_codon as A
import evaluate_full_grammar_original_train as E
from src.m28.full_grammar_codon import centered_codon_scores

ROOT = I.ROOT
FIT = ROOT/'outputs/M28-FULL-GRAMMAR-CODON-FIT'
VERIFY = ROOT/'outputs/M28-FULL-GRAMMAR-CODON-VERIFICATION'
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-CODON-TRAIN-RAW'
EVAL = ROOT/'outputs/M28-FULL-GRAMMAR-CODON-TRAIN-DIAGNOSTIC'
CHECKPOINT = FIT/'step_009216.pt'
BACKEND = 'exact_codon_exon_budget_skyline'
ACCUMULATION = 'FP64_per_base_forward_order_role_plus_completed_codon'
SCORE_KIND = 'unit_fifteen_role_logit_minus_background_plus_complete_spliced_CDS_codon_energy'


def source(fit_job):
    fitted = json.loads((FIT/'summary.json').read_text())
    verified = json.loads((VERIFY/'summary.json').read_text())
    if ((FIT/'STATUS').read_text().strip() != 'COMPLETED' or not fitted['completed']
            or fitted['job_id'] != fit_job or fitted['cumulative_head_optimizer_steps'] != 9216
            or fitted['codon_optimizer_steps'] != 4608 or fitted['optimizer_steps_this_run'] != 4608
            or fitted['fresh_draws'] != 4608 or fitted['fresh_unique_windows'] != 4219
            or fitted['DEV_used'] or fitted['sealed_data_used']
            or any(fitted['configuration'].get(k) != v for k, v in A.IDENTITY.items())
            or fitted['final_checkpoint'] != str(CHECKPOINT.relative_to(ROOT))):
        raise ValueError('Completed fixed codon fresh/joint final9216 fit required')
    if ((VERIFY/'STATUS').read_text().strip() != 'COMPLETED' or not verified['completed']
            or verified['source_fit_job'] != fit_job or verified['configuration'] != fitted['configuration']
            or not verified['fixed_draws_original_observations_and_objective_arithmetic_match']
            or not verified['all_checkpoint_parameters_and_moments_finite_and_correct_age']
            or not verified['final_checkpoint_strict_load'] or verified['verified_final_head_steps'] != 9216
            or verified['verified_codon_steps'] != 4608 or verified['verified_additional_optimizer_steps'] != 4608
            or verified['final_checkpoint'] != fitted['final_checkpoint']
            or verified['DEV_used'] or verified['sealed_data_used']):
        raise ValueError('Completed codon exposure and checkpoint verification required')


def validate_raw_report(raw, seqid, length, expected_scores):
    expected_windows = [[a, min(a+I.R.M.WIDTH, length)] for a in I.R.M.starts(length)]
    if (not raw['inference_complete'] or raw['checkpoint'] != str(CHECKPOINT)
            or raw['lengths'] != {seqid:length} or raw['reference_used'] or raw['labels_or_feature_cache_used']
            or raw['observation_model_used_at_inference'] or raw['optimizer_steps_this_run'] != 0
            or raw['backend'] != BACKEND or raw['score_kind'] != SCORE_KIND
            or raw['raw_codon_scores'] != expected_scores.tolist()
            or raw['centered_codon_scores'] != centered_codon_scores(expected_scores).tolist()
            or any(raw['configuration'].get(k) != v for k, v in A.IDENTITY.items())
            or len(raw['scopes']) != 2 or {s['strand'] for s in raw['scopes']} != {'+', '-'}):
        raise ValueError('Strict final codon checkpoint, coefficients or raw identity differs')
    for scope in raw['scopes']:
        if (scope['seqid'] != seqid or scope['windows'] != expected_windows or scope['reference_used']
                or not scope['joint_whole_oriented_contig_decode'] or scope['context_policy'] != I.R.CONTEXT
                or scope['pruning'] != 'exact_same_stage_actual_carry_lower_exon_budget_dominance_max_only'
                or scope['score_accumulation'] != ACCUMULATION or scope['backend'] != BACKEND):
            raise ValueError('Whole-contig coverage, codon scoring or original structural grammar differs')


def reconcile_inputs():
    state = torch.load(CHECKPOINT, map_location='cpu', weights_only=True)
    A.validate_checkpoint(state)
    scores = state['raw_codon_scores']
    runs = {}
    for task, (species, seqid, length) in enumerate(I.SCOPE):
        folder = OUT/seqid
        if (folder/'STATUS').read_text().strip() != 'COMPLETED':
            raise ValueError('All six codon TRAIN inference tasks must complete before reference access')
        outer = json.loads((folder/'summary.json').read_text())
        raw_report = json.loads((folder/'prediction/summary.json').read_text())
        if (not outer['inference_complete'] or (outer['task'], outer['species'], outer['seqid'], outer['length']) != (task, species, seqid, length)
                or outer['checkpoint'] != str(CHECKPOINT) or outer['DEV_used'] or outer['sealed_data_used']
                or outer['reference_used'] or outer['labels_or_feature_cache_used'] or outer['optimizer_steps_this_run'] != 0):
            raise ValueError('Original TRAIN-only codon inference identity differs')
        validate_raw_report(raw_report, seqid, length, scores)
        gff = folder/'prediction/full_grammar.gff3'
        if not gff.is_file() or not gff.stat().st_size:
            raise ValueError('Raw codon GFF3 missing')
        runs[seqid] = raw_report
    return runs


def raw(task, fit_job):
    if not 0 <= task < 6:
        raise ValueError('One original-six TRAIN chromosome required')
    source(fit_job)
    began = time.perf_counter()
    species, seqid, length = I.SCOPE[task]
    out = OUT/seqid
    with (out/'started.json').open('x') as handle:
        json.dump(dict(job_id=os.environ['SLURM_JOB_ID'], task=task, species=species, seqid=seqid,
            length=length, checkpoint=str(CHECKPOINT), source_fit_job=fit_job,
            reference_used=False, DEV_used=False), handle, indent=2)
    sequence = I.R.M.read_allowed(ROOT/'data/m1_screen'/species/'genome.fa', {seqid:length})[seqid]
    sequence = re.sub('[^ACGT]', 'N', sequence.upper())
    fasta = out/'input.fasta'
    with fasta.open('x') as handle:
        handle.write('>'+seqid+'\n')
        for start in range(0, length, 80):
            handle.write(sequence[start:start+80]+'\n')
    del sequence
    A.main(fasta, CHECKPOINT, out/'prediction')
    with (out/'summary.json').open('x') as handle:
        json.dump(dict(inference_complete=True, task=task, species=species, seqid=seqid, length=length,
            job_id=os.environ['SLURM_JOB_ID'], source_fit_job=fit_job, checkpoint=str(CHECKPOINT),
            prediction_summary='prediction/summary.json', prediction_gff3='prediction/full_grammar.gff3',
            input_fasta='input.fasta', raw_scope='whole_original_TRAIN_contig_both_strands',
            reference_used=False, labels_or_feature_cache_used=False, optimizer_steps_this_run=0,
            DEV_used=False, sealed_data_used=False, seconds=time.perf_counter()-began),
            handle, indent=2, allow_nan=False)
        handle.write('\n')


def evaluate(fit_job, baseline_fit_job):
    began = time.perf_counter()
    source(fit_job)
    B.source(baseline_fit_job)
    runs = reconcile_inputs()
    E.reconcile_inputs(B.OUT, B.CHECKPOINT, B.A.IDENTITY)
    for species, seqid, length in I.SCOPE:
        original = I.R.M.read_allowed(ROOT/'data/m1_screen'/species/'genome.fa', {seqid:length})[seqid]
        expected = re.sub('[^ACGT]', 'N', original.upper())
        for root in (OUT, B.OUT):
            if I.R.M.read_allowed(root/seqid/'input.fasta', {seqid:length})[seqid] != expected:
                raise ValueError('Both caller inputs must directly equal normalized original TRAIN DNA')
        del original, expected
    with (EVAL/'paired_input_match.json').open('x') as handle:
        json.dump(dict(job_id=os.environ['SLURM_JOB_ID'], protocol=A.PROTOCOL,
            codon_fit_job=fit_job, baseline_fit_job=baseline_fit_job,
            original_six_TRAIN_inputs_directly_equal=True, comparison_before_reference_access=True,
            DNA_comparison_policy='original_uppercase_nonACGT_to_N_preserving_length',
            DEV_used=False, sealed_data_used=False), handle, indent=2, allow_nan=False)
        handle.write('\n')
    # Same reference selection, roundtrip checks and count-based ruler as E.
    C = E.C
    report = dict(completed=True, status='original_TRAIN_end_to_end_diagnostic_not_DEV_or_generalization',
        job_id=os.environ['SLURM_JOB_ID'], checkpoint=str(CHECKPOINT),
        reference_policy='M26_CDS_assessable_longest_eligible_CDS_tie_transcript_ID',
        reference_used_in_inference=False, training_eligibility_does_not_filter_reference=True,
        unsupported_reference_strata_retained=True, DEV_used=False, sealed_data_used=False,
        optimizer_steps_this_run=0, chromosomes={}, per_species={})
    by_species = {}
    for species, seqid, length in I.SCOPE:
        folder = OUT/seqid
        sequences = C.M.read_allowed(folder/'input.fasta', {seqid:length})
        data = ROOT/'data/m1_screen'/species
        coding = C.S.E.parse_annotation(data/'reference.gff3', {seqid:length}, protein_coding_only=True)
        annotation = C.S.E.parse_annotation(data/'reference.gff3', {seqid:length})
        with (data/'reference.gff3').open() as handle:
            _, _, raw_cds, _ = C.P.raw_features(handle, seqid)
        refs, isoforms = C.select_reference(coding, raw_cds, sequences)
        if any(not C.R7.phase_chain_consistent(t) for t in refs):
            raise ValueError('Assessable primary reference phase inconsistent; do not silently filter')
        txs = C.S.E.primary_transcripts(C.S.E.parse_annotation(folder/'prediction/full_grammar.gff3', {seqid:length}), complete_only=False)
        if len(txs) != runs[seqid]['primary_CDS_chain_predictions']:
            raise ValueError('GFF3 roundtrip prediction count differs')
        for tx in txs:
            failures, _ = C.M.chain_checks(tx, sequences[seqid])
            if failures or len(tx['CDS']) > 128:
                raise ValueError('Exported codon prediction DNA/phase/cap differs')
        for strand in ('+', '-'):
            ordered = sorted((t['CDS'][0][0], t['CDS'][-1][1]) for t in txs if t['strand'] == strand)
            if any(b > c for (_, b), (c, _) in zip(ordered, ordered[1:])):
                raise ValueError('Exported selected same-strand genes overlap')
        overlaps, unsupported = C.M.C.overlapping_ids(refs), C.Counter()
        for tx in refs:
            failures, _ = C.M.chain_checks(tx, sequences[seqid])
            if tx['id'] in overlaps:
                failures.append('same_strand_primary_overlap')
            if len(tx['CDS']) > 128:
                failures.append('over_128_exons')
            if tx['CDS'][-1][1]-tx['CDS'][0][0] > C.M.WIDTH:
                failures.append('longer_than_native_window')
            unsupported.update(set(failures))
        row = C.score(txs, refs, isoforms, length, annotation)
        row['prediction_nonexact_primary'] = row.pop('reference_nonexact_primary')
        row.update(species=species, seqid=seqid, length=length, reference_primary_records=len(refs),
            reference_assessable_isoform_records=len(isoforms), retained_unsupported_reference_strata=dict(unsupported),
            raw_prediction_path=str((folder/'prediction/full_grammar.gff3').relative_to(ROOT)),
            phase_checked_prediction_records=len(txs), raw_inference_loop_seconds=runs[seqid]['wall_loop_seconds'],
            native_encoding_seconds=sum(s['encoding_seconds'] for s in runs[seqid]['scopes']),
            joint_decoding_seconds=sum(s['decode_seconds'] for s in runs[seqid]['scopes']))
        report['chromosomes'][seqid] = row
        by_species.setdefault(species, []).append(row)
        print(json.dumps(dict(completed_TRAIN_diagnostic=seqid, reference_primary_records=len(refs))), flush=True)
    report['per_species'] = {species:E.aggregate(rows) for species, rows in by_species.items()}
    report['pooled'] = E.aggregate(list(report['chromosomes'].values()))
    f1s = [r['exact_chain']['f1'] for r in report['per_species'].values()]
    report['macro_species_exact_chain_F1'] = sum(f1s)/len(f1s)
    report['worst_species_exact_chain_F1'] = min(f1s)
    report['evaluation_seconds'] = time.perf_counter()-began
    with (EVAL/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('raw', 'evaluate'))
    parser.add_argument('--task', type=int)
    parser.add_argument('--fit-job', required=True)
    parser.add_argument('--baseline-fit-job')
    args = parser.parse_args()
    torch.set_num_threads(2)
    if args.mode == 'raw':
        if args.task is None:
            parser.error('--task required for raw inference')
        raw(args.task, args.fit_job)
    elif args.baseline_fit_job is None:
        parser.error('--baseline-fit-job required for evaluate')
    else:
        evaluate(args.fit_job, args.baseline_fit_job)
