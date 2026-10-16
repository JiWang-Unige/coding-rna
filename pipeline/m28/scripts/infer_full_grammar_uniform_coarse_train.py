#!/usr/bin/env python3
"""Uniform-coarse final9216, original-six raw TRAIN and paired-input ruler."""
import argparse
import json
import os
import re
import time
import torch
import infer_full_grammar_original_train as I
import infer_full_grammar_fresh_joint24_train as B
import annotate_full_grammar_uniform_coarse as A

ROOT = I.ROOT
FIT = ROOT/'outputs/M28-FULL-GRAMMAR-UNIFORM-COARSE-FIT'
VERIFY = ROOT/'outputs/M28-FULL-GRAMMAR-UNIFORM-COARSE-VERIFICATION'
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-UNIFORM-COARSE-TRAIN-RAW'
EVAL = ROOT/'outputs/M28-FULL-GRAMMAR-UNIFORM-COARSE-TRAIN-DIAGNOSTIC'
CHECKPOINT = FIT/'step_009216.pt'


def source(fit_job):
    fitted = json.loads((FIT/'summary.json').read_text())
    verified = json.loads((VERIFY/'summary.json').read_text())
    if ((FIT/'STATUS').read_text().strip() != 'COMPLETED' or not fitted['completed']
            or fitted['job_id'] != fit_job or fitted['cumulative_optimizer_steps'] != 9216
            or fitted['optimizer_steps_this_run'] != 4608 or fitted['fresh_draws'] != 4608
            or fitted['fresh_unique_windows'] != 4219 or fitted['DEV_used'] or fitted['sealed_data_used']
            or any(fitted['configuration'].get(k) != v for k,v in A.IDENTITY.items())
            or fitted['final_checkpoint'] != str(CHECKPOINT.relative_to(ROOT))):
        raise ValueError('Completed fixed uniform-coarse fresh/joint final9216 fit required')
    if ((VERIFY/'STATUS').read_text().strip() != 'COMPLETED' or not verified['completed']
            or verified['source_fit_job'] != fit_job or verified['configuration'] != fitted['configuration']
            or not verified['fixed_fresh_draw_order_observations_and_importance_match']
            or not verified['logged_objective_arithmetic_match'] or not verified['final_checkpoint_strict_load']
            or verified['final_checkpoint'] != fitted['final_checkpoint']
            or verified['DEV_used'] or verified['sealed_data_used']):
        raise ValueError('Completed uniform-coarse exposure and checkpoint verification required')


def raw(task, fit_job):
    if not 0 <= task < 6:
        raise ValueError('One original-six TRAIN chromosome required')
    source(fit_job)
    began = time.perf_counter()
    species,seqid,length = I.SCOPE[task]
    out = OUT/seqid
    with (out/'started.json').open('x') as handle:
        json.dump(dict(job_id=os.environ['SLURM_JOB_ID'],task=task,species=species,
            seqid=seqid,length=length,checkpoint=str(CHECKPOINT),source_fit_job=fit_job,
            reference_used=False,DEV_used=False),handle,indent=2)
    sequence = I.R.M.read_allowed(ROOT/'data/m1_screen'/species/'genome.fa',{seqid:length})[seqid]
    sequence = re.sub('[^ACGT]','N',sequence.upper())
    fasta = out/'input.fasta'
    with fasta.open('x') as handle:
        handle.write('>'+seqid+'\n')
        for start in range(0,length,80):
            handle.write(sequence[start:start+80]+'\n')
    del sequence
    A.main(fasta,CHECKPOINT,out/'prediction')
    with (out/'summary.json').open('x') as handle:
        json.dump(dict(inference_complete=True,task=task,species=species,seqid=seqid,length=length,
            job_id=os.environ['SLURM_JOB_ID'],source_fit_job=fit_job,checkpoint=str(CHECKPOINT),
            prediction_summary='prediction/summary.json',prediction_gff3='prediction/full_grammar.gff3',
            input_fasta='input.fasta',raw_scope='whole_original_TRAIN_contig_both_strands',
            reference_used=False,labels_or_feature_cache_used=False,optimizer_steps_this_run=0,
            DEV_used=False,sealed_data_used=False,seconds=time.perf_counter()-began),
            handle,indent=2,allow_nan=False)
        handle.write('\n')


def evaluate(fit_job, informative_fit_job):
    source(fit_job)
    B.source(informative_fit_job)
    import evaluate_full_grammar_original_train as E
    E.reconcile_inputs(OUT,CHECKPOINT,A.IDENTITY)
    E.reconcile_inputs(B.OUT,B.CHECKPOINT,B.A.IDENTITY)
    for species,seqid,length in I.SCOPE:
        original = I.R.M.read_allowed(ROOT/'data/m1_screen'/species/'genome.fa',{seqid:length})[seqid]
        expected = re.sub('[^ACGT]','N',original.upper())
        for root in (OUT,B.OUT):
            if I.R.M.read_allowed(root/seqid/'input.fasta',{seqid:length})[seqid] != expected:
                raise ValueError('Both caller inputs must directly equal normalized original TRAIN DNA')
        del original,expected
    with (EVAL/'paired_input_match.json').open('x') as handle:
        json.dump(dict(job_id=os.environ['SLURM_JOB_ID'],protocol=A.PROTOCOL,
            uniform_fit_job=fit_job,informative_fit_job=informative_fit_job,
            original_six_TRAIN_inputs_directly_equal=True,comparison_before_reference_access=True,
            DNA_comparison_policy='original_uppercase_nonACGT_to_N_preserving_length',
            DEV_used=False,sealed_data_used=False),handle,indent=2,allow_nan=False)
        handle.write('\n')
    E.main(OUT,CHECKPOINT,A.IDENTITY,EVAL)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('raw','evaluate'))
    parser.add_argument('--task',type=int)
    parser.add_argument('--fit-job',required=True)
    parser.add_argument('--informative-fit-job')
    args = parser.parse_args()
    torch.set_num_threads(2)
    if args.mode == 'raw':
        if args.task is None:
            parser.error('--task required for raw inference')
        raw(args.task,args.fit_job)
    elif args.informative_fit_job is None:
        parser.error('--informative-fit-job required for evaluate')
    else:
        evaluate(args.fit_job,args.informative_fit_job)
