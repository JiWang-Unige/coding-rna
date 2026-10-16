#!/usr/bin/env python3
"""Joint24576 strict raw original-six TRAIN annotation and unchanged evaluation."""
import argparse
import json
import os
import re
import time
import torch
import infer_full_grammar_original_train as I
import annotate_full_grammar_joint24 as A

ROOT = I.ROOT
FIT = ROOT/'outputs/M28-FULL-GRAMMAR-JOINT24-FIT'
VERIFY = ROOT/'outputs/M28-FULL-GRAMMAR-JOINT24-VERIFICATION'
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-JOINT24-TRAIN-RAW'
EVAL = ROOT/'outputs/M28-FULL-GRAMMAR-JOINT24-TRAIN-DIAGNOSTIC'
CHECKPOINT = FIT/'step_009216.pt'


def source(fit_job):
    fitted = json.loads((FIT/'summary.json').read_text())
    verified = json.loads((VERIFY/'summary.json').read_text())
    if ((FIT/'STATUS').read_text().strip() != 'COMPLETED' or not fitted['completed']
            or fitted['job_id'] != fit_job or fitted['cumulative_optimizer_steps'] != 9216
            or fitted['optimizer_steps_this_run'] != 4608 or fitted['original_draws_per_pass'] != 1536
            or fitted['original_unique_windows'] != 1491 or fitted['DEV_used'] or fitted['sealed_data_used']
            or any(fitted['configuration'].get(k) != v for k,v in A.IDENTITY.items())
            or fitted['final_checkpoint'] != str(CHECKPOINT.relative_to(ROOT))):
        raise ValueError('Completed fixed joint24576 final9216 fit required')
    if ((VERIFY/'STATUS').read_text().strip() != 'COMPLETED' or not verified['completed']
            or verified['source_fit_job'] != fit_job or verified['configuration'] != fitted['configuration']
            or not verified['original_draw_order_observations_and_importance_match']
            or not verified['logged_objective_arithmetic_match'] or not verified['final_checkpoint_strict_load']
            or verified['final_checkpoint'] != fitted['final_checkpoint']
            or verified['DEV_used'] or verified['sealed_data_used']):
        raise ValueError('Completed original-exposure and checkpoint verification required')


def raw(task,fit_job):
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('raw','evaluate'))
    parser.add_argument('--task',type=int)
    parser.add_argument('--fit-job',required=True)
    args = parser.parse_args()
    torch.set_num_threads(2)
    if args.mode == 'raw':
        if args.task is None:
            parser.error('--task required for raw inference')
        raw(args.task,args.fit_job)
    else:
        source(args.fit_job)
        import evaluate_full_grammar_original_train as E
        E.main(OUT,CHECKPOINT,A.IDENTITY,EVAL)
