#!/usr/bin/env python3
"""Two prospectively fixed DEV candidates, only after all six TRAIN closures."""
import argparse
import json
import os
import time
import torch
import infer_full_grammar_fixed_dev as I
import infer_full_grammar_native_cds_prior_dev as N
import infer_full_grammar_fresh_draws_train as F
import infer_full_grammar_joint24_train as J

ROOT = I.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-SAMPLING-CONTEXT-DEV-RAW'
COMPARISON = ROOT/'outputs/M28-FULL-GRAMMAR-SAMPLING-CONTEXT-DEV-COMPARISON'
PLAN = 'reports/M28-METHOD-RESTART/full_grammar_sampling_context_dev_plan.md'
BRANCHES = ('FRESH_DRAWS','JOINT24576')
METHODS = ('M28_FRESH_DRAWS_final9216','M28_JOINT24576_final9216')


def candidate(branch):
    if branch == 'FRESH_DRAWS':
        return F,METHODS[0]
    if branch == 'JOINT24576':
        return J,METHODS[1]
    raise ValueError('One of the two prospectively fixed candidates required')


def require_train_closure():
    N.require_train_closure()  # Original three objectives and dual-CDS-prior, without score selection.
    for branch,fit_job,eval_job in (('FRESH_DRAWS','13453573','13453583'),
                                  ('JOINT24576','13453777','13453780')):
        train,_ = candidate(branch)
        train.source(fit_job)  # Existing full fit/exposure/checkpoint verification identity.
        evaluated = json.loads((train.EVAL/'summary.json').read_text())
        if ((train.EVAL/'STATUS').read_text().strip() != 'COMPLETED' or not evaluated['completed']
                or evaluated['job_id'] != eval_job or evaluated['checkpoint'] != str(train.CHECKPOINT)
                or evaluated['status'] != 'original_TRAIN_end_to_end_diagnostic_not_DEV_or_generalization'
                or set(evaluated['chromosomes']) != {r[1] for r in train.I.SCOPE}
                or evaluated['DEV_used'] or evaluated['sealed_data_used']
                or evaluated['optimizer_steps_this_run'] != 0):
            raise ValueError('Both new full original-TRAIN raw diagnostics must close before DEV')


def main(branch,task):
    train,_ = candidate(branch)
    if task not in range(len(I.SCOPE)):
        raise ValueError('Two original fixed DEV tasks only')
    require_train_closure()
    species,seqid,length = I.SCOPE[task]
    out = OUT/branch/seqid
    began = time.perf_counter()
    with (out/'started.json').open('x') as handle:
        json.dump(dict(job_id=os.environ['SLURM_JOB_ID'],branch=branch,task=task,
            species=species,seqid=seqid,length=length,checkpoint=str(train.CHECKPOINT),plan=PLAN,
            all_six_TRAIN_branches_closed=True,reference_used=False,DEV_used=True,sealed_data_used=False),
            handle,indent=2)
    sequence = I.R.M.read_allowed(ROOT/'data/m1_screen'/species/'genome.fa',{seqid:length})[seqid]
    fasta = out/'input.fasta'
    with fasta.open('x') as handle:
        handle.write('>'+seqid+'\n')
        for start in range(0,length,80):
            handle.write(sequence[start:start+80].upper()+'\n')
    del sequence
    train.A.main(fasta,train.CHECKPOINT,out/'prediction')
    with (out/'summary.json').open('x') as handle:
        json.dump(dict(inference_complete=True,branch=branch,task=task,species=species,seqid=seqid,
            length=length,job_id=os.environ['SLURM_JOB_ID'],checkpoint=str(train.CHECKPOINT),plan=PLAN,
            raw_scope=I.RAW_SCOPE,all_six_TRAIN_branches_closed=True,DEV_used=True,sealed_data_used=False,
            reference_used=False,labels_or_feature_cache_used=False,optimizer_steps_this_run=0,
            input_case_policy='uppercase_original_record_caller_non_ACGT_to_N_before_model',
            seconds=time.perf_counter()-began),handle,indent=2,allow_nan=False)
        handle.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--branch',choices=BRANCHES,required=True)
    parser.add_argument('--task',type=int,required=True)
    args = parser.parse_args()
    torch.set_num_threads(2)
    main(args.branch,args.task)
