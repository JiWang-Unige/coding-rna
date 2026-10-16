#!/usr/bin/env python3
"""Prepared fixed DEV entry for all three final9216 arms, after TRAIN closure."""
import argparse
import json
import os
import time
import torch
import infer_full_grammar_fixed_dev as I
import annotate_full_grammar_objective_contrast as A

ROOT = I.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-DEV-RAW'
TRAIN_EVAL = ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-TRAIN-DIAGNOSTIC'
VERIFICATION = ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-VERIFICATION'
TRAIN_JOBS = dict(STRUCTURED_CONTINUE='13451625', LOCAL_CHANNEL='13451623', HALF_STRUCTURED_HALF_LOCAL='13451624')


def final_checkpoint(arm):
    if arm not in A.ARMS:
        raise ValueError('One of all three predeclared arms required')
    return ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST'/arm/'step_009216.pt'


def require_train_closure():
    # Read completion identities, not losses or TRAIN ranking. No DEV input here.
    for arm in A.ARMS:
        expected = final_checkpoint(arm)
        verify, evaluated = VERIFICATION/arm, TRAIN_EVAL/arm
        v = json.loads((verify/'summary.json').read_text())
        t = json.loads((evaluated/'summary.json').read_text())
        if ((verify/'STATUS').read_text().strip() != 'COMPLETED'
                or (evaluated/'STATUS').read_text().strip() != 'COMPLETED'
                or not v['completed'] or not t['completed'] or v['arm'] != arm
                or v['verified_cumulative_optimizer_steps'] != 9216
                or not v['final_checkpoint_strict_load']
                or ROOT/v['final_checkpoint'] != expected or t['checkpoint'] != str(expected)
                or t['job_id'] != TRAIN_JOBS[arm]
                or t['status'] != 'original_TRAIN_end_to_end_diagnostic_not_DEV_or_generalization'
                or v['DEV_used'] or t['DEV_used'] or v['sealed_data_used'] or t['sealed_data_used']):
            raise ValueError('All three fixed final TRAIN continuations must close before DEV')


def main(arm, task):
    if arm not in A.ARMS or task not in range(len(I.SCOPE)):
        raise ValueError('Declared arm and two fixed DEV tasks only')
    require_train_closure()
    species, seqid, length = I.SCOPE[task]
    out, checkpoint = OUT/arm/seqid, final_checkpoint(arm)
    began = time.perf_counter()
    with (out/'started.json').open('x') as handle:
        json.dump(dict(job_id=os.environ['SLURM_JOB_ID'], arm=arm, task=task,
            species=species, seqid=seqid, length=length, checkpoint=str(checkpoint),
            all_three_TRAIN_branches_closed=True, reference_used=False), handle, indent=2)
    sequence = I.R.M.read_allowed(ROOT/'data/m1_screen'/species/'genome.fa', {seqid:length})[seqid]
    fasta = out/'input.fasta'
    with fasta.open('x') as handle:
        handle.write('>'+seqid+'\n')
        for start in range(0, length, 80):
            handle.write(sequence[start:start+80].upper()+'\n')
    del sequence
    A.main(fasta, checkpoint, out/'prediction', arm)
    with (out/'summary.json').open('x') as handle:
        json.dump(dict(inference_complete=True, arm=arm, task=task, species=species, seqid=seqid, length=length,
            job_id=os.environ['SLURM_JOB_ID'], checkpoint=str(checkpoint), raw_scope=I.RAW_SCOPE,
            all_three_TRAIN_branches_closed=True, DEV_used=True, sealed_data_used=False,
            reference_used=False, labels_or_feature_cache_used=False, optimizer_steps_this_run=0,
            input_case_policy='uppercase_original_record_caller_non_ACGT_to_N_before_model',
            seconds=time.perf_counter()-began), handle, indent=2, allow_nan=False)
        handle.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arm', choices=A.ARMS, required=True)
    parser.add_argument('--task', type=int, required=True)
    args = parser.parse_args()
    torch.set_num_threads(2)
    main(args.arm, args.task)
