#!/usr/bin/env python3
"""Prepared fixed final native-prior DEV caller, only after all TRAIN closure."""
import argparse
import json
import os
import time
import torch
import infer_full_grammar_fixed_dev as I
import infer_full_grammar_objective_contrast_dev as C
import infer_full_grammar_native_cds_prior_train as T
import annotate_full_grammar_native_cds_prior as A

ROOT = I.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CDS-PRIOR-DEV-RAW'
VERIFY = T.VERIFY
TRAIN_EVAL = ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CDS-PRIOR-TRAIN-DIAGNOSTIC'
CHECKPOINT = T.CHECKPOINT
METHOD = 'M28_NATIVE_CDS_PRIOR_final9216'
PLAN = 'reports/M28-METHOD-RESTART/full_grammar_native_cds_prior_dev_plan.md'


def require_train_closure():
    C.require_train_closure()  # All original objectives, not a TRAIN-score winner.
    v = json.loads((VERIFY/'summary.json').read_text())
    t = json.loads((TRAIN_EVAL/'summary.json').read_text())
    if ((VERIFY/'STATUS').read_text().strip() != 'COMPLETED'
            or (TRAIN_EVAL/'STATUS').read_text().strip() != 'COMPLETED'
            or not v['completed'] or not t['completed']
            or v['job_id'] != '13453180' or v['source_fit_job'] != '13453177'
            or t['job_id'] != '13453184'
            or v['verified_old_C0_optimizer_steps'] != 9216 or v['verified_new_prior_optimizer_steps'] != 4608
            or not v['final_checkpoint_strict_load'] or ROOT/v['final_checkpoint'] != CHECKPOINT
            or any(v['configuration'].get(k) != value for k, value in A.IDENTITY.items())
            or t['checkpoint'] != str(CHECKPOINT)
            or t['status'] != 'original_TRAIN_end_to_end_diagnostic_not_DEV_or_generalization'
            or any(row[k] for row in (v, t) for k in ('DEV_used', 'sealed_data_used'))):
        raise ValueError('Native-prior final TRAIN verification and full raw diagnostic must close before DEV')


def main(task):
    if task not in range(len(I.SCOPE)): raise ValueError('Two original fixed DEV tasks only')
    require_train_closure()
    species, seqid, length = I.SCOPE[task]; out = OUT/seqid; began = time.perf_counter()
    with (out/'started.json').open('x') as handle:
        json.dump(dict(job_id=os.environ['SLURM_JOB_ID'], task=task, species=species, seqid=seqid,
            length=length, checkpoint=str(CHECKPOINT), plan=PLAN, all_four_TRAIN_branches_closed=True,
            reference_used=False, DEV_used=True, sealed_data_used=False), handle, indent=2)
    sequence = I.R.M.read_allowed(ROOT/'data/m1_screen'/species/'genome.fa', {seqid:length})[seqid]
    fasta = out/'input.fasta'
    with fasta.open('x') as handle:
        handle.write('>'+seqid+'\n')
        for start in range(0, length, 80): handle.write(sequence[start:start+80].upper()+'\n')
    del sequence
    A.main(fasta, CHECKPOINT, out/'prediction')
    with (out/'summary.json').open('x') as handle:
        json.dump(dict(inference_complete=True, task=task, species=species, seqid=seqid, length=length,
            job_id=os.environ['SLURM_JOB_ID'], checkpoint=str(CHECKPOINT), raw_scope=I.RAW_SCOPE,
            all_four_TRAIN_branches_closed=True, DEV_used=True, sealed_data_used=False,
            reference_used=False, labels_or_feature_cache_used=False, optimizer_steps_this_run=0,
            input_case_policy='uppercase_original_record_caller_non_ACGT_to_N_before_model',
            seconds=time.perf_counter()-began), handle, indent=2, allow_nan=False)
        handle.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--task', type=int, required=True)
    args = parser.parse_args(); torch.set_num_threads(2); main(args.task)
