#!/usr/bin/env python3
"""Prepared three-final-weight DEV comparison; no TRAIN-score selection."""
import argparse
import json
import os
import subprocess
import time
import torch
import infer_full_grammar_fixed_dev as I
import infer_full_grammar_fresh_joint24_train as B
import infer_full_grammar_uniform_coarse_train as U
import infer_full_grammar_codon_train as D
import evaluate_full_grammar_dev_baselines as E

ROOT = I.ROOT
PLAN = 'reports/M28-METHOD-RESTART/full_grammar_joint_controls_dev_plan.md'
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-JOINT-CONTROLS-DEV-RAW-R1'
COMPARISON = ROOT/'outputs/M28-FULL-GRAMMAR-JOINT-CONTROLS-DEV-COMPARISON-R1'
BRANCHES = ('FRESH_JOINT24576', 'UNIFORM_COARSE', 'CODON_ENERGY')
METHODS = tuple('M28_'+branch+'_final9216' for branch in BRANCHES)
# Actual already-submitted TRAIN jobs, not guessed future DEV jobs.
SPECS = {
    'FRESH_JOINT24576': (B, '13494895', '13494896', '13494898', '13494899',
        ROOT/'outputs/M28-FULL-GRAMMAR-FRESH-JOINT24-FOUR-CELL'),
    'UNIFORM_COARSE': (U, '13496011', '13496012', '13496014', '13496015',
        ROOT/'outputs/M28-FULL-GRAMMAR-UNIFORM-COARSE-PAIRED'),
    'CODON_ENERGY': (D, '13497914', '13497943', '13497946', '13497947',
        ROOT/'outputs/M28-FULL-GRAMMAR-CODON-PAIRED'),
}


def candidate(branch):
    if branch not in BRANCHES:
        raise ValueError('One of the three fixed final9216 branches required')
    return SPECS[branch][0], METHODS[BRANCHES.index(branch)]


def completed(folder, job):
    value = json.loads((folder/'summary.json').read_text())
    if ((folder/'STATUS').read_text().strip() != 'COMPLETED'
            or (folder/'JOBID').read_text().strip() != job or not value['completed']
            or value['job_id'] != job or value['DEV_used'] or value['sealed_data_used']):
        raise ValueError('Fixed complete TRAIN identity required: '+job)
    return value


def require_train_closure():
    jobs = set()
    for branch in BRANCHES:
        train, fit, verify, evaluation, paired, pair_root = SPECS[branch]
        train.source(fit)
        completed(train.FIT, fit)
        completed(train.VERIFY, verify)
        report = completed(train.EVAL, evaluation)
        if (report['checkpoint'] != str(train.CHECKPOINT)
                or report['status'] != 'original_TRAIN_end_to_end_diagnostic_not_DEV_or_generalization'
                or set(report['chromosomes']) != {r[1] for r in B.I.SCOPE}
                or report['optimizer_steps_this_run'] != 0):
            raise ValueError('All three six-chromosome TRAIN evaluations must close before DEV')
        pair = completed(pair_root, paired)
        if branch == BRANCHES[0]:
            valid = (pair['all_four_caller_inputs_and_normalized_original_DNA_directly_equal']
                and pair['counts_and_species_pooled_PRF_reconciled']
                and pair['source_jobs'] == dict(F00_original_split='13451625', F10_fresh_split='13453583',
                    F01_original_joint='13453780', F11_fresh_joint=evaluation))
        else:
            arm = 'uniform' if branch == 'UNIFORM_COARSE' else 'codon'
            control = 'informative' if arm == 'uniform' else 'baseline'
            valid = (pair['all_paired_caller_inputs_and_normalized_original_DNA_directly_equal']
                and pair['counts_and_PRF_reconciled']
                and pair['DNA_comparison_before_'+arm+'_reference_access']
                and pair['source_jobs'] == {control:'13494898', arm:evaluation}
                and pair['source_fit_jobs'] == {control:'13494895', arm:fit})
        if not valid or pair['optimizer_steps_this_run'] != 0:
            raise ValueError('Complete fixed TRAIN comparison required before DEV')
        jobs.update((fit, verify, evaluation, paired))
        for task, (species, seqid, length) in enumerate(B.I.SCOPE):
            folder = train.OUT/seqid
            row = json.loads((folder/'summary.json').read_text())
            job = (folder/'JOBID').read_text().strip()
            if ((folder/'STATUS').read_text().strip() != 'COMPLETED'
                    or not row['inference_complete'] or row['job_id'] != job
                    or (row['task'],row['species'],row['seqid'],row['length']) != (task,species,seqid,length)
                    or row['source_fit_job'] != fit or row['checkpoint'] != str(train.CHECKPOINT)
                    or row['DEV_used'] or row['sealed_data_used'] or row['reference_used']
                    or row['labels_or_feature_cache_used'] or row['optimizer_steps_this_run'] != 0):
                raise ValueError('Every fixed TRAIN raw task must close before DEV')
            jobs.add(job)
    if len(jobs) != 30 or any(not job.isdigit() for job in jobs):
        raise ValueError('Thirty distinct actual numeric TRAIN Slurm roots required')
    result = subprocess.run(['sacct','-j',','.join(sorted(jobs)),'-n','-P','-o','JobIDRaw,State,ExitCode'],
        check=True, text=True, stdout=subprocess.PIPE)
    E.check_terminal_rows([line.split('|') for line in result.stdout.splitlines() if line.strip()], jobs)


def main(branch, task):
    train, _ = candidate(branch)
    if task not in range(len(I.SCOPE)):
        raise ValueError('Only the two original fixed DEV chromosomes are supported')
    require_train_closure()
    species, seqid, length = I.SCOPE[task]
    out = OUT/branch/seqid
    began = time.perf_counter()
    receipt = dict(branch=branch, task=task, species=species, seqid=seqid, length=length,
        job_id=os.environ['SLURM_JOB_ID'], checkpoint=str(train.CHECKPOINT),
        source_fit_job=SPECS[branch][1], plan=PLAN, raw_scope=I.RAW_SCOPE,
        all_three_TRAIN_flows_closed=True, DEV_used=True, sealed_data_used=False,
        reference_used=False, labels_or_feature_cache_used=False, optimizer_steps_this_run=0,
        input_case_policy='uppercase_original_record_caller_non_ACGT_to_N_before_model')
    with (out/'started.json').open('x') as handle:
        json.dump(receipt, handle, indent=2)
    sequence = I.R.M.read_allowed(ROOT/'data/m1_screen'/species/'genome.fa', {seqid:length})[seqid]
    fasta = out/'input.fasta'
    with fasta.open('x') as handle:
        handle.write('>'+seqid+'\n')
        for start in range(0,length,80):
            handle.write(sequence[start:start+80].upper()+'\n')
    del sequence
    train.A.main(fasta, train.CHECKPOINT, out/'prediction')
    with (out/'summary.json').open('x') as handle:
        json.dump(dict(receipt, inference_complete=True, seconds=time.perf_counter()-began),
            handle, indent=2, allow_nan=False)
        handle.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--branch', choices=BRANCHES, required=True)
    parser.add_argument('--task', type=int, required=True)
    args = parser.parse_args()
    torch.set_num_threads(2)
    main(args.branch, args.task)
