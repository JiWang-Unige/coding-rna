#!/usr/bin/env python3
"""Original six TRAIN chromosomes only, one predeclared arm per array."""
import argparse
import json
import os
import re
import time
import torch
import infer_full_grammar_original_train as I
import annotate_full_grammar_objective_contrast as A

ROOT = I.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-TRAIN-RAW'


def main(arm, task, fit_job):
    if arm not in A.ARMS or not 0 <= task < 6:
        raise ValueError('Declared arm and original TRAIN chromosome required')
    fit = ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST'/arm
    checkpoint = fit/'step_009216.pt'
    fitted = json.loads((fit/'summary.json').read_text())
    if ((fit/'STATUS').read_text().strip() != 'COMPLETED' or fitted['job_id'] != fit_job
            or not fitted['completed'] or fitted['configuration'].get('arm') != arm
            or fitted['cumulative_optimizer_steps'] != 9216 or fitted['optimizer_steps_this_run'] != 4608
            or fitted['DEV_used'] or fitted['sealed_data_used']
            or fitted['final_checkpoint'] != str(checkpoint.relative_to(ROOT))):
        raise ValueError('Completed declared final arm fit required')
    verification = ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-VERIFICATION'/arm
    verified = json.loads((verification/'summary.json').read_text())
    if ((verification/'STATUS').read_text().strip() != 'COMPLETED' or not verified['completed']
            or verified['source_fit_job'] != fit_job or not verified['final_checkpoint_strict_load']
            or verified['configuration'] != fitted['configuration']
            or verified['final_checkpoint'] != fitted['final_checkpoint']):
        raise ValueError('Complete exposure/arithmetic/checkpoint reconciliation required')
    species, seqid, length = I.SCOPE[task]
    out = OUT/arm/seqid
    began = time.perf_counter()
    with (out/'started.json').open('x') as handle:
        json.dump(dict(job_id=os.environ['SLURM_JOB_ID'], arm=arm, task=task,
            species=species, seqid=seqid, length=length, checkpoint=str(checkpoint),
            source_fit_job=fit_job, reference_used=False, DEV_used=False), handle, indent=2)
    sequence = I.R.M.read_allowed(ROOT/'data/m1_screen'/species/'genome.fa', {seqid:length})[seqid]
    sequence = re.sub('[^ACGT]', 'N', sequence.upper())
    fasta = out/'input.fasta'
    with fasta.open('x') as handle:
        handle.write('>'+seqid+'\n')
        for start in range(0, length, 80):
            handle.write(sequence[start:start+80]+'\n')
    del sequence
    A.main(fasta, checkpoint, out/'prediction', arm)
    with (out/'summary.json').open('x') as handle:
        json.dump(dict(inference_complete=True, task=task, arm=arm, species=species, seqid=seqid,
            length=length, job_id=os.environ['SLURM_JOB_ID'], source_fit_job=fit_job,
            checkpoint=str(checkpoint), prediction_summary='prediction/summary.json',
            prediction_gff3='prediction/full_grammar.gff3', input_fasta='input.fasta',
            raw_scope='whole_original_TRAIN_contig_both_strands', reference_used=False,
            labels_or_feature_cache_used=False, optimizer_steps_this_run=0,
            DEV_used=False, sealed_data_used=False, seconds=time.perf_counter()-began),
            handle, indent=2, allow_nan=False)
        handle.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arm', choices=A.ARMS, required=True)
    parser.add_argument('--task', type=int, required=True)
    parser.add_argument('--fit-job', required=True)
    args = parser.parse_args()
    torch.set_num_threads(2)
    main(args.arm, args.task, args.fit_job)
