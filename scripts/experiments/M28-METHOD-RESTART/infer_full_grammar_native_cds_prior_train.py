#!/usr/bin/env python3
"""Original six TRAIN chromosomes, strict final native-prior raw inference."""
import argparse
import json
import os
import re
import time
import torch
import infer_full_grammar_original_train as I
import annotate_full_grammar_native_cds_prior as A

ROOT = I.ROOT
FIT = ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CDS-PRIOR-FIT'
VERIFY = ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CDS-PRIOR-VERIFICATION'
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CDS-PRIOR-TRAIN-RAW'
CHECKPOINT = FIT/'step_009216.pt'


def main(task, fit_job):
    if not 0 <= task < 6: raise ValueError('Original six TRAIN tasks required')
    fitted = json.loads((FIT/'summary.json').read_text()); verified = json.loads((VERIFY/'summary.json').read_text())
    if ((FIT/'STATUS').read_text().strip() != 'COMPLETED' or not fitted['completed'] or fitted['job_id'] != fit_job
            or fitted['cumulative_optimizer_steps'] != 9216 or fitted['optimizer_steps_this_run'] != 4608
            or fitted['DEV_used'] or fitted['sealed_data_used']
            or fitted['final_checkpoint'] != str(CHECKPOINT.relative_to(ROOT))
            or (VERIFY/'STATUS').read_text().strip() != 'COMPLETED' or not verified['completed']
            or verified['source_fit_job'] != fit_job or not verified['final_checkpoint_strict_load']
            or verified['configuration'] != fitted['configuration'] or verified['final_checkpoint'] != fitted['final_checkpoint']):
        raise ValueError('Completed final fit and full exposure/checkpoint reconciliation required')
    species, seqid, length = I.SCOPE[task]; out = OUT/seqid; began = time.perf_counter()
    with (out/'started.json').open('x') as handle:
        json.dump(dict(job_id=os.environ['SLURM_JOB_ID'], task=task, species=species, seqid=seqid,
            length=length, checkpoint=str(CHECKPOINT), source_fit_job=fit_job,
            reference_used=False, DEV_used=False), handle, indent=2)
    sequence = I.R.M.read_allowed(ROOT/'data/m1_screen'/species/'genome.fa', {seqid:length})[seqid]
    sequence = re.sub('[^ACGT]', 'N', sequence.upper()); fasta = out/'input.fasta'
    with fasta.open('x') as handle:
        handle.write('>'+seqid+'\n')
        for start in range(0, length, 80): handle.write(sequence[start:start+80]+'\n')
    del sequence
    A.main(fasta, CHECKPOINT, out/'prediction')
    with (out/'summary.json').open('x') as handle:
        json.dump(dict(inference_complete=True, task=task, species=species, seqid=seqid, length=length,
            job_id=os.environ['SLURM_JOB_ID'], source_fit_job=fit_job, checkpoint=str(CHECKPOINT),
            prediction_summary='prediction/summary.json', prediction_gff3='prediction/full_grammar.gff3',
            input_fasta='input.fasta', raw_scope='whole_original_TRAIN_contig_both_strands',
            reference_used=False, labels_or_feature_cache_used=False, optimizer_steps_this_run=0,
            DEV_used=False, sealed_data_used=False, seconds=time.perf_counter()-began), handle, indent=2, allow_nan=False)
        handle.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=int, required=True); parser.add_argument('--fit-job', required=True)
    args = parser.parse_args(); torch.set_num_threads(2); main(args.task, args.fit_job)
