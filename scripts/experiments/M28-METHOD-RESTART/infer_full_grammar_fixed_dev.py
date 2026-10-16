#!/usr/bin/env python3
"""Prepared fixed-DEV entry; call only after DEV execution is authorized."""
import argparse
import json
import os
import time
import torch
import annotate_full_grammar_coarse_channel as R

ROOT = R.M.C.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-FIXED-DEV-RAW'
CHECKPOINT = ROOT/'outputs/M28-FULL-GRAMMAR-COARSE-CHANNEL-FIRST-FIT/step_004608.pt'
SCOPE = (('arabidopsis_thaliana','NC_003074.8',23459830),
         ('oryza_sativa','NC_089041.1',29936421))
RAW_SCOPE = 'whole_fixed_DEV_contig_both_strands'


def main(task):
    if task not in range(len(SCOPE)):
        raise ValueError('Only the two fixed DEV tasks are supported')
    species,seqid,length = SCOPE[task]
    out = OUT/seqid
    began = time.perf_counter()
    with (out/'started.json').open('x') as handle:
        json.dump(dict(job_id=os.environ['SLURM_JOB_ID'],task=task,species=species,
            seqid=seqid,length=length,checkpoint=str(CHECKPOINT),raw_scope=RAW_SCOPE,
            DEV_used=True,sealed_data_used=False,reference_used=False,
            optimizer_steps_this_run=0),handle,indent=2)
    sequence = R.M.read_allowed(ROOT/'data/m1_screen'/species/'genome.fa',{seqid:length})[seqid]
    # Preserve the uppercase original record for direct comparator-input equality.
    # The unchanged production caller normalizes non-ACGT to N before the model.
    fasta = out/'input.fasta'
    with fasta.open('x') as handle:
        handle.write('>'+seqid+'\n')
        for a in range(0,length,80):
            handle.write(sequence[a:a+80].upper()+'\n')
    del sequence
    R.main(fasta,CHECKPOINT,out/'prediction',.125)
    report = dict(inference_complete=True,task=task,species=species,seqid=seqid,length=length,
        job_id=os.environ['SLURM_JOB_ID'],checkpoint=str(CHECKPOINT),
        prediction_summary='prediction/summary.json',prediction_gff3='prediction/full_grammar.gff3',
        input_fasta='input.fasta',raw_scope=RAW_SCOPE,DEV_used=True,sealed_data_used=False,
        reference_used=False,labels_or_feature_cache_used=False,optimizer_steps_this_run=0,
        input_case_policy='uppercase_original_record_caller_non_ACGT_to_N_before_model',
        seconds=time.perf_counter()-began)
    with (out/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps(report),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task',type=int,required=True)
    args = parser.parse_args(); torch.set_num_threads(2); main(args.task)
