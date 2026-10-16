#!/usr/bin/env python3
"""One whole original TRAIN contig per task; no reference or label access."""
import argparse
import json
import os
from pathlib import Path
import re
import time
import torch
import annotate_full_grammar_coarse_channel as R

ROOT=R.M.C.ROOT
OUT=ROOT/'outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-RAW'
CHECKPOINT=ROOT/'outputs/M28-FULL-GRAMMAR-COARSE-CHANNEL-FIRST-FIT/step_004608.pt'
SCOPE=(('arabidopsis_thaliana','NC_003070.9',30427671),
       ('arabidopsis_thaliana','NC_003071.7',19698289),
       ('arabidopsis_thaliana','NC_003076.8',26975502),
       ('oryza_sativa','NC_089043.1',27474823),
       ('oryza_sativa','NC_089044.1',23931887),
       ('oryza_sativa','NC_089046.1',28271460))


def main(task):
    if not 0<=task<len(SCOPE): raise ValueError('Original six TRAIN tasks only')
    species,seqid,length=SCOPE[task]; out=OUT/seqid
    began=time.perf_counter()
    with (out/'started.json').open('x') as handle:
        json.dump(dict(job_id=os.environ['SLURM_JOB_ID'],task=task,species=species,
            seqid=seqid,length=length,checkpoint=str(CHECKPOINT),
            raw_scope='whole_original_TRAIN_contig_both_strands',
            DEV_used=False,reference_used=False,optimizer_steps_this_run=0),handle,indent=2)
    # Filter before retaining DNA: other records are discarded, not decoded.
    # Lengths/IDs are the original frozen TRAIN geometry, not inferred from labels.
    sequence=R.M.read_allowed(ROOT/'data/m1_screen'/species/'genome.fa',{seqid:length})[seqid]
    sequence=re.sub('[^ACGT]','N',sequence.upper())
    fasta=out/'input.fasta'
    with fasta.open('x') as handle:
        handle.write('>'+seqid+'\n')
        for a in range(0,length,80): handle.write(sequence[a:a+80]+'\n')
    del sequence
    R.main(fasta,CHECKPOINT,out/'prediction',0.125)
    report=dict(inference_complete=True,task=task,species=species,seqid=seqid,length=length,
        job_id=os.environ['SLURM_JOB_ID'],checkpoint=str(CHECKPOINT),
        prediction_summary='prediction/summary.json',prediction_gff3='prediction/full_grammar.gff3',
        input_fasta='input.fasta',raw_scope='whole_original_TRAIN_contig_both_strands',
        DEV_used=False,sealed_data_used=False,reference_used=False,
        labels_or_feature_cache_used=False,optimizer_steps_this_run=0,
        seconds=time.perf_counter()-began)
    with (out/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task',type=int,required=True)
    args=parser.parse_args(); torch.set_num_threads(2); main(args.task)
