#!/usr/bin/env python3
"""Pinned native comparators on the six existing, shared original TRAIN FASTAs."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time
import run_native_baseline as N

ROOT = N.ROOT
OUT = ROOT/'outputs/M28-NATIVE-TRAIN-BASELINES'
CHECKS = ROOT/'outputs/M28-NATIVE-TRAIN-BASELINE-CHECKS-R2'
ORIGINAL_TRAIN = ROOT/'outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-RAW'
METHODS = N.METHODS
SCOPE = (('arabidopsis_thaliana','NC_003070.9',30427671),
         ('arabidopsis_thaliana','NC_003071.7',19698289),
         ('arabidopsis_thaliana','NC_003076.8',26975502),
         ('oryza_sativa','NC_089043.1',27474823),
         ('oryza_sativa','NC_089044.1',23931887),
         ('oryza_sativa','NC_089046.1',28271460))
INPUT_MODE = 'fixed_original_TRAIN_shared_raw_FASTA'


def task_scope(task):
    if not 0 <= task < len(METHODS)*len(SCOPE):
        raise ValueError('Only the eighteen fixed native TRAIN tasks are allowed')
    return METHODS[task//len(SCOPE)], SCOPE[task%len(SCOPE)]


def input_fasta(seqid):
    return ORIGINAL_TRAIN/seqid/'input.fasta'


def runtime_receipt():
    prep = json.loads((N.PREP/'summary.json').read_text())
    if ((N.PREP/'STATUS').read_text().strip() != 'COMPLETED'
            or not all(prep[k] for k in ('completed','inputs_complete','runtime_assets_complete'))
            or prep['reused_asset_root'] != str(N.ASSET_PREP.relative_to(ROOT))
            or any(prep[k] for k in ('reference_used','test_setaria_used','native_inference_executed'))):
        raise ValueError('Existing isolated native runtime is incomplete or differs')
    baselines = prep['baselines']
    if (baselines['ANNEVO']['source_commit'] != '42c920f70f52f96c7c06c665a639ec47d9304fff'
            or baselines['Helixer']['native_metadata']['source_commit'] != 'd17bb496b8842542590644b2437798711c0125d5'
            or baselines['Helixer']['selected_model']['model_file_name'] != 'land_plant_v0.3_a_0080.h5'):
        raise ValueError('Fixed native source or weight selection differs')
    tib = baselines['Tiberius']
    if (tib['native_metadata']['version'] != '2.0.8'
            or tib['native_metadata']['source_commit'] != '871d4979568ae09e56de7eca5d7489046ee833ea'
            or tib['base_container_actual_metadata']['version'] != '2.0.7'
            or not tib['native_metadata']['CPU_model_load_completed']):
        raise ValueError('Isolated Tiberius runtime identity differs')
    return prep


def reconcile_train_input(index):
    species, seqid, length = SCOPE[index]
    folder = ORIGINAL_TRAIN/seqid
    row = json.loads((folder/'summary.json').read_text())
    if ((folder/'STATUS').read_text().strip() != 'COMPLETED' or not row['inference_complete']
            or (row['task'],row['species'],row['seqid'],row['length']) != (index,species,seqid,length)
            or row['input_fasta'] != 'input.fasta' or row['optimizer_steps_this_run'] != 0
            or any(row[k] for k in ('DEV_used','sealed_data_used','reference_used','labels_or_feature_cache_used'))
            or not input_fasta(seqid).is_file()):
        raise ValueError('Completed original TRAIN input identity differs')
    return row


def check_fasta(path, seqid, length):
    bp = n_bp = 0
    with path.open() as handle:
        if handle.readline().rstrip() != '>'+seqid:
            raise ValueError('Shared TRAIN FASTA record differs')
        for line in handle:
            sequence = line.strip()
            if set(sequence)-set('ACGTN'):
                raise ValueError('Shared TRAIN FASTA is not original uppercase ACGTN')
            bp += len(sequence)
            n_bp += sequence.count('N')
    if bp != length:
        raise ValueError('Shared TRAIN FASTA length differs')
    return {'seqid':seqid,'length':bp,'N_bases':n_bp,
            'path':str(path.relative_to(ROOT)),
            'case_policy':'reuse_original_TRAIN_uppercase_non_ACGT_to_N'}


def run(task):
    began = time.perf_counter()
    method, (species,seqid,length) = task_scope(task)
    checked = json.loads((CHECKS/'summary.json').read_text())
    if ((CHECKS/'STATUS').read_text().strip() != 'COMPLETED' or not checked['completed']
            or checked['scope'] != [list(s) for s in SCOPE] or checked['input_mode'] != INPUT_MODE):
        raise ValueError('Native TRAIN preflight has not passed')
    prep = runtime_receipt()
    source = reconcile_train_input(task%len(SCOPE))
    visible = os.environ['CUDA_VISIBLE_DEVICES']
    if not visible or visible in ('-1','NoDevFiles'):
        raise ValueError('Slurm GPU assignment missing')
    folder = OUT/method/seqid
    folder.parent.mkdir(parents=True,exist_ok=True)
    folder.mkdir()
    status = folder/'STATUS'
    status.write_text('RUNNING\n')
    try:
        (folder/'JOBID').write_text(os.environ['SLURM_JOB_ID']+'\n')
        (folder/'temporary').mkdir()
        native, probe = N.commands(method,species,folder,visible,fasta=input_fasta(seqid))
        with (folder/'commands.json').open('x') as handle:
            json.dump({'native':native,'GPU_probe':probe,'CUDA_VISIBLE_DEVICES':visible},handle,indent=2)
        with (folder/'gpu_probe.stderr').open('x') as err:
            result = subprocess.run(probe,check=True,text=True,stdout=subprocess.PIPE,stderr=err)
        (folder/'gpu_probe.stdout').write_text(result.stdout)
        gpu = json.loads(result.stdout.strip().splitlines()[-1])
        with (folder/'command.log').open('x') as handle:
            subprocess.run(native,check=True,stdout=handle,stderr=subprocess.STDOUT)
        gff = N.check_output(folder/'prediction.gff3',folder/'command.log',method,seqid,length)
        summary = dict(completed=True,task=task,method=method,species=species,seqid=seqid,length=length,
            job_id=os.environ['SLURM_JOB_ID'],raw_GFF='prediction.gff3',output_check=gff,GPU_runtime=gpu,
            baseline_asset_receipt=prep['baselines'][method],input_mode=INPUT_MODE,
            input_fasta=str(input_fasta(seqid)),original_input_receipt_job_id=source['job_id'],
            reference_used=False,labels_or_feature_cache_used=False,optimizer_steps=0,
            DEV_used=False,sealed_data_used=False,test_setaria_used=False,native_inference_executed=True,
            seconds=time.perf_counter()-began,
            limitations='original TRAIN diagnostic, not independent generalization or a fair speed benchmark')
        with (folder/'summary.json').open('x') as handle:
            json.dump(summary,handle,indent=2,allow_nan=False); handle.write('\n')
        status.write_text('COMPLETED\n')
        print(json.dumps(dict(completed=True,method=method,seqid=seqid,seconds=summary['seconds'])),flush=True)
    except BaseException:
        status.write_text('FAILED\n')
        raise


def completed_runs():
    runs = {}
    prep = runtime_receipt()
    for task in range(len(METHODS)*len(SCOPE)):
        method, (species,seqid,length) = task_scope(task)
        folder = OUT/method/seqid
        if (folder/'STATUS').read_text().strip() != 'COMPLETED':
            raise ValueError('All eighteen native TRAIN runs must complete before reference access')
        row = json.loads((folder/'summary.json').read_text())
        if (not row['completed'] or (row['task'],row['method'],row['species'],row['seqid'],row['length'])
                != (task,method,species,seqid,length)
                or row['input_mode'] != INPUT_MODE or row['input_fasta'] != str(input_fasta(seqid))
                or row['optimizer_steps'] != 0 or not row['native_inference_executed']
                or any(row[k] for k in ('reference_used','labels_or_feature_cache_used','DEV_used','sealed_data_used','test_setaria_used'))
                or row['baseline_asset_receipt'] != prep['baselines'][method]
                or row['raw_GFF'] != 'prediction.gff3' or not (folder/'prediction.gff3').is_file()):
            raise ValueError('Native TRAIN completed-run identity differs')
        runs[(method,seqid)] = row
    return runs


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task',type=int,required=True)
    run(parser.parse_args().task)
