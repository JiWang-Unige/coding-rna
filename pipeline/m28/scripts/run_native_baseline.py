#!/usr/bin/env python3
"""Six fixed raw-FASTA native runs; no reference, training or selection."""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
PREP = ROOT/'outputs/M28-NATIVE-RUNTIME208'
ASSET_PREP = ROOT/'outputs/M28-NATIVE-BASELINE-PREP-R2'
OUT = ROOT/'outputs/M28-NATIVE-BASELINES'
METHODS = ('ANNEVO','Tiberius','Helixer')
SCOPE = (('arabidopsis_thaliana','NC_003074.8',23459830),
         ('oryza_sativa','NC_089041.1',29936421))
FINISHED = {'ANNEVO':'The gene annotation took ',
            'Tiberius':'Tiberius took ',
            'Helixer':'Helixer successfully finished the annotation'}

def task_scope(task):
    if task not in range(6):
        raise ValueError('Task must identify one of the six fixed native runs')
    return METHODS[task//2], SCOPE[task%2]

def commands(method, species, run, visible, cpus=8, *, fasta=None):
    assets = ASSET_PREP/'assets'
    supplied_fasta = fasta
    fasta = Path(fasta) if fasta is not None else ASSET_PREP/'inputs'/(species+'.fa')
    if method == 'ANNEVO':
        native = [sys.executable,str(assets/'ANNEVO/annotation.py'),
                  '-g',str(fasta),'-m',str(assets/'ANNEVO/saved_model/ANNEVO_Magnoliopsida.pt'),
                  '-l','Magnoliopsida','-o',str(run/'prediction.gff3'),
                  '--batch_size','32','-t',str(cpus),'--num_workers',str(cpus),
                  '--tmp_path',str(run/'temporary'),'--show_log']
        probe = [sys.executable,'-c',
                 "import json,torch; assert torch.cuda.is_available(); "
                 "assert torch.cuda.device_count()==1; "
                 "print(json.dumps({'backend':'torch','version':torch.__version__,"
                 "'cuda':torch.version.cuda,'gpu':torch.cuda.get_device_name(0),'visible_gpus':1}))"]
        return native, probe
    image = assets/('Tiberius-2.0.8.sif' if method=='Tiberius' else 'Helixer-0.3.7.sif')
    prefix = ['apptainer','exec','--nv','--cleanenv',
              '--env','CUDA_VISIBLE_DEVICES='+visible,
              '--env','OMP_NUM_THREADS='+str(cpus),
              '--env','TF_NUM_INTRAOP_THREADS='+str(cpus),
              '--env','TF_NUM_INTEROP_THREADS=1',
              '--bind',str(ASSET_PREP)+':/m28-prep:ro',
              '--bind',str(run)+':/m28-run',
              '--pwd','/m28-run' if method=='Tiberius' else '/home/helixer_user/Helixer']
    container_fasta = '/m28-prep/inputs/'+species+'.fa'
    if supplied_fasta is not None:
        prefix += ['--bind',str(fasta)+':/m28-input.fa:ro']
        container_fasta = '/m28-input.fa'
    if method == 'Tiberius':
        prefix += ['--bind',str(PREP)+':/m28-runtime:ro',
                   '--env','PYTHONPATH=/m28-runtime/site',
                   '--bind',str(assets/'Tiberius-model-weights')+':/m28-runtime/site/model_weights:ro']
        native = prefix+[str(image),'python3','-m','tiberius.main',
                         '--genome',container_fasta,
                         '--model_cfg','/m28-prep/assets/angiosperms.yaml',
                         '--out','/m28-run/prediction.gff3',
                         '--batch_size','4','--seq_len','200070','--no_softmasking']
    elif method == 'Helixer':
        native = prefix+[str(image),'python3','/home/helixer_user/Helixer/Helixer.py',
                         '--fasta-path',container_fasta,
                         '--gff-output-path','/m28-run/prediction.gff3',
                         '--species',species,
                         '--model-filepath','/m28-prep/assets/land_plant_v0.3_a_0080.h5',
                         '--subsequence-length','64152','--batch-size','16',
                         '--temporary-dir','/m28-run/temporary','--deterministic',
                         '--no-multiprocess']
    else:
        raise ValueError('Unknown native method')
    identity = ("import importlib.metadata as m,tiberius; "
                "assert m.version('tiberius')=='2.0.8'; "
                "assert tiberius.__file__.startswith('/m28-runtime/site/tiberius/'); "
                if method=='Tiberius' else '')
    probe = prefix+[str(image),'python3','-c',
                   identity+"import json,tensorflow as tf; g=tf.config.list_physical_devices('GPU'); "
                   "assert len(g)==1,g; "
                   "print(json.dumps({'backend':'tensorflow','version':tf.__version__,"
                   "'gpu':tf.config.experimental.get_device_details(g[0]),"
                   "'build':tf.sysconfig.get_build_info(),'visible_gpus':len(g)}))"]
    return native, probe

def check_output(path, logfile, method, seqid, length):
    # Helixer 0.3.7 does not propagate its postprocessor exit code to main.
    # A completion marker plus raw output checks distinguish failure from zero predictions.
    with logfile.open() as handle:
        finished = any(FINISHED[method] in line for line in handle)
    if not finished:
        raise RuntimeError(method+' native completion marker missing; inspect command.log')
    counts = Counter()
    with path.open() as handle:
        for lineno,line in enumerate(handle,1):
            if not line.strip() or line.startswith('#'):
                continue
            row = line.rstrip('\n').split('\t')
            if len(row)!=9 or row[0]!=seqid:
                raise ValueError('Unexpected GFF scope or fields at line '+str(lineno))
            start,end = int(row[3]),int(row[4])
            if not 1<=start<=end<=length:
                raise ValueError('GFF coordinates outside fixed FASTA at line '+str(lineno))
            if row[2]=='CDS' and (row[6] not in ('+','-') or row[7] not in ('0','1','2') or
                                 not any(a.startswith('Parent=') and a[7:] for a in row[8].split(';'))):
                raise ValueError('CDS strand/phase/parent invalid at line '+str(lineno))
            counts[row[2]] += 1
    return {'bytes':path.stat().st_size,'feature_counts':dict(counts),
            'native_completion_marker':FINISHED[method],
            'empty_prediction_is_not_filtered_or_rescued':not counts.get('CDS',0)}

def run(task):
    began = time.perf_counter()
    method,(species,seqid,length) = task_scope(task)
    prep = json.loads((PREP/'summary.json').read_text())
    if (PREP/'STATUS').read_text().strip()!='COMPLETED' or not all(
            prep[k] for k in ('completed','inputs_complete','runtime_assets_complete')):
        raise ValueError('Native assets are not complete')
    if {(r['species'],r['seqid'],r['length']) for r in prep['inputs']}!=set(SCOPE):
        raise ValueError('Prepared FASTA scope differs')
    tib = prep['baselines']['Tiberius']
    if (prep['reused_asset_root']!=str(ASSET_PREP.relative_to(ROOT)) or
        tib['native_metadata']['version']!='2.0.8' or
        tib['native_metadata']['source_commit']!='871d4979568ae09e56de7eca5d7489046ee833ea' or
        tib['base_container_actual_metadata']['version']!='2.0.7' or
        not tib['native_metadata']['CPU_model_load_completed']):
        raise ValueError('Corrected isolated Tiberius runtime identity differs')
    if any(prep[k] for k in ('reference_used','test_setaria_used','native_inference_executed')):
        raise ValueError('Preparation boundary differs')
    visible = os.environ['CUDA_VISIBLE_DEVICES']
    if not visible or visible in ('-1','NoDevFiles'):
        raise ValueError('Slurm GPU assignment missing')
    run_dir = OUT/method/species
    run_dir.parent.mkdir(parents=True,exist_ok=True)
    run_dir.mkdir()
    status = run_dir/'STATUS'
    status.write_text('RUNNING\n')
    try:
        (run_dir/'JOBID').write_text(os.environ['SLURM_JOB_ID']+'\n')
        (run_dir/'temporary').mkdir()
        native,probe = commands(method,species,run_dir,visible)
        with (run_dir/'commands.json').open('x') as handle:
            json.dump({'native':native,'GPU_probe':probe,'CUDA_VISIBLE_DEVICES':visible},handle,indent=2)
        with (run_dir/'gpu_probe.stderr').open('x') as err:
            result = subprocess.run(probe,check=True,text=True,stdout=subprocess.PIPE,stderr=err)
        (run_dir/'gpu_probe.stdout').write_text(result.stdout)
        gpu = json.loads(result.stdout.strip().splitlines()[-1])
        with (run_dir/'command.log').open('x') as handle:
            subprocess.run(native,check=True,stdout=handle,stderr=subprocess.STDOUT)
        gff = check_output(run_dir/'prediction.gff3',run_dir/'command.log',method,seqid,length)
        summary = {'completed':True,'method':method,'species':species,'seqid':seqid,'length':length,
                   'raw_GFF':'prediction.gff3','output_check':gff,'GPU_runtime':gpu,
                   'baseline_asset_receipt':prep['baselines'][method],
                   'input_mode':'fixed_DEV_raw_FASTA','case_policy':'uppercase_same_as_M28_raw_caller',
                   'reference_used':False,'labels_or_feature_cache_used':False,
                   'optimizer_steps':0,'test_setaria_used':False,
                   'native_inference_executed':True,'seconds':time.perf_counter()-began,
                   'limitations':'native output on two DEV contigs; accuracy not yet scored; not independent generalization'}
        with (run_dir/'summary.json').open('x') as handle:
            json.dump(summary,handle,indent=2,allow_nan=False); handle.write('\n')
        status.write_text('COMPLETED\n')
        print(json.dumps({'completed':True,'method':method,'species':species,'seconds':summary['seconds'],
                          'GFF_feature_counts':gff['feature_counts']}),flush=True)
    except BaseException:
        status.write_text('FAILED\n')
        raise

if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--task',type=int,required=True)
    run(p.parse_args().task)
