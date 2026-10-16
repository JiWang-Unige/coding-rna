#!/usr/bin/env python3
"""Correct the 2.0.7 OCI base with an isolated, pinned 2.0.8 package only."""
import json
import os
from pathlib import Path
import subprocess
import time
import yaml
import prepare_native_baselines as P

ROOT = P.ROOT
BASE = ROOT/'outputs/M28-NATIVE-BASELINE-PREP-R2'
OUT = ROOT/'outputs/M28-NATIVE-RUNTIME208'
SOURCE = P.TIBERIUS_CONFIG_COMMIT
BASE_SOURCE = 'a1472d26829769c22c9678c9f042830bbb75e544'


def dump(path, row):
    with path.open('x') as handle:
        json.dump(row, handle, indent=2, allow_nan=False)
        handle.write('\n')


def container(image, weights=False):
    command = ['apptainer','exec','--cleanenv',
               '--env','CUDA_VISIBLE_DEVICES=',
               '--env','PYTHONPATH=/m28-runtime/site',
               '--env','PYTHONDONTWRITEBYTECODE=1',
               '--env','OMP_NUM_THREADS=2',
               '--env','TF_NUM_INTRAOP_THREADS=2',
               '--env','TF_NUM_INTEROP_THREADS=1',
               '--bind',str(BASE)+':/m28-prep:ro',
               '--bind',str(OUT)+':/m28-runtime',
               '--pwd','/m28-runtime']
    if weights:
        command += ['--bind',str(BASE/'assets/Tiberius-model-weights')+
                    ':/m28-runtime/site/model_weights:ro']
    return command+[str(image)]


def main():
    began = time.perf_counter()
    inputs = json.loads((BASE/'inputs.json').read_text())
    if (not inputs['complete'] or inputs['reference_used'] or inputs['test_setaria_used'] or
        {(r['species'],r['seqid'],r['length']) for r in inputs['inputs']} != set(P.SCOPE)):
        raise ValueError('Previously staged input boundary differs')
    assets = BASE/'assets'; receipts = OUT/'receipts'; receipts.mkdir()
    image = assets/'Tiberius-2.0.8.sif'  # Historical filename, NOT its actual version.
    base_identity = P.container_metadata(image,'tiberius','/opt/Tiberius')
    if base_identity['version'] != '2.0.7' or base_identity['source_commit'] != BASE_SOURCE:
        raise ValueError('Corrected runtime requires the recorded immutable 2.0.7 base')
    repo = OUT/'Tiberius-source'
    subprocess.run(['git','init',str(repo)],check=True)
    subprocess.run(['git','-C',str(repo),'remote','add','origin',
                    'https://github.com/Gaius-Augustus/Tiberius.git'],check=True)
    subprocess.run(['git','-C',str(repo),'fetch','--depth','1','origin',SOURCE],check=True)
    subprocess.run(['git','-C',str(repo),'checkout','--detach','FETCH_HEAD'],check=True)
    if P.capture(['git','-C',str(repo),'rev-parse','HEAD']) != SOURCE:
        raise ValueError('Tiberius source pin differs')
    site = OUT/'site'; site.mkdir()
    # Build isolation may install hatchling in a temporary build environment only.
    # --no-deps leaves the base TensorFlow/CUDA/bricks2marble/hidten stack untouched.
    install = container(image)+['python3','-m','pip','install','--no-deps',
                               '--no-cache-dir','--no-compile',
                               '--target','/m28-runtime/site',
                               '/m28-runtime/Tiberius-source']
    dump(receipts/'install-command.json',{'command':install,'dependency_upgrade':False})
    with (OUT/'install.log').open('x') as log:
        subprocess.run(install,check=True,stdout=log,stderr=subprocess.STDOUT)
    (site/'model_weights').mkdir()
    source_python = sorted((repo/'tiberius').rglob('*.py'))
    for path in source_python:
        installed = site/'tiberius'/path.relative_to(repo/'tiberius')
        if not installed.is_file() or installed.read_bytes() != path.read_bytes():
            raise ValueError('Installed package differs from pinned source: '+str(path))
    cfg = yaml.safe_load((assets/'angiosperms.yaml').read_text())
    if (cfg['weights_url']!=P.WEIGHTS_URL or cfg['softmasking'] is not False or
        cfg['clamsa'] is not False or cfg['default_seq_len']!=200070):
        raise ValueError('Fixed Tiberius model config differs')
    archive = assets/'Tiberius-model-weights/angiosperms.tar.gz'
    if archive.stat().st_size != 29825066:
        raise ValueError('Previously downloaded Tiberius weight archive incomplete')
    # No DNA is opened and no forward or optimizer is run. Loading the unchanged
    # model catches incompatible package/dependency/weight APIs before native jobs.
    code = """
import importlib.metadata as m, json, platform
import tiberius
import tiberius.main as native
import tensorflow as tf
from tiberius.eval_model_class import PredictionGTF
assert m.version('tiberius') == '2.0.8'
assert tiberius.__version__ == '2.0.8'
assert tiberius.__file__.startswith('/m28-runtime/site/tiberius/')
assert native.__file__.startswith('/m28-runtime/site/tiberius/')
assert str(m.distribution('tiberius').locate_file('')).rstrip('/') == '/m28-runtime/site'
assert not tf.config.list_physical_devices('GPU')
cfg = native.load_model_config('/m28-prep/assets/angiosperms.yaml')
weights = native.resolve_weight_download(cfg)
assert weights == '/m28-runtime/site/tiberius/../model_weights/angiosperms', weights
prediction = PredictionGTF(model_path=weights, seq_len=200070, batch_size=4,
    hmm=True, hmm_emitter_epsilon=0.01, hmm_initial_exon_len=200,
    hmm_initial_intron_len=4500, hmm_initial_ir_len=10000, temp_dir=None,
    num_hmm=1, hmm_factor=1, genome=None, softmask=False,
    parallel_factor=native.compute_parallel_factor(200070))
prediction.load_model(summary=0)
print(json.dumps({'package':'tiberius','version':m.version('tiberius'),
    'package_version':tiberius.__version__, 'package_file':tiberius.__file__,
    'main_file':native.__file__, 'python':platform.python_version(),
    'packages':{n:m.version(n) for n in ('tensorflow','bricks2marble','hidten','numpy','pyyaml','rich')},
    'tensorflow':tf.__version__, 'tensorflow_build':tf.sysconfig.get_build_info(),
    'weights_directory':weights, 'CPU_model_load_completed':True,
    'GPU_forward_executed':False, 'DNA_opened':False, 'optimizer_steps':0}))
"""
    probe = container(image,weights=True)+['python3','-c',code]
    dump(receipts/'model-load-command.json',{'command':probe})
    with (OUT/'model-load.stderr').open('x') as err:
        result = subprocess.run(probe,check=True,text=True,stdout=subprocess.PIPE,stderr=err)
    (OUT/'model-load.stdout').write_text(result.stdout)
    runtime = json.loads(result.stdout.strip().splitlines()[-1])
    runtime['source_commit'] = SOURCE
    runtime['installed_source_python_files_exact'] = len(source_python)
    cli = container(image,weights=True)+['python3','-m','tiberius.main','--help']
    with (OUT/'cli-help.txt').open('x') as handle:
        subprocess.run(cli,check=True,stdout=handle,stderr=subprocess.STDOUT)
    annevo = json.loads((BASE/'receipts/ANNEVO.json').read_text())
    if (annevo['source_commit']!=P.ANNEVO_COMMIT or
        (assets/annevo['weight'].removeprefix('assets/')).stat().st_size!=annevo['weight_bytes']):
        raise ValueError('Previously staged ANNEVO source/weight differs')
    hel = P.container_metadata(assets/'Helixer-0.3.7.sif','helixer','/home/helixer_user/Helixer')
    if hel['version']!='0.3.7' or hel['source_commit']!=P.HELIXER_COMMIT:
        raise ValueError('Previously staged Helixer identity differs')
    hel['HelixerPost_source_commit'] = P.capture(
        ['apptainer','exec','--cleanenv',str(assets/'Helixer-0.3.7.sif'),
         'git','-C','/home/helixer_user/HelixerPost','rev-parse','HEAD'])
    selected = P.best_plant_model((assets/'Helixer-model-list.csv').read_text())
    if selected['model_file_name']!=P.HELIXER_WEIGHT or int(selected['priority'])!=80:
        raise ValueError('Fixed Helixer official priority model differs')
    tiberius = {'native_metadata':runtime,'configuration_commit':SOURCE,
        'configuration':cfg,'base_container_uri':P.TIBERIUS_IMAGE,
        'base_container_path':str(image.relative_to(ROOT)),
        'base_container_actual_metadata':base_identity,
        'runtime_kind':'pinned_2.0.8_isolated_package_over_immutable_2.0.7_dependency_SIF',
        'entry_point':'python3 -m tiberius.main',
        'package_directory':str(site.relative_to(ROOT)),
        'weights_directory':str((assets/'Tiberius-model-weights/angiosperms').relative_to(ROOT)),
        'weights_archive_bytes':archive.stat().st_size,
        'weights_requested_url':P.WEIGHTS_URL,
        'download_HTTP_receipt':'not_persisted_before_PREP_R2_failure; not reconstructed',
        'training_exposure':'author_metadata_includes_rice_Kitaake_not_species_unseen'}
    helixer = {'container_uri':P.HELIXER_IMAGE,'native_metadata':hel,
        'container_path':str((assets/'Helixer-0.3.7.sif').relative_to(ROOT)),
        'model_source_commit':P.HELIXER_COMMIT,'selected_model':selected,
        'model_selection':'lowest_official_land_plant_priority_not_DEV_score',
        'weight':str((assets/P.HELIXER_WEIGHT).relative_to(ROOT)),
        'weight_bytes':(assets/P.HELIXER_WEIGHT).stat().st_size,
        'download_HTTP_receipt':'not_persisted_before_PREP_R2_failure; not reconstructed'}
    if helixer['weight_bytes']<=0:
        raise ValueError('Existing Helixer weight empty')
    for name,row in (('ANNEVO',annevo),('Tiberius',tiberius),('Helixer',helixer)):
        dump(receipts/(name+'.json'),row)
    report = {'experiment':'M28-NATIVE-RUNTIME208','completed':True,
        'inputs_complete':True,'runtime_assets_complete':True,'inputs':inputs['inputs'],
        'reused_asset_root':str(BASE.relative_to(ROOT)),
        'original_failed_preparation_preserved':True,
        'baselines':{'ANNEVO':annevo,'Tiberius':tiberius,'Helixer':helixer},
        'native_inference_executed':False,'GPU_forward_executed':False,
        'reference_used':False,'optimizer_steps':0,'test_setaria_used':False,
        'old_assets_modified':False,'seconds':time.perf_counter()-began,
        'limitations':'runtime/CPU model load only; native GPU predictions and accuracy unverified'}
    dump(OUT/'summary.json',report)
    print(json.dumps({'completed':True,'version':runtime['version'],
        'CPU_model_load_completed':True,'seconds':report['seconds']}),flush=True)


if __name__=='__main__': main()
