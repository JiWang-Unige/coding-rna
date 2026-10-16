#!/usr/bin/env python3
"""Stage pinned native comparator runtimes and only the two existing DEV sequences."""
import argparse
import csv
import importlib.metadata
import io
import json
import os
import platform
import shutil
import subprocess
import time
import urllib.request
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT/'outputs/M28-NATIVE-BASELINE-PREP'
SCOPE = (('arabidopsis_thaliana','NC_003074.8',23459830),
         ('oryza_sativa','NC_089041.1',29936421))
ANNEVO_COMMIT = '42c920f70f52f96c7c06c665a639ec47d9304fff'
TIBERIUS_CONFIG_COMMIT = '871d4979568ae09e56de7eca5d7489046ee833ea'
HELIXER_COMMIT = 'd17bb496b8842542590644b2437798711c0125d5'
TIBERIUS_IMAGE = 'docker://larsgabriel23/tiberius@sha256:6c853ae0664d6739d4226c3c0bc08d39d89aefa88d1bf0d20941a44d7ba2b72c'
HELIXER_IMAGE = 'docker://gglyptodon/helixer-docker@sha256:e2294eb2c282c35b919933daa0d6c145b635bfac8b717dbff88e88accbde4303'
WEIGHTS_URL = 'https://bioinf.uni-greifswald.de/bioinf/tiberius/models/angiosperms.tar.gz'
HELIXER_WEIGHT = 'land_plant_v0.3_a_0080.h5'


def extract_fasta(source, target, expected_length, destination):
    """Stream past excluded records without materializing them; uppercase like M28."""
    found = False
    chunks = []
    with source.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if found:
                    break
                found = line[1:].split()[0] == target
            elif found:
                chunks.append(line.strip())
    if not found:
        raise ValueError('Missing permitted DEV seqid '+target)
    original = ''.join(chunks)
    if len(original) != expected_length:
        raise ValueError('Permitted DEV sequence length mismatch')
    sequence = original.upper()
    with destination.open('x') as handle:
        handle.write('>'+target+'\n')
        for offset in range(0,len(sequence),80):
            handle.write(sequence[offset:offset+80]+'\n')
    return {'seqid':target,'length':len(sequence),'records':1,
            'lowercase_source_bases':sum(c.islower() for c in original),
            'case_policy':'uppercase_same_as_M28_raw_caller',
            'path':str(destination.relative_to(OUT))}


def download(url, destination):
    print(json.dumps({'download_start':url,'destination':str(destination)}),flush=True)
    with urllib.request.urlopen(url,timeout=60) as response, destination.open('xb') as handle:
        receipt = {'requested_url':url,'effective_url':response.geturl(),
                   'ETag':response.headers.get('ETag'),
                   'Last-Modified':response.headers.get('Last-Modified'),
                   'Content-Length':response.headers.get('Content-Length')}
        shutil.copyfileobj(response,handle,length=1024*1024)
    receipt['bytes'] = destination.stat().st_size
    if receipt['bytes'] == 0 or (receipt['Content-Length'] is not None and
                               receipt['bytes'] != int(receipt['Content-Length'])):
        raise ValueError('Incomplete downloaded asset '+str(destination))
    return receipt


def best_plant_model(text):
    lines = [line for line in text.splitlines() if line and not line.startswith('#')]
    models = list(csv.DictReader(lines))
    return min((m for m in models if m['lineage']=='land_plant'),
               key=lambda m:int(m['priority']))


def capture(command, cwd=None):
    return subprocess.check_output(command,cwd=cwd,text=True).strip()


def container_metadata(image, package, repository):
    code = ("import importlib.metadata as m,json,platform,subprocess;"
            "from pathlib import Path;"
            "p=Path("+repr(repository)+");"
            "revision=subprocess.check_output(['git','-C',str(p),'rev-parse','HEAD'],text=True).strip() "
            "if (p/'.git').exists() else None;"
            "print(json.dumps({'package':"+repr(package)+",'version':m.version("+repr(package)+"),"
            "'python':platform.python_version(),'source_commit':revision}))")
    return json.loads(capture(['apptainer','exec','--cleanenv',str(image),'python3','-c',code]))


def main():
    began = time.perf_counter()
    inputs = OUT/'inputs'; assets = OUT/'assets'; receipts = OUT/'receipts'
    for path in (inputs,assets,receipts):
        path.mkdir()
    input_rows = []
    for species,seqid,length in SCOPE:
        row = extract_fasta(ROOT/'data/m1_screen'/species/'genome.fa',
                            seqid,length,inputs/(species+'.fa'))
        row['species'] = species
        input_rows.append(row)
    with (OUT/'inputs.json').open('x') as handle:
        json.dump({'complete':True,'inputs':input_rows,'reference_used':False,
                   'test_setaria_used':False,'genomic_scope':'two_existing_DEV_contigs_only'},handle,indent=2)
        handle.write('\n')
    print(json.dumps({'input_complete':True,'inputs':input_rows}),flush=True)

    # New checkout only: never update the shared ANNEVO source used by active C0.
    repo = assets/'ANNEVO'
    subprocess.run(['git','init',str(repo)],check=True)
    subprocess.run(['git','-C',str(repo),'remote','add','origin','https://github.com/xjtu-omics/ANNEVO.git'],check=True)
    subprocess.run(['git','-C',str(repo),'fetch','--depth','1','origin',ANNEVO_COMMIT],check=True)
    subprocess.run(['git','-C',str(repo),'checkout','--detach','FETCH_HEAD'],check=True)
    if capture(['git','-C',str(repo),'rev-parse','HEAD']) != ANNEVO_COMMIT:
        raise ValueError('ANNEVO source pin differs')
    weight = repo/'saved_model/ANNEVO_Magnoliopsida.pt'
    if weight.stat().st_size != 60520287:
        raise ValueError('ANNEVO pinned weight missing or incomplete')
    annevo = {'source_commit':ANNEVO_COMMIT,'version':'2.3.3',
              'weight':str(weight.relative_to(OUT)),
              'weight_git_blob':capture(['git','-C',str(repo),'rev-parse','HEAD:saved_model/ANNEVO_Magnoliopsida.pt']),
              'weight_bytes':weight.stat().st_size,'lineage':'Magnoliopsida',
              'existing_conda_env':'annevo','python':platform.python_version(),
              'packages':{name:importlib.metadata.version(name) for name in
                          ('torch','numpy','h5py','biopython','torchmetrics','numba')}}
    with (receipts/'ANNEVO.json').open('x') as handle:
        json.dump(annevo,handle,indent=2); handle.write('\n')
    print(json.dumps({'ANNEVO_staged':True,'source_commit':ANNEVO_COMMIT}),flush=True)

    tib_cfg = assets/'angiosperms.yaml'
    tib_cfg_receipt = download(
        'https://raw.githubusercontent.com/Gaius-Augustus/Tiberius/'+TIBERIUS_CONFIG_COMMIT+'/model_cfg/angiosperms.yaml',
        tib_cfg)
    cfg = yaml.safe_load(tib_cfg.read_text())
    if (cfg['weights_url'] != WEIGHTS_URL or cfg['softmasking'] is not False or
        cfg['clamsa'] is not False or cfg['default_seq_len'] != 200070):
        raise ValueError('Fixed sequence-only angiosperm model config differs')
    weights = assets/'Tiberius-model-weights'; weights.mkdir()
    archive = weights/'angiosperms.tar.gz'
    tib_weight_receipt = download(WEIGHTS_URL,archive)
    subprocess.run(['tar','-xzf',str(archive),'-C',str(weights)],check=True)
    if not (weights/'angiosperms').is_dir() or not any((weights/'angiosperms').iterdir()):
        raise ValueError('Tiberius archive did not provide expected model directory')
    print(json.dumps({'Tiberius_weights_staged':True,'archive_bytes':archive.stat().st_size}),flush=True)

    model_list = assets/'Helixer-model-list.csv'
    helixer_list_receipt = download(
        'https://raw.githubusercontent.com/usadellab/Helixer/'+HELIXER_COMMIT+'/resources/model_list.csv',
        model_list)
    selected = best_plant_model(model_list.read_text())
    if selected['model_file_name'] != HELIXER_WEIGHT or int(selected['priority']) != 80:
        raise ValueError('Pinned official default land-plant selector changed')
    helixer_weight = assets/HELIXER_WEIGHT
    helixer_weight_receipt = download(selected['download_link'],helixer_weight)

    tib_image,hel_image = assets/'Tiberius-2.0.8.sif',assets/'Helixer-0.3.7.sif'
    cache = assets/'apptainer-cache'; cache.mkdir()
    temporary = assets/'apptainer-tmp'; temporary.mkdir()
    env = dict(os.environ,APPTAINER_CACHEDIR=str(cache),APPTAINER_TMPDIR=str(temporary))
    for image,uri in ((tib_image,TIBERIUS_IMAGE),(hel_image,HELIXER_IMAGE)):
        # Original pull exceeded its 8 GiB allocation at mksquashfs. Bound the
        # compressor directly; the pinned OCI contents and runtime are unchanged.
        command = ['apptainer','build','--disable-cache','--mksquashfs-args',
                   '-processors 2 -mem 2G',str(image),uri]
        print(json.dumps({'container_build_start':uri,'command':command}),flush=True)
        subprocess.run(command,check=True,env=env)
        labels = json.loads(capture(['apptainer','inspect','--json',str(image)]))
        with (receipts/(image.stem+'-labels.json')).open('x') as handle:
            json.dump(labels,handle,indent=2); handle.write('\n')
    tib_native = container_metadata(tib_image,'tiberius','/opt/Tiberius')
    hel_native = container_metadata(hel_image,'helixer','/home/helixer_user/Helixer')
    if tib_native['version'] != '2.0.8' or hel_native['version'] != '0.3.7':
        raise ValueError('Actual container package version differs from frozen comparator')
    hel_native['HelixerPost_source_commit'] = capture(
        ['apptainer','exec','--cleanenv',str(hel_image),'git','-C',
         '/home/helixer_user/HelixerPost','rev-parse','HEAD'])
    if hel_native['source_commit'] != HELIXER_COMMIT:
        raise ValueError('Helixer actual source differs from v0.3.7')
    tib_container_cfg = capture(
        ['apptainer','exec','--cleanenv',str(tib_image),'cat','/opt/Tiberius/model_cfg/angiosperms.yaml'])
    with (receipts/'Tiberius-container-angiosperms.yaml').open('x') as handle:
        handle.write(tib_container_cfg+'\n')
    tiberius = {'container_uri':TIBERIUS_IMAGE,'container_path':str(tib_image.relative_to(OUT)),
                'native_metadata':tib_native,'configuration_commit':TIBERIUS_CONFIG_COMMIT,
                'explicit_operating_config':str(tib_cfg.relative_to(OUT)),
                'configuration_receipt':tib_cfg_receipt,'weights_receipt':tib_weight_receipt,
                'weights_directory':str((weights/'angiosperms').relative_to(OUT)),
                'configuration':cfg,
                'training_exposure':'author_metadata_includes_rice_Kitaake_not_species_unseen'}
    helixer = {'container_uri':HELIXER_IMAGE,'container_path':str(hel_image.relative_to(OUT)),
               'native_metadata':hel_native,'model_source_commit':HELIXER_COMMIT,
               'model_selection':'lowest_official_land_plant_priority_not_DEV_score',
               'selected_model':selected,'model_list_receipt':helixer_list_receipt,
               'weight':str(helixer_weight.relative_to(OUT)),'weight_receipt':helixer_weight_receipt}
    for name,record in (('Tiberius',tiberius),('Helixer',helixer)):
        with (receipts/(name+'.json')).open('x') as handle:
            json.dump(record,handle,indent=2); handle.write('\n')
    report = {'experiment':'M28-NATIVE-BASELINE-PREP','completed':True,
              'inputs_complete':True,'runtime_assets_complete':True,'inputs':input_rows,
              'baselines':{'ANNEVO':annevo,'Tiberius':tiberius,'Helixer':helixer},
              'native_inference_executed':False,'GPU_forward_executed':False,
              'reference_used':False,'optimizer_steps':0,'test_setaria_used':False,
              'old_assets_or_jobs_modified':False,'seconds':time.perf_counter()-began,
              'limitations':'staged assets and versions only; GPU compatibility, native predictions, common-ruler accuracy and independent generalization unverified'}
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps({'completed':True,'seconds':report['seconds'],
                      'native_inference_executed':False}),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=OUT)
    OUT = parser.parse_args().output_dir.resolve()
    main()
