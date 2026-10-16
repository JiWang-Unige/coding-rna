#!/usr/bin/env python3
"""Acquire the fixed public Carbon-A pipeline without reading research inputs."""
import ast
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from huggingface_hub import HfApi, hf_hub_download, snapshot_download
from packaging.specifiers import SpecifierSet
import torch

ROOT = Path(__file__).resolve().parents[3]
REPOSITORY = 'HuggingFaceBio/Carbon-A-1.2B'
REVISION = 'f6789ee79e6ef397de333577ffec1cb734496620'
SOURCE = 'GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview'
SOURCE_REVISION = 'b0483c23b6b63787b61a6d3a204a9b517d6ba345'
ASSETS = ROOT/'refs/carbon-a'/REVISION
OUT = ROOT/'outputs/M28-CARBON-A-ASSET-PREPARATION'/os.environ['SLURM_JOB_ID']


def main():
    began = time.perf_counter()
    info = HfApi(token=False).model_info(REPOSITORY, revision=REVISION, files_metadata=True)
    if info.sha != REVISION or info.private or info.gated:
        raise ValueError('Matching ungated public Carbon-A revision required')
    files = {item.rfilename: item for item in info.siblings}
    required = {'annotation_pipeline/run_annotation.sh', 'annotation_pipeline/src/model_inference.py',
        'annotation_pipeline/src/post_processing.py', 'annotation_pipeline/src/viterbi.py',
        'annotation_pipeline/weights/transition_network.safetensors',
        'annotation_pipeline/weights/confidence_network.safetensors', 'model.safetensors',
        'provenance.json', 'requirements.txt'}
    if len(files) != 35 or not required <= set(files):
        raise ValueError('Fixed complete release file inventory differs')
    weight = files['model.safetensors']
    source = Path(hf_hub_download(SOURCE, 'model.safetensors', revision=SOURCE_REVISION,
        local_files_only=True, token=False)).resolve()
    if (source.name != weight.lfs.sha256 or source.stat().st_size != weight.size
            or weight.size != 4614654504):
        raise ValueError('Original local HF cache object does not match the published shared weight identity')
    if ASSETS.exists():
        raise FileExistsError(ASSETS)
    snapshot_download(REPOSITORY, revision=REVISION, token=False, local_dir=ASSETS,
        ignore_patterns=['model.safetensors'], max_workers=2)
    (ASSETS/'model.safetensors').symlink_to(source)
    for name, item in files.items():
        if not (ASSETS/name).is_file() or (ASSETS/name).stat().st_size != item.size:
            raise ValueError('Incomplete fixed asset: '+name)
    provenance = json.loads((ASSETS/'provenance.json').read_text())
    if (provenance['source_repository'] != SOURCE or provenance['source_revision'] != SOURCE_REVISION
            or provenance['preserved_source_files']['model.safetensors']['sha256'] != weight.lfs.sha256):
        raise ValueError('Author provenance and published weight identity disagree')
    python_files = sorted(name for name in files if name.endswith('.py'))
    shell_files = sorted(name for name in files if name.endswith('.sh'))
    for name in python_files:
        ast.parse((ASSETS/name).read_text(), filename=name)
    for name in shell_files:
        subprocess.run(['bash', '-n', str(ASSETS/name)], check=True)
    sys.path.insert(0, str(ASSETS/'annotation_pipeline/src'))
    from transition_network import load_transition_network
    from confidence_network import load_confidence_network
    networks = {}
    for name, loader in [('transition_network', load_transition_network),
                         ('confidence_network', load_confidence_network)]:
        network = loader(ASSETS/'annotation_pipeline/weights'/(name+'.safetensors'), torch.device('cpu'))
        state = network.state_dict()
        if not all(value.dtype == torch.float32 and bool(torch.isfinite(value).all()) for value in state.values()):
            raise ValueError('Nonfinite or wrong-dtype auxiliary network weights')
        networks[name] = dict(strict_CPU_load=True, tensors=len(state),
            parameters=sum(value.numel() for value in state.values()),
            bytes=(ASSETS/'annotation_pipeline/weights'/(name+'.safetensors')).stat().st_size,
            finite_FP32=True, forward_executed=False)
    dependencies = {}
    for name, required_version in [('torch', '>=2.7,<3'), ('transformers', '>=4.56,<5')]:
        version = importlib.metadata.version(name)
        dependencies[name] = dict(installed=version, required=required_version,
            satisfied=SpecifierSet(required_version).contains(version))
    report = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'],
        protocol='reports/M28-METHOD-RESTART/carbon_a_pipeline_readiness.md',
        repository=REPOSITORY, revision=REVISION, public_file_count=len(files),
        asset_directory=str(ASSETS.relative_to(ROOT)),
        downloaded_non_model_files=len(files)-1,
        downloaded_non_model_bytes=sum(item.size for name, item in files.items() if name != 'model.safetensors'),
        all_published_file_sizes_match=True, python_AST_files=len(python_files), bash_syntax_files=len(shell_files),
        source_repository=SOURCE, source_revision=SOURCE_REVISION,
        shared_model_weight_reused_from_existing_HF_cache=True, model_weight_target=str(source),
        model_weight_bytes=weight.size, published_HF_LFS_object_identity=weight.lfs.sha256,
        model_weight_binary_read=False, model_weight_checksum_recomputed=False,
        auxiliary_networks=networks, current_environment_dependencies=dependencies,
        supported_native_runtime_ready=all(row['satisfied'] for row in dependencies.values()),
        environment_changed=False, actual_native_GLM_forward=False,
        research_DNA_or_reference_read=False, raw_annotation_executed=False,
        DEV_used=False, sealed_data_used=False, optimizer_steps_this_run=0,
        seconds=time.perf_counter()-began,
        interpretation='fixed_assets_and_auxiliary_CPU_load_not_raw_inference_or_accuracy')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(report, allow_nan=False), flush=True)


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
