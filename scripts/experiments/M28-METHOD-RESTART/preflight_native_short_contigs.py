#!/usr/bin/env python3
"""Unchanged native baselines on two artificial short contigs, no accuracy labels."""
import argparse
from collections import Counter
import json
import os
import subprocess
import time
import run_native_baseline as N

OUT = N.ROOT/'outputs/M28-NATIVE-SHORT-CONTIG-PREFLIGHT'
SOURCE = N.ROOT/'outputs/M28-TRANSFER-SHORT-RAW-PREFLIGHT'
PROTOCOL = 'reports/M28-METHOD-RESTART/native_short_contig_preflight.md'
FIXTURES = {'artificial_332': 'ATG'+'CAA'*108+'TAA'+'CC',
            'artificial_7536': 'C'*6140+'ATG'+'CAA'*463+'TAA'+'C'}
LENGTHS = {k: len(v) for k, v in FIXTURES.items()}


def commands(method, run, visible):
    native, probe = N.commands(method, 'artificial_short_preflight', run, visible, cpus=8)
    option = {'ANNEVO': '-g', 'Tiberius': '--genome', 'Helixer': '--fasta-path'}[method]
    native[native.index(option)+1] = str(run/'artificial.fa') if method == 'ANNEVO' else '/m28-run/artificial.fa'
    return native, probe


def check_output(path, logfile, method, lengths):
    if N.FINISHED[method] not in logfile.read_text():
        raise RuntimeError(method+' native completion marker missing; inspect command.log')
    counts = Counter(); by_contig = {name: Counter() for name in lengths}
    with path.open() as handle:
        for lineno, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            row = line.rstrip('\n').split('\t')
            if len(row) != 9 or row[0] not in lengths:
                raise ValueError('Unexpected GFF fields/scope at line '+str(lineno))
            a, b = int(row[3]), int(row[4])
            if not 1 <= a <= b <= lengths[row[0]]:
                raise ValueError('Prediction crosses the true contig boundary at line '+str(lineno))
            if row[2] == 'CDS' and (row[6] not in ('+', '-') or row[7] not in ('0', '1', '2') or
                    not any(a.startswith('Parent=') and a[7:] for a in row[8].split(';'))):
                raise ValueError('Actual CDS strand/phase/Parent invalid at line '+str(lineno))
            counts[row[2]] += 1; by_contig[row[0]][row[2]] += 1
    return dict(bytes=path.stat().st_size, feature_counts=dict(counts),
        feature_counts_per_contig={k: dict(v) for k, v in by_contig.items()},
        native_completion_marker=N.FINISHED[method],
        empty_prediction_is_not_filtered_or_rescued=not counts.get('CDS', 0),
        nonempty_CDS_structure_observed=bool(counts.get('CDS', 0)))


def source_identity():
    if (SOURCE/'STATUS').read_text().strip() != 'COMPLETED' or (SOURCE/'JOBID').read_text().strip() != '13401642':
        raise ValueError('Closed original artificial FASTA preflight required')
    source = json.loads((SOURCE/'summary.json').read_text())
    if (not source['completed'] or source['input_kind'] != 'artificial_DNA_only' or
            source['input_lengths'] != LENGTHS or source['biological_genome_or_annotation_opened']):
        raise ValueError('Artificial input source boundary differs')
    records = {}; current = None
    for line in (SOURCE/'artificial.fa').read_text().splitlines():
        if line.startswith('>'):
            current = line[1:].split()[0]
            if current in records:
                raise ValueError('Repeated artificial input record')
            records[current] = ''
        elif line.strip():
            records[current] += line.strip().upper()
    if records != FIXTURES:
        raise ValueError('Original deterministic artificial DNA changed')
    prep = json.loads((N.PREP/'summary.json').read_text())
    if ((N.PREP/'STATUS').read_text().strip() != 'COMPLETED' or
            not all(prep[k] for k in ('completed', 'inputs_complete', 'runtime_assets_complete')) or
            prep['reused_asset_root'] != str(N.ASSET_PREP.relative_to(N.ROOT)) or
            any(prep[k] for k in ('reference_used', 'test_setaria_used', 'native_inference_executed'))):
        raise ValueError('Closed original native runtime identity required')
    baselines = prep['baselines']; tib = baselines['Tiberius']
    if (baselines['ANNEVO']['source_commit'] != '42c920f70f52f96c7c06c665a639ec47d9304fff' or
            tib['native_metadata']['source_commit'] != '871d4979568ae09e56de7eca5d7489046ee833ea' or
            tib['native_metadata']['version'] != '2.0.8' or
            tib['base_container_actual_metadata']['version'] != '2.0.7' or
            not tib['native_metadata']['CPU_model_load_completed'] or
            baselines['Helixer']['native_metadata']['version'] != '0.3.7' or
            baselines['Helixer']['native_metadata']['source_commit'] != 'd17bb496b8842542590644b2437798711c0125d5'):
        raise ValueError('Fixed native versions differ')
    return baselines


def main(task):
    began = time.perf_counter(); method = N.METHODS[task]; out = OUT/method
    baselines = source_identity(); visible = os.environ['CUDA_VISIBLE_DEVICES']
    if not visible or visible in ('-1', 'NoDevFiles'):
        raise ValueError('Actual Slurm GPU assignment required')
    with (out/'artificial.fa').open('x') as handle:
        for name, sequence in FIXTURES.items():
            handle.write('>'+name+' artificial DNA, not an organism\n'+sequence+'\n')
    (out/'temporary').mkdir()
    native, probe = commands(method, out, visible)
    with (out/'commands.json').open('x') as handle:
        json.dump(dict(native=native, GPU_probe=probe, CUDA_VISIBLE_DEVICES=visible), handle, indent=2)
    with (out/'gpu_probe.stderr').open('x') as err:
        result = subprocess.run(probe, check=True, text=True, stdout=subprocess.PIPE, stderr=err)
    (out/'gpu_probe.stdout').write_text(result.stdout)
    gpu = json.loads(result.stdout.strip().splitlines()[-1])
    with (out/'command.log').open('x') as handle:
        subprocess.run(native, check=True, stdout=handle, stderr=subprocess.STDOUT)
    gff = check_output(out/'prediction.gff3', out/'command.log', method, LENGTHS)
    summary = dict(completed=True, experiment='M28-NATIVE-SHORT-CONTIG-PREFLIGHT', method=method,
        protocol=PROTOCOL, input_kind='artificial_DNA_only', input_lengths=LENGTHS,
        original_artificial_source_identity_verified=True, original_native_runtime_identity_verified=True,
        baseline_asset_receipt=baselines[method], GPU_runtime=gpu, output_check=gff,
        raw_GFF='prediction.gff3', requested_native_operation_points_unchanged=True,
        native_inference_executed=True, reference_used=False, labels_or_feature_cache_used=False,
        optimizer_steps=0, new_parameters=0, biological_genome_or_annotation_opened=False,
        DEV_outputs_opened=False, test_setaria_used=False, independent_panel_selected=False,
        seconds=time.perf_counter()-began,
        limitations='Two artificial lengths: native completion and produced GFF structure only. '
            'Not accuracy, all-contig model-consideration proof, fair speed, fragmentation-scale support or generalization.')
    with (out/'summary.json').open('x') as handle:
        json.dump(summary, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(dict(completed=True, method=method, seconds=summary['seconds'], output_check=gff)), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', required=True, type=int, choices=range(3))
    main(parser.parse_args().task)
