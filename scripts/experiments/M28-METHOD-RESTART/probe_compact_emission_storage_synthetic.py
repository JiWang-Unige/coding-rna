#!/usr/bin/env python3
"""Fixed artificial storage measurements and two same-answer contig replays."""
import argparse
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

import torch

import profile_contig_emission_lattice_synthetic as P
import profile_emission_lattice_synthetic as S
from replay_indexed_contig_synthetic import ANSWER_FIELDS, IDENTITY_FIELDS
from replay_phase_block_contig_synthetic import SCOPE_FIELDS
from src.m28.contig_emission_lattice import stitch_contig_emissions
from src.m28.emission_lattice import _Margins
from src.m28.emission_lattice_compact import _CompactMargins, decode_compact_emission_lattice


SOURCE = Path('outputs/M28-A42-CACHED-METADATA-REPLAY')
SOURCE_ARRAY = '13423648'
MARGIN_LENGTH = 1_000_000


def save(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)


def current_rss():
    for line in Path('/proc/self/status').read_text().splitlines():
        if line.startswith('VmRSS:'):
            return int(line.split()[1])
    raise ValueError('Linux VmRSS measurement unavailable')


def margin_worker(kind):
    torch.set_num_threads(2)
    logits = torch.zeros(MARGIN_LENGTH, 15, dtype=torch.float32)
    baseline = current_rss()
    began = time.perf_counter()
    margins = (_Margins if kind == 'original' else _CompactMargins)(logits)
    seconds = time.perf_counter() - began
    numeric_bytes = (sum(getattr(margins, key).nbytes for key in ('gain', 'prefix', 'coding'))
                     if kind == 'compact' else None)
    if kind == 'compact' and numeric_bytes != 264 * MARGIN_LENGTH + 144:
        raise ValueError('Persistent compact FP64 storage differs from its declared arrays')
    if any(margins.interval(0, MARGIN_LENGTH, role) != 0 for role in range(15)):
        raise ValueError('Artificial zero margins changed')
    return dict(kind=kind, length=MARGIN_LENGTH, input='all_zero_CPU_FP32_N_by_15',
        threads=2, baseline_rss_kib=baseline, after_constructor_rss_kib=current_rss(),
        process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        constructor_seconds=seconds, persistent_numeric_bytes=numeric_bytes,
        interpretation='separate_process_observation_includes_allocator_not_full_decoder_or_chromosome')


def replay(case, root, job_id):
    source = SOURCE / ('case_' + str(case))
    expected = SOURCE_ARRAY + '_' + str(case)
    if (source.joinpath('STATUS').read_text().strip() != 'COMPLETED' or
            source.joinpath('JOBID').read_text().strip() != expected):
        raise ValueError('The original cached-metadata artificial source is not closed')
    out = root / ('case_' + str(case))
    out.mkdir()
    dna_generator = torch.Generator().manual_seed(0)
    sequence = ''.join('ACGT'[i] for i in torch.randint(0, 4, (P.LENGTH,), generator=dna_generator).tolist())
    generator = torch.Generator().manual_seed(41)
    started = dict(case_index=case, case=P.CASES[case], length=P.LENGTH, windows=P.WINDOWS,
        strand='+', DNA_seed=0, potential_seed=41 if case else None,
        potential_definition='all_zero' if not case else 'independent_FP32_integers_minus3_to3_all15_columns',
        context_policy='nearest_containing_genomic_window_center_per_base_lower_start_ties',
        max_exons=128, threshold=0, null=0, role_weights=1, score_accumulation='FP64',
        motif_counts=S.motif_counts(sequence), synthetic_only=True, job_id=job_id,
        source_job=expected, started_utc=S.utc_now())
    save(out / 'started.json', started)

    def records():
        for i, (a, b) in enumerate(P.WINDOWS):
            logits = (torch.zeros(b - a, 15) if not case else
                      torch.randint(-3, 4, (b - a, 15), generator=generator).float())
            yield dict(window_index=i, sequence=sequence[a:b], local_potentials=logits)

    joined = stitch_contig_emissions(sequence, P.WINDOWS, '+', records())
    logits = joined['local_potentials']
    upper_bound = float((logits.double() - logits[:, :1].double()).max(dim=1).values.sum())
    print(json.dumps(dict(case=case, stage='new_compact_natural_decode_started')), flush=True)
    began = time.perf_counter()
    result = decode_compact_emission_lattice(sequence, logits, 128)
    seconds = time.perf_counter() - began
    S.verify_selected(sequence, logits, result, 128)
    if not 0 <= result['gain'] <= upper_bound or (not case and (result['chains'] or result['gain'] != 0)):
        raise ValueError('Pointwise upper bound or zero-null contract differs')
    record = dict(decode_seconds=seconds, pointwise_upper_bound=upper_bound,
        selected_chains=result['chains'], chain_gains=result['chain_gains'], gain=result['gain'],
        selected_count=len(result['chains']), counts=result['counts'],
        genomic_context_segments=joined['genomic_context_segments'])
    with (out / 'natural_records.jsonl').open('x') as handle:
        handle.write(json.dumps(record, allow_nan=False) + '\n')
    # Saved prior ANSWERS are read only after the new natural output closes.
    previous = json.loads(source.joinpath('summary.json').read_text())
    identity, actual = json.loads(json.dumps(started)), json.loads(json.dumps(record))
    if (previous['completed'] is not True or previous['job_id'] != expected or
            any(previous[key] != identity[key] for key in IDENTITY_FIELDS) or
            any(previous[key] is not False for key in SCOPE_FIELDS) or previous['optimizer_steps'] != 0):
        raise ValueError('Assigned artificial source configuration or scope changed')
    if any(actual[key] != previous['run'][key] for key in ANSWER_FIELDS):
        raise ValueError('Compact storage changes chains, gains, ownership or any cached-metadata counter')
    summary = dict(**started, completed=True, completed_utc=S.utc_now(), run=record,
        actual_original_configuration_and_saved_source_job_verified=True,
        actual_new_natural_output_closed_before_saved_answers_opened=True,
        actual_all_chain_identities_order_FP64_gains_and_20_counters_equal=True,
        actual_DNA_ORF_motifs_roles_cap_spans_and_integer_FP64_scores_recomputed=True,
        optimizer_steps=0, **{key: False for key in SCOPE_FIELDS})
    save(out / 'summary.json', summary)
    print(json.dumps(dict(case=case, completed=True, decode_seconds=seconds,
                          selected_count=len(result['chains']), gain=result['gain'])), flush=True)
    return summary


def main(out):
    torch.set_num_threads(2)
    job_id = os.environ['SLURM_JOB_ID']
    if (out.joinpath('UNIT_STATUS').read_text().strip() != 'COMPLETED' or
            out.joinpath('JOBID').read_text().strip() != job_id):
        raise ValueError('Actual unit checks in this job must complete first')
    measurements = {}
    for kind in ('original', 'compact'):
        result = subprocess.run([sys.executable, __file__, '--margin-worker', kind],
                                capture_output=True, text=True)
        if result.returncode:
            sys.stderr.write(result.stderr)
            result.check_returncode()
        measurement = json.loads(result.stdout)
        if measurement['kind'] != kind or measurement['length'] != MARGIN_LENGTH:
            raise ValueError('Assigned isolated storage worker identity changed')
        save(out / ('margin_' + kind + '.json'), measurement)
        measurements[kind] = measurement
        print(json.dumps(measurement), flush=True)
    cases = [replay(case, out, job_id) for case in (0, 1)]
    summary = dict(job_id=job_id, completed=True, completed_utc=S.utc_now(), synthetic_only=True,
        protocol='a42_compact_storage_preparation.md', margin_measurements=measurements,
        all_two_fixed_contig_replays_equal=True, case_counts=[row['run']['selected_count'] for row in cases],
        case_gains=[row['run']['gain'] for row in cases], optimizer_steps=0,
        **{key: False for key in SCOPE_FIELDS},
        frozen_or_research_decoder_replaced=False,
        interpretation='engineering_storage_and_same_answer_evidence_not_speed_or_biological_performance')
    save(out / 'summary.json', summary)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--margin-worker', choices=('original', 'compact'))
    args = parser.parse_args()
    if args.margin_worker:
        print(json.dumps(margin_worker(args.margin_worker), allow_nan=False))
    elif args.out:
        main(args.out)
    else:
        parser.error('an output directory or the isolated margin worker is required')
