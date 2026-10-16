#!/usr/bin/env python3
"""Bounded cProfile of the closed artificial phase-block replay, not a benchmark."""
import argparse
import cProfile
import json
import os
from pathlib import Path

import replay_phase_block_contig_synthetic as R
from src.m28.candidates import Budget
from src.m28.emission_lattice import _Live, _Solution


UNPROFILED = Path('outputs/M28-A42-PHASE-BLOCK-REPLAY')
UNPROFILED_ARRAY = '13423292'


def main(case, out):
    prior_root = UNPROFILED / ('case_' + str(case))
    prior_job = UNPROFILED_ARRAY + '_' + str(case)
    if (prior_root.joinpath('STATUS').read_text().strip() != 'COMPLETED' or
            prior_root.joinpath('JOBID').read_text().strip() != prior_job):
        raise ValueError('Assigned unprofiled phase-block case must be closed first')
    profiler = cProfile.Profile()
    # Reuse the unchanged replay, including independent DNA/role checks and
    # saving/closing the new natural answer before opening original answers.
    # This covers reconstruction, decode, verification and JSON, not imports.
    # Function cumulative times distinguish decode from that surrounding work.
    profiler.enable()
    try:
        R.main(case, out)
    finally:
        profiler.disable()
    current = json.loads(out.joinpath('summary.json').read_text())
    prior = json.loads(prior_root.joinpath('summary.json').read_text())
    if (prior['completed'] is not True or prior['job_id'] != prior_job or
            any(prior[key] != current[key] for key in R.IDENTITY_FIELDS) or
            any(prior[key] is not False for key in R.SCOPE_FIELDS) or prior['optimizer_steps'] != 0):
        raise ValueError('Closed unprofiled phase-block configuration or scope differs')
    if any(current['run'][key] != prior['run'][key] for key in R.ANSWER_FIELDS):
        raise ValueError('Profiling changes chain identity, FP64 gains or any phase-block counter')
    # Standard pstats keys use file/line/co_name, so the generated dataclass
    # __init__ methods all labelled <string>:2 can overwrite one another.
    # Use distinct raw code entries, including their original self times.
    entries = profiler.getstats()
    constructor_names = {kind.__init__.__code__: kind.__name__ + '.__init__'
                         for kind in (_Live, _Solution, Budget)}
    rows = []
    standard_labels = []
    for entry in entries:
        code = entry.code
        if isinstance(code, str):
            filename, line, function = '~', 0, code
        else:
            filename, line, function = code.co_filename, code.co_firstlineno, code.co_name
        standard_labels.append((filename, line, function))
        function = constructor_names.get(code, function)
        rows.append(dict(file=filename, line=line, function=function,
            primitive_calls=entry.callcount - entry.reccallcount, total_calls=entry.callcount,
            self_seconds=entry.inlinetime, cumulative_seconds=entry.totaltime))
    by_self = sorted(rows, key=lambda row: row['self_seconds'], reverse=True)
    by_cumulative = sorted(rows, key=lambda row: row['cumulative_seconds'], reverse=True)
    decoder_functions = [row for row in rows if row['file'].endswith((
        '/src/m28/emission_lattice_phase_block.py', '/src/m28/emission_lattice.py',
        '/src/m28/exon_end_bound.py', '/src/m28/long_exon_phase.py', '/src/m28/candidates.py'))
        or row['function'] in constructor_names.values()]
    summary = dict(case_index=case,
        job_id=os.environ['SLURM_ARRAY_JOB_ID'] + '_' + os.environ['SLURM_ARRAY_TASK_ID'],
        unprofiled_phase_block_job=prior_job, completed=True,
        actual_unprofiled_configuration_scope_and_job_identity_verified=True,
        actual_new_natural_search_closed_before_prior_phase_block_answer_opened=True,
        actual_all_chains_order_FP64_gains_and_all_phase_block_counts_equal=True,
        profile_scope='unchanged_replay_main_reconstruction_decode_verification_JSON_not_imports',
        profiled_decode_wall_seconds=current['run']['decode_seconds'],
        unprofiled_decode_wall_seconds=prior['run']['decode_seconds'],
        actual_raw_code_entries_used_without_standard_label_collapse=True,
        standard_label_collision_surplus=len(entries) - len(set(standard_labels)),
        total_profiled_self_seconds=sum(row['self_seconds'] for row in rows),
        total_calls=sum(row['total_calls'] for row in rows),
        primitive_calls=sum(row['primitive_calls'] for row in rows),
        top_self=by_self[:25], top_cumulative=by_cumulative[:25],
        all_raw_function_entries=by_self,
        decoder_functions=sorted(decoder_functions, key=lambda row: row['self_seconds'], reverse=True),
        optimizer_steps=0, **{key: False for key in R.SCOPE_FIELDS},
        interpretation='instrumented_function_cost_diagnostic_not_unprofiled_timing_or_causal_speed_benchmark')
    with (out / 'profile_summary.json').open('x') as handle:
        json.dump(summary, handle, indent=2, allow_nan=False)
    print(json.dumps({key: value for key, value in summary.items()
        if key not in ('top_self', 'top_cumulative', 'decoder_functions', 'all_raw_function_entries')}, allow_nan=False), flush=True)
    print(json.dumps(dict(top_self=by_self[:10], decoder_functions=summary['decoder_functions']),
        allow_nan=False), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=int, choices=(0, 1), required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    main(args.case, args.out)
