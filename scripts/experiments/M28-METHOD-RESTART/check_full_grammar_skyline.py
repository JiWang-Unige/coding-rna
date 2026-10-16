#!/usr/bin/env python3
"""Actual graph proof and artificial math for the new sparse max skyline."""
import json
import os
from pathlib import Path
import resource
import sys
import time
import unittest

import torch
from src.m28.full_grammar_skyline import _forward, _trace
from test_m28_full_grammar_skyline import FullGrammarSkylineTests

OUT = Path('outputs/M28-FULL-GRAMMAR-SKYLINE-CHECKS')
CONFIG = dict(scope='six_new_graph_dominance_and_artificial_max_math_checks',
    grammar='unchanged_complete_boundary_graph', max_exons_up_to=128,
    pruning='same_stage_actual_carry_lower_count_no_worse_score_max_only',
    execution='sorted_reachable_frontier_not_dense_state_scan', checkpoint_default=8192,
    actual_research_DNA_read=False, actual_native_forward=False, optimizer_steps=0,
    raw_caller_backend_changed=False, stopped_dense_chromosome_job_rerun=False)


def main():
    began = time.perf_counter()
    job = os.environ['SLURM_JOB_ID']
    with (OUT/'started.json').open('x') as handle:
        json.dump(dict(job_id=job, configuration=CONFIG), handle, indent=2)
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(FullGrammarSkylineTests)
    result = unittest.TextTestRunner(stream=sys.stdout, verbosity=2).run(suite)
    summary = dict(completed=result.wasSuccessful() and result.testsRun == 6, job_id=job,
        configuration=CONFIG, tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        forward_nopython_signatures=len(_forward.nopython_signatures), trace_nopython_signatures=len(_trace.nopython_signatures),
        seconds=time.perf_counter()-began, peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        interpretation='exact_max_dominance_graph_and_artificial_math_not_scale_or_model_performance')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(summary, handle, indent=2, allow_nan=False)
    print(json.dumps(summary, allow_nan=False), flush=True)
    if not summary['completed'] or not _forward.nopython_signatures or not _trace.nopython_signatures:
        raise SystemExit(1)


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
