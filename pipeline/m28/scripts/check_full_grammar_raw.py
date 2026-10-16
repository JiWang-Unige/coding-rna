#!/usr/bin/env python3
"""Four artificial tests for the new raw caller; no research files or updates."""
import json
import os
from pathlib import Path
import resource
import sys
import time
import unittest
import torch
from test_m28_full_grammar_raw import FullGrammarRawTests

OUT=Path('outputs/M28-FULL-GRAMMAR-RAW-CHECKS')
CONFIG=dict(scope='four_new_raw_caller_artificial_integration_tests',
    native_encoder='stub_no_released_GLM_weights_or_forward',head='stub_integer_FP32_potentials_no_parameters',
    decoder='actual_nopython_compiled_complete_gene_set_cap128_FP64_no_pruning',
    tests=['long_gene_both_strands','carry_minus_phase_intronic_N','shifted_last_grid_owner_rows','fixed_checkpoint_identity'],
    actual_trained_checkpoint_read=False,actual_research_DNA_read=False,optimizer_steps=0)


def main():
    began=time.perf_counter(); job=os.environ['SLURM_JOB_ID']
    with (OUT/'started.json').open('x') as handle:
        json.dump(dict(job_id=job,configuration=CONFIG),handle,indent=2)
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(FullGrammarRawTests)
    result=unittest.TextTestRunner(stream=sys.stdout,verbosity=2).run(suite)
    summary=dict(completed=result.wasSuccessful(),job_id=job,configuration=CONFIG,
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        seconds=time.perf_counter()-began,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        research_inputs_opened=False,actual_native_GLM_forward=False,actual_trained_checkpoint_read=False,
        reference_or_labels_used=False,optimizer_steps_this_run=0,DEV_used=False,independent_species_used=False,
        A33_outputs_opened=False,A44_execution_chain_reopened=False,
        interpretation='new_raw_interface_artificial_checks_not_actual_model_or_chromosome_performance')
    with (OUT/'summary.json').open('x') as handle: json.dump(summary,handle,indent=2,allow_nan=False)
    print(json.dumps(summary,allow_nan=False),flush=True)
    if not result.wasSuccessful() or result.testsRun!=4: raise SystemExit(1)


if __name__=='__main__':
    torch.set_num_threads(2); main()
