#!/usr/bin/env python3
"""Prospective skyline wiring through the unchanged raw caller, synthetic only.

Only the two existing cases that execute the changed decoder are run. A
test-local binding leaves the production raw backend and scientific identity
unchanged while the independent chromosome resource job is still pending.
"""
import json
import os
from pathlib import Path
import resource
import sys
import time
import unittest
from unittest.mock import patch

import torch
import annotate_full_grammar_coarse_channel as R
from test_m28_full_grammar_raw import FullGrammarRawTests
from src.m28.full_grammar_skyline import decode_full_grammar_skyline, _forward, _trace

OUT = Path('outputs/M28-FULL-GRAMMAR-SKYLINE-RAW-CHECKS')
CONFIG = dict(scope='two_existing_raw_cases_with_actual_new_skyline_decoder',
    cases=['30000bp_complete_gene_crossing_native_windows_both_strands',
           'actual_codon_carry_intronic_N_minus_GFF3_phase'],
    emissions='artificial_FP32_native_encoder_and_head_stubs', max_exons=128,
    backend='sorted_sparse_frontier_exact_exon_budget_skyline',
    binding='temporary_test_local_production_raw_caller_unchanged',
    unchanged_owner_and_checkpoint_identity_cases_rerun=False,
    research_inputs_opened=False, actual_native_GLM_forward=False,
    actual_trained_checkpoint_read=False, optimizer_steps=0)


def main():
    began = time.perf_counter()
    job = os.environ['SLURM_JOB_ID']
    with (OUT/'started.json').open('x') as handle:
        json.dump(dict(job_id=job, configuration=CONFIG), handle, indent=2)
    calls = []

    def prospective_decode(sequence, logits, cap):
        if cap != 128:
            raise ValueError('Raw caller must retain its original cap128')
        result = decode_full_grammar_skyline(sequence, logits, cap)
        if result['pruning'] != 'exact_same_stage_actual_carry_lower_exon_budget_dominance_max_only':
            raise ValueError('Prospective binding did not execute the new exact skyline')
        calls.append(dict(length=len(sequence), gain=result['gain'],
            selected_chains=result['chains'], chain_gains=result['chain_gains'], counts=result['counts'],
            tie_policy=result['tie_policy'], score_accumulation=result['score_accumulation']))
        return result

    names = ['test_long_gene_crosses_native_window_and_exports_both_strands',
             'test_actual_codon_carry_minus_GFF_phase_and_intronic_N']
    suite = unittest.TestSuite(FullGrammarRawTests(name) for name in names)
    with patch.object(R, 'decode_full_grammar_viterbi', prospective_decode):
        result = unittest.TextTestRunner(stream=sys.stdout, verbosity=2).run(suite)
    completed = (result.wasSuccessful() and result.testsRun == 2 and len(calls) == 6
        and bool(_forward.nopython_signatures) and bool(_trace.nopython_signatures))
    summary = dict(completed=completed, job_id=job, configuration=CONFIG,
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        actual_skyline_decoder_invocations=len(calls), invocations=calls,
        forward_nopython_signatures=len(_forward.nopython_signatures), trace_nopython_signatures=len(_trace.nopython_signatures),
        seconds=time.perf_counter()-began, peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        production_raw_backend_changed=False, actual_chromosome_resource_result_claimed=False,
        DEV_used=False, independent_species_used=False, A33_outputs_opened=False,
        stopped_A43_endpoint_prototype_used=False, dense_decoder_executed=False, A44_execution_chain_reopened=False,
        interpretation='prospective_raw_interface_GFF3_with_new_exact_decoder_not_real_checkpoint_or_model_performance')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(summary, handle, indent=2, allow_nan=False)
    print(json.dumps(summary, allow_nan=False), flush=True)
    if not completed:
        raise SystemExit(1)


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
