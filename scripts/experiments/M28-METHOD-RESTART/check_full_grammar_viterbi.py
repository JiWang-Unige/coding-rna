#!/usr/bin/env python3
"""New checkpoint max-semiring and raw wiring; not the stopped A43 prototype."""
import json
import os
from pathlib import Path
import resource
import sys
import time
import unittest
import torch
from src.m28.full_grammar_viterbi import _forward,_trace
from test_m28_full_grammar_viterbi import FullGrammarViterbiTests
from test_m28_full_grammar_raw import FullGrammarRawTests

OUT=Path('outputs/M28-FULL-GRAMMAR-VITERBI-CHECKS')
CONFIG=dict(scope='four_new_max_automaton_math_plus_four_existing_raw_interface_tests',
    graph='validated_complete_boundary_state_graph_not_stopped_A43_endpoint_arena',
    score='FP64_per_base_role_minus_background',ties='background_continuation_then_state_table_arrival',
    predecessor='uint32_edge_tokens_single_recomputed_block',default_checkpoint_bp=8192,
    actual_trained_checkpoint_read=False,actual_research_DNA_read=False,optimizer_steps=0)


def main():
    began=time.perf_counter(); job=os.environ['SLURM_JOB_ID']
    with (OUT/'started.json').open('x') as handle: json.dump(dict(job_id=job,configuration=CONFIG),handle,indent=2)
    loader=unittest.defaultTestLoader
    suite=unittest.TestSuite([loader.loadTestsFromTestCase(FullGrammarViterbiTests),loader.loadTestsFromTestCase(FullGrammarRawTests)])
    result=unittest.TextTestRunner(stream=sys.stdout,verbosity=2).run(suite)
    summary=dict(completed=result.wasSuccessful() and result.testsRun==8,job_id=job,configuration=CONFIG,
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        forward_nopython_signatures=len(_forward.nopython_signatures),trace_nopython_signatures=len(_trace.nopython_signatures),
        seconds=time.perf_counter()-began,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        research_inputs_opened=False,actual_native_GLM_forward=False,actual_trained_checkpoint_read=False,
        optimizer_steps_this_run=0,DEV_used=False,independent_species_used=False,A33_outputs_opened=False,
        stopped_A43_endpoint_prototype_used=False,A44_execution_chain_reopened=False,
        interpretation='new_exact_max_representation_artificial_math_and_raw_wiring_not_chromosome_or_model_performance')
    with (OUT/'summary.json').open('x') as handle: json.dump(summary,handle,indent=2,allow_nan=False)
    print(json.dumps(summary,allow_nan=False),flush=True)
    if not summary['completed']: raise SystemExit(1)


if __name__=='__main__':
    torch.set_num_threads(2); main()
