#!/usr/bin/env python3
"""Numerical candidate qualification only; production posterior stays unchanged."""
import os
import time
import unittest
from unittest.mock import patch

import torch
import src.m28.full_grammar_gene_posterior as G
from src.m28.full_grammar_scaled_mass import _scaled_prefix_queries, _scaled_suffix_queries
import test_m28_full_grammar_gene_posterior as T
import test_m28_full_grammar_scaled_mass as S
import readout_full_grammar_gene_posterior_train as P

OUT = P.ROOT/'outputs/M28-FULL-GRAMMAR-SCALED-MASS-CHECKS'/os.environ['SLURM_JOB_ID']


def main():
    began = time.perf_counter()
    P.write(OUT/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],protocol=P.PROTOCOL,
        candidate_only=True,production_algorithm_modified=False,research_DNA_read=False,reference_used=False))
    suite = unittest.TestSuite([unittest.defaultTestLoader.loadTestsFromModule(S),
        unittest.defaultTestLoader.loadTestsFromModule(T)])
    # Only this test process routes the four established public-API tests to the candidate.
    with patch.multiple(G,_prefix_queries=_scaled_prefix_queries,_suffix_queries=_scaled_suffix_queries):
        result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful(): raise RuntimeError('Scaled mass candidate failed qualification')
    P.write(OUT/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],protocol=P.PROTOCOL,
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        candidate='all_state_background_gauge_shift_every256bp_compensated_offset_restored_at_queries',
        original_probability_and_mass_tolerances_unchanged=True,production_algorithm_modified=False,
        research_DNA_read=False,reference_used=False,native_forward=False,checkpoint_read=False,
        optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,seconds=time.perf_counter()-began))


if __name__=='__main__':
    torch.set_num_threads(2)
    main()
