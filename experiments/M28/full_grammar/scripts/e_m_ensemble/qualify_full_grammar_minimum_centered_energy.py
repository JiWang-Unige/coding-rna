#!/usr/bin/env python3
"""Three synthetic local-agreement tests and immutable source binding only."""
import json
import os
from pathlib import Path
import time
import unittest

import torch
import test_m28_minimum_centered_energy_pool as T
import run_full_grammar_e_m_energy_pool_train as Q


if __name__=='__main__':
    torch.set_num_threads(2); began=time.perf_counter()
    job=os.environ['SLURM_JOB_ID']
    out=Q.ROOT/'outputs/M28-FULL-GRAMMAR-MINIMUM-CENTERED-CHECKS'/job
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(T))
    if not result.wasSuccessful(): raise RuntimeError('Fixed minimum-centered energy arithmetic failed')
    sources={}
    for name in Q.P.MODELS:
        _,_,reports=Q.P.source(name)
        sources[name]={seqid:str(Q.emissions(name,task,reports[seqid]))
            for task,(_,seqid,_) in enumerate(Q.P.E.I.SCOPE)}
    report=dict(completed=True,job_id=job,
        protocol='reports/M28-METHOD-RESTART/full_grammar_minimum_centered_energy_train_plan.md',
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        all_twelve_completed_emission_sources_verified=True,source_folders=sources,
        energy='minimum_of_FP64_per_branch_role_minus_own_background',
        equivalent_fixed_penalty='mean_centered_energy_minus_half_absolute_disagreement',
        whole_gene_minimum_used=False,research_DNA_read=False,research_reference_read=False,
        role_arrays_read=False,checkpoint_read=False,native_forward=False,
        optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,seconds=time.perf_counter()-began)
    with (out/'summary.json').open('x') as handle: json.dump(report,handle,indent=2,allow_nan=False)
    print(json.dumps(report),flush=True)
