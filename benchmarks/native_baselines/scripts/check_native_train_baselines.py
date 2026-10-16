#!/usr/bin/env python3
"""Synthetic runner/ruler checks, then the actual six permitted TRAIN input files."""
import json
import os
from pathlib import Path
import unittest
import run_native_train_baseline as N
import evaluate_full_grammar_original_train as E


def main():
    if N.SCOPE != E.I.SCOPE:
        raise ValueError('Native and grammar original TRAIN scopes differ')
    suite = unittest.TestSuite()
    for pattern in ('test_m28_native_baseline_runner.py','test_m28_native_train_baselines.py'):
        suite.addTests(unittest.defaultTestLoader.discover(str(N.ROOT/'tests'),pattern=pattern))
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise SystemExit(1)
    prep = N.runtime_receipt()
    required = [N.N.ASSET_PREP/'assets/ANNEVO/annotation.py',
                N.N.ASSET_PREP/'assets/ANNEVO/saved_model/ANNEVO_Magnoliopsida.pt',
                N.N.ASSET_PREP/'assets/Tiberius-2.0.8.sif',N.N.ASSET_PREP/'assets/angiosperms.yaml',
                N.N.PREP/'site/tiberius/main.py',N.N.ASSET_PREP/'assets/Helixer-0.3.7.sif',
                N.ROOT/prep['baselines']['Helixer']['weight']]
    if any(not p.is_file() for p in required) or not (N.ROOT/prep['baselines']['Tiberius']['weights_directory']).is_dir():
        raise ValueError('Fixed native runtime file missing')
    inputs = []
    for index,(_,seqid,length) in enumerate(N.SCOPE):
        N.reconcile_train_input(index)
        inputs.append(N.check_fasta(N.input_fasta(seqid),seqid,length))
    report = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],scope=[list(s) for s in N.SCOPE],
        input_mode=N.INPUT_MODE,inputs=inputs,synthetic_tests_passed=result.testsRun,
        actual_original_TRAIN_DNA_checked=True,reference_used=False,DEV_used=False,sealed_data_used=False,
        GPU_forward_executed=False,optimizer_steps=0,baseline_assets_reused=True)
    with (N.CHECKS/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps({'completed':True,'synthetic_tests_passed':result.testsRun,'TRAIN_inputs_checked':len(inputs)}),flush=True)


if __name__ == '__main__': main()
