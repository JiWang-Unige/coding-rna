#!/usr/bin/env python3
"""Synthetic pairing/count checks; no research records or predictions."""
import json
import os
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT/'outputs/M28-NATIVE-TRAIN-REFERENCE-CHECKS'


def main():
    tests = unittest.defaultTestLoader.discover(str(ROOT/'tests'),pattern='test_m28_native_train_reference.py')
    result = unittest.TextTestRunner(verbosity=2).run(tests)
    if not result.wasSuccessful(): raise RuntimeError('Synthetic native TRAIN pairing checks failed')
    report = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],tests_run=result.testsRun,
        fixture_only=True,DNA_read=False,real_reference_GFF_read=False,prediction_GFF_read=False,
        checkpoint_read=False,new_model_inference=False,DEV_used=False,sealed_data_used=False,
        optimizer_steps_this_run=0)
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')


if __name__ == '__main__': main()
