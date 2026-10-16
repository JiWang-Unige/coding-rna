#!/usr/bin/env python3
"""CPU artificial checks for the changed native boundary and raw caller only."""
import os
import sys
import time
import unittest
import torch
import fit_full_grammar_joint24 as J


def main():
    began = time.perf_counter()
    sys.path.insert(0,str(J.ROOT/'tests'))
    import test_m28_full_grammar_joint24_raw as T
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(T.JointRawTests))
    J.C.write(J.CHECKS/'summary.json',dict(completed=result.wasSuccessful(),job_id=os.environ['SLURM_JOB_ID'],
        configuration=J.CONFIG,protocol=J.PROTOCOL,tests_run=result.testsRun,failures=len(result.failures),
        errors=len(result.errors),research_data_read=False,models_or_weights_read=False,
        native_forward=False,optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,
        seconds=time.perf_counter()-began))
    if not result.wasSuccessful():
        raise SystemExit(1)


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
