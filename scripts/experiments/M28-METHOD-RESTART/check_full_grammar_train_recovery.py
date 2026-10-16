#!/usr/bin/env python3
"""Synthetic tests plus a real, completed LOCAL reference replay; no fit or inference."""
import json
import os
from pathlib import Path
import time
import unittest
import compare_full_grammar_train_recovery as R


def main():
    began = time.perf_counter()
    tests = unittest.defaultTestLoader.discover(str(R.ROOT/'tests'),pattern='test_m28_full_grammar_train_recovery.py')
    result = unittest.TextTestRunner(verbosity=2).run(tests)
    if not result.wasSuccessful(): raise RuntimeError('TRAIN recovery synthetic checks failed')
    C,_,reports = R.load_reports((R.ANCHOR,R.LOCAL))
    rows,source = R.source_records(reports[R.ANCHOR],C)
    R.replay(rows,R.LOCAL,reports[R.LOCAL],C)
    value = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],tests_run=result.testsRun,
        plan=R.PLAN,completed_LOCAL_replay_verified=True,source_LOCAL_evaluation_job=reports[R.LOCAL]['job_id'],
        source_paired_reference_job=source['job_id'],reference_primary_chains=len(rows),
        LOCAL=dict(pooled=R.describe(rows,R.LOCAL),
            per_species={s:R.describe([r for r in rows if r['species']==s],R.LOCAL) for s in C.SPECIES},
            measurements={k:reports[R.LOCAL][k] for k in ('pooled','per_species',
                'macro_species_exact_chain_F1','worst_species_exact_chain_F1')}),
        fixture_only=False,saved_paired_reference_profile_read=True,prediction_GFF_read=True,
        DNA_read=False,real_reference_GFF_read=False,checkpoint_read=False,new_model_inference=False,
        predictions_modified=False,optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,
        full_eight_model_comparison_completed=False,scientific_winner_or_DEV_candidate_selected=False,
        seconds=time.perf_counter()-began)
    with (R.CHECKS/'summary.json').open('x') as handle:
        json.dump(value,handle,indent=2,allow_nan=False); handle.write('\n')


if __name__ == '__main__': main()
