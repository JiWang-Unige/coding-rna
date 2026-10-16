#!/usr/bin/env python3
"""Final-checkpoint load and actual production skyline raw/GFF3 binding."""
import json
import os
from pathlib import Path
import sys
import time
import unittest
from unittest.mock import patch
import torch
import annotate_full_grammar_coarse_channel as R
from test_m28_full_grammar_raw import FullGrammarRawTests
from src.m28.full_grammar_skyline import decode_full_grammar_skyline, _forward, _trace

OUT=Path('outputs/M28-FULL-GRAMMAR-PRODUCTION-RAW-CHECKS')
CHECKPOINT=Path('outputs/M28-FULL-GRAMMAR-COARSE-CHANNEL-FIRST-FIT/step_004608.pt')


def main():
    began=time.perf_counter(); calls=[]
    if R.decode_full_grammar_skyline is not decode_full_grammar_skyline:
        raise ValueError('Production caller must bind the actual exact skyline')
    state=torch.load(CHECKPOINT,map_location='cpu',weights_only=True)
    configuration=R.validate_checkpoint(state)
    model=R.M28Core('C0',2048,128)
    model.load_state_dict(state['model'],strict=True); del state,model
    def instrument(sequence,logits,cap):
        if cap!=128: raise ValueError('Unchanged cap128 required')
        result=decode_full_grammar_skyline(sequence,logits,cap)
        calls.append(dict(length=len(sequence),gain=result['gain'],chains=result['chains'],
            pruning=result['pruning'],score_accumulation=result['score_accumulation']))
        return result
    suite=unittest.TestSuite(FullGrammarRawTests(name) for name in (
        'test_long_gene_crosses_native_window_and_exports_both_strands',
        'test_actual_codon_carry_minus_GFF_phase_and_intronic_N'))
    with patch.object(R,'decode_full_grammar_skyline',instrument):
        result=unittest.TextTestRunner(stream=sys.stdout,verbosity=2).run(suite)
    completed=(result.wasSuccessful() and result.testsRun==2 and len(calls)==6
        and bool(_forward.nopython_signatures) and bool(_trace.nopython_signatures))
    summary=dict(completed=completed,job_id=os.environ['SLURM_JOB_ID'],
        checkpoint=str(CHECKPOINT),configuration=configuration,final_checkpoint_strict_load=True,
        production_backend='exact_exon_budget_skyline',tests_run=result.testsRun,
        actual_decoder_calls=len(calls),invocations=calls,seconds=time.perf_counter()-began,
        actual_native_GLM_forward=False,research_DNA_or_labels_opened=False,
        DEV_used=False,optimizer_steps_this_run=0,
        interpretation='production_binding_and_GFF3_correctness_not_model_accuracy')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(summary,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps(summary,allow_nan=False),flush=True)
    if not completed: raise SystemExit(1)


if __name__=='__main__':
    torch.set_num_threads(2); main()
