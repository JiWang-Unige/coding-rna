#!/usr/bin/env python3
"""Artificial-only joint encoding, two-head slot alignment and checkpoint identity."""
import json
import os
import time
import unittest
from unittest.mock import patch
import torch
from torch import nn
import annotate_full_grammar_joint24_native_cds_prior as A
from native_cds_role_prior import NativeCDSRolePrior,native_cds_log_odds

OUT = A.R.M.C.ROOT/'outputs/M28-FULL-GRAMMAR-JOINT24-NATIVE-PRIOR-CHECKS'


def heads():
    values = nn.ModuleList([nn.Linear(2048,12,bias=False) for _ in range(2)]).to(torch.bfloat16)
    with torch.no_grad():
        for index,head in enumerate(values):
            head.weight.zero_()
            for slot in range(6):
                head.weight[2*slot+1,0] = (index+1)*(slot+1)
    return values.eval().requires_grad_(False)


class JointPriorChecks(unittest.TestCase):
    def test_only_new_combined_final_identity_accepted(self):
        state = dict(configuration=dict(A.IDENTITY),protocol=A.PROTOCOL,step=9216,epoch=6,
            encoder_revision=A.R.F.REVISION)
        self.assertEqual(A.validate_checkpoint(state),A.IDENTITY)
        for configuration in (A.N.IDENTITY,A.J.IDENTITY):
            bad = dict(state,configuration=configuration)
            with self.assertRaises(ValueError): A.validate_checkpoint(bad)
        for field,value in (('step',4609),('epoch',5),('protocol',A.N.PROTOCOL)):
            bad = dict(state,**{field:value})
            with self.assertRaises(ValueError): A.validate_checkpoint(bad)

    def test_both_six_position_heads_align_to_all_24576_bases(self):
        features = torch.zeros(1,4096,2048)
        tokens = torch.arange(4096)%8
        features[0,:,0] = tokens
        actual = native_cds_log_odds(heads(),features)
        expected = (tokens[:,None]*torch.arange(1,7)[None,:]).reshape(-1).float()
        self.assertEqual(tuple(actual.shape),(1,24576,2))
        self.assertTrue(torch.equal(actual[0,:,0],expected))
        self.assertTrue(torch.equal(actual[0,:,1],expected*2))

    def test_zero_prior_joint_window_and_short_valid_rows(self):
        torch.manual_seed(0)
        core = NativeCDSRolePrior().eval()
        features = torch.zeros(1,4096,2048)
        features[0,:,0] = torch.arange(4096)%8
        adapter = A.N.NativeCDSHeadAdapter(core,heads()).eval()
        for sequence in ('N','ATGCGTNAGTTTAA'):
            dna = A.R.one_hot(sequence,A.R.M.WIDTH)[None]
            with torch.inference_mode():
                expected = A.R.M28Core.forward(core,features,dna,[len(sequence)])['segmentation_logits'][0,:len(sequence)]
                with patch.object(A.J,'encode_window',return_value=features):
                    actual = A.J.window_logits(adapter,None,None,sequence,'cpu')
            self.assertEqual(tuple(actual.shape),(len(sequence),15))
            self.assertTrue(torch.equal(actual,expected))


if __name__ == '__main__':
    torch.set_num_threads(2)
    began = time.perf_counter()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(JointPriorChecks))
    with (OUT/'summary.json').open('x') as handle:
        json.dump(dict(completed=result.wasSuccessful(),job_id=os.environ['SLURM_JOB_ID'],identity=A.IDENTITY,
            tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
            research_data_read=False,models_or_weights_read=False,native_forward_run=False,
            optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,seconds=time.perf_counter()-began),
            handle,indent=2,allow_nan=False); handle.write('\n')
    if not result.wasSuccessful():
        raise SystemExit(1)
