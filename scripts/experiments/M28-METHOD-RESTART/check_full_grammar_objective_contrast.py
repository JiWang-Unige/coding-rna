#!/usr/bin/env python3
"""Independent artificial observation likelihood and mixture gradient checks."""
import json
import math
import os
import time
import unittest
import torch
import fit_full_grammar_objective_contrast as C
from src.m28.local_coarse_channel import local_coarse_channel_nll
from src.m28.labels import REGION_COLUMNS
from check_full_grammar_coarse_channel import reference, targets


def local_reference(logits, target, epsilon):
    terms = []
    for i, row in enumerate(logits):
        fine, coarse = int(target['seg'][i]), int(target['region'][i])
        if fine < 0 and coarse < 0:
            continue
        # Enumerate all fifteen role outcomes, independently of masks/log-softmax.
        p = row.exp()/row.exp().sum()
        factors = [float(j == fine) if fine >= 0 else
                   (1-epsilon if j in REGION_COLUMNS[coarse] else epsilon/2) for j in range(15)]
        terms.append(-sum(p[j]*factors[j] for j in range(15)).log())
    return sum(terms, logits.sum()*0)


class ObjectiveChecks(unittest.TestCase):
    def test_local_probability_gradient_and_fine_priority(self):
        torch.manual_seed(18)
        x = torch.randn(8, 15, dtype=torch.float64, requires_grad=True)
        t = targets(8)
        t['seg'][:3] = torch.tensor([0, 5, 14])
        t['region'][:3] = torch.tensor([2, 0, 1])  # Deliberately conflicts; fine must win.
        t['region'][3:6] = torch.tensor([0, 1, 2])
        a = local_coarse_channel_nll(x, t)['loss']
        b = local_reference(x, t, .1)
        self.assertAlmostEqual(float(a), float(b), delta=1e-10)
        ga = torch.autograd.grad(a, x)[0]
        gb = torch.autograd.grad(b, x)[0]
        self.assertTrue(torch.allclose(ga, gb, atol=1e-10, rtol=0))
        self.assertEqual(float(ga[6:].abs().sum()), 0.)
        self.assertLess(float(ga.sum(1).abs().max()), 1e-10)
        shifted = x.detach()+torch.randn(8, 1, dtype=torch.float64)
        self.assertAlmostEqual(float(a), float(local_coarse_channel_nll(shifted, t)['loss']), delta=1e-10)
        result = local_coarse_channel_nll(x, t)
        self.assertEqual((result['fine_positions'], result['coarse_positions'], result['fully_unknown_positions']), (3, 3, 2))

    def test_unknown_has_zero_likelihood_loss_and_vjp(self):
        x = torch.randn(5, 15, dtype=torch.float64, requires_grad=True)
        loss = local_coarse_channel_nll(x, targets(5))['loss']
        self.assertEqual(float(loss), 0.)
        self.assertEqual(float(torch.autograd.grad(loss, x)[0].abs().sum()), 0.)

    def test_three_arms_independent_values_and_gradients(self):
        sequence = 'CCCATGAAATAACCCC'
        t = targets(10)
        t['seg'][0] = 13
        t['region'][4] = 2  # No compatible hard intron; channel must stay finite.
        torch.manual_seed(21)
        initial = torch.randn(len(sequence), 15, dtype=torch.float64)*.1
        for arm in C.ARMS:
            x = initial.clone().requires_grad_()
            actual, _ = C.objective(arm, sequence, x, 3, 10, t)
            structured = reference(x, [(3, 12)], 3, 13, t, torch.tensor(.1, dtype=torch.float64))
            local = local_reference(x[3:13], t, .1)
            expected = structured if arm == C.ARMS[0] else local if arm == C.ARMS[1] else .5*(structured+local)
            self.assertAlmostEqual(float(actual), float(expected), delta=1e-9)
            ga = torch.autograd.grad(actual, x)[0]
            gb = torch.autograd.grad(expected, x)[0]
            self.assertTrue(torch.allclose(ga, gb, atol=1e-9, rtol=0))
            if arm == 'LOCAL_CHANNEL':
                self.assertEqual(float(ga[:3].abs().sum()+ga[13:].abs().sum()), 0.)


if __name__ == '__main__':
    torch.set_num_threads(2)
    began = time.perf_counter()
    C.write(C.CHECKS/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], configuration=C.CONFIG))
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ObjectiveChecks))
    C.write(C.CHECKS/'summary.json', dict(completed=result.wasSuccessful(), job_id=os.environ['SLURM_JOB_ID'],
        configuration=C.CONFIG, tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        research_data_read=False, models_or_weights_read=False, optimizer_steps_this_run=0,
        seconds=time.perf_counter()-began))
    if not result.wasSuccessful():
        raise SystemExit(1)
