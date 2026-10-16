"""Artificial dense q/v casting, gradient, state and native-geometry checks."""
from types import SimpleNamespace
import argparse
import json
import os
from pathlib import Path
import time
import unittest
import torch
from src.m28 import joint_encoder_dense_qv as D


def small_projection():
    return torch.nn.Linear(5, 3, bias=False).to(torch.bfloat16)


class DenseQVTests(unittest.TestCase):
    def test_initial_native_output_exact_and_master_gradient(self):
        base = small_projection()
        original = base.weight.detach().clone()
        module = D.FullRankLinear(base)
        x = (torch.arange(10).reshape(2, 5).float()/10).to(torch.bfloat16).requires_grad_()
        self.assertTrue(torch.equal(module(x), base(x)))
        self.assertEqual(module.weight.dtype, torch.float32)
        module(x).float().square().sum().backward()
        self.assertGreater(float(module.weight.grad.norm()), 0)
        self.assertTrue(torch.isfinite(module.weight.grad).all())
        self.assertGreater(float(x.grad.float().norm()), 0)
        self.assertIsNone(base.weight.grad)
        self.assertTrue(torch.equal(base.weight, original))

    def test_changed_forward_and_fp32_adamw_state(self):
        base = small_projection()
        module = D.FullRankLinear(base)
        with torch.no_grad():
            module.weight.add_(.01)
        x = torch.ones(2, 5, dtype=torch.bfloat16)
        self.assertTrue(torch.equal(module(x), torch.nn.functional.linear(x, module.weight.to(x.dtype))))
        optimizer = torch.optim.AdamW([module.weight], lr=D.NATIVE_LR)
        before = module.weight.detach().clone()
        module(x).float().square().sum().backward()
        optimizer.step()
        self.assertFalse(torch.equal(before, module.weight))
        self.assertEqual(int(optimizer.state[module.weight]['step']), 1)
        self.assertEqual(optimizer.state[module.weight]['exp_avg'].dtype, torch.float32)

    def test_state_roundtrip_and_wrong_names_shapes_dtype_nonfinite(self):
        modules = {'fixture': D.FullRankLinear(small_projection())}
        state = D.dense_state(modules)
        with torch.no_grad():
            modules['fixture'].weight.zero_()
        D.restore_dense(modules, state)
        self.assertTrue(torch.equal(modules['fixture'].weight, state['fixture.weight']))
        for bad in ({}, {'fixture.weight': torch.zeros(3, 4)},
                    {'fixture.weight': state['fixture.weight'].double()},
                    {'fixture.weight': torch.full((3, 5), float('nan'))}):
            with self.assertRaises(ValueError):
                D.restore_dense(modules, bad)

    def test_only_eight_pinned_qv_weights_are_trainable(self):
        backbone = torch.nn.Module()
        backbone.config = SimpleNamespace(hidden_size=2048, num_attention_heads=32,
            num_key_value_heads=4, pretraining_tp=1, _attn_implementation='sdpa')
        backbone.layers = torch.nn.ModuleList([torch.nn.Module() for _ in range(26)])
        backbone.other = torch.nn.Linear(2, 2)
        originals = {}
        for index in D.LAYERS:
            attention = torch.nn.Module()
            attention.q_proj = torch.nn.Linear(2048, 2048, bias=False).to(torch.bfloat16)
            attention.v_proj = torch.nn.Linear(2048, 256, bias=False).to(torch.bfloat16)
            backbone.layers[index].self_attn = attention
            for name in D.PROJECTIONS:
                originals[f'layers.{index}.self_attn.{name}.weight'] = getattr(attention, name).weight.detach().float().clone()
        projections = D.attach_dense_qv(backbone)
        parameters = D.dense_parameters(projections)
        self.assertEqual(len(parameters), 8)
        self.assertEqual(sum(p.numel() for p in parameters), D.PARAMETERS)
        self.assertEqual({id(p) for p in parameters},
                         {id(p) for p in backbone.parameters() if p.requires_grad})
        self.assertTrue(all(p.dtype == torch.float32 for p in parameters))
        self.assertFalse(backbone.training)
        self.assertTrue(all(torch.equal(v, originals[k]) for k, v in D.dense_state(projections).items()))


if __name__ == '__main__':
    torch.set_num_threads(2)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    began = time.perf_counter()
    result = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(DenseQVTests))
    with (args.output_dir/'summary.json').open('x') as handle:
        json.dump(dict(completed=result.wasSuccessful(), job_id=os.environ['SLURM_JOB_ID'],
            tests_run=result.testsRun, errors=len(result.errors), failures=len(result.failures),
            research_DNA_read=False, research_reference_read=False, checkpoint_read=False,
            native_model_read=False, optimizer_steps_this_run=0, DEV_used=False,
            sealed_data_used=False, seconds=time.perf_counter()-began),
            handle, indent=2, allow_nan=False)
        handle.write('\n')
    raise SystemExit(0 if result.wasSuccessful() else 1)
