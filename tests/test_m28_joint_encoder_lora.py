"""Artificial adapter arithmetic, state identity, and joint-window gradient path."""
from types import SimpleNamespace
import argparse
import json
import os
from pathlib import Path
import time
import unittest
import torch
from src.m28 import joint_encoder_lora as L


def small_adapter():
    base = torch.nn.Linear(5, 3, bias=False).to(torch.bfloat16)
    return L.LowRankLinear(base, torch.Generator().manual_seed(0))


class AdapterTests(unittest.TestCase):
    def test_zero_nesting_and_first_gradient_without_input_gradient(self):
        module = small_adapter()
        x = torch.arange(10).reshape(2, 5).to(torch.bfloat16)/10
        original = module.base(x)
        actual = module(x)
        self.assertTrue(torch.equal(original, actual))
        self.assertEqual(actual.dtype, torch.bfloat16)
        actual.float().square().sum().backward()
        self.assertIsNone(module.base.weight.grad)
        self.assertFalse(module.base.weight.requires_grad)
        self.assertTrue(torch.equal(module.A.grad, torch.zeros_like(module.A)))
        self.assertGreater(float(module.B.grad.norm()), 0)
        self.assertTrue(torch.isfinite(module.B.grad).all())

    def test_nonzero_formula_and_upstream_gradient(self):
        module = small_adapter()
        with torch.no_grad():
            module.B.copy_(torch.arange(12).reshape(3, 4).float()/100)
        x = (torch.arange(10).reshape(2, 5).float()/10).to(torch.bfloat16).requires_grad_()
        expected = module.base(x) + ((x.float() @ module.A.T) @ module.B.T).to(torch.bfloat16)
        self.assertTrue(torch.equal(module(x), expected))
        module(x).float().square().sum().backward()
        self.assertGreater(float(module.A.grad.norm()), 0)
        self.assertGreater(float(module.B.grad.norm()), 0)
        self.assertGreater(float(x.grad.float().norm()), 0)
        self.assertIsNone(module.base.weight.grad)

    def test_state_roundtrip_and_mismatches(self):
        modules = {'fixture': small_adapter()}
        with torch.no_grad():
            modules['fixture'].B.fill_(.1)
        state = L.adapter_state(modules)
        with torch.no_grad():
            modules['fixture'].A.zero_(); modules['fixture'].B.zero_()
        L.restore_adapters(modules, state)
        self.assertTrue(all(torch.equal(v, state[k]) for k, v in L.adapter_state(modules).items()))
        for bad in ({'fixture.A': state['fixture.A']},
                    dict(state, **{'fixture.B': torch.zeros(3, 3)}),
                    dict(state, **{'fixture.B': state['fixture.B'].double()}),
                    dict(state, **{'fixture.B': torch.full((3, 4), float('nan'))})):
            with self.assertRaises(ValueError):
                L.restore_adapters(modules, bad)

    def test_only_fixed_native_projections_are_trainable(self):
        backbone = torch.nn.Module()
        backbone.config = SimpleNamespace(hidden_size=2048, num_attention_heads=32,
            num_key_value_heads=4, pretraining_tp=1, _attn_implementation='sdpa')
        backbone.layers = torch.nn.ModuleList([torch.nn.Module() for _ in range(26)])
        backbone.other = torch.nn.Linear(2, 2)
        for index in L.LAYERS:
            attention = torch.nn.Module()
            attention.q_proj = torch.nn.Linear(2048, 2048, bias=False).to(torch.bfloat16)
            attention.v_proj = torch.nn.Linear(2048, 256, bias=False).to(torch.bfloat16)
            backbone.layers[index].self_attn = attention
        adapters = L.attach_adapters(backbone)
        parameters = L.adapter_parameters(adapters)
        self.assertEqual(len(adapters), 8)
        self.assertEqual(sum(p.numel() for p in parameters), 102400)
        self.assertEqual({id(p) for p in parameters}, {id(p) for p in backbone.parameters() if p.requires_grad})
        self.assertFalse(backbone.training)
        self.assertTrue(all(p.dtype == torch.float32 for p in parameters))
        initial = L.adapter_state(adapters)
        self.assertTrue(all(torch.count_nonzero(v) == 0 for k, v in initial.items() if k.endswith('.B')))
        with self.assertRaises(ValueError):
            L.attach_adapters(backbone)

    def test_joint_features_retain_gradient_padding_and_literal_tokenization(self):
        class Tokenizer:
            def __call__(self, sequence, **kwargs):
                self.sequence, self.kwargs = sequence, kwargs
                return {'input_ids': torch.zeros(1, len(sequence)//6, dtype=torch.long)}
        class Backbone(torch.nn.Module):
            def __init__(self):
                super().__init__(); self.scale = torch.nn.Parameter(torch.tensor(2.))
            def forward(self, input_ids, attention_mask):
                self.ids, self.mask = input_ids, attention_mask
                values = self.scale.to(torch.bfloat16).expand(1, 4096, 2048)
                return SimpleNamespace(last_hidden_state=values)
        tokenizer, backbone = Tokenizer(), Backbone()
        values = L.encode_joint_window(backbone, tokenizer, 'ATGCGTNAGTTTAA', 'cpu')
        self.assertEqual(tokenizer.sequence, 'ATGCGTNAGTTTAA'+'N'*24562)
        self.assertEqual(tokenizer.kwargs, dict(add_special_tokens=False, return_tensors='pt', truncation=False))
        self.assertEqual(values.dtype, torch.float32)
        self.assertTrue(torch.equal(backbone.mask, torch.ones_like(backbone.ids)))
        values.mean().backward()
        self.assertEqual(float(backbone.scale.grad), 1.)
        with self.assertRaises(ValueError):
            L.encode_joint_window(backbone, tokenizer, '', 'cpu')


if __name__ == '__main__':
    torch.set_num_threads(2)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args()
    began = time.perf_counter()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(AdapterTests))
    if args.output_dir:
        with (args.output_dir/'summary.json').open('x') as handle:
            json.dump(dict(completed=result.wasSuccessful(), job_id=os.environ['SLURM_JOB_ID'],
                tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
                research_DNA_read=False, research_reference_read=False, checkpoint_read=False,
                native_model_read=False, cached_features_used=False, DEV_used=False,
                sealed_data_used=False, optimizer_steps_this_run=0,
                seconds=time.perf_counter()-began), handle, indent=2, allow_nan=False)
            handle.write('\n')
    raise SystemExit(0 if result.wasSuccessful() else 1)
