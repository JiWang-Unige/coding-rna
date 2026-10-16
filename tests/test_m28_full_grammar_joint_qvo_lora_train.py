"""Four new q/v/o wiring/restore/registration checks; no research data or weights."""
import copy
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import torch
from src.m28 import joint_encoder_lora as QV
from src.m28 import joint_encoder_qvo_lora as L
import annotate_full_grammar_joint_qvo_lora as A
import fit_full_grammar_joint_qvo_lora as F
import infer_full_grammar_joint_qvo_lora_train as D


def backbone_fixture():
    backbone = torch.nn.Module()
    backbone.config = SimpleNamespace(hidden_size=2048, num_attention_heads=32,
        num_key_value_heads=4, pretraining_tp=1, _attn_implementation='sdpa')
    backbone.layers = torch.nn.ModuleList([torch.nn.Module() for _ in range(26)])
    backbone.other = torch.nn.Linear(2, 2)
    for index in L.LAYERS:
        attention = torch.nn.Module()
        attention.q_proj = torch.nn.Linear(2048, 2048, bias=False).to(torch.bfloat16)
        attention.v_proj = torch.nn.Linear(2048, 256, bias=False).to(torch.bfloat16)
        attention.o_proj = torch.nn.Linear(2048, 2048, bias=False).to(torch.bfloat16)
        backbone.layers[index].self_attn = attention
    return backbone


def checkpoint_fixture():
    return dict(configuration=dict(F.CONFIG), step=9216, epoch=6, adapter_optimizer_step=4608,
        protocol=A.PROTOCOL, encoder_revision=A.R.F.REVISION,
        model=A.R.M28Core('C0', 2048, 128).state_dict(),
        adapters={k:torch.full(shape, .001, dtype=torch.float32) for k, shape in A.adapter_shapes().items()})


class QVOTests(unittest.TestCase):
    def test_new_target_nesting_preserves_original_qv_initial_state(self):
        parent, extended = backbone_fixture(), backbone_fixture()
        qv = QV.attach_adapters(parent)
        modules = L.attach_adapters(extended)
        state = L.adapter_state(modules)
        self.assertEqual(len(modules), 12)
        self.assertEqual(len(state), 24)
        self.assertEqual(sum(p.numel() for p in L.adapter_parameters(modules)), 167936)
        self.assertEqual(list(state), list(A.adapter_shapes()))
        self.assertEqual({id(p) for p in L.adapter_parameters(modules)},
            {id(p) for p in extended.parameters() if p.requires_grad})
        self.assertTrue(all(torch.equal(state[k], v) for k, v in QV.adapter_state(qv).items()))
        x = torch.ones(2, 2048, dtype=torch.bfloat16)
        for index in L.LAYERS:
            module = modules[f'layers.{index}.self_attn.o_proj']
            self.assertTrue(torch.equal(module(x), module.base(x)))
            self.assertFalse(module.base.weight.requires_grad)
        self.assertFalse(extended.other.weight.requires_grad)
        bad = backbone_fixture()
        bad.layers[24].self_attn.o_proj = torch.nn.Linear(2048, 256, bias=False).to(torch.bfloat16)
        with self.assertRaises(ValueError):
            L.attach_adapters(bad)
        self.assertIsInstance(bad.layers[22].self_attn.q_proj, torch.nn.Linear)

    def test_all_four_output_adapters_receive_gradients_and_restore_named_values(self):
        backbone = backbone_fixture()
        modules = L.attach_adapters(backbone)
        x = torch.ones(2, 2048, dtype=torch.bfloat16)
        for index in L.LAYERS:
            module = modules[f'layers.{index}.self_attn.o_proj']
            module(x).float().square().mean().backward()
            self.assertIsNone(module.base.weight.grad)
            self.assertEqual(float(module.A.grad.norm()), 0.)
            self.assertGreater(float(module.B.grad.norm()), 0.)
            with torch.no_grad():
                module.B.fill_(.01)
        state = L.adapter_state(modules)
        for module in modules.values():
            with torch.no_grad():
                module.A.zero_(); module.B.zero_()
        L.restore_adapters(modules, state)
        self.assertTrue(all(torch.equal(v, state[k]) for k, v in L.adapter_state(modules).items()))
        o_name = 'layers.22.self_attn.o_proj.B'
        for bad in ({k:v for k, v in state.items() if 'o_proj' not in k},
                    dict(state, **{o_name:torch.zeros(256, 4)})):
            with self.assertRaises(ValueError):
                L.restore_adapters(modules, bad)

    def test_new_checkpoint_and_actual_loader_reject_old_sixteen_tensor_state(self):
        state = checkpoint_fixture()
        A.model_from_state(state, 'cpu')
        for key, value in (('step', 7680), ('adapter_optimizer_step', 9216)):
            with self.assertRaises(ValueError):
                A.validate_checkpoint(dict(state, **{key:value}))
        bad = dict(state, adapters={k:v for k, v in state['adapters'].items() if 'o_proj' not in k})
        with self.assertRaises(ValueError):
            A.validate_checkpoint(bad)
        with self.assertRaises(ValueError):
            A.validate_checkpoint(dict(state, configuration=dict(state['configuration'], adapter_parameters=102400)))
        with patch.object(A.R, 'load_native', return_value=(backbone_fixture(), object())):
            backbone, _ = A.load_native_adapted(state)
        self.assertTrue(all(not p.requires_grad for p in backbone.parameters()))
        self.assertTrue(all(torch.equal(dict(backbone.named_parameters())[k], v)
            for k, v in state['adapters'].items()))
        good = dict(encoder_condition=A.CONDITION, native_adapter_tensors_restored=24,
            native_adapter_parameters=167936, encoder_parameters_requires_grad=0,
            score_kind='unit_fifteen_role_logit_minus_background_joint_oriented_contig_gain',
            scopes=[dict(encoder_condition=A.CONDITION) for _ in range(2)])
        D.validate_raw_report(good)
        with self.assertRaises(ValueError):
            D.validate_raw_report(dict(good, native_adapter_tensors_restored=16))

    def test_twenty_four_optimizer_parameters_have_correct_names_shapes_and_ages(self):
        model = A.R.M28Core('C0', 2048, 128)
        adapters = {k:torch.zeros(shape) for k, shape in A.adapter_shapes().items()}
        parameters = list(model.parameters()) + list(adapters.values())
        n = len(list(model.parameters()))
        group = {k:F.CONFIG[k] for k in ('lr', 'weight_decay', 'eps')}
        group['betas'] = tuple(F.CONFIG['betas'])
        saved = dict(param_groups=[dict(group, params=list(range(n))),
            dict(group, params=list(range(n, len(parameters))), param_names=list(A.adapter_shapes()))],
            state={i:dict(step=torch.tensor(9216 if i<n else 4608), exp_avg=torch.zeros_like(p),
                exp_avg_sq=torch.zeros_like(p)) for i, p in enumerate(parameters)})
        F.verify_optimizer(saved, model, adapters, 9216, 4608)
        for change in ('names', 'output_age', 'native_lr'):
            bad = copy.deepcopy(saved)
            if change == 'names':
                bad['param_groups'][1]['param_names'] = [k for k in A.adapter_shapes() if 'o_proj' not in k]
            elif change == 'output_age':
                output = n + list(A.adapter_shapes()).index('layers.22.self_attn.o_proj.B')
                bad['state'][output]['step'] += 1
            else:
                bad['param_groups'][1]['lr'] = 1e-5
            with self.assertRaises(ValueError):
                F.verify_optimizer(bad, model, adapters, 9216, 4608)


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
