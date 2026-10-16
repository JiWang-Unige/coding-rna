"""Artificial strict native_qv endpoints and optimizer registration; no research data."""
import copy
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import torch
import annotate_full_grammar_dense_qv as A
import fit_full_grammar_dense_qv as F
import infer_full_grammar_dense_qv_train as D
from compare_full_grammar_dense_qv_train import paired_scopes, SPECIES
from test_m28_full_grammar_codon_train import paired_fixture as codon_fixture


def checkpoint_fixture():
    return dict(configuration=dict(F.CONFIG), step=9216, epoch=6, native_qv_optimizer_step=4608,
        protocol=A.PROTOCOL, encoder_revision=A.R.F.REVISION,
        model=A.R.M28Core('C0', 2048, 128).state_dict(),
        qv_weights={k:torch.full(shape, .001, dtype=torch.float32) for k, shape in A.native_qv_shapes().items()})


def paired_fixture():
    value = codon_fixture()
    return {'adapted': value['codon'], 'baseline': value['baseline']}


class AdapterTrainTests(unittest.TestCase):
    def test_checkpoint_final_ages_identity_and_all_named_tensors(self):
        state = checkpoint_fixture()
        A.model_from_state(state, 'cpu')
        for key, value in (('step', 7680), ('epoch', 5), ('native_qv_optimizer_step', 9216),
                           ('protocol', A.B.PROTOCOL)):
            bad = dict(state, **{key:value})
            with self.assertRaises(ValueError):
                A.validate_checkpoint(bad)
        for key, value in (('native_qv_parameters', 64), ('native_qv_lr', .0003), ('encoder_frozen', True),
                           ('sampling_seed', 2), ('encoder_blocks_bp', 6144)):
            bad = dict(state, configuration=dict(state['configuration'], **{key:value}))
            with self.assertRaises(ValueError):
                A.validate_checkpoint(bad)
        name = next(iter(state['qv_weights']))
        for changed in ({k:v for k, v in state['qv_weights'].items() if k != name},
                        dict(state['qv_weights'], **{name:torch.zeros(2048, 2047)}),
                        dict(state['qv_weights'], **{name:state['qv_weights'][name].double()}),
                        dict(state['qv_weights'], **{name:torch.full(state['qv_weights'][name].shape, float('nan'))})):
            with self.assertRaises(ValueError):
                A.validate_checkpoint(dict(state, qv_weights=changed))

    def test_native_loader_restores_actual_values_without_unfreezing_base(self):
        state = checkpoint_fixture()
        backbone = torch.nn.Module()
        backbone.config = SimpleNamespace(hidden_size=2048, num_attention_heads=32,
            num_key_value_heads=4, pretraining_tp=1, _attn_implementation='sdpa')
        backbone.layers = torch.nn.ModuleList([torch.nn.Module() for _ in range(26)])
        for index in A.D.LAYERS:
            attention = torch.nn.Module()
            attention.q_proj = torch.nn.Linear(2048, 2048, bias=False).to(torch.bfloat16)
            attention.v_proj = torch.nn.Linear(2048, 256, bias=False).to(torch.bfloat16)
            backbone.layers[index].self_attn = attention
        original = backbone.layers[22].self_attn.q_proj.weight.detach().clone()
        tokenizer = object()
        with patch.object(A.R, 'load_native', return_value=(backbone, tokenizer)):
            restored, actual_tokenizer = A.load_native_adapted(state)
        self.assertIs(actual_tokenizer, tokenizer)
        self.assertIs(restored, backbone)
        self.assertTrue(all(not p.requires_grad for p in restored.parameters()))
        self.assertTrue(all(not module.training for module in restored.modules()))
        for name, value in state['qv_weights'].items():
            self.assertTrue(torch.equal(dict(restored.named_parameters())[name], value))
        module = restored.layers[22].self_attn.q_proj
        x = torch.ones(1, 2048, dtype=torch.bfloat16)
        expected = torch.nn.functional.linear(x,
            state['qv_weights']['layers.22.self_attn.q_proj.weight'].to(torch.bfloat16))
        self.assertTrue(torch.equal(module(x), expected))
        self.assertFalse(torch.equal(module.weight.detach().to(torch.bfloat16), original))

    def test_adamw_named_registration_ages_shapes_and_hyperparameters(self):
        model = A.R.M28Core('C0', 2048, 128)
        weights = {k:torch.zeros(shape) for k, shape in A.native_qv_shapes().items()}
        parameters = list(model.parameters()) + list(weights.values())
        n = len(list(model.parameters()))
        group = {k:F.CONFIG[k] for k in ('lr', 'weight_decay', 'eps')}
        group['betas'] = tuple(F.CONFIG['betas'])
        saved = dict(param_groups=[dict(group, params=list(range(n))),
            dict(group, lr=A.D.NATIVE_LR, params=list(range(n, len(parameters))), param_names=list(A.native_qv_shapes()))],
            state={i:dict(step=torch.tensor(9216 if i < n else 4608), exp_avg=torch.zeros_like(parameter),
                exp_avg_sq=torch.zeros_like(parameter)) for i, parameter in enumerate(parameters)})
        F.verify_optimizer(saved, model, weights, 9216, 4608)
        for target, change in ((0, 'head_age'), (n, 'native_qv_age'), (n, 'moment_shape'),
                               (0, 'lr'), (n, 'native_lr'), (n, 'names')):
            bad = copy.deepcopy(saved)
            if change.endswith('age'):
                bad['state'][target]['step'] += 1
            elif change == 'moment_shape':
                bad['state'][target]['exp_avg'] = torch.zeros(1)
            elif change in ('lr', 'native_lr'):
                bad['param_groups'][0 if change == 'lr' else 1]['lr'] = .001
            else:
                bad['param_groups'][1]['param_names'].reverse()
            with self.assertRaises(ValueError):
                F.verify_optimizer(bad, model, weights, 9216, 4608)

    def test_raw_report_rejects_frozen_model_or_missing_native_qv_restore(self):
        good = dict(encoder_condition=A.CONDITION, native_qv_tensors_restored=8,
            native_qv_parameters=18874368, encoder_parameters_requires_grad=0,
            score_kind='unit_fifteen_role_logit_minus_background_joint_oriented_contig_gain',
            scopes=[dict(encoder_condition=A.CONDITION) for _ in range(2)])
        D.validate_raw_report(good)
        for key, value in (('encoder_condition', 'JOINT24576'), ('native_qv_tensors_restored', 0),
                           ('native_qv_parameters', 64), ('encoder_parameters_requires_grad', 18874368),
                           ('score_kind', 'codon_energy')):
            with self.assertRaises(ValueError):
                D.validate_raw_report(dict(good, **{key:value}))
        changed = copy.deepcopy(good); changed['scopes'][1]['encoder_condition'] = 'JOINT24576'
        with self.assertRaises(ValueError):
            D.validate_raw_report(changed)

    def test_paired_counts_macro_and_worst_use_each_arm_species(self):
        scopes = paired_scopes(paired_fixture())
        self.assertEqual(scopes['pooled']['arms']['adapted']['exact_chain']['reference'], 22688)
        self.assertEqual(scopes['macro_species']['adapted_minus_baseline_PRF_percentage_points']['exact_chain']['f1'], 0.)
        self.assertEqual(scopes['worst_species']['adapted_minus_baseline_PRF_percentage_points']['exact_chain']['f1'], 12.5)
        self.assertEqual(min(scopes[s]['adapted_minus_baseline']['exact_chain']['PRF_percentage_points']['f1']
            for s in SPECIES), -12.5)

    def test_changed_paired_reference_or_aggregate_counts_are_rejected(self):
        reports = paired_fixture()
        changed = copy.deepcopy(reports)
        changed['adapted']['chromosomes'][D.I.SCOPE[0][1]]['reference_primary_records'] -= 1
        with self.assertRaises(ValueError):
            paired_scopes(changed)
        changed = copy.deepcopy(reports)
        changed['adapted']['pooled']['exact_chain']['tp'] += 1
        with self.assertRaises(ValueError):
            paired_scopes(changed)


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()

