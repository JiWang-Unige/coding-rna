"""Only changed learning-rate routing and checkpoint identities; artificial data."""
import copy
import unittest
import torch
import annotate_full_grammar_joint_lora_lr_low as LA
import fit_full_grammar_joint_lora_lr_low as LF
import annotate_full_grammar_dense_qv_lr_high as DA
import fit_full_grammar_dense_qv_lr_high as DF
from test_m28_full_grammar_joint_lora_train import checkpoint_fixture as low_fixture
from test_m28_full_grammar_dense_qv_train import checkpoint_fixture as dense_fixture


def checkpoint_case(case, annotation, fitted, fixture, key, rate, old_rate):
    state = fixture()
    old = dict(state)
    state = dict(state, configuration=dict(fitted.CONFIG), protocol=annotation.PROTOCOL)
    annotation.model_from_state(state, 'cpu')
    with case.assertRaises(ValueError):
        annotation.validate_checkpoint(old)
    for configuration in (dict(state['configuration'], **{key:old_rate}),
                          {k:v for k,v in state['configuration'].items() if k != key},
                          dict(state['configuration'], lr=1e-5)):
        with case.assertRaises(ValueError):
            annotation.validate_checkpoint(dict(state, configuration=configuration))
    case.assertEqual(state['configuration'][key], rate)


def optimizer_case(case, annotation, fitted, shapes, key, rate, old_rate):
    model = torch.nn.Linear(3, 2)
    names = list(shapes())
    native = {name:torch.nn.Parameter(torch.ones(2, 2)) for name in names}
    optimizer = torch.optim.AdamW(model.parameters(), lr=fitted.CONFIG['lr'],
        weight_decay=fitted.CONFIG['weight_decay'], betas=tuple(fitted.CONFIG['betas']),
        eps=fitted.CONFIG['eps'])
    for parameter in model.parameters():
        optimizer.state[parameter] = dict(step=torch.tensor(4608.),
            exp_avg=torch.zeros_like(parameter), exp_avg_sq=torch.zeros_like(parameter))
    fitted.add_native_group(optimizer, list(native.values()))
    for parameter in [*model.parameters(), *native.values()]:
        parameter.grad = torch.ones_like(parameter)
    optimizer.step()
    fitted.verify_optimizer(optimizer.state_dict(), model, native, 4609, 1)
    case.assertEqual([g['lr'] for g in optimizer.param_groups], [3e-4, rate])
    saved = optimizer.state_dict()
    for group_index, changed_rate in ((0, 1e-5), (1, old_rate)):
        wrong = copy.deepcopy(saved)
        wrong['param_groups'][group_index]['lr'] = changed_rate
        with case.assertRaises(ValueError):
            fitted.verify_optimizer(wrong, model, native, 4609, 1)
    wrong = copy.deepcopy(saved)
    native_id = wrong['param_groups'][1]['params'][0]
    wrong['state'][native_id]['step'] += 1
    with case.assertRaises(ValueError):
        fitted.verify_optimizer(wrong, model, native, 4609, 1)
    case.assertEqual(fitted.CONFIG[key], rate)


class LowLearningRateTests(unittest.TestCase):
    def test_checkpoint_rejects_original_or_missing_native_lr(self):
        checkpoint_case(self, LA, LF, low_fixture, 'native_adapter_lr', 1e-5, 3e-4)

    def test_actual_adamw_groups_and_wrong_rate_or_age_rejection(self):
        optimizer_case(self, LA, LF, LA.adapter_shapes, 'native_adapter_lr', 1e-5, 3e-4)


class DenseLearningRateTests(unittest.TestCase):
    def test_checkpoint_rejects_original_or_missing_native_lr(self):
        checkpoint_case(self, DA, DF, dense_fixture, 'native_qv_lr', 3e-4, 1e-5)

    def test_actual_adamw_groups_and_wrong_rate_or_age_rejection(self):
        optimizer_case(self, DA, DF, DA.native_qv_shapes, 'native_qv_lr', 3e-4, 1e-5)
