import copy
import unittest
import torch
from src.m28.chain_set import chain_set_loss
from src.m28.core import select_nonoverlapping
from src.m28.natural_chain_set import natural_record


def row(chain, label, origin='free', gain=1.):
    return dict(chain=chain, label=label, origin=origin, gain=gain)


def score(local):
    return dict(logits=torch.tensor([r['gain'] for r in local], dtype=torch.float64,
                                   requires_grad=True),
                null=torch.tensor(0., dtype=torch.float64, requires_grad=True))


class NaturalChainSetTests(unittest.TestCase):
    def test_projection_preserves_free_order_geometry_and_all_original_labels(self):
        rec = dict(window={'id': 'train'}, counts={'old': 4}, candidates=[
            row([[0, 9]], 1), row([[3, 12]], 0),
            row([[15, 24]], -1), row([[30, 39]], 1, 'reference_injected')])
        original = copy.deepcopy(rec); natural = natural_record(rec)
        self.assertEqual(rec, original)
        self.assertEqual(natural['candidates'], original['candidates'][:3])
        self.assertEqual(natural['counts'], dict(windows=1, free_candidates=3,
            reference_injected_candidates=0, zero_known_windows=0,
            **{'label_1': 1, 'label_0': 1, 'label_-1': 1}))

    def test_losing_injected_gold_is_empty_gold_not_a_relabelled_candidate(self):
        rec = dict(candidates=[row([[0, 18]], 1, 'reference_injected', 20.),
                               row([[0, 9]], 0), row([[9, 18]], 0)])
        local = natural_record(rec)['candidates']; s = score(local)
        result = chain_set_loss(s, [r['chain'] for r in local], [r['label'] for r in local])
        self.assertEqual(result['known_positives'], 0)
        self.assertEqual(float(result['gold_gain']), 0.)
        result['loss'].backward()
        self.assertTrue(bool((s['logits'].grad > 0).all()))

    def test_empty_or_unknown_natural_domain_keeps_differentiable_zero(self):
        for rows in ([row([[0, 9]], 1, 'reference_injected')],
                     [row([[0, 9]], 1, 'reference_injected'), row([[0, 18]], -1, gain=100.)]):
            natural = natural_record(dict(candidates=rows)); local = natural['candidates']
            s = score(local)
            result = chain_set_loss(s, [r['chain'] for r in local], [r['label'] for r in local])
            result['loss'].backward()
            self.assertEqual(natural['counts']['windows'], 1)
            self.assertEqual(natural['counts']['zero_known_windows'], 1)
            self.assertEqual(float(result['loss']), 0.)
            self.assertTrue(torch.equal(s['logits'].grad, torch.zeros_like(s['logits'])))

    def test_selection_preserves_unknown_and_never_uses_injected_gold(self):
        rec = dict(candidates=[row([[0, 9]], 0, gain=2.), row([[9, 18]], -1, gain=3.),
                               row([[0, 18]], 1, 'reference_injected', 100.)])
        local = natural_record(rec)['candidates']
        chosen = select_nonoverlapping([r['chain'] for r in local], [r['gain'] for r in local])
        self.assertEqual([local[i]['label'] for i in chosen], [0, -1])


if __name__ == '__main__':
    unittest.main()
