"""A25 checks that change readiness decisions, not a performance gate."""
import unittest
import torch
from src.m28.query_local_context import LocalContextGeneQueryModel, local_context_attention
from src.m28.query_start_weighted import (StartWeightedGeneQueryModel,
    start_weighted_attention, initiation_components)
from src.m28.query_decode import decode_query_chains
from src.m28.query_beam import decode_query_beam


def models():
    torch.manual_seed(0)
    old = LocalContextGeneQueryModel(8, 4, 8)
    torch.manual_seed(0)
    new = StartWeightedGeneQueryModel(8, 4, 8)
    return old, new


class StartWeightedTests(unittest.TestCase):
    def test_fresh_common_parameters_and_zero_prior_equal_A17(self):
        old, new = models()
        for key, value in old.state_dict().items():
            self.assertTrue(torch.equal(value, new.state_dict()[key]), key)
        self.assertEqual(sum(p.numel() for p in new.parameters())-sum(p.numel() for p in old.parameters()), 8)
        self.assertTrue((new.query.initiation.weight == 0).all())
        hidden = torch.randn(11, 8)
        a, b = old.query_window(hidden, 11), new.query_window(hidden, 11)
        self.assertTrue(torch.equal(a.initial, b.initial))
        self.assertTrue(torch.equal(a.query_logits, b.query_logits))
        for stage, lower in (('start', 0), ('end', 2), ('next_start', 5)):
            torch.testing.assert_close(a.pointer_logits(a.initial, stage, lower),
                                       b.pointer_logits(b.initial, stage, lower), rtol=0, atol=0)

    def test_nonzero_prior_reweights_only_context_and_short_contigs(self):
        _, model = models()
        for length in (1, 3, 11):
            hidden = torch.randn(length, 8)
            positions = model.query.coordinate(torch.arange(length)[:, None].float()/length)
            prior = torch.zeros(length)
            baseline = local_context_attention(model.query, hidden, positions)
            prior[0] = 4
            attention = start_weighted_attention(model.query, hidden, positions, prior)
            self.assertTrue(torch.isfinite(attention).all())
            torch.testing.assert_close(attention.sum(-1), torch.ones(4))
            self.assertTrue(torch.equal(attention == 0, baseline == 0))
            if length == 11:
                self.assertGreater(float(attention[0, 0]), float(baseline[0, 0]))
        with torch.no_grad():
            model.query.initiation.weight.fill_(.5)
        window = model.query_window(torch.randn(11, 8), 11)
        # Query0's local context does not crop its full-chain entry domain.
        self.assertTrue(torch.isfinite(window.pointer_logits(window.initial[0], 'start', 0)[9]))

    def test_auxiliary_masks_unknown_triplet_and_no_positive(self):
        logits = torch.zeros(8, requires_grad=True)
        targets = dict(seg=[0, 13, 13, 13, -100, -100, 1, -100],
                       positive_chains=[((1, 4),), ((4, 7),), ((1, 4),)])
        result = initiation_components(logits, targets)
        self.assertEqual(result['starts'], (1, 4))
        self.assertEqual(result['support'].tolist(), [True, True, False, False, True, False, True, False])
        result['loss'].backward()
        self.assertTrue(torch.isfinite(logits.grad).all())
        self.assertTrue((logits.grad[[2, 3, 5, 7]] == 0).all())
        self.assertLess(float(logits.grad[1]), 0.)
        other = torch.randn(8, requires_grad=True)
        zero = initiation_components(other, dict(seg=targets['seg'], positive_chains=[]))['loss']
        zero.backward()
        self.assertEqual(float(zero), 0.)
        self.assertTrue((other.grad == 0).all())

    def test_teacher_free_padding_old_components_and_new_gradient(self):
        old, model = models()
        hidden = torch.randn(15, 8, requires_grad=True)
        teacher = model.query_window(hidden, 11)
        free = model.query.for_window(hidden, 11)
        short = model.query_window(hidden[:11], 11)
        self.assertTrue(torch.equal(teacher.initial, free.initial))
        self.assertTrue(torch.equal(teacher.initial, short.initial))
        labels, chains = [1, -1, -1, -1], {0: ((1, 4), (7, 10))}
        parts = teacher.teacher_components(labels, chains)
        baseline = old.query_window(hidden, 11).teacher_components(labels, chains)
        for key in ('existence_bce', 'boundary_ce', 'termination_bce'):
            torch.testing.assert_close(parts[key], baseline[key], rtol=0, atol=0)
        aux = initiation_components(teacher.initiation_logits,
            dict(seg=[0]*15, positive_chains=list(chains.values())))['loss']
        (sum(parts[k] for k in ('existence_bce', 'boundary_ce', 'termination_bce'))+aux).backward()
        self.assertGreater(float(model.query.initiation.weight.grad.abs().sum()), 0.)
        self.assertTrue(torch.isfinite(hidden.grad).all())
        self.assertEqual(float(hidden.grad[11:].abs().sum()), 0.)

    def test_reference_free_original_greedy_and_beam(self):
        _, model = models()
        sequence = 'CCCATGAAATAACC'
        hidden = torch.randn(len(sequence), 8)
        with torch.no_grad():
            model.query.initiation.weight.fill_(.2)
        greedy = decode_query_chains(model.query, hidden, sequence, 8, len(sequence))
        beam = decode_query_beam(model.query, hidden, sequence, 8, 4, 4, len(sequence))
        self.assertEqual(len(greedy), 4)
        self.assertTrue(all(r['terminated'] and r['chain'] == ((3, 12),) for r in greedy))
        self.assertTrue(all(len(r['completed']) == 1 and r['completed'][0]['chain'] == ((3, 12),) for r in beam))


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
