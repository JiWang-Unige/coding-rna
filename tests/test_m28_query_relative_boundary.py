"""A23 initialization, actual distance dependence and original raw interfaces."""
import unittest
import torch
from src.m28.query_collaboration import CollaborativeGeneQueryModel, CollaborativeQueryWindow
from src.m28.query_relative_boundary import RelativeBoundaryGeneQueryModel
from src.m28.query_decode import decode_query_chains
from src.m28.query_beam import decode_query_beam


def models():
    torch.manual_seed(0)
    old = CollaborativeGeneQueryModel(8, 4, 8)
    torch.manual_seed(0)
    new = RelativeBoundaryGeneQueryModel(8, 4, 8)
    return old, new


class RelativeBoundaryTests(unittest.TestCase):
    def test_common_parameters_and_zero_delta_scores_match_A3(self):
        old, new = models()
        for key, value in old.state_dict().items():
            self.assertTrue(torch.equal(value, new.state_dict()[key]), key)
        self.assertEqual(sum(p.numel() for p in new.parameters())-sum(p.numel() for p in old.parameters()), 160)
        hidden = torch.randn(14, 8)
        a, b = old.query_window(hidden, 14), new.query_window(hidden, 14)
        self.assertTrue(torch.equal(a.initial, b.initial))
        self.assertTrue(torch.equal(a.query_logits, b.query_logits))
        for stage, lower in (('start', 0), ('end', 2), ('next_start', 5)):
            for index in (0, None):
                sa = a.initial if index is None else a.initial[index]
                sb = b.initial if index is None else b.initial[index]
                torch.testing.assert_close(a.pointer_logits(sa, stage, lower), b.pointer_logits(sb, stage, lower), rtol=0, atol=0)

    def test_bias_depends_on_lower_without_changing_first_start(self):
        _, model = models()
        with torch.no_grad():
            model.query.relative_boundary[0].weight.zero_()
            model.query.relative_boundary[0].bias.zero_()
            model.query.relative_boundary[0].weight[0, 0] = 1
            model.query.relative_boundary[-1].weight.zero_()
            model.query.relative_boundary[-1].weight[:, 0] = 1
        w = model.query_window(torch.randn(14, 8), 14)
        state = w.initial[0]
        delta = []
        for lower in (2, 5):
            plain = CollaborativeQueryWindow.pointer_logits(w, state, 'end', lower)
            delta.append(float((w.pointer_logits(state, 'end', lower)-plain)[9]))
        self.assertGreater(delta[0], delta[1])
        torch.testing.assert_close(w.pointer_logits(state, 'start', 0),
                                   CollaborativeQueryWindow.pointer_logits(w, state, 'start', 0), rtol=0, atol=0)
        self.assertTrue(torch.isneginf(w.pointer_logits(state, 'end', 5)[:6]).all())

    def test_padding_teacher_free_and_real_teacher_gradients(self):
        old, new = models()
        hidden = torch.randn(20, 8, requires_grad=True)
        teacher, free = new.query_window(hidden, 14), new.query.for_window(hidden, 14)
        torch.testing.assert_close(teacher.initial, free.initial, rtol=0, atol=0)
        short = new.query_window(hidden[:14], 14)
        torch.testing.assert_close(teacher.pointer_logits(teacher.initial, 'end', 2),
                                   short.pointer_logits(short.initial, 'end', 2), rtol=0, atol=0)
        chain = {0: ((1, 4), (9, 12))}
        parts = teacher.teacher_components([1, -1, -1, -1], chain)
        prior = old.query_window(hidden, 14).teacher_components([1, -1, -1, -1], chain)
        for key in ('existence_bce', 'boundary_ce', 'termination_bce'):
            torch.testing.assert_close(parts[key], prior[key], rtol=0, atol=0)
        sum(parts[k] for k in ('existence_bce', 'boundary_ce', 'termination_bce')).backward()
        final = new.query.relative_boundary[-1].weight.grad
        self.assertTrue(torch.isfinite(final).all())
        self.assertTrue((final.abs().sum(dim=1) > 0).all())
        self.assertTrue((new.query.relative_boundary[0].weight.grad == 0).all())
        self.assertGreater(float(hidden.grad[:14].abs().sum()), 0.)
        self.assertEqual(float(hidden.grad[14:].abs().sum()), 0.)

    def test_original_greedy_and_beam_use_the_new_pointer_without_reference(self):
        _, model = models()
        sequence = 'CCCATGAAATAACC'
        hidden = torch.randn(len(sequence), 8)
        with torch.no_grad():
            model.query.relative_boundary[-1].weight.fill_(.2)
        greedy = decode_query_chains(model.query, hidden, sequence, 8, len(sequence))
        beam = decode_query_beam(model.query, hidden, sequence, 8, 4, 4, len(sequence))
        self.assertEqual(len(greedy), 4)
        self.assertTrue(all(r['terminated'] and r['chain'] == ((3, 12),) for r in greedy))
        self.assertTrue(all(len(r['completed']) == 1 and r['completed'][0]['chain'] == ((3, 12),) for r in beam))


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
