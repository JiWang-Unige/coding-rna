import unittest
import torch
from src.m28.query_collaboration import CollaborativeGeneQueryModel
from src.m28.query_local_context import (LocalContextGeneQueryModel, LocalContextQueryPointer,
                                        context_intervals, local_context_attention)


class LocalContextTests(unittest.TestCase):
    def test_exact_A3_initial_parameters_without_extra_parameters(self):
        torch.manual_seed(0); original = CollaborativeGeneQueryModel(8, 4, 8)
        torch.manual_seed(0); local = LocalContextGeneQueryModel(8, 4, 8)
        self.assertEqual(original.state_dict().keys(), local.state_dict().keys())
        for name, value in original.state_dict().items():
            self.assertTrue(torch.equal(value, local.state_dict()[name]), name)

    def test_context_support_uses_valid_geometry_and_handles_short_contigs(self):
        pointer = LocalContextQueryPointer(4, 8)
        for length in (1, 3, 8, 11):
            features = torch.randn(length, 8)
            positions = pointer.coordinate(torch.arange(length)[:, None].float()/length)
            attention = local_context_attention(pointer, features, positions)
            self.assertTrue(torch.isfinite(attention).all())
            self.assertTrue(torch.allclose(attention.sum(-1), torch.ones(4)))
            for i, (a, b) in enumerate(context_intervals(length, 4)):
                self.assertTrue((attention[i, :a] == 0).all())
                self.assertTrue((attention[i, b:] == 0).all())
                self.assertTrue((attention[i, a:b] > 0).all())
            if length >= 4:
                self.assertEqual(sum(b-a for a, b in context_intervals(length, 4)), length)

    def test_shared_teacher_and_free_state_global_pointers_and_gradient(self):
        torch.manual_seed(0); model = LocalContextGeneQueryModel(8, 4, 8)
        hidden = torch.randn(12, 8, requires_grad=True)
        teacher = model.query_window(hidden, 8)
        free = model.query.for_window(hidden, 8)
        self.assertTrue(torch.equal(teacher.initial, free.initial))
        self.assertTrue(torch.equal(teacher.query_logits, free.query_logits))
        self.assertEqual(len(teacher.features), 8)
        # Query 0 pools bases 0..2 but may predict an entire exon elsewhere.
        self.assertTrue(torch.isfinite(teacher.pointer_logits(teacher.initial[0], 'start', 0)[6]))
        result = teacher.teacher_components([1, -1, -1, -1], {0: ((5, 8),)})
        sum(result[k] for k in ('existence_bce', 'boundary_ce', 'termination_bce')).backward()
        self.assertGreater(float(hidden.grad[:8].abs().sum()), 0.)
        self.assertEqual(float(hidden.grad[8:].abs().sum()), 0.)
        for name in ('context_key', 'context_value', 'query_attention', 'start_key', 'end_key', 'recurrent'):
            gradients = [p.grad for p in getattr(model.query, name).parameters()]
            self.assertTrue(all(g is not None and torch.isfinite(g).all() for g in gradients), name)
            self.assertGreater(sum(float(g.abs().sum()) for g in gradients), 0., name)


if __name__ == '__main__':
    torch.set_num_threads(2); unittest.main()
