import unittest
import torch
from src.m28.query_pointer import GeneQueryPointer


def scripted_geometry(length=16):
    """Real head parameters encode two synthetic exon boundaries, not an oracle API."""
    model = GeneQueryPointer(num_queries=2, hidden=4)
    with torch.no_grad():
        for parameter in model.parameters():
            parameter.zero_()
        model.start_key.weight.copy_(torch.eye(4))
        model.end_key.weight.copy_(torch.eye(4))
        model.pointer_query.weight.copy_(torch.eye(4))
        model.event.weight[0, 0] = model.event.weight[2, 0] = 1
        model.event.weight[1, 1] = 1
        model.exists.bias.fill_(1)
        # GRU update z~0; hidden component 2 reads the last boundary feature.
        model.recurrent.bias_ih[4:8].fill_(-20)
        model.recurrent.weight_ih[10, 2] = 20
        model.terminate.weight[0, 2] = 10
    features = torch.zeros(length, 4)
    features[1, 0], features[8, 0] = 3, 2
    features[3, 1], features[-1, 1] = 3, 2
    features[3, 2], features[-1, 2] = -1, 1
    return model, features


class QueryPointerTests(unittest.TestCase):
    def test_half_open_terminal_boundary_and_multiexon_roundtrip(self):
        model, features = scripted_geometry()
        outputs = model.greedy_queries(features, max_exons=3)
        self.assertEqual(len(outputs), 2)
        for row in outputs:
            self.assertEqual(row['chain'], ((1,4),(8,16)))
            self.assertTrue(row['terminated'])
            self.assertEqual(row['stop_reason'], 'predicted_termination')
            self.assertGreater(row['query_logit'], 0)

    def test_cap_does_not_manufacture_an_end_marker(self):
        model, features = scripted_geometry()
        row = model.greedy_queries(features, max_exons=1)[0]
        self.assertEqual(row['chain'], ((1,4),))
        self.assertFalse(row['terminated'])
        self.assertEqual(row['stop_reason'], 'exon_limit')

    def test_continue_at_window_end_is_not_a_forced_stop(self):
        model, features = scripted_geometry()
        with torch.no_grad():
            model.terminate.weight.zero_()
            model.terminate.bias.fill_(-1)
        row = model.greedy_queries(features, max_exons=4)[0]
        self.assertEqual(row['chain'], ((1,4),(8,16)))
        self.assertFalse(row['terminated'])
        self.assertEqual(row['stop_reason'], 'no_next_start')

    def test_padding_cannot_change_context_or_pointers(self):
        model, features = scripted_geometry()
        padded = torch.cat([features, torch.full((7,4), 10000.)])
        self.assertEqual(model.greedy_queries(features, max_exons=3),
                         model.greedy_queries(padded, max_exons=3, valid_bases=16))

    def test_new_head_has_finite_gradients_but_input_features_are_frozen(self):
        torch.manual_seed(0)
        model = GeneQueryPointer(num_queries=2, hidden=8)
        features = torch.randn(19,8,requires_grad=True)
        parts = model.teacher_components(features, [1,-1], {0:[(1,4),(8,16)]}, valid_bases=16)
        loss = sum(parts[k] for k in ('existence_bce','boundary_ce','termination_bce'))
        self.assertTrue(torch.isfinite(loss))
        loss.backward()
        self.assertIsNone(features.grad)
        gradients = [p.grad for p in model.parameters() if p.grad is not None]
        self.assertTrue(all(torch.isfinite(g).all() for g in gradients))
        self.assertGreater(sum(float(g.abs().sum()) for g in gradients), 0)
        self.assertEqual(float(model.queries.grad[1].abs().sum()), 0)

    def test_all_unknown_queries_are_protected_differentiable_zero(self):
        model = GeneQueryPointer(num_queries=2, hidden=8)
        features = torch.randn(16,8,requires_grad=True)
        parts = model.teacher_components(features, [-1,-1], {})
        loss = sum(parts[k] for k in ('existence_bce','boundary_ce','termination_bce'))
        self.assertEqual(float(loss), 0)
        loss.backward()
        self.assertEqual(float(model.queries.grad.abs().sum()), 0)
        self.assertIsNone(features.grad)

    def test_empty_reference_background_has_no_fabricated_pointer_target(self):
        model = GeneQueryPointer(num_queries=2, hidden=8)
        parts = model.teacher_components(torch.randn(16,8), [0,0], {})
        self.assertEqual((parts['positive_queries'], parts['known_queries']), (0,2))
        self.assertTrue(torch.isfinite(parts['existence_bce']))
        self.assertEqual(float(parts['boundary_ce']), 0)
        self.assertEqual(float(parts['termination_bce']), 0)

    def test_labels_assignment_and_coordinate_errors_are_not_silently_repaired(self):
        model = GeneQueryPointer(num_queries=2, hidden=8)
        features = torch.randn(16,8)
        for labels, chains in [
                ([1,0], {}),
                ([0,-1], {1:[(1,4)]}),
                ([1,0], {0:[(8,12),(1,4)]}),
                ([1,0], {0:[(1,17)]}),
                ([1,0], {0:[(1,4),(4,8)]})]:
            with self.assertRaises(ValueError):
                model.teacher_components(features, labels, chains)

    def test_queries_are_not_a_single_mutually_exclusive_gene_class(self):
        model, features = scripted_geometry()
        parts = model.teacher_components(features, [1,1], {0:[(1,4)],1:[(8,16)]})
        self.assertEqual(parts['positive_queries'], 2)
        expected = torch.nn.functional.softplus(torch.tensor(-1.))
        self.assertAlmostEqual(float(parts['existence_bce']), float(expected), places=6)


if __name__ == '__main__':
    unittest.main()
