import itertools
import unittest
import torch
from src.m28.query_pointer import GeneQueryPointer
from src.m28.query_training import matched_query_components


COEFFICIENTS = (0.7, 1.3, 0.2)  # Synthetic fixture, not real-fit coefficients.
CHAINS = (((0,4),(8,14)), ((3,11),))


def fixture():
    torch.manual_seed(13)
    model = GeneQueryPointer(num_queries=3, hidden=8)
    with torch.no_grad():
        model.queries.copy_(2*torch.randn_like(model.queries))
        model.exists.weight.mul_(4)
    features = torch.randn(18,8,requires_grad=True)
    return model, features


def explicit_total(model, features, chains, gene_to_query, unmatched_label):
    labels = [unmatched_label]*model.num_queries
    assigned = {}
    for gene, query in enumerate(gene_to_query):
        labels[query] = 1
        assigned[query] = chains[gene]
    parts = model.teacher_components(features, labels, assigned)
    return sum(c*parts[name] for c,name in zip(
        COEFFICIENTS, ('existence_bce','boundary_ce','termination_bce')))


class QueryMatchingTests(unittest.TestCase):
    def test_assignment_minimizes_actual_whole_window_objective(self):
        model, features = fixture()
        for unmatched in (-1,0):
            result = matched_query_components(model, features, CHAINS, COEFFICIENTS, unmatched)
            alternatives = [float(explicit_total(model, features, CHAINS, assignment, unmatched))
                            for assignment in itertools.permutations(range(3),2)]
            self.assertAlmostEqual(float(result['total']), min(alternatives), places=6)
            self.assertEqual(len(set(result['gene_to_query'])), 2)
            self.assertEqual(sorted(result['positive_chains'].values()), sorted(CHAINS))

    def test_existence_null_baseline_and_component_normalizers_match_loss_differences(self):
        model, features = fixture()
        # Single-gene cost differences must equal full-objective differences:
        # catches omission of unmatched negative BCE or inconsistent means.
        for unmatched in (-1,0):
            result = matched_query_components(model, features, CHAINS[:1], COEFFICIENTS, unmatched)
            costs = result['assignment_costs'][0]
            totals = [float(explicit_total(model, features, CHAINS[:1], (q,), unmatched))
                      for q in range(3)]
            for q in (1,2):
                self.assertAlmostEqual(float(costs[q]-costs[0]), totals[q]-totals[0], places=6)

    def test_reference_order_does_not_change_minimum_training_loss(self):
        model, features = fixture()
        forward = matched_query_components(model, features, CHAINS, COEFFICIENTS, -1)
        reverse = matched_query_components(model, features, CHAINS[::-1], COEFFICIENTS, -1)
        self.assertAlmostEqual(float(forward['total']), float(reverse['total']), places=6)
        self.assertEqual(sorted(forward['positive_chains'].values()),
                         sorted(reverse['positive_chains'].values()))

    def test_capacity_failure_cannot_silently_discard_genes(self):
        model, features = fixture()
        four = [((a,a+3),) for a in (0,4,8,12)]
        with self.assertRaisesRegex(ValueError, 'capacity'):
            matched_query_components(model, features, four, COEFFICIENTS, -1)

    def test_unknown_queries_protected_and_frozen_input_not_backpropagated(self):
        model, features = fixture()
        result = matched_query_components(model, features, CHAINS[:1], COEFFICIENTS, -1)
        self.assertFalse(result['assignment_costs'].requires_grad)
        result['total'].backward()
        self.assertIsNone(features.grad)
        chosen = result['gene_to_query'][0]
        self.assertGreater(float(model.queries.grad[chosen].abs().sum()), 0)
        for q in set(range(3)) - {chosen}:
            self.assertEqual(result['query_labels'][q], -1)
            self.assertEqual(float(model.queries.grad[q].abs().sum()), 0)
        self.assertTrue(all(torch.isfinite(p.grad).all()
                            for p in model.parameters() if p.grad is not None))

    def test_empty_reference_has_explicit_negative_or_unknown_policy_and_no_fake_chain(self):
        for unmatched in (-1,0):
            model, features = fixture()
            result = matched_query_components(model, features, [], COEFFICIENTS, unmatched)
            self.assertEqual(result['gene_to_query'], ())
            self.assertEqual(result['positive_chains'], {})
            self.assertEqual(result['query_labels'], (unmatched,)*3)
            self.assertEqual(result['components']['positive_queries'], 0)
            self.assertEqual(float(result['components']['boundary_ce']), 0)
            self.assertEqual(float(result['components']['termination_bce']), 0)
            if unmatched == -1:
                self.assertEqual(float(result['total']), 0)
                result['total'].backward()
                self.assertEqual(float(model.queries.grad.abs().sum()), 0)
            else:
                self.assertGreater(float(result['total']), 0)
                self.assertEqual(result['components']['known_queries'], 3)

    def test_padding_ignored_but_out_of_window_reference_not_repaired(self):
        model, features = fixture()
        baseline = matched_query_components(model, features, CHAINS, COEFFICIENTS, -1)
        padded = torch.cat([features, torch.full((4,8), 10000.)])
        other = matched_query_components(model, padded, CHAINS, COEFFICIENTS, -1, valid_bases=18)
        self.assertEqual(baseline['gene_to_query'], other['gene_to_query'])
        self.assertTrue(torch.equal(baseline['assignment_costs'], other['assignment_costs']))
        with self.assertRaisesRegex(ValueError, 'ordered oriented'):
            matched_query_components(model, padded, [((0,19),)], COEFFICIENTS, -1, valid_bases=18)

    def test_nonfinite_matching_cost_stops_instead_of_replacing_target(self):
        model, features = fixture()
        with torch.no_grad():
            model.exists.weight.fill_(float('nan'))
        with self.assertRaisesRegex(ValueError, 'Nonfinite query assignment'):
            matched_query_components(model, features, CHAINS, COEFFICIENTS, -1)


if __name__ == '__main__':
    unittest.main()
