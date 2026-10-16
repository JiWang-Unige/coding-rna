import unittest
import torch
from src.m28.query_pointer import GeneQueryPointer
from src.m28.query_training import teacher_chain_components


# Fixed FP32 absolute tolerances, specified before the Slurm run.
LOSS_ATOL = 1e-6
STATE_ATOL = 2e-6


def fixture(q=4):
    torch.manual_seed(17)
    model = GeneQueryPointer(num_queries=q, hidden=8)
    features = torch.randn(64, 8)
    return model, features


class BatchedTeacherTests(unittest.TestCase):
    def test_batched_pointer_and_recurrence_match_each_scalar_query(self):
        for q in (1, 4):
            model, features = fixture(q)
            window = model.for_window(features)
            for stage, lower, coordinate in (
                    ('start', 0, 0), ('end', 0, 64),
                    ('next_start', 12, 13), ('end', 62, 64)):
                batched = window.pointer_logits(window.initial, stage, lower)
                scalar = torch.stack([
                    window.pointer_logits(state, stage, lower)
                    for state in window.initial])
                torch.testing.assert_close(batched, scalar, atol=LOSS_ATOL, rtol=0)
                advanced = window.advance(window.initial, stage, coordinate)
                scalar_states = torch.stack([
                    window.advance(state, stage, coordinate)
                    for state in window.initial])
                torch.testing.assert_close(advanced, scalar_states, atol=STATE_ATOL, rtol=0)

    def test_chain_event_means_match_independent_scalar_teachers_without_grad(self):
        model, features = fixture()
        window = model.for_window(features)
        chains = (((0, 64),), ((0, 3), (8, 17), (23, 64)),
                  tuple((6*j, 6*j+3) for j in range(10)) + ((60, 64),))
        for chain in chains:
            batched = teacher_chain_components(window, chain)
            scalar = []
            for i in range(model.num_queries):
                labels = [-1]*model.num_queries
                labels[i] = 1
                scalar.append(window.teacher_components(labels, {i: chain}))
            for name in ('boundary_ce', 'termination_bce'):
                expected = torch.stack([part[name] for part in scalar])
                self.assertFalse(batched[name].requires_grad)
                torch.testing.assert_close(batched[name], expected, atol=LOSS_ATOL, rtol=0)

    def test_padding_does_not_enter_batched_costs(self):
        model, features = fixture()
        chain = ((0, 3), (8, 17), (23, 64))
        baseline = teacher_chain_components(model.for_window(features), chain)
        padded = torch.cat([features, torch.full((4, 8), 10000.)])
        other = teacher_chain_components(model.for_window(padded, valid_bases=64), chain)
        for name in ('boundary_ce', 'termination_bce'):
            self.assertTrue(torch.equal(baseline[name], other[name]))


if __name__ == '__main__':
    unittest.main()
