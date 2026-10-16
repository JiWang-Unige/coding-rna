"""A31 fixed emission arithmetic, original fine roles and null selection."""
import itertools
import unittest
import torch
from src.m28.labels import _render
from src.m28.local_path_score import local_path_scores
from src.m28.core import select_nonoverlapping


class LocalPathScoreTests(unittest.TestCase):
    def test_original_reference_render_and_log_probability_identity_both_strands(self):
        torch.manual_seed(0)
        logits = torch.randn(24, 15, dtype=torch.float64)
        chain = ((2, 6), (11, 22))
        actual = local_path_scores(logits, [chain])[0]
        for strand, cds in (('+', [(102, 106, 0), (111, 122, 2)]),
                            ('-', [(102, 113, 2), (118, 122, 0)])):
            labels = _render(dict(CDS=cds, CDS_complete_supported=True),
                dict(start=100, end=124, width=24, valid_bases=24, strand=strand))['seg']
            logp = logits.log_softmax(-1)
            expected = (logp[torch.arange(24), torch.as_tensor(labels)]-logp[:, 0]).sum()
            torch.testing.assert_close(actual, expected, rtol=0, atol=1e-12)

    def test_offset_padding_and_outside_span_do_not_change_gain(self):
        torch.manual_seed(1)
        logits = torch.randn(24, 15, dtype=torch.float64)
        chains = [((2, 6), (11, 22)), ((5, 14),)]
        actual = local_path_scores(logits, chains)
        torch.testing.assert_close(actual, local_path_scores(logits+torch.arange(24)[:, None]*8, chains),
                                   rtol=0, atol=1e-12)
        padded = torch.cat([logits, torch.full((4, 15), 1e9, dtype=torch.float64)])
        torch.testing.assert_close(actual, local_path_scores(padded, chains, 24), rtol=0, atol=0)
        changed = logits.clone(); changed[[0, 1, 22, 23]] = 1e9
        torch.testing.assert_close(actual, local_path_scores(changed, chains), rtol=0, atol=0)

    def test_zero_emissions_select_null_and_positive_interval_objective(self):
        chains = [((0, 6),), ((10, 16),), ((0, 16),)]
        zero = local_path_scores(torch.zeros(18, 15), chains).tolist()
        self.assertEqual(zero, [0., 0., 0.])
        self.assertEqual(select_nonoverlapping(chains, zero), [])
        logits = torch.zeros(18, 15, dtype=torch.float64)
        logits[:, 13:] = 1
        gains = local_path_scores(logits, chains).tolist()
        selected = select_nonoverlapping(chains, gains)
        feasible = []
        for mask in itertools.product((False, True), repeat=3):
            indices = [i for i, use in enumerate(mask) if use]
            if any(chains[i][0][0] < chains[j][-1][1] and chains[j][0][0] < chains[i][-1][1]
                   for i, j in itertools.combinations(indices, 2)):
                continue
            feasible.append(sum(gains[i] for i in indices))
        self.assertEqual(sum(gains[i] for i in selected), max(feasible))
        self.assertEqual(selected, [0, 1])

    def test_empty_pool_no_gradient_and_invalid_emissions(self):
        logits = torch.randn(24, 15, requires_grad=True)
        self.assertEqual(tuple(local_path_scores(logits, []).shape), (0,))
        self.assertFalse(local_path_scores(logits, [((2, 20),)]).requires_grad)
        with self.assertRaisesRegex(ValueError, '15-class'):
            local_path_scores(torch.zeros(24, 14), [])
        bad = logits.detach().clone(); bad[4, 2] = float('nan')
        with self.assertRaisesRegex(ValueError, 'Nonfinite'):
            local_path_scores(bad, [((2, 20),)])


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
