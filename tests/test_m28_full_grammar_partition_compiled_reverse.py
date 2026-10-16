"""Compiled reverse versus independent oracle, Torch and original reverse."""
import torch
import unittest

import test_m28_full_grammar_partition as P
from src.m28.full_grammar_partition import full_grammar_log_partition
from src.m28.full_grammar_partition_compiled_reverse import (compiled_reverse_full_grammar_log_partition,
    _accumulate_reverse_nodes)
from src.m28.full_grammar_partition_reverse import reverse_full_grammar_log_partition as plain_reverse


def compare_compiled_reverse(test, sequence, logits, cap, scale=1.):
    old_input = logits.detach().requires_grad_(True)
    new_input = logits.detach().requires_grad_(True)
    old = full_grammar_log_partition(sequence, old_input, cap)
    new = compiled_reverse_full_grammar_log_partition(sequence, new_input, cap)
    test.assertEqual(new['log_partition'].dtype, torch.float64)
    test.assertAlmostEqual(float(new['log_partition'].detach()),
                           float(old['log_partition'].detach()), places=10)
    grad, = torch.autograd.grad(scale * new['log_partition'], new_input)
    previous, = torch.autograd.grad(scale * old['log_partition'], old_input)
    relative, absolute = (1e-9, 1e-10) if logits.dtype == torch.float64 else (1e-6, 1e-6)
    torch.testing.assert_close(grad, previous, rtol=relative, atol=absolute)
    for key, value in old['counts'].items():
        test.assertEqual(new['counts'][key], value)
    for key in ('reference_used', 'enumerated_candidate_pool_used', 'Viterbi_winner_pruning'):
        test.assertIs(new[key], False)
    return new, grad


class CompiledReverseFullGrammarPartitionTests(P.FullGrammarPartitionTests):
    def compare(self, sequence, logits, cap):
        _, chains, sets, _ = super().compare(sequence, logits, cap)
        new, gradient = compare_compiled_reverse(self, sequence, logits, cap)
        return new, chains, sets, gradient

    def test_FP32_common_logit_shift_and_large_finite_gain(self):
        sequence = 'ATGCGTTAGTTTAA'
        logits = torch.arange(len(sequence)*15).reshape(len(sequence), 15).float()%7-3
        first = compiled_reverse_full_grammar_log_partition(sequence, logits.requires_grad_(True), 2)['log_partition']
        shifted = (logits.detach()+torch.arange(len(sequence))[:, None]*8).requires_grad_(True)
        second = compiled_reverse_full_grammar_log_partition(sequence.lower(), shifted, 2)['log_partition']
        self.assertEqual(float(first.detach()), float(second.detach()))
        self.assertTrue(torch.equal(torch.autograd.grad(first, logits)[0],
                                   torch.autograd.grad(second, shifted)[0]))
        for scale in (-1000., 1000.):
            values = torch.full((len(sequence), 15), scale, dtype=torch.float64)
            values[:, 0] = 0.
            result, _, _, gradient = self.compare(sequence, values, 2)
            self.assertTrue(bool(torch.isfinite(result['log_partition'])))
            self.assertTrue(bool(torch.isfinite(gradient).all()))

    def test_shape_nonfinite_and_explicit_cap(self):
        sequence = 'ATGAAATAA'
        for values in (torch.zeros(8, 15), torch.zeros(9, 14),
                       torch.full((9, 15), float('nan')), torch.zeros(9, 15, dtype=torch.int64)):
            with self.assertRaises(ValueError): compiled_reverse_full_grammar_log_partition(sequence, values, 2)
        with self.assertRaises(ValueError): compiled_reverse_full_grammar_log_partition(sequence, torch.zeros(9, 15), 0)

    def test_dense_cap128_same_logZ_gradients_and_original_counts(self):
        generator = torch.Generator().manual_seed(89)
        sequence = 'ATG'+''.join('ACGT'[i] for i in torch.randint(0, 4, (384,), generator=generator).tolist())+'TAA'
        logits = torch.randn(len(sequence), 15, generator=generator)
        compare_compiled_reverse(self, sequence, logits, 128)

    def test_nonunit_grad_output_FP32_FP64_and_first_order_graph_attachment(self):
        sequence = 'ATGCGTTAGAGTTAGAGCTAGTAA'
        generator = torch.Generator().manual_seed(91)
        for dtype in (torch.float32, torch.float64):
            logits = torch.randn(len(sequence), 15, dtype=dtype, generator=generator)
            for scale in (-0.7, 2.3):
                result, _ = compare_compiled_reverse(self, sequence, logits, 4, scale)
                self.assertTrue(result['log_partition'].requires_grad)
                self.assertIsNotNone(result['log_partition'].grad_fn)


    def test_same_forward_full_counts_all_gradients_and_actual_nopython_execution(self):
        sequence = 'ATGCGTTAGAGTTAGAGCTAGTAA' * 4
        generator = torch.Generator().manual_seed(93)
        for dtype in (torch.float32, torch.float64):
            values = torch.randn(len(sequence), 15, dtype=dtype, generator=generator)
            original = values.detach().requires_grad_(True)
            compiled = values.detach().requires_grad_(True)
            old = plain_reverse(sequence, original, 128)
            new = compiled_reverse_full_grammar_log_partition(sequence, compiled, 128)
            self.assertEqual(float(new['log_partition'].detach()), float(old['log_partition'].detach()))
            self.assertEqual(new['counts'], old['counts'])
            for scale in (-0.7, 2.3):
                before, = torch.autograd.grad(scale * old['log_partition'], original, retain_graph=True)
                after, = torch.autograd.grad(scale * new['log_partition'], compiled, retain_graph=True)
                relative, absolute = (1e-9, 1e-10) if dtype == torch.float64 else (1e-6, 1e-6)
                torch.testing.assert_close(after, before, rtol=relative, atol=absolute)
        self.assertGreater(len(_accumulate_reverse_nodes.nopython_signatures), 0)


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
