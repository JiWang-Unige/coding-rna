"""Indexed sum-product: original answers/gradients and independent DNA oracle."""
import torch
import unittest

import test_m28_full_grammar_partition as P
from src.m28.full_grammar_partition import full_grammar_log_partition
from src.m28.full_grammar_partition_indexed import indexed_full_grammar_log_partition


def compare_indexed(test, sequence, logits, cap):
    old_input = logits.detach().requires_grad_(True)
    new_input = logits.detach().requires_grad_(True)
    old = full_grammar_log_partition(sequence, old_input, cap)
    new = indexed_full_grammar_log_partition(sequence, new_input, cap)
    test.assertEqual(float(new['log_partition'].detach()), float(old['log_partition'].detach()))
    grad, = torch.autograd.grad(new['log_partition'], new_input)
    previous, = torch.autograd.grad(old['log_partition'], old_input)
    test.assertTrue(torch.equal(grad, previous))
    for key, value in old['counts'].items(): test.assertEqual(new['counts'][key], value)
    for key in ('reference_used', 'enumerated_candidate_pool_used', 'Viterbi_winner_pruning'):
        test.assertIs(new[key], False)
    test.assertLessEqual(new['counts']['indexed_events_inspected'],
                         new['counts']['original_all_event_loop_iterations'])
    return new


class IndexedFullGrammarPartitionTests(P.FullGrammarPartitionTests):
    def compare(self, sequence, logits, cap):
        independent = super().compare(sequence, logits, cap)
        compare_indexed(self, sequence, logits, cap)
        return independent

    def test_FP32_common_logit_shift_and_large_finite_gain(self):
        sequence = 'ATGCGTTAGTTTAA'
        logits = torch.arange(len(sequence)*15).reshape(len(sequence), 15).float()%7-3
        first = indexed_full_grammar_log_partition(sequence, logits.requires_grad_(True), 2)['log_partition']
        shifted = (logits.detach()+torch.arange(len(sequence))[:, None]*8).requires_grad_(True)
        second = indexed_full_grammar_log_partition(sequence.lower(), shifted, 2)['log_partition']
        self.assertEqual(float(first.detach()), float(second.detach()))
        self.assertTrue(torch.equal(torch.autograd.grad(first, logits)[0],
                                   torch.autograd.grad(second, shifted)[0]))
        for scale in (-1000., 1000.):
            values = torch.full((len(sequence), 15), scale, dtype=torch.float64)
            values[:, 0] = 0.
            self.compare(sequence, values, 2)

    def test_shape_nonfinite_and_explicit_cap(self):
        sequence = 'ATGAAATAA'
        for values in (torch.zeros(8, 15), torch.zeros(9, 14),
                       torch.full((9, 15), float('nan')), torch.zeros(9, 15, dtype=torch.int64)):
            with self.assertRaises(ValueError): indexed_full_grammar_log_partition(sequence, values, 2)
        with self.assertRaises(ValueError): indexed_full_grammar_log_partition(sequence, torch.zeros(9, 15), 0)

    def test_dense_cap128_same_answers_gradients_counts_and_fewer_event_inspections(self):
        generator = torch.Generator().manual_seed(89)
        sequence = 'ATG'+''.join('ACGT'[i] for i in torch.randint(0, 4, (384,), generator=generator).tolist())+'TAA'
        logits = torch.randn(len(sequence), 15, generator=generator)
        actual = compare_indexed(self, sequence, logits, 128)
        self.assertLess(actual['counts']['indexed_events_inspected'],
                        actual['counts']['original_all_event_loop_iterations'])


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
