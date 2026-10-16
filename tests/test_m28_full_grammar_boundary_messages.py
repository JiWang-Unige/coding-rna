"""Artificial chunk messages versus independent complete-chain-set enumeration."""
import math
import unittest

import torch

from src.m28.full_grammar_boundary_messages import (
    BACKGROUND, BoundaryMessage, boundary_forward,
    boundary_full_grammar_log_partition, boundary_core_observation_nll)
from src.m28.full_grammar_observed_partition import (
    InfeasibleObservation, observed_full_grammar_log_partition, partial_observation_nll)
from src.m28.phase_evidence_chain import candidate_roles
from test_m28_full_grammar_observed_partition import exhaustive_observed, exact_roles


class BoundaryMessageTests(unittest.TestCase):
    def compare(self, sequence, values, cap, allowed=None, cuts=(), scale=1.):
        logits = values.detach().requires_grad_(True)
        allowed = torch.ones_like(logits, dtype=torch.bool) if allowed is None else allowed
        expected, sets = exhaustive_observed(sequence, logits, cap, allowed)
        if expected is None:
            with self.assertRaises(InfeasibleObservation):
                boundary_full_grammar_log_partition(sequence, logits, cap, allowed, cuts)
            return None, sets, None
        result = boundary_full_grammar_log_partition(sequence, logits, cap, allowed, cuts)
        old = observed_full_grammar_log_partition(sequence, logits, cap, allowed)
        self.assertEqual(result['log_partition'].dtype, torch.float64)
        for value in (expected, old['log_partition']):
            self.assertAlmostEqual(float(result['log_partition']), float(value), places=10)
        actual_grad, = torch.autograd.grad(result['log_partition'] * scale, logits, retain_graph=True)
        tolerance = (1e-9, 1e-10) if logits.dtype == torch.float64 else (1e-6, 1e-6)
        for value in (expected, old['log_partition']):
            gradient, = torch.autograd.grad(value * scale, logits, retain_graph=True)
            torch.testing.assert_close(actual_grad, gradient, rtol=tolerance[0], atol=tolerance[1])
        self.assertFalse(result['interior_edges_closed'])
        self.assertTrue(result['global_complete_gene_grammar'])
        return result, sets, actual_grad

    def test_global_gene_set_mass_and_full_gradients_match_two_oracles(self):
        sequences = ('ATGAAATAAATGCCCTAA', 'ATGCGTTAGTTTAA', 'ATGGCTAGAAATAA',
            'ATGCGTTAGAGTTAGAGCTAGTAA', 'ATGTGTTAGATGCGTTAGAATAA',
            'ATGCGTNAGTTTAA', 'ATGNAATAA', 'ATGAAAGTAGTCCCCAGAAATAA')
        generator = torch.Generator().manual_seed(103)
        for sequence in sequences:
            for cap in (1, 2, 4):
                for dtype, scale in ((torch.float32, -.7), (torch.float64, 2.3)):
                    with self.subTest(sequence=sequence, cap=cap, dtype=dtype):
                        values = torch.randn(len(sequence), 15, dtype=dtype, generator=generator)
                        self.compare(sequence, values, cap, scale=scale)

    def test_every_cut_and_single_base_chunks_preserve_mass_and_gradients(self):
        for sequence in ('ATGCGTTAGAGTTAGAGCTAGTAA', 'ATGGCTAGAAATAA', 'ATGAAATAAATGCCCTAA'):
            logits = (torch.arange(len(sequence)*15).reshape(-1, 15).double()/37).requires_grad_(True)
            whole = boundary_full_grammar_log_partition(sequence, logits, 4)['log_partition']
            whole_grad, = torch.autograd.grad(whole, logits, retain_graph=True)
            for cuts in [(i,) for i in range(1, len(sequence))] + [tuple(range(1, len(sequence)))]:
                split = boundary_full_grammar_log_partition(sequence, logits, 4, cuts=cuts)['log_partition']
                self.assertEqual(float(split), float(whole))
                gradient, = torch.autograd.grad(split, logits, retain_graph=True)
                torch.testing.assert_close(gradient, whole_grad, rtol=0, atol=0)

    def test_literal_partial_motif_carry_and_maturity_states(self):
        sequence = 'ATGCGTTAGTTTAA'
        roles = torch.tensor([13,13,13,7,4,4,4,4,4,11,3,14,14,14])
        logits = torch.zeros(14, 15, dtype=torch.float64)
        allowed = exact_roles(roles)
        message = BoundaryMessage(2, {BACKGROUND: logits.sum()*0})
        expected = {1:('S1',1,''), 2:('S2',1,''), 3:('F',1,''), 4:('D0',1,'C'),
            5:('D1',1,'C'), 6:('I2',1,'C'), 7:('I3',1,'C'), 9:('X',2,'C'),
            10:('E',2,'CT'), 11:('E',2,''), 12:('T',2,''), 13:('TA',2,''), 14:BACKGROUND}
        for i in range(14):
            message = boundary_forward(sequence[i:i+1], logits[i:i+1], message, allowed[i:i+1])
            if i+1 in expected:
                self.assertIn(expected[i+1], message.log_masses)
        self.assertEqual(set(message.log_masses), {BACKGROUND})
        self.assertEqual(float(message.log_masses[BACKGROUND]), 0.)

    def test_core_observation_has_real_prefix_suffix_support_and_NLL_gradients(self):
        sequence = 'ATGCGTTAGTTTAA'
        roles = torch.tensor([13,13,13,7,4,4,4,4,4,11,3,14,14,14])
        a, b = 3, 12
        generator = torch.Generator().manual_seed(104)
        for dtype, scale in ((torch.float32, -.7), (torch.float64, 2.3)):
            logits = torch.randn(14, 15, dtype=dtype, generator=generator).requires_grad_(True)
            core = exact_roles(roles[a:b])
            with self.assertRaises(InfeasibleObservation):
                observed_full_grammar_log_partition(sequence[a:b], logits[a:b], 2, core)
            loss = boundary_core_observation_nll(sequence, logits, 2, a, b, core, cuts=(a,b))
            allowed = torch.ones(14, 15, dtype=torch.bool); allowed[a:b] = core
            expected_obs, sets = exhaustive_observed(sequence, logits, 2, allowed)
            expected_all, _ = exhaustive_observed(sequence, logits, 2, torch.ones_like(allowed))
            expected = expected_all - expected_obs
            old = partial_observation_nll(sequence, logits, 2, allowed)['loss']
            self.assertEqual(sets, 1)
            for value in (expected, old):
                self.assertAlmostEqual(float(loss['loss']), float(value), places=10)
            actual_grad, = torch.autograd.grad(loss['loss']*scale, logits, retain_graph=True)
            tolerance = (1e-9, 1e-10) if dtype == torch.float64 else (1e-6, 1e-6)
            for value in (expected, old):
                gradient, = torch.autograd.grad(value*scale, logits, retain_graph=True)
                torch.testing.assert_close(actual_grad, gradient, rtol=tolerance[0], atol=tolerance[1])
            self.assertGreater(float(actual_grad[:a].abs().sum()), 0.)
            self.assertGreater(float(actual_grad[b:].abs().sum()), 0.)
            self.assertFalse(loss['loss_normalized'])
            self.assertFalse(loss['context_extent_selected'])

    def test_outside_and_unresolved_genes_stay_marginalized(self):
        sequence = 'ATGAAATAA'*6
        logits = torch.zeros(len(sequence), 15, dtype=torch.float64)
        allowed = torch.ones(27, 15, dtype=torch.bool)
        allowed[:3] = False; allowed[:3,13] = True
        allowed[18:27] = False; allowed[18:27,0] = True
        loss = boundary_core_observation_nll(sequence, logits, 1, 0, 27, allowed, cuts=(9,18,27,36,45))
        self.assertAlmostEqual(float(loss['log_partition_all']), math.log(64), places=12)
        self.assertAlmostEqual(float(loss['log_partition_observed']), math.log(16), places=12)
        self.assertAlmostEqual(float(loss['loss']), math.log(4), places=12)

    def test_exon_cap_and_short_exon_priority_are_not_reset_at_cuts(self):
        sequence = 'ATGCGTTAGAGTTAGAGCTAGTAA'
        chain = ((0,4),(9,10),(15,16),(21,24))
        allowed = exact_roles(candidate_roles(chain, len(sequence)))
        logits = torch.zeros(len(sequence), 15, dtype=torch.float64)
        cuts = (2,4,8,9,10,14,15,16,20,21,22,23)
        _, sets, _ = self.compare(sequence, logits, 4, allowed, cuts)
        self.assertEqual(sets, 1)
        _, sets, _ = self.compare(sequence, logits, 3, allowed, cuts)
        self.assertEqual(sets, 0)

    def test_actual_carry_rejects_internal_STOP_across_chunks_and_exons(self):
        sequence = 'ATGTGTTAGATGCGTTAGAATAA'
        bad_chain = ((0,4),(18,23))
        allowed = exact_roles(candidate_roles(bad_chain, len(sequence)))
        _, sets, _ = self.compare(sequence, torch.zeros(len(sequence),15,dtype=torch.float64),
            3, allowed, cuts=(4,8,18,19,20))
        self.assertEqual(sets, 0)

    def test_all_unknown_zero_NLL_and_common_row_shift(self):
        sequence = 'ATGCGTTAGTTTAA'
        logits = (torch.arange(14*15).reshape(14,15).float()%7-3).requires_grad_(True)
        unknown = torch.ones(14,15,dtype=torch.bool)
        loss = boundary_core_observation_nll(sequence, logits, 128, 3, 12, unknown[3:12], cuts=(1,2,3,12,13))
        self.assertEqual(float(loss['loss']), 0.)
        self.assertTrue(torch.equal(torch.autograd.grad(loss['loss'], logits)[0], torch.zeros_like(logits)))
        values = logits.detach().requires_grad_(True)
        shifted = (values.detach()+torch.arange(14)[:,None]*8).requires_grad_(True)
        first = boundary_full_grammar_log_partition(sequence, values, 128, cuts=(4,7,12))['log_partition']
        second = boundary_full_grammar_log_partition(sequence.lower(), shifted, 128, cuts=(4,7,12))['log_partition']
        self.assertEqual(float(first), float(second))
        torch.testing.assert_close(torch.autograd.grad(first,values)[0], torch.autograd.grad(second,shifted)[0], rtol=0, atol=0)

    def test_GC_donor_minimum_intron_and_background_observations(self):
        for sequence, chain, sets in (('ATGGCTAGAAATAA',((0,3),(8,14)),1),
                ('ATGGTAGAAATAA',((0,3),(7,13)),0)):
            allowed = exact_roles(candidate_roles(chain,len(sequence)))
            _, actual, _ = self.compare(sequence, torch.zeros(len(sequence),15,dtype=torch.float64),
                2, allowed, cuts=tuple(range(1,len(sequence))))
            self.assertEqual(actual,sets)
        sequence = 'CCCATGAAATAACCCATGCCCTAACCC'
        logits = torch.zeros(len(sequence),15,dtype=torch.float64)
        for p in (0,13,len(sequence)-1):
            allowed = torch.ones(len(sequence),15,dtype=torch.bool)
            allowed[p] = False; allowed[p,1] = True
            _, sets, _ = self.compare(sequence,logits,1,allowed,cuts=(3,12,15,24))
            self.assertEqual(sets,0)

    def test_global_completion_empty_identity_and_invalid_inputs(self):
        sequence = 'ATGAAATA'
        logits = torch.zeros(8,15,dtype=torch.float64)
        allowed = torch.ones(8,15,dtype=torch.bool)
        allowed[:3] = False; allowed[:3,13] = True
        with self.assertRaises(InfeasibleObservation):
            boundary_full_grammar_log_partition(sequence,logits,2,allowed,cuts=(1,2,7))
        message = BoundaryMessage(2,{BACKGROUND:logits.sum()*0})
        self.assertEqual(boundary_forward('',logits[:0],message).log_masses,message.log_masses)
        for values,cap,cuts in ((torch.zeros(7,15),2,()),(torch.full((8,15),float('nan')),2,()),
                (logits,0,()),(logits,2,(3,3)),(logits,2,(0,))):
            with self.assertRaises(ValueError):
                boundary_full_grammar_log_partition(sequence,values,cap,cuts=cuts)
        with self.assertRaises(ValueError):
            boundary_core_observation_nll(sequence,logits,2,0,8,torch.ones(8,15))


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
