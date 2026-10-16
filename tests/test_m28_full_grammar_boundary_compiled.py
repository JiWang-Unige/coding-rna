"""Compiled messages: independent set oracle and native incoming-message gradients."""
import math
import unittest

import torch

from src.m28.full_grammar_boundary_compiled import (
    CompiledBoundaryMessage, _state_table, _forward_messages, _reverse_messages,
    compiled_boundary_forward, compiled_boundary_full_grammar_log_partition,
    compiled_boundary_core_observation_nll)
from src.m28.full_grammar_boundary_messages import (
    BACKGROUND, BoundaryMessage, boundary_forward, boundary_full_grammar_log_partition,
    boundary_core_observation_nll)
from src.m28.full_grammar_observed_partition import (
    InfeasibleObservation, observed_full_grammar_log_partition)
from src.m28.phase_evidence_chain import candidate_roles
from test_m28_full_grammar_observed_partition import exhaustive_observed, exact_roles


class CompiledBoundaryTests(unittest.TestCase):
    def close_gradients(self, actual, expected, dtype):
        rtol, atol = (1e-9, 1e-10) if dtype == torch.float64 else (1e-6, 1e-6)
        torch.testing.assert_close(actual, expected, rtol=rtol, atol=atol)

    def compare(self, sequence, values, cap, allowed=None, cuts=(), scale=1., span=3):
        logits = values.detach().requires_grad_(True)
        mask = torch.ones_like(logits, dtype=torch.bool) if allowed is None else allowed
        expected, sets = exhaustive_observed(sequence, logits, cap, mask)
        if expected is None:
            with self.assertRaises(InfeasibleObservation):
                compiled_boundary_full_grammar_log_partition(sequence, logits, cap, mask, cuts, span)
            return None, sets, None
        result = compiled_boundary_full_grammar_log_partition(sequence, logits, cap, mask, cuts, span)
        native = boundary_full_grammar_log_partition(sequence, logits, cap, mask, cuts)['log_partition']
        original = observed_full_grammar_log_partition(sequence, logits, cap, mask)['log_partition']
        value = result['log_partition']
        self.assertEqual(value.dtype, torch.float64)
        actual, = torch.autograd.grad(value*scale, logits, retain_graph=True)
        for oracle in (expected, native, original):
            self.assertAlmostEqual(float(value), float(oracle), places=10)
            gradient, = torch.autograd.grad(oracle*scale, logits, retain_graph=True)
            self.close_gradients(actual, gradient, logits.dtype)
        return result, sets, actual

    def initial(self, cap):
        states, _, _, bg = _state_table(cap)
        values = torch.full((len(states),), -float('inf'), dtype=torch.float64)
        values[bg] = 0.
        return CompiledBoundaryMessage(cap, values)

    def test_closed_mass_and_all_gradients_match_three_oracles(self):
        generator = torch.Generator().manual_seed(105)
        sequences = ('ATGAAATAAATGCCCTAA', 'ATGCGTTAGTTTAA', 'ATGGCTAGAAATAA',
            'ATGCGTTAGAGTTAGAGCTAGTAA', 'ATGTGTTAGATGCGTTAGAATAA',
            'ATGCGTNAGTTTAA', 'ATGNAATAA', 'ATGAAAGTAGTCCCCAGAAATAA')
        for sequence in sequences:
            for cap in (1, 2, 4):
                for dtype, scale in ((torch.float32, -.7), (torch.float64, 2.3)):
                    with self.subTest(sequence=sequence, cap=cap, dtype=dtype):
                        self.compare(sequence, torch.randn(len(sequence),15,dtype=dtype,generator=generator),
                            cap, scale=scale)

    def test_every_supported_output_message_and_vector_gradient_match_native(self):
        generator = torch.Generator().manual_seed(106)
        for sequence in ('ATGCGTTAG', 'ATGCGTTAGTT', 'ATGCGTRAG', 'ATGAAATAAAT'):
            logits = torch.randn(len(sequence),15,dtype=torch.float64,generator=generator).requires_grad_(True)
            native = boundary_forward(sequence, logits, BoundaryMessage(4,{BACKGROUND:logits.sum()*0}))
            compiled = compiled_boundary_forward(sequence, logits, self.initial(4), checkpoint_bp=2)
            states = _state_table(4)[0]
            indices = [i for i,s in enumerate(states) if s in native.log_masses]
            self.assertEqual([states[i] for i in torch.where(torch.isfinite(compiled.log_masses))[0].tolist()],
                sorted(native.log_masses))
            oracle = torch.stack([native.log_masses[states[i]] for i in indices])
            torch.testing.assert_close(compiled.log_masses[indices],oracle,rtol=1e-9,atol=1e-10)
            weights = torch.randn(len(indices),dtype=torch.float64,generator=generator)
            actual, = torch.autograd.grad((compiled.log_masses[indices]*weights).sum(),logits,retain_graph=True)
            expected, = torch.autograd.grad((oracle*weights).sum(),logits)
            self.close_gradients(actual,expected,logits.dtype)

    def test_incoming_message_gradients_survive_motif_carry_and_splice_cuts(self):
        sequence = 'ATGCGTTAGTTTAA'
        generator = torch.Generator().manual_seed(107)
        values = torch.randn(14,15,dtype=torch.float64,generator=generator)
        states = _state_table(2)[0]
        for cut in (1,2,4,6,7,9,10,12,13):
            prefix = boundary_forward(sequence[:cut],values[:cut],BoundaryMessage(2,{BACKGROUND:values.sum()*0}))
            incoming = torch.full((len(states),),-float('inf'),dtype=torch.float64)
            for i,state in enumerate(states):
                if state in prefix.log_masses: incoming[i] = prefix.log_masses[state].detach()
            incoming.requires_grad_(True)
            logits = values[cut:].detach().requires_grad_(True)
            sparse = {state:incoming[i] for i,state in enumerate(states) if bool(torch.isfinite(incoming[i]))}
            native = boundary_forward(sequence[cut:],logits,BoundaryMessage(2,sparse))
            compiled = compiled_boundary_forward(sequence[cut:],logits,CompiledBoundaryMessage(2,incoming),checkpoint_bp=3)
            indices = [i for i,s in enumerate(states) if s in native.log_masses]
            weights = torch.randn(len(indices),dtype=torch.float64,generator=generator)
            oracle = (torch.stack([native.log_masses[states[i]] for i in indices])*weights).sum()
            actual = (compiled.log_masses[indices]*weights).sum()
            self.assertAlmostEqual(float(actual),float(oracle),places=10)
            a = torch.autograd.grad(actual,(logits,incoming),retain_graph=True)
            b = torch.autograd.grad(oracle,(logits,incoming))
            for x,y in zip(a,b): self.close_gradients(x,y,torch.float64)
            self.assertGreater(float(a[1].abs().sum()),0.)

    def test_all_interior_cuts_and_one_base_chunks_preserve_full_history(self):
        for sequence in ('ATGCGTTAGAGTTAGAGCTAGTAA','ATGGCTAGAAATAA','ATGAAATAAATGCCCTAA'):
            logits = (torch.arange(len(sequence)*15).reshape(-1,15).double()/37).requires_grad_(True)
            whole = compiled_boundary_full_grammar_log_partition(sequence,logits,4,checkpoint_bp=3)['log_partition']
            gradient, = torch.autograd.grad(whole,logits,retain_graph=True)
            for cuts in [(i,) for i in range(1,len(sequence))]+[tuple(range(1,len(sequence)))]:
                split = compiled_boundary_full_grammar_log_partition(sequence,logits,4,cuts=cuts,checkpoint_bp=2)
                self.assertEqual(float(split['log_partition']),float(whole))
                self.close_gradients(torch.autograd.grad(split['log_partition'],logits,retain_graph=True)[0],gradient,logits.dtype)
                self.assertFalse(split['interior_edges_closed'])

    def test_core_NLL_keeps_original_rows_and_actual_outside_gradients(self):
        sequence = 'ATGCGTTAGTTTAA'; a,b = 3,12
        roles = torch.tensor([13,13,13,7,4,4,4,4,4,11,3,14,14,14])
        core = exact_roles(roles[a:b]); mask = torch.ones(14,15,dtype=torch.bool); mask[a:b] = core
        generator = torch.Generator().manual_seed(108)
        for dtype,scale in ((torch.float32,-.7),(torch.float64,2.3)):
            logits = torch.randn(14,15,dtype=dtype,generator=generator).requires_grad_(True)
            with self.assertRaises(InfeasibleObservation):
                observed_full_grammar_log_partition(sequence[a:b],logits[a:b],2,core)
            result = compiled_boundary_core_observation_nll(sequence,logits,2,a,b,core,cuts=(a,b),checkpoint_bp=2)
            observed,sets = exhaustive_observed(sequence,logits,2,mask)
            all_mass,_ = exhaustive_observed(sequence,logits,2,torch.ones_like(mask))
            native = boundary_core_observation_nll(sequence,logits,2,a,b,core,cuts=(a,b))['loss']
            self.assertEqual(sets,1)
            actual, = torch.autograd.grad(result['loss']*scale,logits,retain_graph=True)
            for expected in (all_mass-observed,native):
                self.assertAlmostEqual(float(result['loss']),float(expected),places=10)
                self.close_gradients(actual,torch.autograd.grad(expected*scale,logits,retain_graph=True)[0],dtype)
            self.assertGreater(float(actual[:a].abs().sum()),0.)
            self.assertGreater(float(actual[b:].abs().sum()),0.)
            self.assertFalse(result['loss_normalized']); self.assertFalse(result['context_extent_selected'])

    def test_partial_role_masks_match_independent_sets(self):
        generator = torch.Generator().manual_seed(109)
        for sequence in ('ATGAAATAAATGCCCTAA','ATGCGTTAGTTTAA','ATGGCTAGAAATAA',
                'ATGCGTTAGAGTTAGAGCTAGTAA','ATGTGTTAGATGCGTTAGAATAA','ATGCGTNAGTTTAA'):
            for cap in (1,2,4):
                mask = torch.rand(len(sequence),15,generator=generator) > .55
                mask[:,0] = True
                for dtype,scale in ((torch.float32,-.7),(torch.float64,2.3)):
                    self.compare(sequence,torch.randn(len(sequence),15,dtype=dtype,generator=generator),
                        cap,mask,cuts=(2,len(sequence)-2),scale=scale,span=2)

    def test_cap_short_exons_internal_STOP_GC_and_minimum_intron(self):
        fixtures = (('ATGCGTTAGAGTTAGAGCTAGTAA',((0,4),(9,10),(15,16),(21,24)),4,1),
            ('ATGCGTTAGAGTTAGAGCTAGTAA',((0,4),(9,10),(15,16),(21,24)),3,0),
            ('ATGTGTTAGATGCGTTAGAATAA',((0,4),(18,23)),3,0),
            ('ATGGCTAGAAATAA',((0,3),(8,14)),2,1),
            ('ATGGTAGAAATAA',((0,3),(7,13)),2,0))
        for sequence,chain,cap,sets in fixtures:
            _,actual,_ = self.compare(sequence,torch.zeros(len(sequence),15,dtype=torch.float64),cap,
                exact_roles(candidate_roles(chain,len(sequence))),tuple(range(1,len(sequence))))
            self.assertEqual(actual,sets)
        sequence = 'ATGCGTTAGTTTAA'
        roles = exact_roles(torch.tensor([13,13,13,7,4,4,4,4,4,11,3,14,14,14]))
        message = self.initial(2); states = _state_table(2)[0]
        expected = {4:('D0',1,'C'),6:('I2',1,'C'),7:('I3',1,'C'),9:('X',2,'C'),10:('E',2,'CT')}
        for i in range(14):
            message = compiled_boundary_forward(sequence[i:i+1],torch.zeros(1,15,dtype=torch.float64),message,roles[i:i+1])
            if i+1 in expected: self.assertTrue(bool(torch.isfinite(message.log_masses[states.index(expected[i+1])])) )

    def test_unknown_and_outside_genes_are_marginalized_and_row_shift_invariant(self):
        sequence = 'ATGAAATAA'*6
        logits = torch.zeros(len(sequence),15,dtype=torch.float64)
        core = torch.ones(27,15,dtype=torch.bool); core[:3] = False; core[:3,13] = True
        core[18:27] = False; core[18:27,0] = True
        result = compiled_boundary_core_observation_nll(sequence,logits,1,0,27,core,cuts=(9,18,27,36,45))
        for key,value in (('log_partition_all',64),('log_partition_observed',16),('loss',4)):
            self.assertAlmostEqual(float(result[key]),math.log(value),places=12)
        sequence = 'ATGCGTTAGTTTAA'
        values = (torch.arange(14*15).reshape(14,15).float()%7-3).requires_grad_(True)
        unknown = torch.ones(9,15,dtype=torch.bool)
        loss = compiled_boundary_core_observation_nll(sequence,values,128,3,12,unknown,cuts=(3,12))['loss']
        self.assertEqual(float(loss),0.)
        self.assertTrue(torch.equal(torch.autograd.grad(loss,values)[0],torch.zeros_like(values)))
        shifted = (values.detach()+torch.arange(14)[:,None]*8).requires_grad_(True)
        first = compiled_boundary_full_grammar_log_partition(sequence,values,128)['log_partition']
        second = compiled_boundary_full_grammar_log_partition(sequence.lower(),shifted,128)['log_partition']
        self.assertEqual(float(first),float(second))
        torch.testing.assert_close(torch.autograd.grad(first,values)[0],torch.autograd.grad(second,shifted)[0],rtol=0,atol=0)

    def test_checkpoint_recomputation_is_invariant_and_kernels_are_nopython(self):
        sequence = 'ATGCGTTAGAGTTAGAGCTAGTAA'
        logits = (torch.arange(len(sequence)*15).reshape(-1,15).double()/43).requires_grad_(True)
        reference = None
        for span in (1,2,3,5,256):
            result = compiled_boundary_full_grammar_log_partition(sequence,logits,4,checkpoint_bp=span)
            gradient, = torch.autograd.grad(result['log_partition'],logits,retain_graph=True)
            if reference is None: reference = (float(result['log_partition']),gradient)
            self.assertEqual(float(result['log_partition']),reference[0])
            torch.testing.assert_close(gradient,reference[1],rtol=0,atol=0)
            stats = result['statistics'][0]
            rows = len([0]+list(range(span,len(sequence),span))+[len(sequence)])
            self.assertEqual(stats['checkpoint_bytes'],rows*stats['state_count']*8)
            self.assertEqual(stats['maximum_recomputed_alpha_bytes'],(min(span,len(sequence))+1)*stats['state_count']*8)
        self.assertGreater(len(_forward_messages.nopython_signatures),0)
        self.assertGreater(len(_reverse_messages.nopython_signatures),0)

    def test_global_completion_empty_identity_and_unsupported_readout(self):
        sequence = 'ATGAAATA'; logits = torch.zeros(8,15,dtype=torch.float64)
        mask = torch.ones(8,15,dtype=torch.bool); mask[:3] = False; mask[:3,13] = True
        with self.assertRaises(InfeasibleObservation):
            compiled_boundary_full_grammar_log_partition(sequence,logits,2,mask,cuts=(1,2,7))
        initial = self.initial(2)
        self.assertIs(compiled_boundary_forward('',logits[:0],initial).log_masses,initial.log_masses)
        for values,cap,cuts,span in ((torch.zeros(7,15),2,(),2),(torch.full((8,15),float('nan')),2,(),2),
                (logits,0,(),2),(logits,2,(3,3),2),(logits,2,(0,),2),(logits,2,(),0)):
            with self.assertRaises(ValueError):
                compiled_boundary_full_grammar_log_partition(sequence,values,cap,cuts=cuts,checkpoint_bp=span)
        with self.assertRaises(ValueError):
            compiled_boundary_core_observation_nll(sequence,logits,2,0,8,torch.ones(8,15))
        values = torch.zeros(1,15,dtype=torch.float64,requires_grad=True)
        message = compiled_boundary_forward('A',values,self.initial(2))
        with self.assertRaises(ValueError):
            message.log_masses[_state_table(2)[0].index(('TG',1,''))].backward()


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
