"""Focused A2 synthetic wiring and frozen evaluation tests; no real input CLI."""
import copy
import itertools
import sys
import unittest
from pathlib import Path
from unittest.mock import patch
import torch
from torch.nn import functional as F
from src.m28.query_model import GeneQueryModel, balanced_query_components
from src.m28.query_pointer import QueryPointerWindow
from src.m28.query_collaboration import CollaborativeGeneQueryModel, CollaborativeQueryWindow
from src.m28.query_decode import decode_query_chains
from src.m28.query_output import query_window_candidates
from src.m28.training import one_hot

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import run_query_collaboration as Q
import evaluate_query_collaboration as E


def fixture():
    torch.manual_seed(13)
    model = CollaborativeGeneQueryModel(12,3,8)
    tokens = torch.randn(1,4,12,requires_grad=True)
    dna = one_hot('ATGAAATAA',24)[None].requires_grad_()
    outputs = model(tokens,dna,[9])
    return model,tokens,dna,outputs


class CollaborationTests(unittest.TestCase):
    def test_full_frozen_initialization_preserves_every_A_parameter(self):
        torch.manual_seed(0)
        old = GeneQueryModel(2048,16,128)
        torch.manual_seed(0)
        new = CollaborativeGeneQueryModel(2048,16,128)
        old_state,new_state = old.state_dict(),new.state_dict()
        self.assertTrue(all(torch.equal(value,new_state[key]) for key,value in old_state.items()))
        self.assertEqual(set(new_state)-set(old_state),{
            'query.query_norm.weight','query.query_norm.bias',
            'query.query_attention.in_proj_weight','query.query_attention.in_proj_bias',
            'query.query_attention.out_proj.weight','query.query_attention.out_proj.bias'})
        self.assertEqual(sum(p.numel() for p in new.parameters()),649457)
        self.assertEqual(new.query.query_attention.num_heads,4)
        self.assertEqual(new.query.query_attention.dropout,0.)

    def test_residual_pre_norm_and_zero_attention_recover_original_state(self):
        model,_,_,outputs = fixture()
        h = outputs['features'][0]
        old = QueryPointerWindow(model.query,h,9)
        new = model.query_window(h,9)
        q = model.query.query_norm(old.initial)[None]
        mixed,_ = model.query.query_attention(q,q,q,need_weights=False)
        torch.testing.assert_close(new.initial,old.initial+mixed[0],rtol=0,atol=0)
        torch.testing.assert_close(new.query_logits,model.query.exists(new.initial).squeeze(-1))
        self.assertFalse(torch.equal(new.initial,old.initial))
        with torch.no_grad():
            model.query.query_attention.out_proj.weight.zero_()
            model.query.query_attention.out_proj.bias.zero_()
        recovered = model.query_window(h,9)
        torch.testing.assert_close(recovered.initial,old.initial,rtol=0,atol=0)
        torch.testing.assert_close(recovered.query_logits,old.query_logits,rtol=0,atol=0)

    def test_off_diagonal_query_communication_exists_without_reference_masks(self):
        model,_,_,outputs = fixture()
        h = outputs['features'][0]
        old = QueryPointerWindow(model.query,h,9).initial.detach().clone()
        before = model.query_window(h,9).initial.detach().clone()
        with torch.no_grad():
            model.query.queries[2,0] += 1
        after = model.query_window(h,9).initial
        independent = QueryPointerWindow(model.query,h,9).initial
        torch.testing.assert_close(old[:2],independent[:2],rtol=0,atol=0)
        self.assertGreater(float((before[:2]-after[:2]).abs().sum()),0)

    def test_unknown_direct_loss_absent_but_indirect_gradients_reach_attention_and_shared(self):
        model,tokens,dna,outputs = fixture()
        window = model.query_window(outputs['features'][0],9)
        window.query_logits.retain_grad()
        result = balanced_query_components(window,(((0,9),),),(1.,1.,1.),-1)
        result['total'].backward()
        self.assertIsNone(tokens.grad); self.assertIsNone(dna.grad)
        self.assertIsNone(model.segmentation.weight.grad)
        chosen = result['gene_to_query'][0]
        unknown = sorted(set(range(3))-{chosen})
        self.assertTrue(all(result['query_labels'][i] == -1 for i in unknown))
        self.assertTrue(all(float(window.query_logits.grad[i]) == 0 for i in unknown))
        self.assertGreater(float(model.query.queries.grad[unknown].abs().sum()),0)
        for module in (model.shared,model.query.query_norm,model.query.query_attention):
            self.assertGreater(sum(float(p.grad.abs().sum()) for p in module.parameters() if p.grad is not None),0)
        self.assertTrue(all(torch.isfinite(p.grad).all() for p in model.parameters() if p.grad is not None))

    def test_balanced_matching_still_minimizes_whole_window_objective(self):
        model,_,_,outputs = fixture()
        window = model.query_window(outputs['features'][0],9)
        chains = (((0,4),),((3,9),))
        with torch.no_grad():
            for unmatched in (0,-1):
                result = balanced_query_components(window,chains,(1.,1.,1.),unmatched)
                alternatives = []
                for assignment in itertools.permutations(range(3),2):
                    labels = [unmatched]*3
                    assigned = {}
                    for chain,index in zip(chains,assignment):
                        labels[index] = 1; assigned[index] = chain
                    parts = window.teacher_components(labels,assigned)
                    means = [F.softplus(-window.query_logits[list(assignment)]).mean()]
                    negative = [i for i,v in enumerate(labels) if v == 0]
                    if negative: means.append(F.softplus(window.query_logits[negative]).mean())
                    alternatives.append(float(torch.stack(means).mean()+parts['boundary_ce']+parts['termination_bce']))
                self.assertAlmostEqual(float(result['total']),min(alternatives),delta=1e-6)

    def test_teacher_and_free_decoder_use_same_collaborative_initial_states_and_logits(self):
        model,_,_,outputs = fixture()
        h = outputs['features'][0]
        train = model.query_window(h,9)
        self.assertIsInstance(train,CollaborativeQueryWindow)
        self.assertIsInstance(model.query.for_window(h,9),CollaborativeQueryWindow)
        def capture(window,index,*args):
            self.assertIsInstance(window,CollaborativeQueryWindow)
            torch.testing.assert_close(window.initial[index],train.initial[index])
            torch.testing.assert_close(window.query_logits[index],train.query_logits[index])
            return {'query_index':index}
        with patch('src.m28.query_decode._decode_query',side_effect=capture):
            self.assertEqual(len(decode_query_chains(model.query,h,'ATGAAATAA',128,9)),3)
        free = query_window_candidates(model.query,h,'ATGAAATAA',0,'synthetic',128,0.,9)
        self.assertEqual(free['counts']['queries'],3)
        self.assertFalse(free['reference_used'])
        for row in free['raw_queries']:
            self.assertAlmostEqual(row['query_logit'],float(train.query_logits[row['query_index']]),places=6)

    def test_A2_identity_rejects_A_or_wrong_last_step(self):
        valid = {'arm':'A2','completed':True,'configuration':copy.deepcopy(Q.CONFIG),
                 'feature_revision':Q.F.REVISION,'optimizer_steps':4608}
        Q.check_last_identity(valid)
        for key,value in (('arm','A'),('optimizer_steps',3072),('completed',False)):
            invalid = {**valid,key:value}
            with self.assertRaises(ValueError): Q.check_last_identity(invalid)

    def test_full_A2_inference_identity_requires_the_complete_frozen_grid(self):
        fitted = {'arm':'A2','completed':True,'configuration':copy.deepcopy(Q.CONFIG),
                  'feature_revision':Q.F.REVISION,'optimizer_steps':4608,
                  'final_checkpoint':'outputs/M28-A2-FIRST-FIT/step_004608.pt'}
        decoded = {**fitted,'inference_complete':True,'windows':8690,
                   'reference_used':False,'test_setaria_used':False,
                   'checkpoint':fitted['final_checkpoint'],'scopes':[]}
        for species,seqid,length,_ in E.E.S.SCOPE:
            for strand in ('+','-'):
                decoded['scopes'].append({'species':species,'seqid':seqid,'strand':strand,
                    'windows':[[a,min(a+Q.M.WIDTH,length)] for a in Q.M.starts(length)]})
        E.validate_complete(fitted,decoded)
        broken = copy.deepcopy(decoded)
        broken['scopes'][0]['windows'].pop()
        with self.assertRaises(ValueError): E.validate_complete(fitted,broken)
        with self.assertRaises(ValueError): E.validate_complete(fitted,{**decoded,'arm':'A'})

    def test_gate_uses_A_for_F1_and_C0_for_background_not_rounded_metrics(self):
        def rows(tp,pred,span):
            return {s:{'exact_chain':{'tp':tp,'predicted':pred,'reference':n,'f1':0.5},
                       'background':{'bp':1000,'predicted_span_bp':span}} for s,n in E.REF_COUNTS.items()}
        old,c0 = rows(10,20,30),rows(5,20,10)
        self.assertTrue(E.continuation_gate(rows(11,20,10),old,c0)['passed'])
        self.assertFalse(E.continuation_gate(rows(10,20,10),old,c0)['passed'])
        self.assertFalse(E.continuation_gate(rows(11,20,20),old,c0)['passed'])
        wrong = rows(11,20,10)
        wrong['oryza_sativa']['background']['bp'] = 999
        with self.assertRaises(ValueError): E.continuation_gate(wrong,old,c0)


if __name__ == '__main__':
    torch.set_num_threads(1)
    unittest.main()
