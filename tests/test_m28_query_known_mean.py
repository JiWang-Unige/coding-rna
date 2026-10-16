"""A3 objective/matching/identity checks only; no real TRAIN or DEV inputs."""
import copy
import itertools
import sys
import unittest
from pathlib import Path
from unittest.mock import patch
import torch
from torch.nn import functional as F
from src.m28.query_collaboration import CollaborativeGeneQueryModel
from src.m28.query_model import balanced_query_components
from src.m28.query_known_mean import known_query_components, joint_known_query_losses
from src.m28.training import one_hot

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import run_query_known_mean as R
import run_query_collaboration as A2
import evaluate_query_known_mean as E


def fixture():
    torch.manual_seed(13)
    model = CollaborativeGeneQueryModel(12,4,8)
    tokens = torch.randn(1,4,12,requires_grad=True)
    dna = one_hot('ATGAAATAA',24)[None].requires_grad_()
    outputs = model(tokens,dna,[9])
    return model,tokens,dna,outputs


class KnownMeanTests(unittest.TestCase):
    def test_assignment_minimizes_exact_whole_objective_for_unequal_and_full_classes(self):
        model,_,_,outputs = fixture()
        window = model.query_window(outputs['features'][0],9)
        all_chains = (((0,3),),((1,6),),((3,9),),((0,4),(5,9)))
        coefficients = (.7,1.2,.4)
        with torch.no_grad():
            for n in (1,2,4):
                for unmatched in (0,-1):
                    chains = all_chains[:n]
                    result = known_query_components(window,chains,coefficients,unmatched)
                    choices = []
                    for assignment in itertools.permutations(range(4),n):
                        labels = [unmatched]*4
                        assigned = {}
                        for chain,index in zip(chains,assignment):
                            labels[index] = 1; assigned[index] = chain
                        parts = window.teacher_components(labels,assigned)
                        choices.append(float(sum(c*parts[k] for c,k in zip(
                            coefficients,('existence_bce','boundary_ce','termination_bce')))))
                    self.assertAlmostEqual(float(result['total']),min(choices),delta=1e-6)
                    self.assertEqual(len(set(result['gene_to_query'])),n)

    def test_known_mean_changes_direct_class_weights_not_operating_point(self):
        model,_,_,outputs = fixture()
        window = model.query_window(outputs['features'][0],9)
        window.query_logits = torch.zeros(4,requires_grad=True)
        result = known_query_components(window,(((0,9),),),(1.,0.,0.),0)
        result['total'].backward()
        chosen = result['gene_to_query'][0]
        expected = torch.full((4,),.125); expected[chosen] = -.125
        torch.testing.assert_close(window.query_logits.grad,expected)
        model,_,_,outputs = fixture()
        window = model.query_window(outputs['features'][0],9)
        window.query_logits = torch.zeros(4,requires_grad=True)
        old = balanced_query_components(window,(((0,9),),),(1.,0.,0.),0)
        old['total'].backward()
        expected = torch.full((4,),1/12); expected[old['gene_to_query'][0]] = -.25
        torch.testing.assert_close(window.query_logits.grad,expected)
        self.assertEqual(R.CONFIG['existence_logit_threshold'],0.)

    def test_unknown_has_no_direct_logit_loss_but_shared_collaboration_remains_trainable(self):
        model,tokens,dna,outputs = fixture()
        window = model.query_window(outputs['features'][0],9)
        window.query_logits.retain_grad()
        result = known_query_components(window,(((0,9),),),(1.,1.,1.),-1)
        result['total'].backward()
        chosen = result['gene_to_query'][0]
        unknown = sorted(set(range(4))-{chosen})
        self.assertTrue(all(result['query_labels'][i] == -1 for i in unknown))
        self.assertTrue(all(float(window.query_logits.grad[i]) == 0 for i in unknown))
        self.assertGreater(float(model.query.queries.grad[unknown].abs().sum()),0)
        self.assertIsNone(tokens.grad); self.assertIsNone(dna.grad)
        for module in (model.shared,model.query.query_attention,model.query.query_norm):
            self.assertGreater(sum(float(p.grad.abs().sum()) for p in module.parameters() if p.grad is not None),0)
        self.assertTrue(all(torch.isfinite(p.grad).all() for p in model.parameters() if p.grad is not None))

    def test_empty_positive_known_negative_and_all_unknown_are_distinct(self):
        model,_,_,outputs = fixture()
        window = model.query_window(outputs['features'][0],9)
        for unmatched in (0,-1):
            result = known_query_components(window,(),(1.,1.,1.),unmatched)
            expected = F.softplus(window.query_logits).mean() if unmatched == 0 else window.query_logits.sum()*0
            torch.testing.assert_close(result['total'],expected)
            self.assertEqual(result['query_labels'],(unmatched,)*4)
            self.assertEqual(result['assignment_costs'].shape,(0,4))
        known_query_components(window,(),(1.,1.,1.),-1)['total'].backward()
        self.assertTrue(all(float(p.grad.abs().sum()) == 0 for p in model.parameters() if p.grad is not None))

    def test_capacity_and_explicit_unknown_policy_are_not_relaxed(self):
        model,_,_,outputs = fixture()
        window = model.query_window(outputs['features'][0],9)
        with self.assertRaises(ValueError):
            known_query_components(window,(((0,9),),)*5,(1.,1.,1.),0)
        with self.assertRaises(ValueError):
            known_query_components(window,(),(1.,1.,1.),1)
        targets = {'positive_chains':(((0,9),),),'protected_chains':(),
                   'unknown_chain_spans':[(1,2)],'callable_bases':torch.ones(9,dtype=torch.bool)}
        with patch('src.m28.query_known_mean.local_loss',return_value=outputs['segmentation_logits'].sum()*0):
            result = joint_known_query_losses(model,outputs,targets,9,2.,(1.,1.,1.,1.))
        self.assertEqual(result['unmatched_label'],-1)
        self.assertEqual(result['query']['query_labels'].count(-1),3)
        torch.testing.assert_close(result['total'],2*result['query']['total'])

    def test_runner_routes_R4_targets_to_new_joint_loss(self):
        model,_,_,outputs = fixture()
        window = {'valid_bases':9}
        targets = object()
        expected = {'total':torch.tensor(1.)}
        with patch.object(R.Q,'build_targets',return_value=targets) as build, \
             patch.object(R,'joint_known_query_losses',return_value=expected) as readout:
            result = R.supervised_readout(model,outputs,window,{},'ATGAAATAA',.7)
        self.assertIs(result,expected)
        build.assert_called_once_with(window,{},'ATGAAATAA')
        readout.assert_called_once_with(model,outputs,targets,9,.7,R.CONFIG['loss_coefficients'])

    def test_complete_A2_initialization_and_parameter_count_are_preserved(self):
        with patch.object(torch.nn.Module,'to',lambda self,*args,**kwargs:self):
            old,new = A2.new_model(),R.new_model()
        self.assertTrue(all(torch.equal(value,new.state_dict()[key]) for key,value in old.state_dict().items()))
        self.assertEqual(set(old.state_dict()),set(new.state_dict()))
        self.assertEqual(sum(p.numel() for p in new.parameters()),649457)
        self.assertEqual({k:v for k,v in R.CONFIG.items() if k != 'existence_reduction'},A2.CONFIG)

    def test_A3_identity_and_complete_grid_reject_old_objective_or_partial_DEV(self):
        fit = {'arm':'A3','completed':True,'configuration':copy.deepcopy(R.CONFIG),
               'feature_revision':R.F.REVISION,'optimizer_steps':4608,
               'final_checkpoint':'outputs/M28-A3-FIRST-FIT/step_004608.pt'}
        R.check_last_identity(fit)
        for changes in ({'arm':'A2'},{'configuration':A2.CONFIG},{'optimizer_steps':3072}):
            with self.assertRaises(ValueError): R.check_last_identity({**fit,**changes})
        infer = {**fit,'inference_complete':True,'windows':8690,'reference_used':False,
                 'test_setaria_used':False,'checkpoint':fit['final_checkpoint'],'scopes':[]}
        for species,seqid,length,_ in E.E.S.SCOPE:
            for strand in ('+','-'):
                infer['scopes'].append({'species':species,'seqid':seqid,'strand':strand,
                    'windows':[[a,min(a+R.M.WIDTH,length)] for a in R.M.starts(length)]})
        E.validate_complete(fit,infer)
        broken = copy.deepcopy(infer); broken['scopes'][0]['windows'].pop()
        with self.assertRaises(ValueError): E.validate_complete(fit,broken)
        with self.assertRaises(ValueError): E.validate_complete(fit,{**infer,'reference_used':True})


if __name__ == '__main__':
    torch.set_num_threads(1)
    unittest.main()
