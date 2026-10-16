#!/usr/bin/env python3
"""Artificial exact-enumeration and derivative checks; no research data."""
import itertools
import json
import os
from pathlib import Path
import time
import unittest
import torch
from src.m28.full_grammar_coarse_channel import coarse_channel_nll
from src.m28.full_grammar_boundary_compiled import compiled_boundary_core_observation_nll
from src.m28.full_grammar_observed_partition import allowed_roles_from_local_targets,InfeasibleObservation
from src.m28.labels import REGION_COLUMNS

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'outputs/M28-FULL-GRAMMAR-COARSE-CHANNEL-CHECKS'
CONFIG=dict(scope='artificial_unique_single_exon_ORFs_only_no_research_inputs',
    channel='three_region_classes_correct_1_minus_epsilon_other_each_epsilon_over_two',
    fine='unchanged_exact_role',unknown='missing_likelihood_one',
    epsilon_examples=[0.1,2/3],epsilon_examples_not_fit_selection=True,
    arithmetic='FP64',enumeration_atol=1e-9,gradient_atol=1e-9,finite_difference_atol=1e-6,
    optimizer_steps=0)


def targets(n):
    return dict(seg=torch.full((n+2,),-100,dtype=torch.long),region=torch.full((n+2,),-100,dtype=torch.long))


def reference(logits,genes,start,end,t,epsilon):
    # These two DNA fixtures have exactly the stated nonoverlapping complete
    # single-exon ORFs. Enumerate null and every subset, independently of DP.
    scores=[]; weighted=[]
    for active in itertools.product((False,True),repeat=len(genes)):
        roles=torch.zeros(len(logits),dtype=torch.long)
        for keep,(a,b) in zip(active,genes):
            if keep:
                roles[a:a+3]=13; roles[a+3:b-3]=torch.arange(b-a-6)%3+1; roles[b-3:b]=14
        score=(logits[torch.arange(len(logits)),roles]-logits[:,0]).sum(); scores.append(score)
        observation=epsilon*0; supported=True
        for i in range(end-start):
            role=int(roles[start+i]); f=int(t['seg'][i]); c=int(t['region'][i])
            if f>=0:
                if role!=f: supported=False; break
            elif c>=0:
                observation=observation+(torch.log1p(-epsilon) if role in REGION_COLUMNS[c] else torch.log(epsilon/2))
        if supported: weighted.append(score+observation)
    return torch.logsumexp(torch.stack(scores),0)-torch.logsumexp(torch.stack(weighted),0)


class ChannelChecks(unittest.TestCase):
    def compare(self,sequence,genes,start,end,t,error):
        torch.manual_seed(7)
        logits=(torch.randn(len(sequence),15,dtype=torch.float64)*0.15).requires_grad_()
        epsilon=torch.tensor(error,dtype=torch.float64,requires_grad=True)
        before={k:v.clone() for k,v in t.items()}
        result=coarse_channel_nll(sequence,logits,2,start,end,t,epsilon,cuts=(5,10),checkpoint_bp=3)
        expected=reference(logits,genes,start,end,t,epsilon)
        self.assertAlmostEqual(float(result['loss']),float(expected),delta=1e-9)
        actual=torch.autograd.grad(result['loss'],(logits,epsilon),allow_unused=True)
        predicted=torch.autograd.grad(expected,(logits,epsilon),allow_unused=True)
        self.assertTrue(torch.allclose(actual[0],predicted[0],rtol=0,atol=1e-9))
        if actual[1] is not None: self.assertAlmostEqual(float(actual[1]),float(predicted[1]),delta=1e-9)
        self.assertTrue(all(torch.equal(t[k],v) for k,v in before.items()))
        self.assertGreaterEqual(float(result['loss']),-1e-9)
        self.assertLessEqual(float(actual[0].sum(1).abs().max()),1e-9)

    def test_single_gene_correct_and_incompatible_coarse(self):
        for label in (1,2):
            t=targets(16); t['seg'][3]=13; t['region'][7]=label
            self.compare('CCCATGAAATAACCCC',[(3,12)],0,16,t,0.1)

    def test_two_genes_core_only_and_row_gauge(self):
        sequence='CCCATGAAATAACCCATGCCCTAGCCC'; t=targets(21)
        t['seg'][0]=13; t['region'][16]=1; t['region'][10]=0
        self.compare(sequence,[(3,12),(15,24)],3,24,t,0.1)
        torch.manual_seed(9); x=torch.randn(len(sequence),15,dtype=torch.float64)
        shift=torch.randn(len(sequence),1,dtype=torch.float64)
        a=coarse_channel_nll(sequence,x,2,3,24,t,0.1)['loss']
        b=coarse_channel_nll(sequence,x+shift,2,3,24,t,0.1)['loss']
        self.assertAlmostEqual(float(a),float(b),delta=1e-9)

    def test_unknown_rows_have_unit_likelihood(self):
        self.compare('CCCATGAAATAACCCC',[(3,12)],0,16,targets(16),0.1)

    def test_uninformative_coarse_limit_is_not_deleted(self):
        t=targets(16); t['seg'][3]=13; t['region'][7]=2; t['region'][8]=1
        self.compare('CCCATGAAATAACCCC',[(3,12)],0,16,t,2/3)
        x=torch.zeros((16,15),dtype=torch.float64,requires_grad=True)
        value=coarse_channel_nll('CCCATGAAATAACCCC',x,2,0,16,t,2/3)['loss']
        fine_only=targets(16); fine_only['seg']=t['seg'].clone()
        old=compiled_boundary_core_observation_nll('CCCATGAAATAACCCC',x,2,0,16,
            allowed_roles_from_local_targets(fine_only,16))['loss']
        self.assertAlmostEqual(float(value-old),2*float(torch.log(torch.tensor(3.,dtype=torch.float64))),delta=1e-9)
        ga=torch.autograd.grad(value,x)[0]; gb=torch.autograd.grad(old,x)[0]
        self.assertTrue(torch.allclose(ga,gb,rtol=0,atol=1e-9))

    def test_incompatible_event_and_epsilon_finite_difference(self):
        sequence='CCCATGAAATAACCCC'; t=targets(16); t['seg'][3]=13; t['region'][7]=2
        x=torch.zeros((16,15),dtype=torch.float64,requires_grad=True)
        with self.assertRaises(InfeasibleObservation):
            compiled_boundary_core_observation_nll(sequence,x,2,0,16,allowed_roles_from_local_targets(t,16))
        e=torch.tensor(0.1,dtype=torch.float64,requires_grad=True)
        value=coarse_channel_nll(sequence,x,2,0,16,t,e)['loss']; gx,ge=torch.autograd.grad(value,(x,e))
        self.assertAlmostEqual(float(ge),-10.,delta=1e-9)
        h=1e-5
        def evaluate(y,error): return float(coarse_channel_nll(sequence,y,2,0,16,t,error)['loss'])
        self.assertAlmostEqual(float(ge),(evaluate(x,0.1+h)-evaluate(x,0.1-h))/(2*h),delta=1e-6)
        for i,j in ((3,13),(7,0),(7,2)):
            plus=x.detach().clone(); minus=x.detach().clone(); plus[i,j]+=h; minus[i,j]-=h
            self.assertAlmostEqual(float(gx[i,j]),(evaluate(plus,0.1)-evaluate(minus,0.1))/(2*h),delta=1e-6)


if __name__=='__main__':
    torch.set_num_threads(2); began=time.perf_counter()
    with (OUT/'started.json').open('x') as handle:
        json.dump(dict(job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG),handle,indent=2)
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ChannelChecks))
    summary=dict(completed=result.wasSuccessful(),job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
        tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
        seconds=time.perf_counter()-began,research_data_read=False,models_or_weights_read=False,
        original_R4_policy_changed=False,fit_error_probability_selected=False,optimizer_steps_this_run=0,
        interpretation='candidate_observation_likelihood_math_only_not_fit_or_model_performance')
    with (OUT/'summary.json').open('x') as handle: json.dump(summary,handle,indent=2,allow_nan=False)
    print(json.dumps(summary,allow_nan=False),flush=True)
    if not result.wasSuccessful(): raise SystemExit(1)
