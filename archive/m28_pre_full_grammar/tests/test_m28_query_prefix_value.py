import unittest
import torch
import numpy as np
from torch import nn
from src.m28.query_prefix_value import (PrefixValueHead,START,NEXT_START,DONOR,STOP,
    root_row,feature_sums,entry_events,end_events)
from src.m28.query_prefix_value_training import prefix_labels,prefix_loss,PrefixSamples
from src.m28.query_prefix_value_decode import decode_prefix_values,complete_query_values
from src.m28.query_beam import decode_query_beam
from test_m28_query_decode import fixture


class GeometryValue(nn.Module):
    def __init__(self,kind):super().__init__();self.kind=kind
    def forward(self,initial,state,mean,geometry,events):
        if self.kind=='later_start':return geometry[:,0]*100
        return torch.where(events==STOP,geometry[:,3],geometry[:,3]*0-10)


def decode(model,features,sequence,head,cap=128,limit=4,valid=None):
    return decode_prefix_values(model,features,sequence,head,cap,4,limit,valid)


class QueryPrefixValueTests(unittest.TestCase):
    def test_joint_readout_is_geometry_slot_deduplicated_and_shared_by_both_arms(self):
        sequence='ATGAAATAA'*2;model,features=fixture(sequence,[(0,9),(9,18)])
        head=GeometryValue('later_start');a=(0,((0,9),));b=(1,((9,18),))
        values=complete_query_values(model,features,sequence,head,[b,a,a])
        self.assertEqual(set(values),{a,b});self.assertGreater(values[b],values[a])
        self.assertEqual(values,complete_query_values(model,features,sequence,head,[a,b]))
        self.assertEqual(complete_query_values(model,features,sequence,head,[]),{})

    def test_fifth_start_is_scored_before_cut_not_raw_pointer_top_four(self):
        sequence='ATGAAATAA'*5
        model,features=fixture(sequence,[(a,a+9) for a in range(0,45,9)])
        old=decode_query_beam(model,features,sequence,1,4,4)
        new=decode(model,features,sequence,GeometryValue('later_start'),cap=1)
        chain=((36,45),)
        for q,previous in zip(new,old):
            self.assertNotIn(chain,[r['chain'] for r in previous['completed']])
            self.assertIn(chain,[r['chain'] for r in q['completed']])
            self.assertEqual(q['counts']['entry_candidates'],5)

    def test_low_live_value_does_not_impose_a_future_bound_or_fake_cap_stop(self):
        sequence='ATGCGTTAGTTTAA';chain=((0,4),(9,14))
        model,features=fixture(sequence,chain)
        head=GeometryValue('more_exons')
        for q in decode(model,features,sequence,head,limit=1):
            self.assertEqual(q['completed'][0]['chain'],chain)
        for q in decode(model,features,sequence,head,cap=1):
            self.assertNotIn(chain,[r['chain'] for r in q['completed']])

    def test_actual_codon_letters_padding_softmask_and_existence(self):
        sequence='ATGC'+'GTNAG'+'TTTAA';chain=((0,4),(9,14))
        model,features=fixture(sequence,chain);head=GeometryValue('more_exons')
        expected=decode(model,features,sequence,head)
        padded=torch.cat([features,torch.full((6,4),10000.)])
        self.assertEqual(expected,decode(model,padded,sequence.lower(),head,valid=len(sequence)))
        with torch.no_grad():model.exists.bias.fill_(-100.)
        self.assertEqual([q['completed'] for q in expected],
                         [q['completed'] for q in decode(model,features,sequence,head)])
        invalid='ATGTGTTAGAAAAATAA';bad=((0,4),(9,17))
        model,features=fixture(invalid,bad)
        for q in decode(model,features,invalid,head):
            self.assertNotIn(bad,[r['chain'] for r in q['completed']])

    def test_content_is_coding_mean_and_prospective_states_match_scalar_replay(self):
        sequence='ATGCGTTAGTTTAA';model,features=fixture(sequence,[(0,4),(9,14)])
        features[:,3]=torch.arange(len(sequence));model.requires_grad_(False)
        window=model.for_window(features);sums=feature_sums(window);row=root_row(window,0)
        entry=entry_events(window,row,[0]);row=dict(row,entry=0,state=entry['states'][0])
        first=end_events(window,sums,row,[(4,'donor','C')])
        self.assertTrue(torch.allclose(entry['states'][0],window.advance(window.initial[0],'start',0)))
        self.assertTrue(torch.allclose(first['states'][0],window.advance(row['state'],'end',4)))
        row=dict(query_index=0,state=first['states'][0],chain=((0,4),),carry='C',
                 previous_end=4,coding_sum=first['coding_sums'][0],coding_bases=4)
        next_entry=entry_events(window,row,[9]);row=dict(row,entry=9,state=next_entry['states'][0])
        last=end_events(window,sums,row,[(14,'stop','')])
        actual=torch.cat([features[:4],features[9:14]]).mean(0)
        self.assertTrue(torch.allclose(last['inputs'][2][0],actual))
        self.assertEqual(int(last['coding_bases'][0]),9)
        self.assertEqual(float(first['inputs'][3][0,4]),.5)

    def test_slot_compatibility_unknown_isoform_and_callable_protection(self):
        primary=((0,9),);other=((12,21),);isoform=((0,6),(12,18))
        def node(q,event,chain=(),entry=None):
            return dict(query_index=q,event=event,chain=chain,entry=entry)
        targets=dict(positive_chains=[primary,other],protected_chains=[primary,other,isoform],
                     unknown_chain_spans=[(24,30)],callable_bases=np.ones(36,dtype=bool))
        nodes=[node(0,STOP,primary),node(0,STOP,other),node(1,STOP,primary),
               node(0,DONOR,((0,6),)),node(0,START,entry=24),node(0,START,entry=31),
               node(0,START,entry=0),node(0,NEXT_START,((0,6),),entry=12)]
        targets['callable_bases'][31]=False
        self.assertEqual(prefix_labels(nodes,{0:primary},[1,-1],targets),[1,0,-1,-1,-1,-1,1,-1])

    def test_equal_event_means_unknown_and_no_class_balancing(self):
        logits=torch.tensor([2.,0.,100.,-1.],requires_grad=True)
        labels=[1,0,-1,1];events=torch.tensor([START,START,START,STOP])
        loss=prefix_loss(logits,labels,events)
        expected=(torch.nn.functional.binary_cross_entropy_with_logits(logits[:2],torch.tensor([1.,0.]))+
                  torch.nn.functional.softplus(-logits[3]))/2
        self.assertTrue(torch.allclose(loss,expected));loss.backward()
        self.assertEqual(float(logits.grad[2]),0.)
        x=torch.randn(2,requires_grad=True);zero=prefix_loss(x,[-1,-1],torch.tensor([START,STOP]))
        zero.backward();self.assertEqual(float(zero),0.);self.assertEqual(float(x.grad.abs().sum()),0.)
        empty=torch.empty(0,requires_grad=True)
        prefix_loss(empty,[],torch.empty(0,dtype=torch.long)).backward()
        with self.assertRaises(ValueError):prefix_loss(torch.zeros(1),[.5],torch.tensor([START]))

    def test_only_new_head_has_gradient_and_samples_deduplicate_origins(self):
        model,features=fixture('ATGAAATAA',[(0,9)]);features.requires_grad_(True)
        model.requires_grad_(False);window=model.for_window(features)
        batch=entry_events(window,root_row(window,0),[0]);samples=PrefixSamples(window)
        samples.add(batch,'free');samples.add(batch,'TRAIN_teacher_or_counterfactual')
        self.assertEqual(len(samples.records()),1)
        self.assertEqual(samples.public_records()[0]['origins'],['TRAIN_teacher_or_counterfactual','free'])
        head=PrefixValueHead(4);before={n:p.detach().clone() for n,p in head.named_parameters()}
        logits=head(*samples.inputs());prefix_loss(logits,[1],torch.tensor([START])).backward()
        self.assertTrue(all(p.grad is not None and torch.isfinite(p.grad).all() for p in head.parameters()))
        self.assertGreater(sum(float(p.grad.abs().sum()) for p in head.parameters()),0)
        self.assertTrue(all(p.grad is None for p in model.parameters()));self.assertIsNone(features.grad)
        self.assertTrue(all(torch.equal(p.detach(),before[n]) for n,p in head.named_parameters()))


if __name__=='__main__':unittest.main()
