import copy
import itertools
import unittest
import torch
from torch.nn import functional as F
from src.m28.core import M28Core
from src.m28.labels import build_targets
from src.m28.training import one_hot
from src.m28.query_model import (GeneQueryModel, balanced_query_components,
                                unmatched_query_policy, joint_query_losses)

COEFFICIENTS = (0.7, 1.3, 0.2)  # Authored fixture, not real-fit choices.
CHAINS = (((0,4),(9,14)), ((3,11),), ((5,17),))


def fixture():
    torch.manual_seed(13)
    model = GeneQueryModel(12, 3, hidden=8)
    tokens = torch.randn(1,4,12,requires_grad=True)
    dna = one_hot('ATGCGTTAGTTTAACCCC',24)[None].requires_grad_()
    outputs = model(tokens,dna,[18])
    return model, tokens, dna, outputs


def explicit_total(window, chains, assignment, unmatched):
    labels = [unmatched]*window.model.num_queries
    assigned = {}
    for gene, query in enumerate(assignment):
        labels[query] = 1
        assigned[query] = chains[gene]
    parts = window.teacher_components(labels, assigned)
    # Independently enumerate full positive/negative BCE terms, not matcher delta.
    positive = [window.query_logits[i] for i,v in enumerate(labels) if v == 1]
    negative = [window.query_logits[i] for i,v in enumerate(labels) if v == 0]
    means = []
    if positive: means.append(F.softplus(-torch.stack(positive)).mean())
    if negative: means.append(F.softplus(torch.stack(negative)).mean())
    existence = torch.stack(means).mean() if means else window.query_logits.sum()*0
    return (COEFFICIENTS[0]*existence + COEFFICIENTS[1]*parts['boundary_ce'] +
            COEFFICIENTS[2]*parts['termination_bce'])


def targets():
    w = {'width':24,'valid_bases':18,'start':0,'end':18,'strand':'+',
         'overlapping_transcript_ids':['tx'],'positive_ids':['tx'],
         'fine_label_unknown_intervals':[],'chain_unknown_intervals':[],
         'no_annotated_gene_either_strand':False}
    tx = {'id':'tx','gene_id':'g','CDS':[(0,4,'0'),(9,14,'2')],
          'training_primary':True,'CDS_complete_supported':True}
    return build_targets(w,{'tx':tx},'ATGCGTTAGTTTAACCCC')


class QueryJointTests(unittest.TestCase):
    def test_common_initialization_matches_c0_and_padding_is_not_supervision(self):
        torch.manual_seed(0)
        control = M28Core('C0',12,8)
        torch.manual_seed(0)
        query = GeneQueryModel(12,3,8)
        for prefix in ('shared','segmentation'):
            actual = getattr(query,prefix).state_dict()
            expected = getattr(control,prefix).state_dict()
            self.assertTrue(all(torch.equal(actual[k],expected[k]) for k in expected))
        _,tokens,dna,_ = fixture()
        with torch.no_grad():
            baseline = query(tokens,dna,[18])
            altered_tokens,altered_dna = tokens.clone(),dna.clone()
            altered_tokens[:,3] = 10000
            altered_dna[:,18:] = 10000
            altered = query(altered_tokens,altered_dna,[18])
            self.assertTrue(torch.equal(baseline['features'],altered['features']))
            self.assertTrue(torch.equal(baseline['segmentation_logits'],altered['segmentation_logits']))
        self.assertEqual(unmatched_query_policy(targets(),18)['unmatched_label'],0)

    def test_query_only_gradient_reaches_shared_but_not_external_cache_or_unknown_rows(self):
        model,tokens,dna,outputs = fixture()
        window = model.query_window(outputs['features'][0],18)
        result = balanced_query_components(window,CHAINS[:1],COEFFICIENTS,-1)
        self.assertFalse(result['assignment_costs'].requires_grad)
        result['total'].backward()
        self.assertIsNone(tokens.grad)
        self.assertIsNone(dna.grad)
        self.assertIsNone(model.segmentation.weight.grad)
        self.assertGreater(sum(float(p.grad.abs().sum()) for p in model.shared.parameters()
                               if p.grad is not None),0)
        self.assertTrue(all(torch.isfinite(p.grad).all() for p in model.parameters()
                            if p.grad is not None))
        chosen = result['gene_to_query'][0]
        for q in set(range(3))-{chosen}:
            self.assertEqual(result['query_labels'][q],-1)
            self.assertEqual(float(model.query.queries.grad[q].abs().sum()),0)

    def test_balanced_assignment_matches_exhaustive_whole_window_minimum(self):
        model,_,_,outputs = fixture()
        with torch.no_grad():
            for n in (1,2,3):
                for unmatched in (-1,0):
                    window = model.query_window(outputs['features'][0],18)
                    result = balanced_query_components(window,CHAINS[:n],COEFFICIENTS,unmatched)
                    alternatives = [float(explicit_total(window,CHAINS[:n],p,unmatched))
                                    for p in itertools.permutations(range(3),n)]
                    self.assertAlmostEqual(float(result['total']),min(alternatives),delta=1e-6)
                    self.assertEqual(len(set(result['gene_to_query'])),n)
                    self.assertEqual(sorted(result['positive_chains'].values()),sorted(CHAINS[:n]))
                    reverse = balanced_query_components(window,CHAINS[:n][::-1],COEFFICIENTS,unmatched)
                    self.assertAlmostEqual(float(result['total']),float(reverse['total']),delta=1e-6)

    def test_cost_plus_null_baseline_equals_each_full_balanced_assignment(self):
        model,_,_,outputs = fixture()
        with torch.no_grad():
            for n in (1,2,3):
                for unmatched in (-1,0):
                    window = model.query_window(outputs['features'][0],18)
                    result = balanced_query_components(window,CHAINS[:n],COEFFICIENTS,unmatched)
                    baseline = (COEFFICIENTS[0]*0.5*F.softplus(window.query_logits).sum()/(3-n)
                                if unmatched == 0 and n < 3 else window.query_logits.sum()*0)
                    for p in itertools.permutations(range(3),n):
                        cost = sum(result['assignment_costs'][gene,q] for gene,q in enumerate(p))
                        actual = explicit_total(window,CHAINS[:n],p,unmatched)
                        self.assertAlmostEqual(float(cost+baseline),float(actual),delta=1e-6)

    def test_empty_window_has_only_known_negative_or_differentiable_unknown_loss(self):
        for unmatched in (-1,0):
            model,tokens,dna,outputs = fixture()
            window = model.query_window(outputs['features'][0],18)
            result = balanced_query_components(window,[],COEFFICIENTS,unmatched)
            self.assertEqual(result['gene_to_query'],())
            self.assertEqual(result['positive_chains'],{})
            self.assertEqual(result['query_labels'],(unmatched,)*3)
            self.assertEqual(float(result['components']['boundary_ce']),0)
            self.assertEqual(float(result['components']['termination_bce']),0)
            if unmatched == 0:
                self.assertTrue(torch.allclose(result['components']['existence_bce'],
                                               F.softplus(window.query_logits).mean(),atol=1e-7,rtol=0))
            else:
                self.assertEqual(float(result['total']),0)
                result['total'].backward()
                self.assertEqual(sum(float(p.grad.abs().sum()) for p in model.parameters()
                                     if p.grad is not None),0)
                self.assertIsNone(tokens.grad)
                self.assertIsNone(dna.grad)

    def test_full_window_protection_keeps_positive_targets_and_ignores_padding(self):
        base = targets()
        positive = copy.deepcopy(base['positive_chains'])
        for kind in ('unknown','alternative','noncallable'):
            value = copy.deepcopy(base)
            if kind == 'unknown': value['unknown_chain_spans'] = [(15,17)]
            if kind == 'alternative': value['protected_chains'].append([(2,10)])
            if kind == 'noncallable': value['callable_bases'][16] = False
            result = unmatched_query_policy(value,18)
            self.assertEqual(result['unmatched_label'],-1)
            self.assertEqual(len(result['protection_reasons']),1)
            self.assertEqual(value['positive_chains'],positive)
        self.assertTrue(not base['callable_bases'][18:].any())
        self.assertEqual(unmatched_query_policy(base,18)['unmatched_label'],0)
        with self.assertRaisesRegex(ValueError,'valid bases'):
            unmatched_query_policy(base,25)

    def test_joint_local_and_query_loss_preserves_importance_weight(self):
        model,_,_,outputs = fixture()
        value = targets()
        value['unknown_chain_spans'] = [(15,17)]
        coefficients = (0.8,0.7,1.3,0.2)
        first = joint_query_losses(model,outputs,value,18,0.25,coefficients)
        second = joint_query_losses(model,outputs,value,18,2.,coefficients)
        self.assertEqual(first['unmatched_label'],-1)
        self.assertEqual(first['query']['components']['positive_queries'],1)
        self.assertTrue(torch.equal(second['total'],8*first['total']))
        self.assertTrue(torch.equal(first['total'],
            0.25*sum(c*first['losses'][k] for c,k in zip(coefficients,first['losses']))))
        second['total'].backward()
        self.assertGreater(float(model.segmentation.weight.grad.abs().sum()),0)
        self.assertGreater(float(model.query.exists.weight.grad.abs().sum()),0)

    def test_capacity_invalid_targets_and_nonfinite_cost_stop_without_filtering(self):
        model,_,_,outputs = fixture()
        window = model.query_window(outputs['features'][0],18)
        with self.assertRaisesRegex(ValueError,'capacity'):
            balanced_query_components(window,CHAINS+(((1,5),),),COEFFICIENTS,-1)
        with self.assertRaisesRegex(ValueError,'ordered oriented'):
            balanced_query_components(window,[((0,19),)],COEFFICIENTS,-1)
        with self.assertRaisesRegex(ValueError,'finite nonnegative'):
            balanced_query_components(window,CHAINS,(1,float('nan'),1),-1)
        with torch.no_grad():
            model.query.exists.weight.fill_(float('nan'))
        window = model.query_window(outputs['features'][0],18)
        with self.assertRaisesRegex(ValueError,'Nonfinite balanced'):
            balanced_query_components(window,CHAINS,COEFFICIENTS,0)


if __name__ == '__main__':
    unittest.main()
