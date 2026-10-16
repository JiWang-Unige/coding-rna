"""Changed joint encoder boundary and raw caller, using artificial inputs only."""
import copy
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import torch
import annotate_full_grammar_joint24 as A
from test_m28_emission_lattice import prefer_chain


class Tokenizer:
    def __init__(self):
        self.sequences = []
    def __call__(self,sequence,**kwargs):
        assert kwargs == dict(add_special_tokens=False,return_tensors='pt',truncation=False)
        self.sequences.append(sequence)
        return dict(input_ids=torch.arange(len(sequence)//6)[None])


class Backbone:
    def __init__(self):
        self.calls = []
    def __call__(self,input_ids,attention_mask):
        assert torch.equal(attention_mask,torch.ones_like(input_ids))
        self.calls.append(tuple(input_ids.shape))
        return SimpleNamespace(last_hidden_state=torch.ones(1,input_ids.shape[1],2048,dtype=torch.bfloat16))


class Head:
    def __init__(self):
        self.weight = torch.nn.Parameter(torch.tensor(2.))
    def __call__(self,features,dna,lengths):
        assert tuple(features.shape) == (1,4096,2048)
        assert tuple(dna.shape) == (1,24576,5)
        assert 0 < lengths[0] <= 24576
        # Backward must save native-derived input outside inference mode.
        rows = (features[:,:,0]*self.weight).repeat_interleave(6,dim=1)
        return dict(segmentation_logits=rows[:,:,None].expand(-1,-1,15))


class JointRawTests(unittest.TestCase):
    def test_complete_and_short_inputs_keep_all_tokens_and_valid_rows(self):
        for sequence in ('A','ATGCGTNAGTTTAA','ACGT'*6144):
            tokenizer,backbone,head = Tokenizer(),Backbone(),Head()
            values = A.window_logits(head,backbone,tokenizer,sequence,'cpu')
            self.assertEqual(tokenizer.sequences,[sequence+'N'*(24576-len(sequence))])
            self.assertEqual(backbone.calls,[(1,4096)])
            self.assertEqual(tuple(values.shape),(len(sequence),15))
            self.assertEqual(values.dtype,torch.float32)
            values.sum().backward()
            self.assertEqual(float(head.weight.grad),len(sequence)*15)

    def test_checkpoint_rejects_old_encoder_and_intermediate_weights(self):
        state = dict(configuration=dict(A.IDENTITY),step=9216,epoch=6,
            protocol=A.PROTOCOL,encoder_revision=A.R.F.REVISION)
        self.assertEqual(A.validate_checkpoint(state),A.IDENTITY)
        for key,value in (('step',6144),('epoch',4),('protocol',A.R.PROTOCOL)):
            bad = copy.deepcopy(state); bad[key] = value
            with self.assertRaises(ValueError):
                A.validate_checkpoint(bad)
        for key,value in (('encoder_blocks_bp',6144),('epsilon',.2),('max_exons',64)):
            bad = copy.deepcopy(state); bad['configuration'][key] = value
            with self.assertRaises(ValueError):
                A.validate_checkpoint(bad)

    def test_short_minus_codon_carry_intronic_N_and_valid_DNA_only(self):
        oriented = 'ATGCGTNAGTTTAA'
        sequence = A.R.M.rc(oriented)
        expected = prefer_chain(oriented,((0,4),(9,14))).float()
        calls = []
        def logits(model,backbone,tokenizer,dna,device):
            calls.append(dna)
            return torch.zeros(14,15) if len(calls) == 1 else expected
        with patch.object(A,'window_logits',logits):
            txs,scopes = A.annotate_contig(None,None,None,'toy14',sequence,0,'cpu')
        self.assertEqual(calls,[sequence,oriented])
        self.assertEqual(len(txs),1)
        self.assertEqual(txs[0]['strand'],'-')
        self.assertEqual(txs[0]['CDS'],[(0,5,'2'),(10,14,'0')])
        self.assertTrue(all(s['windows'] == [(0,14)] and s['encoder_condition'] == 'JOINT24576'
            and not s['reference_used'] and s['joint_whole_oriented_contig_decode'] for s in scopes))

    def test_irregular_final_anchor_ownership_both_strands_no_role_permutation(self):
        sequence = 'N'*40000
        windows = [(a,min(a+24576,len(sequence))) for a in A.R.M.starts(len(sequence))]
        self.assertEqual(windows,[(0,24576),(12288,36864),(15424,40000)])
        for strand in ('+','-'):
            matrices = [torch.arange((b-a)*15).reshape(b-a,15).float()+i*1000000
                for i,(a,b) in enumerate(windows)]
            calls = []
            def logits(model,backbone,tokenizer,dna,device):
                value = matrices[len(calls)]; calls.append(dna); return value
            with patch.object(A,'window_logits',logits):
                joined = A.contig_emissions(None,None,None,sequence,windows,strand,'cpu')
            desired = torch.empty(len(sequence),15)
            for q in range(len(sequence)):
                p = q if strand == '+' else len(sequence)-1-q
                owner = min((i for i,(a,b) in enumerate(windows) if a <= p < b),
                    key=lambda i:(abs(2*p+1-sum(windows[i])),windows[i][0]))
                a,b = windows[owner]
                desired[q] = matrices[owner][p-a if strand == '+' else b-1-p]
            self.assertTrue(torch.equal(joined['local_potentials'],desired))
            self.assertEqual(len(calls),3)


if __name__ == '__main__':
    torch.set_num_threads(2)
    unittest.main()
