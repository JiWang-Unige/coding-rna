import json,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import torch
from torch import nn
from src.m28.query_prefix_pair import paired_query_records,ARMS
from src.m28.query_prefix_value_decode import complete_query_values
from src.m28.query_decode import decode_query_chains
from src.m28.query_beam import decode_query_beam
from src.m28.query_beam_pool import completed_union
from src.m28.assembly import assemble_predictions
from test_m28_query_decode import fixture
from test_m28_query_prefix_value import GeometryValue

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import infer_query_prefix_pair as P

LIMITS=dict(max_exons=128,beam_width=4,completed_limit=4)


class ZeroValue(nn.Module):
    def forward(self,initial,state,mean,geometry,events):return geometry[:,0]*0


class SyntheticGenerator:
    def __init__(self,model):self.query=model;self.calls=0
    def __call__(self,features,dna,valid):
        self.calls+=1
        if len(valid)!=1 or tuple(dna.shape)!=(1,P.I.M.WIDTH,5):raise AssertionError('Raw window contract')
        return dict(features=features)


class QueryPrefixPairTests(unittest.TestCase):
    def test_late_is_literal_original_union_and_shared_readout_occurs_once(self):
        sequence='ATGAAATAA'*5;model,features=fixture(sequence,[(a,a+9) for a in range(0,45,9)])
        head=GeometryValue('later_start')
        with patch('src.m28.query_prefix_pair.complete_query_values',wraps=complete_query_values) as readout:
            pair=paired_query_records(model,features,sequence,head,0,'synthetic',**LIMITS)
        self.assertEqual(readout.call_count,1)
        greedy=decode_query_chains(model,features,sequence,128)
        beam=decode_query_beam(model,features,sequence,128,4,4)
        old,_=completed_union(greedy,beam)
        key=lambda c:tuple(tuple(e) for e in c)
        self.assertEqual({key(c) for c in old},set(pair['arms']['late']['chains']))
        self.assertIn(((36,45),),pair['arms']['early']['chains'])
        self.assertNotIn(((36,45),),pair['arms']['late']['chains'])
        values={(r['query_index'],key(r['chain'])):r['complete_value'] for r in pair['joint_query_chain_values']}
        self.assertEqual(len(values),len(pair['joint_query_chain_values']))
        for arm in ARMS:
            record=pair['arms'][arm]
            for c,v,qs in zip(record['chains'],record['score'],record['query_indices']):
                self.assertEqual(v,max(values[q,c] for q in qs))
            self.assertGreater(record['counts']['unique_query_geometry'],record['counts']['free_candidates'])
        self.assertFalse(pair['reference_used']);self.assertTrue(pair['common_query_chain_values_identical'])

    def test_zero_values_remain_in_pool_but_do_not_emit_genes(self):
        sequence='ATGAAATAA';model,features=fixture(sequence,[(0,9)])
        pair=paired_query_records(model,features,sequence,ZeroValue(),0,'synthetic',**LIMITS)
        for arm in ARMS:
            r=pair['arms'][arm]
            self.assertEqual(r['chains'],[((0,9),)]);self.assertEqual(r['score'],[0.])
            self.assertEqual(r['counts']['nonpositive_candidates'],1)
            self.assertEqual(assemble_predictions([r],[(0,9)],'+')['predictions'],[])

    def test_no_atg_never_calls_complete_head_or_invents_a_gene(self):
        sequence='CCCAAATAA';model,features=fixture(sequence,[(0,9)])
        class Unexpected(nn.Module):
            def forward(self,*args):raise AssertionError('No legal start, no values')
        pair=paired_query_records(model,features,sequence,Unexpected(),0,'synthetic',**LIMITS)
        self.assertEqual(pair['joint_query_chain_values'],[])
        self.assertTrue(all(not pair['arms'][a]['chains'] for a in ARMS))

    def test_padding_and_softmask_do_not_change_either_arm(self):
        sequence='ATGCGTTAGTTTAA';model,features=fixture(sequence,[(0,4),(9,14)])
        head=GeometryValue('more_exons')
        expected=paired_query_records(model,features,sequence,head,0,'synthetic',**LIMITS)
        padded=torch.cat([features,torch.full((7,4),10000.)])
        actual=paired_query_records(model,padded,sequence.lower(),head,0,'synthetic',**LIMITS,valid_bases=len(sequence))
        self.assertEqual(expected,actual)

    def test_raw_helper_reuses_one_forward_and_preserves_minus_coordinates_phases(self):
        oriented='ATGCGTTAGTTTAA';chain=((0,4),(9,14));model,features=fixture(oriented,chain)
        generator=SyntheticGenerator(model);head=GeometryValue('more_exons')
        raw=P.I.M.rc(oriented).lower();seen=[]
        def encode(backbone,tokenizer,sequence,device):seen.append(sequence);return features[None]
        with tempfile.TemporaryDirectory() as directory,patch.object(P.I,'encode_window',side_effect=encode):
            transcripts,scopes=P.annotate_contig(generator,head,None,None,'synthetic','chr',raw,0,Path(directory),'cpu',LIMITS)
            self.assertEqual(generator.calls,2);self.assertEqual(len(scopes),2)
            self.assertEqual(seen,[raw.upper(),oriented])
            self.assertEqual([s['strand'] for s in scopes],['+','-'])
            for arm in ARMS:
                minus=[r for r in transcripts[arm] if r['strand']=='-']
                self.assertEqual(len(minus),1)
                self.assertEqual(minus[0]['CDS'],[(0,5,'2'),(10,14,'0')])
                self.assertIn('M28_A9_'+arm,minus[0]['id'])
            for scope in scopes:
                records=[json.loads(x) for x in (Path(directory)/scope['file']).read_text().splitlines()]
                self.assertEqual(len(records),1)
                self.assertFalse(scope['reference_used']);self.assertFalse(scope['labels_or_feature_cache_used'])
                for arm in ARMS:
                    rebuilt=assemble_predictions(P.arm_records(Path(directory)/scope['file'],arm),scope['windows'],scope['strand'])
                    self.assertEqual(rebuilt['counts'],scope['assembly_counts'][arm])


if __name__=='__main__':unittest.main()
