import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import torch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import annotate_query as A
import evaluate_query as V
import run_query as Q
from src.m28.query_model import GeneQueryModel
from src.m28.export import write_gff3


class Tokenizer:
    def __call__(self, sequence, **kwargs):
        assert kwargs == {'add_special_tokens':False,'return_tensors':'pt','truncation':False}
        values = {'A':1,'C':2,'G':3,'T':4,'N':0}
        return {'input_ids':torch.tensor([[values[sequence[i]] for i in range(0,len(sequence),6)]])}


class Backbone:
    def __init__(self): self.calls = []
    def __call__(self, input_ids, attention_mask):
        self.calls.append(tuple(input_ids.shape))
        assert torch.equal(attention_mask,torch.ones_like(input_ids))
        return SimpleNamespace(last_hidden_state=input_ids.float()[...,None].expand(-1,-1,2048))


def complete_pair():
    fitted = {'arm':'A','configuration':copy.deepcopy(Q.CONFIG),'completed':True,
              'feature_revision':Q.F.REVISION,'optimizer_steps':4608,'final_checkpoint':'synthetic/step_004608.pt'}
    scopes = [{'species':s,'seqid':q,'strand':strand,
               'windows':[[a,min(a+Q.M.WIDTH,length)] for a in Q.M.starts(length)]}
              for s,q,length,_ in V.E.S.SCOPE for strand in ('+','-')]
    decoded = {**fitted,'inference_complete':True,'windows':8690,'reference_used':False,
               'test_setaria_used':False,'checkpoint':fitted['final_checkpoint'],'scopes':scopes}
    return fitted,decoded


def counts():
    return {s:{'exact_chain':{'tp':100,'predicted':1000,'reference':n,'f1':.99},
               'background':{'bp':1000,'predicted_span_bp':20}}
            for s,n in V.REF_COUNTS.items()}


class QueryClosureTests(unittest.TestCase):
    def test_fasta_lowercase_ambiguity_multi_record_and_duplicate_empty_failures(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'artificial.fa'
            path.write_text('>first description\natgaaataary\n>second\nccc\n')
            self.assertEqual(list(A.fasta_records(path)),[('first','ATGAAATAANN'),('second','CCC')])
            for text in ('>q\nACG\n>q\nTTT\n','>empty\n','ACGT\n'):
                path.write_text(text)
                with self.assertRaises(ValueError): list(A.fasta_records(path))

    def test_encoder_uses_exact_four_blocks_no_special_tokens_and_bf16_materialization(self):
        backbone = Backbone()
        features = A.encode_window(backbone,Tokenizer(),'ATGAAATAA','cpu')
        self.assertEqual(tuple(features.shape),(1,4096,2048))
        self.assertEqual(features.dtype,torch.float32)
        self.assertEqual(backbone.calls,[(1,1024)]*4)
        self.assertTrue(torch.equal(features[0,:3,0],torch.tensor([1.,4.,0.])))
        class SpecialTokenizer(Tokenizer):
            def __call__(self, sequence, **kwargs):
                ids = super().__call__(sequence,**kwargs)['input_ids']
                return {'input_ids':torch.cat((ids,ids[:,:1]),dim=1)}
        with self.assertRaisesRegex(ValueError,'alignment'):
            A.encode_window(Backbone(),SpecialTokenizer(),'ATGAAATAA','cpu')
        def invalid_hidden(input_ids,attention_mask):
            return SimpleNamespace(last_hidden_state=torch.zeros(1,1024,12))
        with self.assertRaisesRegex(ValueError,'Invalid frozen'):
            A.encode_window(invalid_hidden,Tokenizer(),'ATGAAATAA','cpu')

    def test_artificial_fasta_to_both_strand_records_owned_gene_and_gff_phase(self):
        torch.manual_seed(0); model = GeneQueryModel(2048,2,8).eval()
        with torch.no_grad():
            model.query.exists.weight.zero_(); model.query.exists.bias.fill_(1.)
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)
            txs,scopes = A.annotate_contig(model,Backbone(),Tokenizer(),'q','ATGAAATAA',0,out,'cpu',1024**2)
            self.assertEqual([s['strand'] for s in scopes],['+','-'])
            self.assertEqual(len(txs),1)
            self.assertEqual(txs[0]['CDS'],[(0,9,'0')]); self.assertEqual(txs[0]['strand'],'+')
            self.assertEqual(scopes[0]['assembly_counts']['duplicate_owner_records'],1)
            self.assertEqual(scopes[1]['query_counts']['unterminated_queries'],2)
            for scope in scopes:
                row = list(A.saved_records(out/scope['file']))[0]
                self.assertFalse(row['reference_used'])
                self.assertEqual(len(row['raw_queries']),2)
            gff = out/'artificial.gff3'; write_gff3(gff,txs,{'q':9},'M28_A_SYNTHETIC')
            parsed = V.E.S.E.primary_transcripts(V.E.S.E.parse_annotation(gff,{'q':9}),complete_only=False)
            V.validate_gff_roundtrip(parsed,txs)
            with self.assertRaisesRegex(RuntimeError,'budget'):
                A.output_bytes(out,1)

    def test_eval_rejects_partial_wrong_checkpoint_or_scope_before_reference_access(self):
        fitted,decoded = complete_pair(); V.validate_complete(fitted,decoded)
        for key,value in (('inference_complete',False),('windows',8689),('reference_used',True),
                          ('checkpoint','synthetic/wrong.pt'),('scopes',decoded['scopes'][:-1])):
            wrong = copy.deepcopy(decoded); wrong[key] = value
            with self.assertRaises(ValueError): V.validate_complete(fitted,wrong)
        wrong = copy.deepcopy(decoded); wrong['scopes'][0]['windows'][0][0] = 1
        with self.assertRaisesRegex(ValueError,'grid'):
            V.validate_complete(fitted,wrong)

    def test_eval_record_operating_point_and_query_provenance_cannot_use_B1(self):
        scope = {'species':'synthetic','seqid':'q','strand':'+','windows':[[0,9]]}
        row = {'window_index':0,'window_id':'synthetic|q|0|+','reference_used':False,
               'score_kind':'query_existence_logit_minus_threshold','existence_logit_threshold':0.,
               'counts':{'queries':16},'raw_queries':[{'terminated':False}]*16,'chains':[],'score':[]}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'artificial.jsonl'
            path.write_text(json.dumps(row)+'\n')
            self.assertEqual(list(V.query_scope_records(path,scope)),[row])
            for key,value in (('score_kind','old_B1'),('existence_logit_threshold',1.),
                              ('counts',{'queries':15}),('raw_queries',row['raw_queries'][:-1])):
                wrong = dict(row); wrong[key] = value
                path.write_text(json.dumps(wrong)+'\n')
                with self.assertRaisesRegex(ValueError,'provenance'):
                    list(V.query_scope_records(path,scope))

    def test_two_species_gate_uses_integer_counts_and_same_full_denominators(self):
        control = counts(); current = copy.deepcopy(control)
        for value in current.values():
            value['exact_chain']['tp'] += 1; value['exact_chain']['f1'] = 0.
        self.assertTrue(V.continuation_gate(current,control)['passed'])  # Ignore misleading rounded field.
        species = next(iter(current))
        for change in ('no_improvement','more_background','reference','background_denominator'):
            wrong = copy.deepcopy(current)
            if change == 'no_improvement': wrong[species]['exact_chain']['tp'] = 100
            if change == 'more_background': wrong[species]['background']['predicted_span_bp'] = 21
            if change == 'reference': wrong[species]['exact_chain']['reference'] -= 1
            if change == 'background_denominator': wrong[species]['background']['bp'] += 1
            if change in ('reference','background_denominator'):
                with self.assertRaisesRegex(ValueError,'denominator'):
                    V.continuation_gate(wrong,control)
            else: self.assertFalse(V.continuation_gate(wrong,control)['passed'])

    def test_reference_loader_preserves_unsupported_rows_and_stops_on_changed_reference(self):
        species = 'synthetic'
        refs = [{'id':'a','gene_id':'ga','seqid':'q','strand':'+','CDS':[(0,9,'0')],
                 'M28_unsupported_strata_not_excluded':['longer_than_window']},
                {'id':'b','gene_id':'gb','seqid':'q','strand':'+','CDS':[(20,29,'0')],
                 'M28_unsupported_strata_not_excluded':['non_ATG_start']}]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); path = root/'data/m1_screen'/species
            path.mkdir(parents=True); (path/'reference.gff3').write_text('')
            with patch.object(V.E,'ROOT',root), patch.object(V.E.S,'SCOPE',[(species,'q',60,2)]), \
                 patch.object(V,'REF_COUNTS',{species:2}), patch.object(Q.M,'read_allowed',return_value={'q':'C'*60}), \
                 patch.object(V.E.S.E,'parse_annotation',return_value={}), \
                 patch.object(V.E.S.E,'primary_transcripts',return_value=refs), \
                 patch.object(V.E.P,'raw_features',return_value=(None,None,{},None)), \
                 patch.object(V.E,'select_reference',return_value=(refs,refs)), \
                 patch.object(V,'rows',return_value=refs) as frozen:
                actual = V.reference_context()[0][species]
                self.assertEqual(actual,refs)
                frozen.return_value = refs[:1]
                with self.assertRaisesRegex(ValueError,'do not filter unsupported'):
                    V.reference_context()

    def test_gff_phase_change_is_not_hidden_by_equal_chain_coordinates(self):
        expected = [{'id':'a.t1','gene_id':'a','seqid':'q','strand':'+','CDS':[(0,4,'0'),(9,14,'2')]}]
        V.validate_gff_roundtrip(expected,expected)
        wrong = copy.deepcopy(expected); wrong[0]['CDS'][1] = (9,14,'0')
        with self.assertRaisesRegex(ValueError,'phase'): V.validate_gff_roundtrip(wrong,expected)


if __name__ == '__main__':
    unittest.main()
