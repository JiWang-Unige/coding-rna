import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock
import torch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import infer_full_train as I
import evaluate_full_train as V


def fit(arm):
    return {'completed':True,'arm':arm,'configuration':copy.deepcopy(I.R.CONFIGS[arm]),
            'optimizer_steps':I.R.STEPS,'feature_revision':I.F.REVISION,'feature_spec':I.R.FEATURE_SPEC,
            'dataset':{'windows':I.R.N},'final_checkpoint':str((I.FIT[arm]/'step_076536.pt').relative_to(I.ROOT)),
            'main_evaluation':'last_step_076536_only','protocol':I.R.PROTOCOL,
            'DEV_used_for_selection':False,'test_setaria_used':False,'epochs':[{'steps':I.R.N}]*3}


def inference(arm):
    fitted = fit(arm)
    scopes = [{'species':s,'seqid':q,'strand':t,
               'windows':[[a,min(a+I.M.WIDTH,n)] for a in I.M.starts(n)]}
              for s,qs in I.DEV_LENGTHS.items() for q,n in qs.items() for t in ('+','-')]
    return {'inference_complete':True,'arm':arm,'configuration':copy.deepcopy(I.R.CONFIGS[arm]),
            'optimizer_steps':I.R.STEPS,'checkpoint':fitted['final_checkpoint'],
            'encoder_revision':I.F.REVISION,'feature_spec':I.R.FEATURE_SPEC,
            'input_mode':'fixed_DEV_raw_FASTA','reference_used':False,'labels_or_feature_cache_used':False,
            'optimizer_steps_this_run':0,'test_setaria_used':False,'DEV_used_for_selection':False,
            'windows':8690,'scopes':scopes,'lengths':{q:n for qs in I.DEV_LENGTHS.values() for q,n in qs.items()}}


class FullInferenceTests(unittest.TestCase):
    def test_only_matching_last_fit_and_checkpoint(self):
        for arm in ('C0','Q'):
            s = fit(arm); I.validate_fit(s,arm)
            state = {'arm':arm,'configuration':s['configuration'],'step':76536,
                     'revision':I.F.REVISION,'protocol':I.R.PROTOCOL}
            I.validate_checkpoint(state,arm)
            for field,value in (('step',51024),('revision','wrong'),('arm','Q' if arm=='C0' else 'C0')):
                with self.assertRaises(ValueError): I.validate_checkpoint({**state,field:value},arm)
            for field,value in (('optimizer_steps',4608),('completed',False),('DEV_used_for_selection',True)):
                with self.assertRaises(ValueError): I.validate_fit({**s,field:value},arm)
            bad = copy.deepcopy(s); bad['configuration']['lr'] = 1e-3
            with self.assertRaises(ValueError): I.validate_fit(bad,arm)

    def test_complete_raw_DEV_grid_and_no_reference_or_cache_provenance(self):
        for arm in ('C0','Q'):
            s = inference(arm); V.validate_complete(fit(arm),s,arm)
            for field,value in (('windows',8689),('input_mode','explicit_raw_FASTA'),
                                ('labels_or_feature_cache_used',True),('reference_used',True),
                                ('optimizer_steps_this_run',1),('test_setaria_used',True)):
                with self.assertRaises(ValueError): V.validate_complete(fit(arm),{**s,field:value},arm)
            for kind in ('partial','duplicate','coordinate'):
                bad = copy.deepcopy(s)
                if kind=='partial': bad['scopes'][0]['windows'].pop()
                elif kind=='duplicate': bad['scopes'][1] = bad['scopes'][0]
                else: bad['scopes'][0]['windows'][0][0] = 1
                with self.assertRaises(ValueError): V.validate_complete(fit(arm),bad,arm)

    def test_DEV_FASTA_materializes_only_permitted_contig_and_no_annotation(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td); path = root/'data/m1_screen/s/genome.fa'; path.parent.mkdir(parents=True)
            path.write_text('>sealed\nNOT_A_VALID_SEQUENCE\n>allowed\nACGTNN\n>train\nTTTT\n')
            with mock.patch.object(I,'ROOT',root), mock.patch.object(I,'DEV_LENGTHS',{'s':{'allowed':6}}), \
                 mock.patch.object(V.A,'reference_context',side_effect=AssertionError('inference cannot read references')):
                self.assertEqual(list(I.input_contigs(None)),[('s','allowed','ACGTNN')])

    def test_C0_bridge_uses_native_gains_and_never_queries(self):
        sequence = 'ATGAAATAA'; outputs = {'segmentation_logits':torch.zeros(1,9,15)}
        native = {'chains':[[(0,9)]],'gains':[2.5],'partial_paths':0}
        with mock.patch.object(I,'decode_c0',return_value=native) as decode, \
             mock.patch.object(I,'query_window_candidates',side_effect=AssertionError('not Q')):
            r = I.window_record('C0',None,outputs,sequence,0,'s|q|0|+')
        self.assertEqual(r['score'],[2.5]); self.assertEqual(r['score_kind'],'native_HMM_chain_removal_gain')
        self.assertEqual(decode.call_args.args[1].shape,(9,15))
        self.assertFalse(r['reference_used'])

    def test_Q_routes_actual_collaborative_model_with_fixed_16_queries_and_tau_zero(self):
        sequence = 'CCC'+'ATGAAATAA'+'C'*24
        model = I.R.new_model('Q','cpu').eval()
        outputs = model(torch.zeros(1,6,2048),I.one_hot(sequence,36)[None],[36])
        with mock.patch.object(I,'decode_c0',side_effect=AssertionError('Q cannot use HMM')):
            r = I.window_record('Q',model,outputs,sequence,0,'s|q|0|+')
        self.assertEqual(r['counts']['queries'],16); self.assertEqual(len(r['raw_queries']),16)
        self.assertEqual(r['existence_logit_threshold'],0.)
        self.assertEqual(len(r['chains']),len(r['score'])); self.assertFalse(r['reference_used'])

    def test_raw_two_strands_tail_assembly_and_exact_GFF_coordinates_phases(self):
        sequence = 'ATGAAATAA'+'C'*30; seen = []
        def encode(backbone,tokenizer,seq,device):
            seen.append(seq); return torch.zeros(1,6,2048)
        def free(arm,model,outputs,seq,index,ident):
            return {'window_index':index,'window_id':ident,'chains':[[[0,4],[7,12]]],
                    'score':[1.],'counts':{'queries':16},'reference_used':False}
        model = lambda features,dna,n: {'features':torch.zeros(1,n[0],128)}
        with tempfile.TemporaryDirectory() as td:
            out = Path(td)
            with mock.patch.object(I.M,'WIDTH',36), mock.patch.object(I.M,'starts',return_value=[0,3]), \
                 mock.patch.object(I,'encode_window',side_effect=encode), mock.patch.object(I,'window_record',side_effect=free):
                txs,scopes = I.annotate_contig('Q',model,None,None,'s','q',sequence,0,out,'cpu')
            self.assertEqual(seen,[sequence[:36],sequence[3:],I.M.rc(sequence[:36]),I.M.rc(sequence[3:])])
            self.assertEqual([s['strand'] for s in scopes],['+','-'])
            self.assertTrue(all(s['windows']==[(0,36),(3,39)] for s in scopes))
            self.assertTrue(all(len(list(I.saved_records(out/s['file'])))==2 for s in scopes))
            self.assertTrue(any(t['strand']=='-' for t in txs))
            gff = out/'Q.gff3'; I.write_gff3(gff,txs,{'q':39},'M28_FULL_Q')
            parsed = V.E.S.E.primary_transcripts(V.E.S.E.parse_annotation(gff,{'q':39}),complete_only=False)
            V.A.validate_gff_roundtrip(parsed,txs)
            bad = copy.deepcopy(parsed); bad[0]['CDS'][0] = (*bad[0]['CDS'][0][:2],'1')
            with self.assertRaises(ValueError): V.A.validate_gff_roundtrip(bad,txs)

    def test_candidate_record_provenance_and_completeness(self):
        scope = {'species':'s','seqid':'q','strand':'+','windows':[[0,9]]}
        row = {'window_index':0,'window_id':'s|q|0|+','reference_used':False,'chains':[],'score':[],
               'score_kind':'query_existence_logit_minus_threshold','existence_logit_threshold':0.,
               'counts':{'queries':16},'raw_queries':[{}]*16}
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'x.jsonl'; path.write_text(json.dumps(row)+'\n')
            self.assertEqual(list(V.scope_records(path,scope,'Q',V.Counter())),[row])
            for field,value in (('reference_used',True),('existence_logit_threshold',1.),('window_id','wrong')):
                path.write_text(json.dumps({**row,field:value})+'\n')
                with self.assertRaises(ValueError): list(V.scope_records(path,scope,'Q',V.Counter()))
            path.write_text('')
            with self.assertRaises(ValueError): list(V.scope_records(path,scope,'Q',V.Counter()))

    def test_historical_comparison_labels_exposure_and_rejects_changed_ruler(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td); path = root/'outputs/M28-A-EVALUATION/summary.json'
            path.parent.mkdir(parents=True)
            scores = {s:{'exact_chain':{'tp':1,'predicted':2,'reference':n}} for s,n in V.A.REF_COUNTS.items()}
            s = {'reference_primary_chains':7728,'test_setaria_used':False,
                 'methods':{'A':{'per_species':scores,'pooled':{'exact_chain':{'tp':2}}}}}
            path.write_text(json.dumps(s))
            with mock.patch.object(I,'ROOT',root):
                r = V.historical_exposure_comparison()
                self.assertEqual(r['A']['optimizer_steps'],4608)
                self.assertIn('different_exposure',r['A']['status'])
                self.assertEqual(r['A3']['status'],'not_available_not_a_result')
                s['reference_primary_chains'] = 6450; path.write_text(json.dumps(s))
                with self.assertRaises(ValueError): V.historical_exposure_comparison()

    def test_predeclared_gate_strict_F1_and_background_both_species(self):
        def scores(tp,span):
            return {s:{'exact_chain':{'tp':tp,'predicted':1000,'reference':n},
                       'background':{'bp':10000,'predicted_span_bp':span}} for s,n in V.A.REF_COUNTS.items()}
        c0 = scores(500,100); q = scores(501,100)
        self.assertTrue(V.A.continuation_gate(q,c0)['passed'])
        self.assertFalse(V.A.continuation_gate(c0,c0)['passed'])
        q['oryza_sativa']['background']['predicted_span_bp'] = 101
        self.assertFalse(V.A.continuation_gate(q,c0)['passed'])
        q = scores(501,100); q['arabidopsis_thaliana']['exact_chain']['reference'] -= 1
        with self.assertRaises(ValueError): V.A.continuation_gate(q,c0)


if __name__=='__main__':
    torch.set_num_threads(1)
    unittest.main()
