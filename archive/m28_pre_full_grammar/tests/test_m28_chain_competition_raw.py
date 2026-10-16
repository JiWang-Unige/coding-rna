import copy
import json
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path
from unittest.mock import patch
import torch
from test_m28_query_decode import fixture

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import infer_chain_competition as I
import evaluate_chain_competition as E
from src.m28.export import write_gff3


class SyntheticGenerator:
    def __init__(self, model):
        self.query = model; self.calls = 0

    def __call__(self, features, dna, valid):
        self.calls += 1
        if tuple(dna.shape) != (1,I.M.WIDTH,5) or len(valid) != 1:
            raise AssertionError('Original raw-window DNA contract changed')
        return dict(features=features)


class SyntheticHead:
    def __init__(self, gain):
        self.gain = gain

    def __call__(self, hidden, chains, valid):
        return dict(logits=hidden.new_full((len(chains),),self.gain))


def fitted_receipt():
    return dict(arm='A13',configuration=I.R.CONFIG,protocol=I.R.PROTOCOL,
        feature_revision=I.C.I.F.REVISION,feature_spec=I.C.R.FEATURE_SPEC,completed=True,
        optimizer_steps=4608,quality_trainable_parameters=165635,generator_trainable_parameters=0,
        final_checkpoint=str((I.R.FIT/'step_004608.pt').relative_to(I.ROOT)),epochs=[dict(steps=1536)]*3,
        source_free_pool='outputs/M28-A6-CHAIN-QUALITY-FIT/free_pool.jsonl',
        source_TRAIN_cohort='outputs/M28-A12-QUALITY-TRAIN-DIAGNOSTIC',
        DEV_used_for_selection=False,test_setaria_used=False)


def preparation_receipt(fitted):
    a = {k:dict(candidates=10,positive_gain=4,positive_candidates_not_above_max_known_overlap=3)
         for k in ('positive','free_positive','known_nonmatch','five_or_more_positive','unknown')}
    b = copy.deepcopy(a)
    for key in ('positive','free_positive','five_or_more_positive'):
        b[key]['positive_gain'] = 5
    pairs = dict(A6=a,A13=b)
    checks = {}
    for species in I.C.R.ALLOW:
        for suffix in ('free_correct_positive_gain_increased','all_correct_positive_gain_increased',
                       'known_nonmatch_positive_gain_not_increased'):
            checks[species+'_'+suffix] = True
    checks['pooled_five_or_more_correct_positive_gain_increased'] = True
    checks['pooled_correct_not_above_known_conflict_not_increased'] = True
    return dict(experiment='M28-A13-CHAIN-COMPETITION-TRAIN-READOUT',completed=True,
        configuration=I.R.CONFIG,protocol=I.R.PROTOCOL,checkpoint=fitted['final_checkpoint'],
        source_generator=I.C.GENERATOR,source_A6='outputs/M28-A12-QUALITY-TRAIN-DIAGNOSTIC',
        windows=1491,original_draws=1536,candidate_window_pairs=16857,optimizer_steps_this_run=0,
        per_species={s:copy.deepcopy(pairs) for s in I.C.R.ALLOW},pooled=copy.deepcopy(pairs),
        candidate_geometry_labels_and_origins_unchanged=True,threshold=0.,
        predeclared_TRAIN_preparation_gate=checks,TRAIN_preparation_supported=True,
        **{k:False for k in ('new_candidate_search','new_backbone_forward','DEV_forward_used',
            'test_setaria_used','threshold_scan','head_weights_updated','checkpoints_written')})


def decoded_receipt(fitted):
    scopes = []
    for ordinal,(species,lengths) in enumerate(I.C.I.DEV_LENGTHS.items()):
        for seqid,length in lengths.items():
            for strand in ('+','-'):
                scopes.append(dict(species=species,seqid=seqid,strand=strand,ordinal=ordinal,
                    windows=[[a,min(a+I.M.WIDTH,length)] for a in I.M.starts(length)],
                    file='contig_'+str(ordinal)+'_'+strand+'.jsonl'))
    return dict(arm='A13',completed=True,inference_complete=True,configuration=I.R.CONFIG,
        protocol=I.DEV_PROTOCOL,TRAIN_protocol=I.R.PROTOCOL,TRAIN_preparation_gate_passed=True,
        checkpoint=fitted['final_checkpoint'],optimizer_steps=4608,generator_checkpoint=I.C.GENERATOR,
        feature_revision=I.C.I.F.REVISION,feature_spec=I.C.R.FEATURE_SPEC,input_mode='fixed_DEV_raw_FASTA',
        reference_used=False,labels_or_feature_cache_used=False,optimizer_steps_this_run=0,
        DEV_used_for_selection=False,test_setaria_used=False,windows=8690,scopes=scopes,
        lengths={q:n for qs in I.C.I.DEV_LENGTHS.values() for q,n in qs.items()})


class CompetitionRawTests(unittest.TestCase):
    def test_failed_or_false_positive_train_gate_rejects_before_generator_load(self):
        fitted = fitted_receipt(); good = preparation_receipt(fitted)
        I.validate_train_preparation(fitted,good)
        failed = copy.deepcopy(good)
        species = next(iter(I.C.R.ALLOW))
        failed['per_species'][species]['A13']['known_nonmatch']['positive_gain'] = 5
        with self.assertRaisesRegex(ValueError,'TRAIN preparation gate'):
            I.validate_train_preparation(fitted,failed)
        for key,value in (('optimizer_steps',1536),('final_checkpoint','early.pt'),
                          ('feature_spec',{'encoder_context':'whole_4096'})):
            with self.subTest(key=key),self.assertRaises(ValueError):
                I.validate_train_preparation({**fitted,key:value},good)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); fit = root/'fit'; readout = root/'readout'
            fit.mkdir(); readout.mkdir()
            with patch.object(I,'ROOT',root),patch.object(I.R,'ROOT',root),\
                 patch.object(I.R,'FIT',fit),patch.object(I.R,'READOUT',readout),\
                 patch.object(I.C,'load_generator') as generator:
                fitted = fitted_receipt(); failed = preparation_receipt(fitted)
                failed['per_species'][species]['A13']['known_nonmatch']['positive_gain'] = 5
                for path,receipt in ((fit,fitted),(readout,failed)):
                    (path/'STATUS').write_text('COMPLETED\n')
                    (path/'summary.json').write_text(json.dumps(receipt))
                with self.assertRaisesRegex(ValueError,'TRAIN preparation gate'):
                    I.main(root/'infer')
                generator.assert_not_called()

    def test_raw_normalization_two_orientations_and_minus_gff_phase(self):
        oriented = 'ATGCGTTAGTTTAA'; model,features = fixture(oriented,[(0,4),(9,14)])
        generator = SyntheticGenerator(model); raw = I.M.rc(oriented).lower()+'r'; seen = []
        def encode(backbone,tokenizer,sequence,device):
            seen.append(sequence)
            values = torch.cat([torch.zeros(1,4),features]) if sequence == 'N'+oriented else torch.zeros(15,4)
            return values[None]
        with tempfile.TemporaryDirectory() as directory,patch.object(I.C.I,'encode_window',side_effect=encode):
            out = Path(directory)
            txs,scopes = I.annotate_contig(generator,SyntheticHead(1.),None,None,
                'synthetic','chr',raw,0,out,'cpu')
            self.assertEqual(generator.calls,2)
            self.assertEqual(seen,[raw[:-1].upper()+'N','N'+oriented])
            minus = [t for t in txs if t['strand'] == '-']
            self.assertEqual(len(minus),1)
            self.assertEqual(minus[0]['CDS'],[(0,5,'2'),(10,14,'0')])
            gff = out/'A13.gff3'; write_gff3(gff,txs,{'chr':15},'M28_A13')
            parsed = E.E.S.E.primary_transcripts(E.E.S.E.parse_annotation(gff,{'chr':15}),complete_only=False)
            E.A.validate_gff_roundtrip(parsed,txs)
            self.assertEqual([s['windows'] for s in scopes],[[[0,15]],[[0,15]]])

    def test_zero_gain_retains_natural_pool_but_selects_null(self):
        sequence = 'ATGAAATAA'; model,features = fixture(sequence,[(0,9)])
        with tempfile.TemporaryDirectory() as directory,patch.object(I.C.I,'encode_window',return_value=features[None]):
            out = Path(directory)
            txs,scopes = I.annotate_contig(SyntheticGenerator(model),SyntheticHead(0.),
                None,None,'synthetic','chr',sequence,0,out,'cpu')
            self.assertEqual(txs,[])
            row = next(I.C.I.saved_records(out/scopes[0]['file']))
            self.assertTrue(row['chains'])
            self.assertTrue(all(v == 0. for v in row['score']))
            self.assertEqual(row['score_kind'],I.SCORE_KIND)
            self.assertEqual(scopes[0]['assembly_counts']['selected_unique'],0)

    def test_partial_grid_explicit_fasta_and_wrong_context_reject(self):
        fitted = fitted_receipt(); good = decoded_receipt(fitted); E.validate_complete(fitted,good)
        for key,value in (('scopes',good['scopes'][:-1]),('input_mode','explicit_raw_FASTA'),
                          ('feature_spec',{'encoder_context':'whole_4096'}),('reference_used',True)):
            with self.subTest(key=key),self.assertRaises(ValueError):
                E.validate_complete(fitted,{**good,key:value})
        truncated = copy.deepcopy(good); truncated['scopes'][0]['windows'].pop()
        with self.assertRaises(ValueError):
            E.validate_complete(fitted,truncated)

    def test_paired_pool_allows_new_scores_but_not_geometry_source_order_drift(self):
        scope = dict(species='synthetic',seqid='chr',strand='+',windows=[[0,9]])
        row = dict(window_index=0,window_id='synthetic|chr|0|+',chains=[[[0,9]]],score=[1.],
            decoder_sources=[['beam4','greedy']],score_kind=I.SCORE_KIND,quality_logit_threshold=0.,
            reference_used=False,counts=dict(free_candidates=1,nonpositive_candidates=0))
        old = {**row,'score':[-1.],'score_kind':'A6_complete_chain_logit_minus_window_null'}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); current = root/'new.jsonl'; previous = root/'old.jsonl'
            previous.write_text(json.dumps(old)+'\n'); current.write_text(json.dumps(row)+'\n')
            self.assertEqual(len(list(E.matched_records(current,scope,previous,scope,Counter()))),1)
            for key,value in (('chains',[[[0,6]]]),('decoder_sources',[['greedy']]),('score',[float('nan')])):
                current.write_text(json.dumps({**row,key:value})+'\n')
                with self.subTest(key=key),self.assertRaises(ValueError):
                    list(E.matched_records(current,scope,previous,scope,Counter()))
            current.write_text('')
            with self.assertRaises(ValueError):
                list(E.matched_records(current,scope,previous,scope,Counter()))


if __name__ == '__main__':
    unittest.main()
