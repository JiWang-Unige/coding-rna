"""A37 complete evaluator: synthetic metadata, components and saved GFF only."""
from collections import Counter
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from test_m28_local_path_caller_evaluation import decoded_receipt as baseline_receipt, raw_row as baseline_row

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import evaluate_context_chain as V
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript, write_gff3


def decoded_receipt():
    old = baseline_receipt('A3')
    return dict(arm='A37',completed=True,inference_complete=True,configuration=V.I.D.CONFIG,
        protocol=V.I.DEV_PROTOCOL,TRAIN_protocol=V.I.D.PROTOCOL,fit_protocol=V.I.D.F.PROTOCOL,
        TRAIN_preparation_gate_passed=True,source_TRAIN_jobs=dict(fit=V.I.D.FIT_JOB,
            readout=V.I.READOUT_JOB,CPU=V.I.CHECK_JOB),generator_checkpoint=V.I.C.GENERATOR,
        scorer_checkpoint=str((V.I.D.F.OUT/'step_004608.pt').relative_to(V.I.ROOT)),
        generator_optimizer_steps=4608,scorer_optimizer_steps=4608,
        model_parameters=dict(generator=649457,scorer=612324),feature_revision=V.I.C.I.F.REVISION,
        feature_spec=V.I.C.R.FEATURE_SPEC,input_mode='fixed_DEV_raw_FASTA',reference_used=False,
        labels_or_feature_cache_used=False,optimizer_steps_this_run=0,DEV_used_for_selection=False,
        test_setaria_used=False,A33_outputs_opened=False,GEN_encoded_once_per_oriented_window=True,
        natural_proposal_generator='unchanged_A3_greedy_plus_beam4',model_weights_unchanged=True,
        all_parameters_frozen_and_gradients_absent=True,score_kind=V.I.SCORE_KIND,
        windows=8690,scopes=old['scopes'],lengths=old['lengths'])


def raw_row():
    return dict(baseline_row(),variant='A37',score_kind=V.I.SCORE_KIND,score=[2.],
                exon_additive=[1.],junction_additive=[0.],additive=[1.],nonadditive=[3.],null=2.)


class ContextChainEvaluationTests(unittest.TestCase):
    def test_complete_grid_last_models_source_identity_and_shared_encoding_required(self):
        good = decoded_receipt(); V.validate_complete(good)
        for key,value in (('scopes',good['scopes'][:-1]),('input_mode','explicit_raw_FASTA'),
                ('feature_spec',{'encoder_context':'whole_4096'}),('scorer_checkpoint','early.pt'),
                ('generator_checkpoint','different.pt'),('scorer_optimizer_steps',1536),
                ('GEN_encoded_once_per_oriented_window',False),('source_TRAIN_jobs',{'fit':'other'}),
                ('model_parameters',dict(generator=649457,scorer=248068))):
            with self.subTest(key=key),self.assertRaises(ValueError):
                V.validate_complete({**good,key:value})
        truncated = copy.deepcopy(good); truncated['scopes'][0]['windows'].pop()
        with self.assertRaises(ValueError):
            V.validate_complete(truncated)

    def test_original_pool_source_and_exact_saved_FP32_components_required(self):
        scope = dict(species='synthetic',seqid='chr',strand='+',windows=[[0,9]])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); new,old = root/'new.jsonl',root/'old.jsonl'
            old.write_text(json.dumps(baseline_row())+'\n'); new.write_text(json.dumps(raw_row())+'\n')
            self.assertEqual(len(list(V.paired_records(new,scope,old,scope,Counter()))),1)
            for key,value in (('chains',[[[0,6]]]),('decoder_sources',[['greedy']]),('score',[float('nan')]),
                    ('score_threshold',.1),('counts',dict(free_candidates=1,nonpositive_candidates=1)),
                    ('score',[3.]),('additive',[2.]),('nonadditive',[])):
                new.write_text(json.dumps({**raw_row(),key:value})+'\n')
                with self.subTest(key=key),self.assertRaises((ValueError,AssertionError)):
                    list(V.paired_records(new,scope,old,scope,Counter()))
            new.write_text(json.dumps(raw_row())+'\n')
            old.write_text(json.dumps({**baseline_row(),'variant':'A31_A23'})+'\n')
            with self.assertRaisesRegex(ValueError,'A31_A3'):
                list(V.paired_records(new,scope,old,scope,Counter()))

    def test_saved_negative_strand_components_selection_and_GFF_reconcile_before_reference(self):
        row = dict(raw_row(),window_id='synthetic|chr|0|-',chains=[[[0,4],[9,14]]])
        old = dict(baseline_row(),window_id=row['window_id'],chains=row['chains'])
        scope = dict(species='synthetic',seqid='chr',strand='-',ordinal=0,windows=[[0,14]],file='scope.jsonl',
            candidate_counts=row['counts'],assembly_counts=assemble_predictions([row],[[0,14]],'-')['counts'])
        decoded = dict(scopes=[scope],lengths={'chr':14},primary_CDS_chain_predictions=1)
        tx = as_transcript([(0,5),(10,14)],'chr','-','M28_A37_0_m_0',2.)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); new,prior = root/'new',root/'old'; new.mkdir(); prior.mkdir()
            (new/'scope.jsonl').write_text(json.dumps(row)+'\n')
            (prior/'scope.jsonl').write_text(json.dumps(old)+'\n')
            write_gff3(new/'A37.gff3',[tx],{'chr':14},'M28_A37')
            with patch.object(V.A,'reference_context') as reference:
                predictions,_,_ = V.reconstruct(decoded,dict(scopes=[scope]),new,prior)
                self.assertEqual(predictions['synthetic'],[tx]); reference.assert_not_called()
                with self.assertRaisesRegex(ValueError,'prediction count'):
                    V.reconstruct({**decoded,'primary_CDS_chain_predictions':2},dict(scopes=[scope]),new,prior)
                reference.assert_not_called()

    def test_full_denominators_bind_all_comparisons_macro_and_worst_species(self):
        species = list(V.H.A3.REF_COUNTS)
        per = {s:dict(exact_chain=V.E.S.prf(i+1,i+3,V.H.A3.REF_COUNTS[s]),background=dict(bp=100))
               for i,s in enumerate(species)}
        current = dict(per_species=per,pooled=dict(exact_chain=V.E.S.prf(3,7,7728)))
        names = ('A37','A31_A3','A31_A23','ANNEVO','Tiberius','Helixer','A3','A6','A','R5_C0')
        methods = {name:copy.deepcopy(current) for name in names}
        table = V.comparison_table(methods,current)
        f1 = [per[s]['exact_chain']['f1'] for s in species]
        self.assertEqual(set(table),set(names))
        self.assertEqual(table['A37']['macro_F1'],sum(f1)/len(f1))
        self.assertEqual(table['A37']['worst_species_F1'],min(f1))
        bad = copy.deepcopy(current); bad['per_species'][species[0]]['exact_chain']['reference'] -= 1
        with self.assertRaisesRegex(ValueError,'denominators'):
            V.comparison_table(dict(A37=current,ANNEVO=bad),current)
        bad = copy.deepcopy(current); bad['pooled']['exact_chain']['reference'] -= 1
        with self.assertRaisesRegex(ValueError,'pooled.*denominator'):
            V.comparison_table(dict(A37=current,Helixer=bad),current)
        bad = copy.deepcopy(current); bad['per_species'][species[0]]['background']['bp'] -= 1
        with self.assertRaisesRegex(ValueError,'denominators'):
            V.comparison_table(dict(A37=current,A31_A3=bad),current)

    def test_empty_candidate_windows_are_required_in_complete_paired_scope(self):
        scope = dict(species='synthetic',seqid='chr',strand='+',windows=[[0,9],[9,18]])
        empty = dict(raw_row(),window_index=1,window_id='synthetic|chr|9|+',chains=[],score=[],
            decoder_sources=[],exon_additive=[],junction_additive=[],additive=[],nonadditive=[],null=.25,
            counts=dict(free_candidates=0,nonpositive_candidates=0))
        old_empty = dict(baseline_row(),window_index=1,window_id=empty['window_id'],chains=[],score=[],
                         decoder_sources=[],counts=empty['counts'])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); new,old = root/'new.jsonl',root/'old.jsonl'
            new_rows = '\n'.join(json.dumps(row) for row in (raw_row(),empty))+'\n'
            old_rows = '\n'.join(json.dumps(row) for row in (baseline_row(),old_empty))+'\n'
            new.write_text(new_rows); old.write_text(old_rows); counts = Counter()
            rows = list(V.paired_records(new,scope,old,scope,counts))
            self.assertEqual(len(rows),2); self.assertEqual(rows[1]['chains'],[])
            self.assertEqual(dict(counts),dict(free_candidates=1,nonpositive_candidates=0))
            new.write_text(json.dumps(raw_row())+'\n')
            with self.assertRaises(ValueError):
                list(V.paired_records(new,scope,old,scope,Counter()))
            new.write_text(new_rows); old.write_text(json.dumps(baseline_row())+'\n')
            with self.assertRaises(ValueError):
                list(V.paired_records(new,scope,old,scope,Counter()))


if __name__ == '__main__':
    unittest.main()
