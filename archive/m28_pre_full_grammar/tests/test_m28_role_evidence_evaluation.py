"""A40 complete-ruler routing, using only artificial metadata, records and GFF."""
from collections import Counter
from copy import deepcopy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from test_m28_context_chain_evaluation import decoded_receipt as old_receipt
from test_m28_local_path_caller_evaluation import raw_row as baseline_row

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import evaluate_role_evidence_context as V
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript, write_gff3


def decoded_receipt():
    value = old_receipt()
    value.update(arm='A40',configuration=V.I.D.CONFIG,protocol=V.I.DEV_PROTOCOL,
        TRAIN_protocol=V.I.D.PROTOCOL,fit_protocol=V.I.D.F.PROTOCOL,
        source_TRAIN_jobs=dict(fit=V.I.FIT_JOB,readout=V.I.READOUT_JOB,CPU=V.I.CHECK_JOB),
        scorer_checkpoint=str((V.I.D.F.OUT/'step_004608.pt').relative_to(V.I.ROOT)),
        frozen_local_source_checkpoint=V.I.D.CONFIG['frozen_local_source'],
        frozen_local_source_matches_actual_fit=True,local_source_updated=False,
        A3_shared_encoded_once_per_oriented_window=True,
        model_parameters=dict(generator=649457,scorer=613348),score_kind=V.I.SCORE_KIND)
    return value


def raw_row():
    return dict(baseline_row(),variant='A40',score_kind=V.I.SCORE_KIND,score=[2.],
        exon_additive=[1.],junction_additive=[0.],additive=[1.],nonadditive=[3.],null=2.,
        frozen_local_source_checkpoint=V.I.D.CONFIG['frozen_local_source'])


class RoleEvidenceEvaluationTests(unittest.TestCase):
    def test_A40_source_grid_and_exposure_metadata_block_before_reference_access(self):
        good = decoded_receipt(); V.validate_complete(good)
        mutations = (('arm','A39'),('configuration',old_receipt()['configuration']),
            ('source_TRAIN_jobs',dict(fit='13415276',readout='13415309',CPU='13415319')),
            ('frozen_local_source_checkpoint','another.pt'),('frozen_local_source_matches_actual_fit',False),
            ('local_source_updated',True),('A3_shared_encoded_once_per_oriented_window',False),
            ('model_parameters',dict(generator=649457,scorer=612324)),
            ('score_kind','A39_FP32_joint_context_ordered_chain_null_relative_gain'),
            ('scopes',good['scopes'][:-1]))
        with tempfile.TemporaryDirectory() as directory,patch.object(V.I,'readiness'),\
             patch.object(Path,'is_file',return_value=True),patch.object(V.I.D.F.A.A,'require_jobs'),\
             patch.object(V.I.D.R,'completed') as completed,patch.object(V.A,'reference_context') as reference:
            for key,value in mutations:
                completed.return_value = {**good,key:value}
                with self.subTest(key=key),self.assertRaisesRegex(ValueError,'8690-window'):
                    V.main(Path(directory)/'fake-infer',Path(directory)/'not-created','synthetic-job')
                self.assertFalse((Path(directory)/'not-created').exists())
            reference.assert_not_called()

    def test_original_pool_local_source_and_FP32_components_include_empty_windows(self):
        scope = dict(species='synthetic',seqid='chr',strand='+',windows=[[0,9],[9,18]])
        empty = dict(raw_row(),window_index=1,window_id='synthetic|chr|9|+',chains=[],score=[],decoder_sources=[],
            exon_additive=[],junction_additive=[],additive=[],nonadditive=[],null=.25,
            counts=dict(free_candidates=0,nonpositive_candidates=0))
        old_empty = dict(baseline_row(),window_index=1,window_id=empty['window_id'],chains=[],score=[],
                         decoder_sources=[],counts=empty['counts'])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); new,old = root/'new.jsonl',root/'old.jsonl'
            old.write_text('\n'.join(json.dumps(r) for r in (baseline_row(),old_empty))+'\n')
            def save(first):
                new.write_text('\n'.join(json.dumps(r) for r in (first,empty))+'\n')
            save(raw_row()); counts = Counter()
            rows = list(V.paired_records(new,scope,old,scope,counts))
            self.assertEqual(len(rows),2); self.assertEqual(rows[1]['chains'],[])
            self.assertEqual(dict(counts),dict(free_candidates=1,nonpositive_candidates=0))
            for key,value in (('frozen_local_source_checkpoint','another.pt'),('variant','A39'),
                    ('chains',[[[0,6]]]),('decoder_sources',[['greedy']]),('score',[3.]),('additive',[2.]),
                    ('score_threshold',.1)):
                save({**raw_row(),key:value})
                with self.subTest(key=key),self.assertRaises((ValueError,AssertionError)):
                    list(V.paired_records(new,scope,old,scope,Counter()))
            new.write_text(json.dumps(raw_row())+'\n')
            with self.assertRaises(ValueError):
                list(V.paired_records(new,scope,old,scope,Counter()))

    def test_A40_negative_strand_GFF_reconstructs_before_reference_access(self):
        row = dict(raw_row(),window_id='synthetic|chr|0|-',chains=[[[0,4],[9,14]]])
        old = dict(baseline_row(),window_id=row['window_id'],chains=row['chains'])
        scope = dict(species='synthetic',seqid='chr',strand='-',ordinal=0,windows=[[0,14]],file='scope.jsonl',
            candidate_counts=row['counts'],assembly_counts=assemble_predictions([row],[[0,14]],'-')['counts'])
        decoded = dict(scopes=[scope],lengths={'chr':14},primary_CDS_chain_predictions=1)
        transcript = as_transcript([(0,5),(10,14)],'chr','-','M28_A40_0_m_0',2.)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); new,prior = root/'new',root/'old'; new.mkdir(); prior.mkdir()
            (new/'scope.jsonl').write_text(json.dumps(row)+'\n'); (prior/'scope.jsonl').write_text(json.dumps(old)+'\n')
            write_gff3(new/'A40.gff3',[transcript],{'chr':14},'M28_A40')
            with patch.object(V.A,'reference_context') as reference:
                predictions,_,_ = V.reconstruct(decoded,dict(scopes=[scope]),new,prior)
                self.assertEqual(predictions['synthetic'],[transcript])
                with self.assertRaisesRegex(ValueError,'prediction count'):
                    V.reconstruct({**decoded,'primary_CDS_chain_predictions':2},dict(scopes=[scope]),new,prior)
                reference.assert_not_called()

    def test_A40_complete_comparator_table_keeps_denominators_macro_and_worst_species(self):
        species = list(V.H.A3.REF_COUNTS)
        per = {s:dict(exact_chain=V.E.S.prf(i+1,i+3,V.H.A3.REF_COUNTS[s]),background=dict(bp=100))
               for i,s in enumerate(species)}
        current = dict(per_species=per,pooled=dict(exact_chain=V.E.S.prf(3,7,7728)))
        names = ('A40','A31_A3','A31_A23','ANNEVO','Tiberius','Helixer','A3','A6','A','R5_C0')
        table = V.comparison_table({name:deepcopy(current) for name in names},current)
        values = [per[s]['exact_chain']['f1'] for s in species]
        self.assertEqual(set(table),set(names)); self.assertEqual(table['A40']['macro_F1'],sum(values)/len(values))
        self.assertEqual(table['A40']['worst_species_F1'],min(values))
        for method,where,field in (('ANNEVO','species','reference'),('Helixer','pooled','reference'),
                ('Tiberius','background','bp')):
            bad = deepcopy(current)
            if where == 'pooled': bad['pooled']['exact_chain'][field] -= 1
            else: bad['per_species'][species[0]]['background' if where == 'background' else 'exact_chain'][field] -= 1
            with self.subTest(method=method),self.assertRaisesRegex(ValueError,'denominator'):
                V.comparison_table(dict(A40=current,**{method:bad}),current)


if __name__ == '__main__':
    unittest.main()
