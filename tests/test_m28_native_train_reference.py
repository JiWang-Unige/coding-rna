import copy
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import diagnose_native_train_reference as D


def reference(name='r',seqid='TRAIN',strand='+',species='plant',hit=False):
    return dict(id=name,gene_id='g'+name,species=species,seqid=seqid,strand=strand,
                CDS=[[10,19,'0'],[30,39,'0']],CDS_segments=2,span_bp=29,
                exact_primary_hit=hit,original_draw_R4_positive_chain=False)


def error(ref):
    return {**{k:ref[k] for k in ('id','species','seqid','strand','CDS_segments')},
            'category':'exact_primary_chain' if ref['exact_primary_hit'] else 'no_same_strand_CDS_overlap'}


def decorated(model,native):
    ref = reference(hit=model)
    ref['native_exact_hits'] = dict(zip(D.METHODS,native))
    ref['original_model_error_category'] = error(ref)['category']
    return ref


class NativeTrainReferenceTests(unittest.TestCase):
    def test_geometry_key_preserves_species_seqid_strand_not_phase_or_id(self):
        ref = reference()
        changed = copy.deepcopy(ref); changed['CDS'][0][2]='2'; changed['id']='other'
        self.assertEqual(D.reference_key(ref),D.reference_key(changed))
        for key,value in (('species','rice'),('seqid','other'),('strand','-')):
            self.assertNotEqual(D.reference_key(ref),D.reference_key({**ref,key:value}))

    def test_same_transcript_id_on_distinct_contigs_and_strands_pairs_correctly(self):
        refs = [reference(),reference(seqid='TRAIN2'),reference(strand='-',hit=True)]
        rows = D.join_profiles(refs,[error(r) for r in reversed(refs)],{('plant','TRAIN'):2,('plant','TRAIN2'):1})
        self.assertEqual([r['original_model_error_category'] for r in rows],
                         ['no_same_strand_CDS_overlap']*2+['exact_primary_chain'])

    def test_missing_duplicated_or_wrong_scope_reference_is_rejected(self):
        ref = reference(); err = error(ref); scope = {('plant','TRAIN'):1}
        for refs,errors,counts in (([ref],[],scope),([ref],[err,err],scope),
                                  ([ref,ref],[err],scope),([ref],[err],{('plant','TRAIN'):2})):
            with self.assertRaises(ValueError): D.join_profiles(refs,errors,counts)

    def test_geometry_duplicate_with_different_id_is_not_new_reference(self):
        refs = [reference('a'),reference('b')]
        with self.assertRaises(ValueError):
            D.join_profiles(refs,[error(r) for r in refs],{('plant','TRAIN'):2})

    def test_error_categories_and_hit_flags_are_bound_to_original_profile(self):
        ref = reference(); err = error(ref)
        for change in ({'category':'exact_primary_chain'},{'CDS_segments':1},{'category':'unknown'}):
            with self.assertRaises(ValueError):
                D.join_profiles([ref],[{**err,**change}],{('plant','TRAIN'):1})

    def test_reference_pairs_patterns_and_consensus_use_the_full_denominator(self):
        rows = [decorated(True,(True,True,True)),decorated(False,(True,True,True)),
                decorated(True,(False,False,False)),decorated(False,(False,False,False)),
                decorated(False,(True,False,False))]
        result = D.summarize(rows)
        self.assertEqual(result['reference_chains'],5)
        self.assertEqual(result['native_vs_model']['ANNEVO']['membership'],
                         dict(both=1,native_only=2,model_only=1,neither=1))
        self.assertEqual(result['any_native_vs_model'],dict(both=1,native_only=2,model_only=1,neither=1))
        self.assertEqual(result['all_three_native_exact_and_model_missed'],1)
        self.assertEqual(result['exact_primary_hits'],dict(original4608=2,ANNEVO=3,Tiberius=2,Helixer=2))
        self.assertEqual(result['reference_recall']['ANNEVO'],3/5)
        self.assertEqual(len(result['membership_patterns']),16)
        self.assertEqual(sum(p['reference_chains'] for p in result['membership_patterns']),5)

    def test_strata_and_empty_denominators_do_not_drop_difficult_references(self):
        ref = decorated(False,(True,False,False)); ref['span_bp']=8000
        result = D.describe([ref])
        for strata in result['strata'].values():
            self.assertEqual(sum(v['reference_chains'] for v in strata.values()),1)
        self.assertEqual(result['strata']['span_bp']['over_6144']['reference_chains'],1)
        empty = result['strata']['CDS_segments']['single']
        self.assertIsNone(empty['reference_recall']['ANNEVO'])
        self.assertEqual(empty['reference_chains'],0)

    def test_unique_prediction_replay_rejects_recount_and_original_hit_drift(self):
        ref = reference(hit=True); key = D.reference_key(ref)
        self.assertEqual(D.exact_counts(set([key,key]),{key}),dict(tp=1,predicted=1,reference=1))
        D.reconcile_hits([decorated(True,(True,False,False))],dict(tp=1,reference=1))
        with self.assertRaises(ValueError):
            D.reconcile_hits([decorated(True,(True,False,False))],dict(tp=0,reference=1))


if __name__ == '__main__': unittest.main()
