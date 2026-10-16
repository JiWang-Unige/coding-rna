import json
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import diagnose_query_context as D
from src.m28.assembly import assemble_predictions
from test_m28_native_reference import raw


def baseline(selected=False):
    row = raw('True','True','True','1','selected') if selected else raw()
    return {**row,**{m+'_exact_primary':'False' for m in D.METHODS}}


class ContextDiagnosticTests(unittest.TestCase):
    def test_native_false_is_false_and_invalid_csv_membership_is_rejected(self):
        row = D.parse_baseline(baseline())
        self.assertTrue(all(row[m+'_exact_primary'] is False for m in D.METHODS))
        with self.assertRaises(ValueError):
            D.parse_baseline({**baseline(),'ANNEVO_exact_primary':'0'})

    def test_gross_new_and_lost_hits_native_gap_and_complete_strata(self):
        rows = []
        for old,new in ((True,True),(False,True),(True,False),(False,False)):
            row = D.parse_baseline(baseline(old)); key = D.chain_key(row['CDS'])
            sets = {s:{key} if new else set() for s in D.STAGES}
            D.add_context(row,sets,{key:1.} if new else {})
            row['ANNEVO_exact_primary'] = not new
            rows.append(row)
        result = D.summarize(rows)
        for stage in D.STAGES:
            self.assertEqual(result['stages'][stage],dict(both=1,A10_only=1,A3_only=1,neither=1))
        self.assertEqual(result['native_vs_A10']['ANNEVO']['native_only_by_A10_state'],{'not_free':2})
        self.assertEqual(sum(r['reference_chains'] for r in result['A3_to_A10_states']),4)
        self.assertTrue(result['not_subset_precision_or_F1'])
        self.assertEqual([D.exon_stratum(n) for n in (1,2,4,5,9,10)],['1','2-4','2-4','5-9','5-9','10+'])
        self.assertEqual([D.span_stratum(n) for n in (6144,6145,24576,24577)],['<=6144','6145-24576','6145-24576','>24576'])

    def test_replay_uses_max_duplicate_owner_gain_and_preserves_nonowner_only_loss(self):
        windows = [(0,20),(10,30)]; chains = [((1,10),),((1,10),),((14,20),)]
        def record(index,chains,scores):
            return dict(window_index=index,window_id='plant|DEV|'+str(windows[index][0])+'|+',
                chains=chains,score=scores,raw_queries=[{}]*16,reference_used=False,
                score_kind='query_existence_logit_minus_threshold',existence_logit_threshold=0.,
                counts=dict(queries=16,completed_queries=len(chains),completed_nonpositive_gain=sum(x<=0 for x in scores)))
        records = [record(0,chains,[-2.,3.,4.]),record(1,[],[])]
        expected = assemble_predictions(records,windows,'+')
        scope = dict(species='plant',seqid='DEV',strand='+',windows=windows,
            query_counts=dict(completed_queries=3,completed_nonpositive_gain=1),assembly_counts=expected['counts'])
        reference = {((1,10),),((14,20),)}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'raw.jsonl'
            path.write_text(''.join(json.dumps(r)+'\n' for r in records))
            assembled,gains = D.replay_scope(path,scope,reference)
        self.assertEqual(gains,{((1,10),):3.})
        sets = {s:{D.chain_key(c) for c in assembled['stages'][s]} for s in D.STAGES}
        row = D.parse_baseline(baseline()); row['CDS'] = [(14,20)]
        D.add_context(row,sets,gains)
        self.assertEqual(row['A10_state'],'lost_at_owner')
        self.assertFalse(row['A10_owner_unique'])
        with self.assertRaisesRegex(ValueError,'recorded gain'):
            missing = {**sets,'selected_unique':set()}
            D.add_context({**row,'CDS':[(1,10)]},missing,{})


if __name__ == '__main__':
    unittest.main()
