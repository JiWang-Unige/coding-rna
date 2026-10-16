import unittest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
from evaluate_window_pool import pool_audit
from src.m28.assembly import assemble_predictions
from src.m28.window_pool import assemble_window_pool
from src.m28.export import as_transcript


def row(i,chains,scores):
    return {'window_index':i,'chains':chains,'score':scores}


class WindowPoolTests(unittest.TestCase):
    def test_nonowner_only_chain_is_preserved_without_reference(self):
        records = [row(0,[[(70,79)]],[1.]),row(1,[],[])]
        windows = [(0,100),(40,140)]
        old = assemble_predictions(records,windows,'+')
        new = assemble_window_pool(records,windows,'+')
        self.assertEqual(old['predictions'],[])
        self.assertEqual(new['predictions'],[{'chain':((70,79),),'score':1.,'window_index':0}])
        self.assertEqual(new['counts']['free_unique'],old['counts']['free_unique'])

    def test_global_chain_dedup_max_gain_and_ties_do_not_depend_on_record_order(self):
        windows = [(0,100),(40,140)]
        records = [row(0,[[(70,79)]],[.5]),row(1,[[(30,39)]],[-.5])]
        self.assertEqual(assemble_predictions(records,windows,'+')['predictions'],[])
        chosen = assemble_window_pool(records,windows,'+')
        self.assertEqual(chosen['predictions'][0]['score'],.5)
        self.assertEqual(chosen['counts']['duplicate_free_records'],1)
        records[1]['score'] = [.5]
        forward = assemble_window_pool(records,windows,'+')
        reverse = assemble_window_pool(records[::-1],windows,'+')
        self.assertEqual(forward,reverse)
        self.assertEqual(forward['predictions'][0]['window_index'],0)

    def test_original_global_nonoverlap_positive_gain_selection_is_retained(self):
        records = [row(0,[[(0,9)],[(9,18)],[(0,18)],[(25,34)]],[2.,2.,3.,0.])]
        chosen = assemble_window_pool(records,[(0,100)],'+')
        self.assertEqual([r['chain'] for r in chosen['predictions']],[((0,9),),((9,18),)])
        self.assertEqual(chosen['counts']['window_pool_unique'],4)
        with self.assertRaisesRegex(ValueError,'gain'):
            assemble_window_pool([row(0,[[(0,9)]],[])],[(0,100)],'+')

    def test_reference_audit_uses_window_pool_not_legacy_owner_stage(self):
        assembled = assemble_window_pool([row(0,[[(0,9)],[(20,29)]],[1.,0.])],[(0,100)],'+')
        refs = [{'CDS':[(0,9,'0')]},{'CDS':[(20,29,'0')]},{'CDS':[(40,49,'0')]}]
        audit = pool_audit(assembled,refs)
        self.assertEqual(audit['exact_reference_hits'],
                         {'free_unique':2,'window_pool_unique':2,'selected_unique':1})
        self.assertEqual((audit['not_in_actual_free_candidates'],audit['lost_at_window_pool'],
                          audit['lost_at_conflict_or_nonpositive_gain']),(1,0,1))

    def test_minus_global_coordinates_and_transcript_phase_are_not_reversed(self):
        chosen = assemble_window_pool([row(0,[[(0,4),(9,14)]],[1.])],[(100,200)],'-')
        self.assertEqual(chosen['predictions'][0]['chain'],((186,191),(196,200)))
        tx = as_transcript(chosen['predictions'][0]['chain'],'q','-','synthetic',1.)
        self.assertEqual(tx['CDS'],[(186,191,'2'),(196,200,'0')])


if __name__ == '__main__':
    unittest.main()
