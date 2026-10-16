import copy
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import diagnose_native_reference as D


def raw(free='False',owner='False',selected='False',gain='',state='not_free'):
    return {'species':'plant','id':'r','gene_id':'g','seqid':'DEV','strand':'-',
            'exon_count':'2','span_bp':'26','coding_bp':'15','CDS_intervals':'[[20,29],[40,46]]',
            'unsupported_strata':'longer_than_window','A3_free_unique':free,
            'A3_owner_unique':owner,'A3_selected_unique':selected,'A3_owner_gain':gain,'A3_state':state}


class NativeReferenceTests(unittest.TestCase):
    def test_csv_false_state_and_nested_membership(self):
        row = D.parse_record(raw())
        self.assertFalse(row['A3_free_unique'])
        self.assertEqual(row['unsupported_strata'],'longer_than_window')
        self.assertEqual((row['exon_count'],row['span_bp'],row['coding_bp']),(2,26,15))
        for changed in ({'A3_free_unique':'0'},{'A3_owner_unique':'True'},
                        {'A3_state':'selected'}):
            with self.assertRaises(ValueError):
                D.parse_record({**raw(),**changed})

    def test_identity_retains_strand_species_and_ignores_phase(self):
        row = D.parse_record(raw())
        phased = copy.deepcopy(row)
        phased['CDS'] = [[20,29,'1'],[40,46,'2']]
        self.assertEqual(D.reference_key(row),D.reference_key(phased))
        for field,value in (('strand','+'),('species','other'),('seqid','other')):
            changed = {**row,field:value}
            self.assertNotEqual(D.reference_key(row),D.reference_key(changed))

    def test_gross_membership_and_native_only_attrition(self):
        records = []
        for a,n in ((True,True),(False,True),(True,False),(False,False)):
            row = D.parse_record(raw('True','True','True','1','selected') if a else raw())
            row.update({m+'_exact_primary':n for m in D.METHODS})
            records.append(row)
        report = D.summarize(records)
        self.assertEqual(report['reference_chains'],4)
        for item in report['native_vs_A3'].values():
            self.assertEqual(item['membership'],{'both':1,'native_only':1,'A3_only':1,'neither':1})
            self.assertEqual(item['native_only_by_A3_state'],{'not_free':1})
        self.assertEqual(sum(p['reference_chains'] for p in report['selected_membership_patterns']),4)
        self.assertEqual(report['A3_selected_and_no_native_exact_hit'],1)
        self.assertEqual(report['any_native_exact_hit_absent_from_A3_free'],1)


if __name__ == '__main__':
    unittest.main()
