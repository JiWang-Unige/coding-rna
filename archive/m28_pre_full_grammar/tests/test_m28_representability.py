import json
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import diagnose_representability as D


def row(chain, unsupported='', selected=False):
    result = {'species':'plant','id':'r','seqid':'DEV','strand':'-',
              'CDS_intervals':json.dumps(chain),'unsupported_strata':unsupported}
    for arm in D.ARMS:
        result[arm+'_state']='selected' if selected else 'not_free'
        for stage in D.STAGES:
            result[arm+'_'+stage]='True' if selected else 'False'
    return result


class RepresentabilityTests(unittest.TestCase):
    def test_fixed_grid_gap_and_tail_containment_do_not_depend_on_strand(self):
        self.assertEqual(D.geometry(((2,12),),[(0,10),(5,15)]),(0,0,False))
        self.assertEqual(D.geometry(((5,15),),[(0,10),(5,15)]),(1,1,True))
        self.assertEqual(D.geometry(((5,10),),[(0,10),(5,15)]),(2,0,True))

    def test_overlap_flag_is_not_grammar_but_min_intron_and_exon_cap_are(self):
        supported=D.classify(row([(1,6),(12,16)],'same_strand_primary_overlap'),[(0,20)],128,5)
        self.assertEqual(supported['representability_group'],D.GROUPS[2])
        self.assertTrue(supported['annotation_overlap_flag'])
        blocked=D.classify(row([(1,6),(10,15)]),[(0,20)],1,5)
        self.assertEqual(blocked['representability_group'],D.GROUPS[1])
        self.assertEqual(set(blocked['decoder_grammar_reasons'].split(';')),
                         {'more_than_decoder_max_exons','intron_shorter_than_decoder_minimum'})

    def test_groups_partition_missing_and_successful_references_without_filtering(self):
        rows=[D.classify(row([(1,9)]),[(0,20)],128,5),
              D.classify(row([(1,9)],'non_ATG_start'),[(0,20)],128,5),
              D.classify(row([(1,9)],selected=True),[(0,20)],128,5)]
        result=D.summarize(rows)
        self.assertEqual(result['reference_chains'],3)
        self.assertEqual(result['methods']['A3']['missing_from_free'],2)
        self.assertEqual(sum(g['reference_chains'] for g in result['groups'].values()),3)
        self.assertEqual(result['groups'][D.GROUPS[2]]['methods']['A3']['missing_from_free'],1)
        self.assertAlmostEqual(result['methods']['A3']['free_coverage'],1/3)

    def test_observed_exact_hit_contradicts_blocked_representation(self):
        with self.assertRaisesRegex(ValueError,'contradicts representability'):
            D.classify(row([(1,9)],'non_ATG_start',selected=True),[(0,20)],128,5)


if __name__=='__main__':
    unittest.main()
