"""New three-arm pairing and unchanged stronger-baseline gate only."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import readout_role_evidence_context_train as D


class Pairing(unittest.TestCase):
    def fixture(self):
        new,old,lab,prior=[],[],[],[]
        for species in ('arabidopsis_thaliana','oryza_sativa'):
            w=dict(id=species+'|toy',species=species,split='train')
            chains=[[[0,6]],[[12,18]],[[24,30]]]; labels=[1,0,-1]
            base=dict(chains=chains,scores=[2.,1.,1.],selected=[0,1,2])
            new.append(dict(window=w,chains=chains,scores=[2.,0.,1.],selected=[0,2],reference_used_in_scoring=False))
            prior.append(dict(window=w,**base,reference_used_in_scoring=False))
            old.append(dict(window=w,variants={'A31_A3':base},reference_used_in_scoring=False))
            lab.append(dict(window=w,labels={'A31_A3':labels},counts={'A31_A3':D.R.selected_counts(base,labels)},
                eligible_references=2,eligible_five_or_more=0))
        return tuple(deepcopy(v) for v in (new,old,lab,prior))

    def test_pairing_preserves_unknown_denominators_and_rejects_bad_sources(self):
        args=self.fixture(); rows,counts,refs=D.paired_readout(*args)
        self.assertEqual(len(rows),2)
        for species in counts:
            self.assertEqual(counts[species]['A40']['selected_unknown'],1)
            self.assertEqual(counts[species]['A40']['selected_known_nonmatch'],0)
            self.assertEqual(refs[species],dict(eligible=2,five_or_more=0))
        for index in (0,1,3):
            changed=deepcopy(args); changed[index][0]['reference_used_in_scoring']=True
            with self.assertRaisesRegex(ValueError,'reference-free'): D.paired_readout(*changed)
        missing=deepcopy(args); missing[3].pop()
        with self.assertRaisesRegex(ValueError,'four-way'): D.paired_readout(*missing)
        changed=deepcopy(args); changed[3][0]['chains'][0][0][1]=7
        with self.assertRaisesRegex(ValueError,'identity'): D.paired_readout(*changed)

    def test_gate_never_substitutes_closed_A39_or_accepts_equal_correct(self):
        counts={}
        for species,correct,nonmatch,long in (('arabidopsis_thaliana',1531,347,548),('oryza_sativa',478,251,152)):
            base=dict(selected_correct=correct,selected_known_nonmatch=nonmatch,selected_five_or_more_correct=long)
            counts[species]=dict(A31_A3=base,A39={**base,'selected_correct':correct-100},A40=dict(base))
        actual=D.gate(counts)
        self.assertEqual(sum(actual.values()),2)
        for values in counts.values():
            values['A40']['selected_correct']+=1; values['A40']['selected_five_or_more_correct']+=1
        self.assertTrue(all(D.gate(counts).values()))
        counts['oryza_sativa']['A40']['selected_known_nonmatch']+=1
        self.assertFalse(all(D.gate(counts).values()))


if __name__=='__main__': unittest.main()
