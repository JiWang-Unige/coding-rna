import sys
import unittest
from pathlib import Path
import numpy as np
import torch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.m28.labels import build_targets,candidate_targets,local_loss

def fixture(partial=False,strand="+"):
    cds=[(3,8,"0"),(14,21,"1")] if strand=="+" else [(3,10,"1"),(16,21,"0")]
    tx={"id":"p","gene_id":"g","CDS":cds,"training_primary":True,
        "CDS_complete_supported":not partial,"strand":strand}
    w={"start":0,"end":24,"valid_bases":24,"width":30,"strand":strand,
       "overlapping_transcript_ids":["p"],"positive_ids":[] if partial else ["p"],
       "fine_label_unknown_intervals":[[3,21]] if partial else [],
       "chain_unknown_intervals":[[3,21]] if partial else [],
       "no_annotated_gene_either_strand":False}
    return w,{"p":tx},"CCCATGAAGTCCAGACCCTAAGGG"

class LabelTests(unittest.TestCase):
    def test_phase_boundaries_and_padding(self):
        w,c,dna=fixture();t=build_targets(w,c,dna)
        self.assertEqual(t["seg"][3:6].tolist(),[13,13,13])
        self.assertEqual(int(t["seg"][7]),8) # donor, cumulative CDS phase1
        self.assertEqual(t["seg"][8:14].tolist(),[5]*6)
        self.assertEqual(int(t["seg"][14]),12) # acceptor, cumulative phase2
        self.assertEqual(t["seg"][18:21].tolist(),[14]*3)
        self.assertTrue((t["seg"][24:]==-100).all())
        self.assertFalse(t["endpoint_known"][24:].any())
        self.assertEqual(candidate_targets([[(3,8),(14,21)],[(3,8)]],t),[1,0])

    def test_reverse_coordinate_equivalence(self):
        w,c,dna=fixture();p=build_targets(w,c,dna)
        w,c,dna=fixture(strand="-");m=build_targets(w,c,dna)
        np.testing.assert_array_equal(p["seg"],m["seg"])
        np.testing.assert_array_equal(p["endpoint"],m["endpoint"])

    def test_partial_preserves_local_evidence_not_false_empty(self):
        w,c,dna=fixture(partial=True);t=build_targets(w,c,dna)
        self.assertTrue((t["seg"][3:21]==-100).all())
        self.assertEqual(int(t["region"][7]),1)
        self.assertEqual(int(t["region"][9]),2)
        self.assertEqual(float(t["endpoint"][7,2]),1.)
        self.assertTrue(t["endpoint_known"][7,2])
        self.assertFalse(t["endpoint_known"][7,0])
        self.assertFalse(t["reference_empty_window"])
        self.assertEqual(candidate_targets([[(3,8)]],t),[-1])
        logits=torch.randn(30,15,requires_grad=True)
        loss=local_loss(logits,t);loss.backward()
        self.assertGreater(float(logits.grad[9].abs().sum()),0.)

    def test_alternative_conflict_and_negative_exclusion(self):
        w,c,dna=fixture()
        c["a"]=dict(c["p"],id="a",training_primary=False,CDS=[(3,21,"0")])
        w["overlapping_transcript_ids"].append("a")
        t=build_targets(w,c,dna)
        self.assertEqual(int(t["region"][9]),-100)
        self.assertEqual(candidate_targets([[(3,21)]],t),[-1])

    def test_reference_background_requires_callable_sequence(self):
        w,c,dna=fixture()
        w.update(overlapping_transcript_ids=[],positive_ids=[],no_annotated_gene_either_strand=True)
        self.assertTrue(build_targets(w,{},dna)["reference_empty_window"])
        t=build_targets(w,{},"N"+dna[1:])
        self.assertFalse(t["reference_empty_window"])
        self.assertEqual(int(t["seg"][0]),-100)
        self.assertEqual(candidate_targets([[(0,9)]],t),[-1])

if __name__=="__main__": unittest.main()
