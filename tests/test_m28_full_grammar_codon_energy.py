"""Score decomposition checks on artificial DNA, never research inputs."""
import copy
import unittest
import numpy as np
import torch
import diagnose_full_grammar_codon_energy as E
from src.m28.full_grammar_codon_skyline import decode_full_grammar_codon


def fixture(hit=False):
    seq, parts = 'ATGCGTNAGTTTAA', ((0,4),(9,14))
    coefficients = E.centered_codon_scores(torch.arange(64).float()/7).numpy()
    values = np.arange(len(seq)*15, dtype=np.float32).reshape(-1,15)/13
    terms = E.chain_terms(seq, values, parts, coefficients)
    row = dict(species='fixture', seqid='fixture', strand='+', reference_CDS=[[0,4],[9,14]],
        single_chain_supported=True, exact_primary_hit=hit, category='exact_primary_chain' if hit else 'scoring_miss',
        reference_gain_vs_background=terms['gain'], **{k:v for k,v in terms.items() if k != 'gain'},
        **E.margin_terms(terms, [terms] if hit else []))
    return row, coefficients


class CodonEnergyTests(unittest.TestCase):
    def test_split_codon_count_and_closure_match_production(self):
        fixtures = (
            ('ATGCGTNAGTTTAA', ((0,4),(9,14)), ('ATG','CTT','TAA')),
            ('ATGCGTTAGAGTTAGAGCTAGTAA', ((0,4),(9,10),(15,16),(21,24)), ('ATG','CAA','TAA')),
            ('ATGGCTAGAAATAA', ((0,3),(8,14)), ('ATG','AAA','TAA')),
            ('ATGAAATAG', ((0,9),), ('ATG','AAA','TAG')),
            ('ATGAAATGA', ((0,9),), ('ATG','AAA','TGA')))
        raw = torch.arange(64).double()/7
        coefficients = E.centered_codon_scores(raw).numpy()
        for sequence, parts, words in fixtures:
            path = E.D.role_path(parts)
            values = torch.zeros(len(sequence),15, dtype=torch.float64)
            values[torch.arange(len(path)), torch.from_numpy(path).long()] = 100.
            self.assertTrue(E.D.supported(sequence, parts, path, len(parts)))
            terms = E.chain_terms(sequence, values.numpy(), parts, coefficients)
            mapped = decode_full_grammar_codon(sequence, values, len(parts), raw, 3)
            self.assertEqual(mapped['chains'], (parts,))
            self.assertAlmostEqual(terms['gain'], mapped['chain_gains'][0], places=9)
            self.assertEqual(sum(terms['codon_counts']), len(words))
            for word in words:
                self.assertEqual(terms['codon_counts'][E.CODON_INDEX[word]], words.count(word))

    def test_zero_coefficients_and_reverse_oriented_geometry(self):
        row, _ = fixture()
        sequence = 'ATGCGTNAGTTTAA'
        genomic = sequence.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
        tx = dict(strand='-', CDS=[(len(sequence)-b,len(sequence)-a,'0') for a,b in reversed(row['reference_CDS'])])
        parts = E.D.oriented_parts(tx, len(sequence))
        dna = genomic.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
        values = np.arange(len(dna)*15, dtype=np.float32).reshape(-1,15)/13
        terms = E.chain_terms(dna, values, parts, np.zeros(64))
        parent, roles, counts = E.D.chain_score(values, parts, E.D.role_path(parts))
        self.assertEqual(terms['gain'], parent)
        self.assertEqual(terms['gain_by_role'], roles)
        self.assertEqual(terms['rows_by_role'], counts)
        self.assertEqual(terms['codon_gain'], 0.)

    def test_entire_multiple_conflicts_and_empty_set_margins(self):
        row, coefficients = fixture()
        reference = dict(gain=row['reference_gain_vs_background'], **{k:row[k] for k in ('gain_by_role','rows_by_role','codon_counts','codon_gain_by_type','codon_gain')})
        margin = E.margin_terms(reference, [reference, reference])
        self.assertAlmostEqual(margin['replacement_gain_difference'], -reference['gain'])
        self.assertEqual(margin['conflicting_selected_codon_counts'], [2*n for n in reference['codon_counts']])
        self.assertEqual(E.margin_terms(reference, [])['replacement_gain_difference'], reference['gain'])
        for intervals, expected in ((((3,18),),2), (((9,12),),0)):
            index = ([dict(parts=[(0,9)]),dict(parts=[(12,21)])],[9,21])
            self.assertEqual(len(E.D.conflicting(intervals, index)), expected)
        row.update(margin)
        E.check_row(row, coefficients)

    def test_wrong_score_codon_count_or_margin_is_rejected(self):
        original, coefficients = fixture()
        for field, index in (('reference_gain_vs_background',None), ('codon_counts', E.CODON_INDEX['ATG']),
                             ('replacement_margin_by_role',4), ('replacement_margin_by_codon', E.CODON_INDEX['CTT'])):
            row = copy.deepcopy(original)
            if index is None:
                row[field] += 1.
            else:
                row[field][index] += 1.
            with self.assertRaises(ValueError):
                E.check_row(row, coefficients)

    def test_four_hit_groups_keep_unsupported_and_same_reference_codons(self):
        control, candidate = {}, {}
        for i, (a,b) in enumerate(((True,True),(True,False),(False,True),(False,False))):
            for mapping, hit in ((control,a),(candidate,b)):
                row, _ = fixture(hit)
                if i == 3:
                    row = dict(single_chain_supported=False, exact_primary_hit=False, category='unsupported',
                               rows_by_role=None, codon_counts=None, reference_gain_vs_background=None)
                mapping[i] = row
        paired = E.paired_summary(control,candidate)
        self.assertEqual((paired['reference'],paired['gained'],paired['lost'],paired['net_hits']), (4,1,1,0))
        self.assertTrue(all(r['reference'] == 1 for r in paired['hit_groups'].values()))
        self.assertEqual(paired['hit_groups']['both_missed']['supported_reference'],0)
        with self.assertRaises(ValueError):
            E.paired_summary(control,{k:v for k,v in candidate.items() if k != 3})
        changed = copy.deepcopy(candidate)
        changed[2]['codon_counts'][E.CODON_INDEX['ATG']] += 1
        with self.assertRaises(ValueError):
            E.paired_summary(control,changed)

    def test_fixed_set_without_codon_score_is_not_a_new_prediction(self):
        row, _ = fixture()
        raw = torch.zeros(64)
        raw[E.CODON_INDEX['CTT']] = -610.
        coefficients = E.centered_codon_scores(raw).numpy()
        sequence, parts = 'ATGCGTNAGTTTAA', ((0,4),(9,14))
        values = np.ones((len(sequence),15), dtype=np.float32)
        values[:,0] = 0.
        terms = E.chain_terms(sequence,values,parts,coefficients)
        row.update({k:v for k,v in terms.items() if k != 'gain'})
        row.update(reference_gain_vs_background=terms['gain'], **E.margin_terms(terms,[]))
        E.check_row(row,coefficients)
        self.assertLess(row['replacement_gain_difference'],0.)
        self.assertGreater(row['role_only_fixed_set_margin'],0.)
        paired = E.paired_summary({0:row}, {0:row})
        group = paired['hit_groups']['both_missed']['per_arm'][E.ARMS[1]]
        self.assertEqual(group['full_margin_nonpositive_but_role_only_margin_positive'],1)
        self.assertEqual(paired['gained'],0)
        self.assertEqual(paired['net_hits'],0)

    def test_incomplete_ambiguous_and_unsupported_scoring_is_rejected(self):
        for sequence in ('ATGNAATAA','ATGAAATA'):
            with self.assertRaises(ValueError):
                E.chain_terms(sequence,np.zeros((len(sequence),15)),((0,len(sequence)),),np.zeros(64))
        row = dict(single_chain_supported=False, exact_primary_hit=False, reference_gain_vs_background=None)
        E.check_row(row,np.zeros(64))
        row['reference_gain_vs_background'] = 1.
        with self.assertRaises(ValueError):
            E.check_row(row,np.zeros(64))
