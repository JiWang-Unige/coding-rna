import unittest
import torch
from src.m28.candidates import ORFIndex
from src.m28.c0_decode import genomic_chain
from src.m28.query_decode import decode_query_chains, _end_events, _legal_ends
from test_m28_query_pointer import scripted_geometry


def fixture(sequence, chain):
    model, _ = scripted_geometry(max(9, len(sequence)))
    features = torch.zeros(len(sequence), 4)
    for j, (a, b) in enumerate(chain):
        features[a, 0] = 3-j
        features[b-1, 1] = 3-j
        features[b-1, 2] = 1 if j == len(chain)-1 else -1
    return model, features


class QueryCDSGrammarTests(unittest.TestCase):
    def test_gt_and_gc_multiexon_complete_cds_without_reference(self):
        for intron in ('GTTAG', 'GCTAG'):
            sequence = 'ATGC' + intron + 'TTTAA'
            model, features = fixture(sequence, [(0,4),(9,14)])
            rows = decode_query_chains(model, features, sequence, max_exons=3)
            for row in rows:
                self.assertEqual(row['chain'], ((0,4),(9,14)))
                self.assertTrue(row['terminated'])
                self.assertEqual(row['stop_reason'], 'canonical_stop')
                self.assertEqual(''.join(sequence[a:b] for a,b in row['chain']), 'ATGCTTTAA')

    def test_same_phase_different_letters_detect_junction_stop(self):
        # Both prefixes carry one base; T+AA is STOP whereas A+AA is not.
        for base, allowed in [('T', False), ('A', True)]:
            sequence = 'ATG' + base + 'GTTAG' + 'AAAAATAA'
            ends = _legal_ends(ORFIndex(sequence), _end_events(sequence), 9, base, False)
            self.assertEqual((17, 'stop', '') in ends, allowed)

    def test_four_base_intron_is_not_silently_accepted(self):
        sequence = 'ATGC' + 'GTAG' + 'TTTAA'
        model, features = fixture(sequence, [(0,4),(8,13)])
        row = decode_query_chains(model, features, sequence, max_exons=3)[0]
        self.assertEqual(row['chain'], ((0,4),))
        self.assertFalse(row['terminated'])
        self.assertEqual(row['stop_reason'], 'no_acceptor')

    def test_exon_cap_cannot_switch_to_an_alternate_stop(self):
        sequence = 'ATGC' + 'GTTAG' + 'TTTAA'
        model, features = fixture(sequence, [(0,4),(9,14)])
        row = decode_query_chains(model, features, sequence, max_exons=1)[0]
        self.assertEqual(row['chain'], ((0,4),))
        self.assertFalse(row['terminated'])
        self.assertEqual(row['stop_reason'], 'exon_limit')

    def test_ambiguous_coding_rejected_but_skipped_intron_is_allowed(self):
        sequence = 'ATGNAATAA'
        model, features = fixture(sequence, [(0,9)])
        row = decode_query_chains(model, features, sequence, max_exons=2)[0]
        self.assertFalse(row['terminated'])
        self.assertEqual(row['stop_reason'], 'no_legal_exon_end')
        sequence = 'ATGC' + 'GTNAG' + 'TTTAA'
        model, features = fixture(sequence, [(0,4),(9,14)])
        row = decode_query_chains(model, features, sequence, max_exons=3)[0]
        self.assertEqual(row['chain'], ((0,4),(9,14)))
        self.assertTrue(row['terminated'])

    def test_split_start_and_stop_are_not_misreported_as_supported(self):
        sequence = 'AT' + 'GTTAG' + 'GAAATAA'
        self.assertNotIn((2, 'donor', 'AT'),
                         _legal_ends(ORFIndex(sequence), _end_events(sequence), 0, '', True))
        sequence = 'ATGT' + 'GTTAG' + 'AA'
        model, features = fixture(sequence, [(0,4),(9,11)])
        row = decode_query_chains(model, features, sequence, max_exons=3)[0]
        self.assertEqual(row['chain'], ((0,4),))
        self.assertFalse(row['terminated'])
        self.assertEqual(row['stop_reason'], 'no_legal_exon_end')

    def test_softmask_padding_and_coordinate_alignment(self):
        sequence = 'ATGC' + 'GTTAG' + 'TTTAA'
        model, features = fixture(sequence, [(0,4),(9,14)])
        expected = decode_query_chains(model, features, sequence, max_exons=3)
        padded = torch.cat([features, torch.full((6,4), 10000.)])
        self.assertEqual(expected, decode_query_chains(model, padded, sequence.lower(), 3, valid_bases=14))
        with self.assertRaisesRegex(ValueError, 'DNA and valid feature'):
            decode_query_chains(model, features, sequence[:-1], 3)

    def test_oriented_chain_maps_both_strands_with_half_open_coordinates(self):
        chain = ((0,4),(9,14))
        self.assertEqual(genomic_chain(chain, 100, 114, '+'), [(100,104),(109,114)])
        self.assertEqual(genomic_chain(chain, 100, 114, '-'), [(100,105),(110,114)])

    def test_all_legal_endpoints_are_available_and_eos_is_joint(self):
        sequence = 'ATG' + 'GTG'*10 + 'TAA'
        model, features = fixture(sequence, [(0,len(sequence))])
        ends = _legal_ends(ORFIndex(sequence), _end_events(sequence), 0, '', True)
        self.assertGreater(len(ends), 8)
        with torch.no_grad():
            # STOP ranks below every donor on boundary scores alone, but
            # donor continuation is disfavored by its prospective EOS state.
            for b, kind, _ in ends:
                features[b-1, 1] = 3 if kind == 'donor' else 1
                features[b-1, 2] = 1
        row = decode_query_chains(model, features, sequence, max_exons=2)[0]
        self.assertEqual(row['chain'], ((0,len(sequence)),))
        self.assertTrue(row['terminated'])

    def test_no_atg_and_low_eos_are_distinct_from_forced_truncation(self):
        sequence = 'CCCAAATAA'
        model, features = fixture(sequence, [(0,9)])
        row = decode_query_chains(model, features, sequence, max_exons=2)[0]
        self.assertEqual(row['chain'], ())
        self.assertFalse(row['terminated'])
        self.assertEqual(row['stop_reason'], 'no_start_motif')
        sequence = 'ATGAAATAA'
        model, features = fixture(sequence, [(0,9)])
        features[-1, 2] = -1
        # The only DNA-legal end event is STOP: conditional joint argmax,
        # not the raw pointer prototype's EOS > 0 decision rule.
        row = decode_query_chains(model, features, sequence, max_exons=2)[0]
        self.assertEqual(row['chain'], ((0,9),))
        self.assertTrue(row['terminated'])


if __name__ == '__main__':
    unittest.main()
