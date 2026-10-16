import sys
import tempfile
import unittest
from pathlib import Path
import torch
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'pipeline/m28/scripts'))
import eval_structure_diagnostic as E
from src.m28.query_output import query_window_candidates
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript, write_gff3
from test_m28_query_decode import fixture


def record(sequence, chain, threshold=0., index=0, max_exons=3):
    model, features = fixture(sequence, chain)
    result = query_window_candidates(model, features, sequence, index,
                                     'synthetic-window-'+str(index),
                                     max_exons, threshold)
    return result


class QueryOutputTests(unittest.TestCase):
    def test_multiexon_query_to_both_strand_gff_roundtrip(self):
        sequence = 'ATGC'+'GTTAG'+'TTTAA'
        row = record(sequence, [(0,4),(9,14)])
        self.assertEqual(row['query_indices'], [0,1])
        self.assertFalse(row['reference_used'])
        self.assertEqual(row['counts']['completed_queries'], 2)
        transcripts = []
        for strand in ('+','-'):
            result = assemble_predictions([row], [(100,114)], strand)
            self.assertEqual(result['counts']['duplicate_owner_records'], 1)
            self.assertEqual(result['counts']['selected_unique'], 1)
            pred = result['predictions'][0]
            tx = as_transcript(pred['chain'], 'q', strand, 'a'+strand, pred['score'])
            expected = ([(100,104,'0'),(109,114,'2')] if strand == '+'
                        else [(100,105,'2'),(110,114,'0')])
            self.assertEqual(tx['CDS'], expected)
            transcripts.append(tx)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'a.gff3'
            write_gff3(path, transcripts, {'q':200}, 'M28A_SYNTHETIC')
            parsed = E.parse_annotation(path, {'q':200}, protein_coding_only=True)
            self.assertEqual(len(parsed['genes']), 2)
            for tx in transcripts:
                self.assertEqual(parsed['transcripts'][tx['id']]['CDS'], tx['CDS'])

    def test_explicit_threshold_keeps_free_records_and_zero_selects_null(self):
        sequence = 'ATGAAATAA'
        for threshold, count in [(0.,1),(1.,0),(2.,0)]:
            row = record(sequence, [(0,9)], threshold)
            self.assertEqual(row['score'], [1.-threshold]*2)
            result = assemble_predictions([row], [(0,9)], '+')
            self.assertEqual(result['counts']['free_records'], 2)
            self.assertEqual(result['counts']['free_unique'], 1)
            self.assertEqual(result['counts']['owner_unique'], 1)
            self.assertEqual(result['counts']['selected_unique'], count)

    def test_unterminated_query_never_becomes_a_gene(self):
        sequence = 'ATGC'+'GTTAG'+'TTTAA'
        row = record(sequence, [(0,4),(9,14)], max_exons=1)
        self.assertEqual(row['chains'], [])
        self.assertEqual(row['score'], [])
        self.assertEqual(row['counts']['unterminated_queries'], 2)
        self.assertTrue(all(x['chain'] == ((0,4),) for x in row['raw_queries']))
        self.assertEqual(assemble_predictions([row], [(0,14)], '+')['predictions'], [])

    def test_no_start_motif_produces_explicit_empty_output(self):
        row = record('CCCAAATAA', [(0,9)])
        self.assertEqual(row['chains'], [])
        self.assertEqual(row['counts']['queries'], 2)
        self.assertTrue(all(x['stop_reason'] == 'no_start_motif' for x in row['raw_queries']))
        self.assertEqual(assemble_predictions([row], [(0,9)], '+')['counts']['free_unique'], 0)

    def test_owner_loss_and_query_duplicates_are_counted_separately(self):
        first = record('C'*15+'ATGAAATAA', [(15,24)])
        windows = [(0,24),(12,36)]
        missing_owner = assemble_predictions([first], windows, '+')
        self.assertEqual(missing_owner['counts']['nonowner_only_unique_lost'], 1)
        self.assertEqual(missing_owner['predictions'], [])
        second = record('CCC'+'ATGAAATAA'+'C'*12, [(3,12)], index=1)
        together = assemble_predictions([first,second], windows, '+')
        self.assertEqual(together['counts']['free_records'], 4)
        self.assertEqual(together['counts']['free_unique'], 1)
        self.assertEqual(together['counts']['owner_records'], 2)
        self.assertEqual(together['counts']['duplicate_owner_records'], 1)
        self.assertEqual(together['counts']['nonowner_only_unique_lost'], 0)
        self.assertEqual(together['predictions'][0]['chain'], ((15,24),))
        self.assertEqual(together['predictions'][0]['window_index'], 1)

    def test_nonfinite_operating_point_or_existence_logit_stops_assembly(self):
        model, features = fixture('ATGAAATAA', [(0,9)])
        for threshold in (float('nan'), float('inf'), float('-inf')):
            with self.assertRaisesRegex(ValueError, 'finite explicit'):
                query_window_candidates(model, features, 'ATGAAATAA', 0, 's', 3, threshold)
        with torch.no_grad():
            model.exists.bias.fill_(float('nan'))
        with self.assertRaisesRegex(ValueError, 'Nonfinite query'):
            query_window_candidates(model, features, 'ATGAAATAA', 0, 's', 3, 0.)


if __name__ == '__main__':
    unittest.main()
