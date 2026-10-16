"""Artificial joins for original local labels and closed finite supervision."""
from copy import deepcopy
import unittest

from src.m28.local_structured_emission import reconcile_local_supervision


def source(strand='+'):
    window = dict(id='toy', species='toy_species', split='train', seqid='toy_chr',
        strand=strand, start=100, end=124, width=30, valid_bases=24,
        label_policy='CDS_assessable_primary_R4', positive_ids=['positive'],
        overlapping_transcript_ids=['positive'], fine_label_unknown_intervals=[(3, 5)])
    chain = [(2, 6), (8, 14)]
    cds = [(102, 106, 0), (108, 114, 2)] if strand == '+' else [(110, 116, 0), (118, 122, 2)]
    catalog = {'positive': dict(id='positive', CDS=cds)}
    record = dict(window={k: window[k] for k in ('id', 'species', 'split', 'seqid',
        'strand', 'start', 'end', 'width', 'valid_bases')}, candidates=[
        dict(chain=chain, label=1), dict(chain=[(15, 21)], label=-1)])
    return record, window, catalog


class LocalStructuredSourceTests(unittest.TestCase):
    def test_both_orientations_match_without_mutating_masks_or_unknowns(self):
        for strand in ('+', '-'):
            values = source(strand); before = deepcopy(values)
            reconcile_local_supervision(*values, 'original.pt', 'original.pt')
            self.assertEqual(values, before)

    def test_window_or_actual_cache_mismatch_rejected(self):
        for field, value in (('end', 125), ('strand', '-'), ('split', 'dev')):
            record, window, catalog = source()
            record['window'][field] = value
            with self.assertRaisesRegex(ValueError, 'geometry'):
                reconcile_local_supervision(record, window, catalog, 'original.pt', 'original.pt')
        with self.assertRaisesRegex(ValueError, 'cache source'):
            reconcile_local_supervision(*source(), 'original.pt', 'other.pt')

    def test_complete_gold_not_outer_span_or_positive_count(self):
        for change in ('different_internal_exon', 'missing_positive', 'extra_positive'):
            record, window, catalog = source()
            if change == 'different_internal_exon':
                record['candidates'][0]['chain'] = [(2, 5), (8, 14)]
            elif change == 'missing_positive':
                record['candidates'][0]['label'] = -1
            else:
                record['candidates'][1]['label'] = 1
            with self.assertRaisesRegex(ValueError, 'gold'):
                reconcile_local_supervision(record, window, catalog, 'original.pt', 'original.pt')


if __name__ == '__main__':
    unittest.main()
