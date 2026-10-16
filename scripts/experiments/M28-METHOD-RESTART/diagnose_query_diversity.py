#!/usr/bin/env python3
"""Read-only within-window A query diversity; no references, scores or model changes."""
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
INFER = ROOT/'outputs/M28-A-DEV-INFER'
OUT = ROOT/'outputs/M28-A-QUERY-DIVERSITY'


def main():
    inference = json.loads((INFER/'summary.json').read_text())
    reports = []
    for scope in inference['scopes']:
        counts = Counter()
        positive_unique_histogram = Counter()
        all_positive_unique_histogram = Counter()
        with (INFER/scope['file']).open() as handle:
            for line in handle:
                row = json.loads(line)
                chains = [tuple(map(tuple, chain)) for chain in row['chains']]
                positive = [chain for chain, score in zip(chains, row['score']) if score > 0]
                unique = set(chains)
                positive_unique = set(positive)
                counts['windows'] += 1
                counts['completed_queries'] += len(chains)
                counts['unique_completed_chains_within_window'] += len(unique)
                counts['within_window_duplicate_completed_queries'] += len(chains)-len(unique)
                counts['positive_completed_queries'] += len(positive)
                counts['unique_positive_chains_within_window'] += len(positive_unique)
                counts['within_window_duplicate_positive_queries'] += len(positive)-len(positive_unique)
                counts['windows_with_positive_duplicates'] += int(len(positive) != len(positive_unique))
                positive_unique_histogram[len(positive_unique)] += 1
                if len(positive) == 16:
                    counts['all_16_queries_positive_completed_windows'] += 1
                    all_positive_unique_histogram[len(positive_unique)] += 1
        q = scope['query_counts']
        if (counts['windows'] != len(scope['windows']) or
            counts['completed_queries'] != q['completed_queries'] or
            counts['positive_completed_queries'] != q['completed_queries']-q['completed_nonpositive_gain'] or
            counts['all_16_queries_positive_completed_windows'] != q['all_queries_positive_completed']):
            raise ValueError('Diversity source/count mismatch with completed inference')
        reports.append({'species':scope['species'],'strand':scope['strand'],
                        'source':scope['file'],'counts':dict(counts),
                        'positive_unique_chain_count_histogram':dict(sorted(positive_unique_histogram.items())),
                        'all_16_positive_unique_chain_count_histogram':dict(sorted(all_positive_unique_histogram.items()))})
    totals = Counter()
    for report in reports:
        totals.update(report['counts'])
    result = {'status':'completed_read_only_existing_A_query_diversity',
              'checkpoint':inference['checkpoint'],'reference_used':False,
              'test_setaria_used':False,'scores_or_candidates_changed':False,
              'positive_definition':'existing_frozen_gain_gt_0_not_biological_truth',
              'duplicate_definition':'identical_ordered_half_open_exon_chains_within_one_oriented_window',
              'scope_reports':reports,'totals':dict(totals)}
    with (OUT/'summary.json').open('x') as handle:
        json.dump(result,handle,indent=2); handle.write('\n')
    print(json.dumps(result),flush=True)


if __name__ == '__main__':
    main()
