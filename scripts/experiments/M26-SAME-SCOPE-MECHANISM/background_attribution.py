#!/usr/bin/env python3
"""Attribute fixed predicted-span burden to disjoint reference-policy regions."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import same_scope as S

LABELS = ['complete_primary_exon_span', 'other_span_of_primary_eligible_coding_gene',
          'partial_flagged_coding_gene', 'other_coding_gene_without_complete_primary',
          'noncoding_or_other_annotated_gene', 'outside_all_annotated_genes']


def regions(length, annotation, coding, primary):
    labels = np.full(length, 5, dtype=np.uint8)
    eligible = {t['gene_id'] for t in primary}
    # Assign in reverse priority so overlaps are accounted for exactly once.
    for g in annotation['genes'].values():
        labels[g['start']:g['end']] = 4
    for priority in (3, 2, 1):
        for gene_id, g in coding['genes'].items():
            category = 1 if gene_id in eligible else 2 if g['partial'] else 3
            if category == priority:
                labels[g['start']:g['end']] = category
    for t in primary:
        parts = t['exon'] or t['CDS']
        labels[parts[0][0]:parts[-1][1]] = 0
    return labels


def burden(txs, labels):
    predicted = np.zeros(len(labels), dtype=bool)
    for t in txs:
        predicted[t['CDS'][0][0]:t['CDS'][-1][1]] = True
    denom = np.bincount(labels, minlength=6)
    counts = np.bincount(labels[predicted], minlength=6)
    return {'genome_bp_by_region': dict(zip(LABELS, map(int, denom))),
            'predicted_union_span_bp_by_region': dict(zip(LABELS, map(int, counts))),
            'frozen_background_FPR': float(counts[1:].sum()/denom[1:].sum()),
            'all_coding_gene_background_FPR': float(counts[4:].sum()/denom[4:].sum()),
            'all_annotated_gene_background_FPR': float(counts[5]/denom[5]),
            'old_FP_burden_inside_partial_coding_fraction': float(counts[2]/counts[1:].sum()) if counts[1:].sum() else None,
            'old_FP_burden_inside_any_annotated_gene_fraction': float(counts[1:5].sum()/counts[1:].sum()) if counts[1:].sum() else None}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    previous = json.loads((S.ROOT/'outputs/M26-SAME-SCOPE-MECHANISM-R1/summary.json').read_text())
    result = {'status': 'posthoc_reference_policy_attribution_not_biological_validation',
              'regions_in_overlap_priority_order': LABELS, 'methods': {}, 'setaria_access': False}
    rows = []
    for species, seqid, length, expected in S.SCOPE:
        reference = S.ROOT/'data/m1_screen'/species/'reference.gff3'
        annotation = S.E.parse_annotation(reference, {seqid: length})
        coding = S.E.parse_annotation(reference, {seqid: length}, protein_coding_only=True)
        primary = S.E.primary_transcripts(coding)
        assert len(primary) == expected
        labels = regions(length, annotation, coding, primary)
        for method, template in S.METHODS.items():
            path = S.ROOT/template.format(species=species)
            txs = S.E.primary_transcripts(S.E.parse_annotation(path, {seqid: length}), complete_only=False)
            stats = burden(txs, labels)
            old = previous['methods'][method]['per_species'][species]
            assert abs(stats['frozen_background_FPR']-old['frozen_metrics']['intergenic_FPR']) < 1e-12
            assert abs(stats['all_annotated_gene_background_FPR']-old['all_annotated_gene_background_FPR']) < 1e-12
            result['methods'].setdefault(method, {})[species] = stats
            for region in LABELS:
                rows.append({'method': method, 'species': species, 'region': region,
                             'genome_bp': stats['genome_bp_by_region'][region],
                             'predicted_union_span_bp': stats['predicted_union_span_bp_by_region'][region]})
    with (out/'background_attribution.json').open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    with (out/'background_attribution.tsv').open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
