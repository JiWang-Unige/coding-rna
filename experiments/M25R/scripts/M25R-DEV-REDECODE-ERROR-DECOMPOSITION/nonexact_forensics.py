#!/usr/bin/env python3
"""Read-only annotation comparison of the nine frozen bypass releases."""
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'scripts'))
from eval_structure_diagnostic import parse_annotation, primary_transcripts

OUT = ROOT / 'outputs/M25R-E1-PHASE-NONEXACT-FORENSICS'
SOURCE = ROOT / 'outputs/M25R-E1-PHASE-GATE-BYPASS-R2'
SCOPE = {'arabidopsis_thaliana': 'NC_003074.8', 'oryza_sativa': 'NC_089041.1'}

def chain(tx):
    return tuple((a, b) for a, b, _ in tx['CDS'])

def overlap(a, b):
    return sum(max(0, min(y, v)-max(x, u)) for x, y in a for u, v in b)

def compare(pred, tx, strand):
    ref = chain(tx)
    introns = lambda c: {(c[i][1], c[i+1][0]) for i in range(len(c)-1)}
    pi, ri = introns(pred), introns(ref)
    start = 0 if strand == '+' else -1
    end = -1 if strand == '+' else 0
    pstart = pred[start][0 if strand == '+' else 1]
    rstart = ref[start][0 if strand == '+' else 1]
    pstop = pred[end][1 if strand == '+' else 0]
    rstop = ref[end][1 if strand == '+' else 0]
    ov = overlap(pred, ref)
    return {'reference_id': tx['id'], 'gene_id': tx['gene_id'], 'reference_CDS': ref,
            'reference_CDS_count': len(ref), 'overlap_bp': ov,
            'predicted_only_CDS_bp': sum(b-a for a,b in pred)-ov,
            'reference_only_CDS_bp': sum(b-a for a,b in ref)-ov,
            'same_intron_chain': pi == ri, 'shared_introns': len(pi & ri),
            'predicted_only_introns': sorted(pi-ri), 'reference_only_introns': sorted(ri-pi),
            'start_boundary_delta_genomic_bp': pstart-rstart,
            'stop_boundary_delta_genomic_bp': pstop-rstop,
            'exact_CDS_chain': pred == ref}

def main():
    OUT.mkdir(exist_ok=False)
    annotations, primaries, lengths = {}, {}, {}
    for species, seqid in SCOPE.items():
        path = ROOT / 'data/m1_screen' / species / 'reference.gff3'
        with path.open() as f:
            for line in f:
                if line.startswith('##sequence-region '):
                    _, name, start, stop = line.split()
                    if name == seqid:
                        assert int(start) == 1
                        lengths[seqid] = int(stop)
                        break
        assert seqid in lengths
        annotation = parse_annotation(str(path), {seqid: lengths[seqid]}, protein_coding_only=True)
        annotations[species] = annotation
        primaries[species] = primary_transcripts(annotation)
    assert sum(map(len, primaries.values())) == 6450
    panel = [json.loads(line) for line in (SOURCE/'paired_results.jsonl').open()]
    results, positive_checks = [], 0
    gff = parse_annotation(str(SOURCE/'B_panel_predictions.gff3'), lengths)
    gff_keys = {(tx['seqid'], tx['strand'], chain(tx)) for tx in gff['transcripts'].values()}
    for row in panel:
        if row['B']['terminal'] != 'emitted':
            continue
        species, strand, seqid = row['species'], row['strand'], row['seqid']
        pred = row['B']['model']['cds']
        n = lengths[seqid]
        pred = tuple(sorted((a,b) if strand == '+' else (n-b,n-a) for a,b in pred))
        assert (seqid, strand, pred) in gff_keys
        refs = [tx for tx in primaries[species] if tx['strand'] == strand]
        exact = [tx['id'] for tx in refs if chain(tx) == pred]
        if row['panel_group'] == 'exact_phase_failure':
            assert row['exact_reference_id'] in exact
            positive_checks += 1
        if row['panel_group'] != 'nonexact_phase_failure':
            continue
        assert not exact
        comparisons = [compare(pred, tx, strand) for tx in refs if overlap(pred, chain(tx)) > 0]
        comparisons.sort(key=lambda x: (-x['overlap_bp'], x['reference_id']))
        all_exact = [tx['id'] for tx in annotations[species]['transcripts'].values()
                     if tx['strand'] == strand and chain(tx) == pred]
        top = comparisons[0] if comparisons else None
        category = ('annotated_nonprimary_exact' if all_exact else
                    'no_same_strand_primary_CDS_overlap' if not top else
                    'same_introns_terminal_boundary_difference' if top['same_intron_chain'] else
                    'intron_chain_difference')
        results.append({'species': species, 'seqid': seqid, 'strand': strand,
                        'lineage_id': row['lineage_id'], 'predicted_CDS': pred,
                        'category': category, 'all_coding_annotation_exact_ids': all_exact,
                        'overlapping_primary_count': len(comparisons),
                        'primary_comparisons': comparisons,
                        'boundary_scores': row['B']['boundary_scores'],
                        'phase_checks': row['B']['phase_checks']})
    assert positive_checks == 64 and len(results) == 9
    with (OUT/'cases.json').open('x') as f:
        json.dump(results, f, indent=2)
    summary = {'n': 9, 'positive_coordinate_reference_checks': positive_checks,
               'categories': dict(Counter(r['category'] for r in results)),
               'primary_reference_total': 6450, 'new_inference': False,
               'setaria_access': False, 'coordinate_system': '0-based half-open genomic',
               'matching_rule': 'same-strand maximum CDS overlap; ties retained and ordered by transcript ID',
               'biological_truth_adjudicated': False}
    with (OUT/'summary.json').open('x') as f:
        json.dump(summary, f, indent=2)
    (OUT/'STATUS').write_text('COMPLETED\n')
    print(json.dumps(summary), flush=True)

if __name__ == '__main__':
    main()
