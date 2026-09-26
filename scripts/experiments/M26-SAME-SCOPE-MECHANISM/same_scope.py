#!/usr/bin/env python3
"""Finite, read-only development comparison of saved M25R and M12B predictions."""
import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT), str(ROOT / 'scripts')]
import eval_structure_diagnostic as E
import eval_m25_structure as V
from src.foundation_probe import train_generanno_structural_heads as M

SCOPE = [('arabidopsis_thaliana', 'NC_003074.8', 23459830, 4151),
         ('oryza_sativa', 'NC_089041.1', 29936421, 2299)]
RECOVERY = 'outputs/M25R-E1-RUN-GRAPH-R1-CPU-RECOVERY'
METHODS = {
    'M25R_A': RECOVERY + '/A_predictions.gff3',
    'M25R_B': RECOVERY + '/B_predictions.gff3',
    'ANNEVO': 'outputs/M12B-SAMEPANEL-BASELINES-ANNEVO/predictions/{species}.gff',
    'Helixer': 'outputs/M12B-SAMEPANEL-BASELINES-HELIXER/predictions/{species}.gff3',
    'Tiberius': 'outputs/M12B-SAMEPANEL-BASELINES-TIBERIUS/predictions/{species}.gtf',
}


def key(tx):
    return (tx['seqid'], tx['strand'], tuple((a, b) for a, b, _ in tx['CDS']))


def models(txs):
    return [{'strand': t['strand'], 'cds': [(a, b) for a, b, _ in t['CDS']]} for t in txs]


def prf(tp, predicted, reference):
    return {'tp': tp, 'predicted': predicted, 'reference': reference,
            'precision': tp / predicted if predicted else 0,
            'recall': tp / reference if reference else 0,
            'f1': 2 * tp / (predicted + reference) if predicted + reference else 0}


def union_stats(left, right):
    a, b = M._merge(left), M._merge(right)
    return sum(y-x for x, y in a), sum(y-x for x, y in b), M._overlap_length(a, b)


def compare(predictions, references, lengths, all_annotations):
    frozen = M._validation_metrics({k: models(v) for k, v in predictions.items()}, references, lengths)
    # The old helper assumes validity by construction; never report it for baselines.
    frozen.pop('structurally_valid_complete_fraction')
    p = {key(t) for rows in predictions.values() for t in rows}
    r = {key(t) for rows in references.values() for t in rows}
    cds_p = cds_r = cds_tp = span_fp = intergenic_all = 0
    for scope, txs in predictions.items():
        a, b, overlap = union_stats([(x, y) for t in txs for x, y, _ in t['CDS']],
                                    [(x, y) for t in references[scope] for x, y, _ in t['CDS']])
        cds_p += a
        cds_r += b
        cds_tp += overlap
        a, b, overlap = union_stats([(t['CDS'][0][0], t['CDS'][-1][1]) for t in txs],
                                    [(g['start'], g['end']) for g in all_annotations[scope]['genes'].values()])
        span_fp += a - overlap
        intergenic_all += lengths[scope] - b
    return {'frozen_metrics': frozen, 'exact_chain': prf(len(p & r), len(p), len(r)),
            'CDS_base': prf(cds_tp, cds_p, cds_r),
            'all_annotated_gene_background_FPR': span_fp / intergenic_all if intergenic_all else None,
            'all_annotated_background_counts': {'predicted_intergenic_span_bp': span_fp,
                                               'intergenic_bp': intergenic_all},
            'prediction_records': sum(map(len, predictions.values())),
            'duplicate_chain_records': sum(map(len, predictions.values())) - len(p)}


def error_rows(predictions, refs, all_coding):
    exact = {key(t) for t in refs}
    isoform = {key(t) for t in all_coding}
    rows = []
    for t in predictions:
        k = key(t)
        a, b = t['CDS'][0][0], t['CDS'][-1][1]
        overlaps = [r for r in refs if r['CDS'][0][0] < b and r['CDS'][-1][1] > a]
        same = [r for r in overlaps if r['strand'] == t['strand']]
        if k in exact:
            category = 'exact_primary'
        elif k in isoform:
            category = 'exact_other_complete_isoform'
        elif same and E.introns(t) and any(E.introns(t) == E.introns(r) for r in same):
            category = 'same_intron_chain_terminal_disagreement'
        elif same:
            category = 'same_strand_reference_overlap'
        elif overlaps:
            category = 'opposite_strand_reference_overlap_only'
        else:
            category = 'no_complete_primary_CDS_span_overlap'
        rows.append({'transcript': t['id'], 'seqid': t['seqid'], 'strand': t['strand'],
                     'category': category, 'cds_bp': sum(y-x for x, y, _ in t['CDS']),
                     'cds_count': len(t['CDS']), 'span_bp': b-a,
                     'same_strand_span_overlap_count': len(same),
                     'opposite_strand_span_overlap_count': len(overlaps)-len(same)})
    return rows


def sequence_convention(txs, seqs):
    counts = Counter()
    for t in txs:
        sequence = V.transcript_sequence(t, seqs)
        counts['transcripts'] += 1
        counts['terminal_stop_in_CDS'] += sequence[-3:] in {'TAA', 'TAG', 'TGA'}
        counts['ATG_start'] += sequence[:3] == 'ATG'
        counts['CDS_length_divisible_by_3'] += len(sequence) % 3 == 0
    return dict(counts)


def oracle_membership(refs, predictions, ledger):
    by_id = {(species, seqid, t['strand'], t['id']): t
             for (species, seqid), txs in refs.items() for t in txs}
    records = {tuple(row[f] for f in ('species', 'seqid', 'strand', 'transcript_id')): row for row in ledger}
    assert len(records) == len(ledger) == 6450 and set(records) == set(by_id)
    reachable = {key(by_id[k]) for k, row in records.items() if row['structurally_valid_chain_in_graph']}
    all_ref = {key(t) for t in by_id.values()}
    exact = {key(t) for txs in predictions.values() for t in txs} & all_ref
    return {'references': len(all_ref), 'reachable_references': len(reachable),
            'B_exact': len(exact), 'B_exact_inside_oracle': len(exact & reachable),
            'B_exact_outside_oracle': len(exact - reachable),
            'reachable_not_recovered': len(reachable - exact),
            'unreachable_references': len(all_ref - reachable),
            'categories': dict(Counter(r['category'] for r in ledger)),
            'B_exact_identity_check_pass': exact <= reachable,
            'oracle_utilization_if_consistent': len(exact & reachable) / len(reachable) if exact <= reachable else None,
            'perfect_filter_current_B_chain_F1': 2 * len(exact) / (len(exact) + len(all_ref)),
            'scope': 'reference-assisted conditional M25R graph diagnostic; not a baseline oracle'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'summary.json').exists():
        raise FileExistsError('Do not overwrite completed diagnostic')
    refs, lengths, annotations, all_coding, sequences = {}, {}, {}, {}, {}
    for species, seqid, length, count in SCOPE:
        scope = (species, seqid)
        fasta = ROOT / 'data/m1_screen' / species / 'genome.fa'
        seq = V.read_fasta(fasta)[seqid]
        assert len(seq) == length
        sequences[seqid] = seq
        lengths[scope] = length
        path = ROOT / 'data/m1_screen' / species / 'reference.gff3'
        annotations[scope] = E.parse_annotation(path, {seqid: length})
        coding = E.parse_annotation(path, {seqid: length}, protein_coding_only=True)
        refs[scope] = E.primary_transcripts(coding)
        assert len(refs[scope]) == count and len({key(t) for t in refs[scope]}) == count
        all_coding[scope] = [t for t in coding['transcripts'].values() if t['CDS'] and not t['partial']
                             and not coding['genes'][t['gene_id']]['partial']]
    summary = {'experiment': 'M26-SAME-SCOPE-MECHANISM-R1', 'scope': SCOPE,
               'status': 'exploratory_development_diagnostic', 'training': False, 'setaria_access': False,
               'reference_policy': 'complete protein-coding primary by longest CDS, then transcript ID',
               'prediction_policy': 'primary by longest CDS; retain partial predictions; no normalization or filtering',
               'frozen_FPR': 'predicted CDS span outside complete-primary reference exon span / complement bases',
               'methods': {}, 'reference_conventions': {s: sequence_convention(refs[(s,q)], sequences)
                                                       for s,q,_n,_c in SCOPE}}
    all_predictions = {}
    table = []
    for method, template in METHODS.items():
        preds, paths, per_species, conventions = {}, {}, {}, {}
        for species, seqid, length, _count in SCOPE:
            scope = (species, seqid)
            path = ROOT / template.format(species=species)
            annotation = E.parse_annotation(path, {seqid: length})
            txs = E.primary_transcripts(annotation, complete_only=False)
            if not txs:
                raise ValueError(f'{method}: missing prediction coverage for {seqid}; do not score as zero')
            preds[scope] = txs
            paths[species] = {'path': str(path.relative_to(ROOT)), 'bytes': path.stat().st_size,
                              'mtime_ns': path.stat().st_mtime_ns, 'retained_records': annotation['retained_records']}
            conventions[species] = sequence_convention(txs, sequences)
            stats = compare({scope: txs}, {scope: refs[scope]}, {scope: length}, {scope: annotations[scope]})
            rows = error_rows(txs, refs[scope], all_coding[scope])
            stats['error_categories'] = dict(Counter(r['category'] for r in rows))
            per_species[species] = stats
            with (out / f'{method}_{species}_chains.tsv').open('x') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
                writer.writeheader()
                writer.writerows(rows)
        result = compare(preds, refs, lengths, annotations)
        result.update(paths=paths, per_species=per_species, sequence_conventions=conventions)
        summary['methods'][method] = result
        all_predictions[method] = preds
        for label, stats in [('pooled', result), *per_species.items()]:
            table.append({'method': method, 'scope': label, **stats['frozen_metrics'],
                          'CDS_base_F1': stats['CDS_base']['f1'],
                          'exact_chains': stats['exact_chain']['tp'],
                          'predicted_chains': stats['exact_chain']['predicted'],
                          'reference_chains': stats['exact_chain']['reference'],
                          'all_annotated_background_FPR': stats['all_annotated_gene_background_FPR']})
        print(method, json.dumps(result['frozen_metrics']), flush=True)
    old = json.loads((ROOT / RECOVERY / 'B_evaluation.json').read_text())['metrics']
    reproduced = summary['methods']['M25R_B']['frozen_metrics']
    summary['B_frozen_reproduction'] = {name: abs(value-old[name]) for name, value in reproduced.items()}
    assert all(error < 1e-12 for error in summary['B_frozen_reproduction'].values())
    with (ROOT / 'outputs/M25R-E1-BLOCK-ASSEMBLY-REACHABILITY/references.jsonl').open() as handle:
        ledger = [json.loads(line) for line in handle]
    summary['oracle_identity'] = oracle_membership(refs, all_predictions['M25R_B'], ledger)
    with (out / 'summary.json').open('x') as handle:
        json.dump(summary, handle, indent=2)
        handle.write('\n')
    with (out / 'comparison.tsv').open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(table)
    print('M26_COMPLETED', flush=True)


if __name__ == '__main__':
    main()
