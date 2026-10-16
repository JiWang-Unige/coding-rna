#!/usr/bin/env python3
"""Feature-level partialness and CDS/span coverage; never rewrites frozen metrics."""
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import same_scope as S
import background_attribution as B

STOPS = {'TAA', 'TAG', 'TGA'}


def raw_features(lines, seqid):
    genes, transcripts, cds = {}, {}, defaultdict(list)
    counts = defaultdict(Counter)
    for line in lines:
        if line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) != 9 or fields[0] != seqid:
            continue
        feature, attrs = fields[2], S.E.parse_attrs(fields[8])
        if feature not in {'gene', 'mRNA', 'CDS'}:
            continue
        explicit = attrs.get('partial', '').lower() == 'true'
        ranged = 'start_range' in attrs or 'end_range' in attrs
        kind = 'partial_true_and_range' if explicit and ranged else 'partial_true_only' if explicit else 'range_only' if ranged else 'unflagged'
        counts[feature][kind] += 1
        if feature == 'gene':
            genes[attrs['ID']] = attrs
        elif feature == 'mRNA':
            transcripts[attrs['ID']] = attrs
        else:
            for parent in attrs['Parent'].split(','):
                cds[parent].append(attrs)
    return genes, transcripts, cds, {k: dict(v) for k, v in counts.items()}


def compatible(tx, sequence):
    parts = tx['CDS']
    ordered = parts if tx['strand'] == '+' else list(reversed(parts))
    return bool(parts and ordered[0][2] == '0' and len(sequence) % 3 == 0
                and set(sequence) <= set('ACGT') and sequence[:3] == 'ATG'
                and sequence[-3:] in STOPS
                and not any(sequence[i:i+3] in STOPS for i in range(0, len(sequence)-3, 3))
                and all(a[1] <= b[0] for a, b in zip(parts, parts[1:])))


def coverage_masks(txs, length):
    span, cds = np.zeros(length, bool), np.zeros(length, bool)
    for t in txs:
        span[t['CDS'][0][0]:t['CDS'][-1][1]] = True
        for a, b, _ in t['CDS']:
            cds[a:b] = True
    assert not np.any(cds & ~span)
    return span, cds


def split_coverage(txs, labels):
    span, cds = coverage_masks(txs, len(labels))
    counts = {}
    for label, mask in [('span', span), ('CDS', cds), ('span_only', span & ~cds)]:
        counts[label] = np.bincount(labels[mask], minlength=6).astype(int).tolist()
    assert counts['span'] == (np.array(counts['CDS']) + counts['span_only']).tolist()
    return counts, span, cds


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', required=True)
    out = Path(parser.parse_args().output_dir)
    out.mkdir(parents=True, exist_ok=True)
    if (out/'partial_semantics.json').exists():
        raise FileExistsError('Do not overwrite completed results')
    old = json.loads((S.ROOT/'outputs/M26-REFERENCE-BACKGROUND-R2/background_attribution.json').read_text())
    result = {'status': 'exploratory_feature_semantics_not_independent_replication',
              'setaria_access': False, 'frozen_metrics_unchanged': True,
              'regions': B.LABELS, 'species': {}}
    ledger = []
    for species, seqid, length, expected in S.SCOPE:
        path = S.ROOT/'data/m1_screen'/species/'reference.gff3'
        annotation = S.E.parse_annotation(path, {seqid: length})
        coding = S.E.parse_annotation(path, {seqid: length}, protein_coding_only=True)
        primary = S.E.primary_transcripts(coding)
        assert len(primary) == expected
        primary_genes = {t['gene_id'] for t in primary}
        raw_g, raw_t, raw_c, tag_counts = raw_features(open(path), seqid)
        seqs = {seqid: S.V.read_fasta(S.ROOT/'data/m1_screen'/species/'genome.fa')[seqid].upper()}
        assert len(seqs[seqid]) == length
        combos, partial_combos = Counter(), Counter()
        eligible = defaultdict(list)
        partial_txs = []
        for tx in coding['transcripts'].values():
            if not tx['CDS']:
                continue
            g = coding['genes'][tx['gene_id']]
            gp = S.E.is_partial(raw_g[tx['gene_id']])
            tp = S.E.is_partial(raw_t.get(tx['id'], {}))
            cp = any(S.E.is_partial(a) for a in raw_c[tx['id']])
            assert gp == g['partial']
            sequence = S.V.transcript_sequence(tx, seqs).upper()
            ok = compatible(tx, sequence)
            exception = any('exception' in a or 'transl_except' in a for a in raw_c[tx['id']])
            combo = f'gene={int(gp)};transcript={int(tp)};CDS={int(cp)};sequence_compatible={int(ok)}'
            combos[combo] += 1
            if gp:
                partial_combos[combo] += 1
                partial_txs.append(tx)
            if not cp and ok:
                eligible[tx['gene_id']].append(tx)
            ledger.append({'species': species, 'seqid': seqid, 'gene_id': tx['gene_id'],
                           'transcript_id': tx['id'], 'gene_partial': gp, 'transcript_partial': tp,
                           'CDS_partial_any_row': cp, 'CDS_sequence_compatible': ok,
                           'CDS_exception_present': exception,
                           'old_primary_eligible_gene': tx['gene_id'] in primary_genes,
                           'CDS_bp': len(sequence)})
        extra = {g: min(ts, key=lambda t: (-sum(b-a for a,b,_ in t['CDS']), t['id']))
                 for g, ts in eligible.items() if g not in primary_genes}
        pg = {g for g, data in coding['genes'].items() if data['partial']}
        partial_compatible = [t for g in pg for t in eligible.get(g, [])]
        labels = B.regions(length, annotation, coding, primary)
        _, ref_partial_cds = coverage_masks(partial_txs, length)
        partial_ok_span, partial_ok_cds = coverage_masks(partial_compatible, length)
        region = labels == 2
        stats = {'raw_feature_row_flags': tag_counts, 'coding_transcript_combinations': dict(combos),
                 'partial_gene_transcript_combinations': dict(partial_combos),
                 'coding_genes': len(coding['genes']), 'frozen_primary_genes': expected,
                 'partial_coding_genes': len(pg),
                 'partial_genes_with_unflagged_compatible_CDS': len(pg & eligible.keys()),
                 'additional_CDS_assessable_primary_candidates': len(extra),
                 'additional_candidates_from_partial_genes': len(pg & extra.keys()),
                 'partial_region_bp': int(region.sum()),
                 'partial_region_reference_CDS_union_bp': int((region & ref_partial_cds).sum()),
                 'partial_region_compatible_CDS_union_bp': int((region & partial_ok_cds).sum()),
                 'partial_region_compatible_CDS_span_union_bp': int((region & partial_ok_span).sum()),
                 'methods': {}}
        extra_keys = {S.key(t) for t in extra.values()}
        extra_any = {S.key(t) for g, ts in eligible.items() if g not in primary_genes for t in ts}
        for method, template in S.METHODS.items():
            txs = S.E.primary_transcripts(S.E.parse_annotation(S.ROOT/template.format(species=species), {seqid: length}), complete_only=False)
            split, span, cds = split_coverage(txs, labels)
            expected_span = list(old['methods'][method][species]['predicted_union_span_bp_by_region'].values())
            assert split['span'] == expected_span
            keys = {S.key(t) for t in txs}
            stats['methods'][method] = {
                'coverage_by_region_bp': split,
                'partial_span_on_reference_CDS_bp': int((span & region & ref_partial_cds).sum()),
                'partial_span_on_compatible_CDS_span_bp': int((span & region & partial_ok_span).sum()),
                'partial_predicted_CDS_on_reference_CDS_bp': int((cds & region & ref_partial_cds).sum()),
                'partial_predicted_CDS_on_compatible_reference_CDS_bp': int((cds & region & partial_ok_cds).sum()),
                'exact_additional_primary_candidate_chains': len(keys & extra_keys),
                'exact_additional_any_compatible_isoform_chains': len(keys & extra_any),
                'excluded_genes_with_exact_compatible_isoform': sum(any(S.key(t) in keys for t in ts) for g, ts in eligible.items() if g not in primary_genes)}
        result['species'][species] = stats
    with (out/'partial_semantics.json').open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    with (out/'transcript_semantics.tsv').open('x') as f:
        writer = csv.DictWriter(f, fieldnames=list(ledger[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(ledger)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
