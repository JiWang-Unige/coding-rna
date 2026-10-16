#!/usr/bin/env python3
"""Reference-only qualification of one TAIR release; never reads predictions."""
import argparse
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SEQID, LENGTH = 'NC_003074.8', 23459830
KEEP = {'ID', 'Parent', 'Name', 'locus_tag', 'locus_type', 'partial',
        'start_range', 'end_range', 'exception', 'transl_except', 'pseudo', 'pseudogene'}


def structural_attrs(text):
    # Free-text descriptions in this release contain unescaped semicolons.
    # Only explicitly delimited structural keys are used; raw input is retained.
    found = defaultdict(list)
    for item in text.rstrip(';').split(';'):
        if '=' not in item:
            continue
        key, value = item.strip().split('=', 1)
        if key in KEEP:
            found[key].append(value)
    if any(len(set(values)) != 1 for values in found.values()):
        raise ValueError('Conflicting repeated structural attributes')
    return {key: values[0] for key, values in found.items()}


def type_table(lines):
    header = next(lines).rstrip('\n').split('\t')
    if header[:2] != ['name', 'gene_model_type']:
        raise ValueError('Unexpected functional-table header')
    mapping, model_types = defaultdict(set), {}
    subfeatures, unknown = [], []
    for number, line in enumerate(lines, 2):
        cells = line.rstrip('\n').split('\t', 2)
        if len(cells) < 2:
            raise ValueError(f'Type table malformed row {number}')
        model, kind = cells[:2]
        if not model.startswith('AT3G'):
            continue
        if kind == 'uORF' and re.fullmatch(r'AT3G\d+\.uORF\d+(?:-\d+)?', model):
            host = model.split('.')[0]
            derives = re.search(r'Derives_from\s+(AT3G\d+)(?:[;\s]|$)', cells[2] if len(cells) > 2 else '')
            if derives and derives.group(1) == host:
                subfeatures.append({'model': model, 'type': kind, 'host': host, 'row': number})
            else:
                unknown.append({'model': model, 'type': kind, 'row': number, 'reason': 'uORF_host_not_confirmed'})
            continue
        if not re.fullmatch(r'AT3G\d+(?:\.\d+)?', model) or not kind:
            unknown.append({'model': model, 'type': kind, 'row': number, 'reason': 'unclassified_table_object'})
            continue
        if model in model_types and model_types[model] != kind:
            raise ValueError(f'Conflicting model type: {model}')
        model_types[model] = kind
        mapping[model.split('.')[0]].add(kind)
    return mapping, model_types, subfeatures, unknown


def decide_type(gene, attrs, mapping):
    stated = attrs.get('locus_type')
    other = mapping.get(gene, set())
    if len(other) > 1:
        return None, 'functional_table_conflict'
    if stated and other and other != {stated}:
        return None, 'GFF_table_disagreement'
    if stated:
        return stated, 'GFF_locus_type'
    if len(other) == 1:
        return next(iter(other)), 'same_release_functional_table'
    return None, 'unresolved'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', required=True)
    out = Path(parser.parse_args().output_dir)
    out.mkdir(parents=True, exist_ok=True)
    data = ROOT/'data/m26_reference_pair_20260926'
    with gzip.open(data/'Araport11_functional_descriptions_20240630.txt.gz', 'rt', encoding='latin-1') as handle:
        mapping, model_types, subfeatures, unknown = type_table(handle)
    records, genes, transcripts, cds = [], {}, {}, []
    issues = [['unknown_table_object', row['model']] for row in unknown]
    feature_counts, resolution_counts, missing_types = Counter(), Counter(), Counter()
    ledger = []
    with gzip.open(data/'Araport11.20240701.gff.gz', 'rt', encoding='latin-1') as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if fields[0] != 'Chr3':
                continue
            if len(fields) != 9:
                raise ValueError(f'GFF malformed row {number}')
            feature, attrs = fields[2], structural_attrs(fields[8])
            start, end = int(fields[3])-1, int(fields[4])
            if not (0 <= start < end <= LENGTH) or fields[6] not in {'+', '-'}:
                issues.append(['coordinates', number])
            feature_counts[feature] += 1
            if feature == 'gene':
                gene = attrs['ID']
                if gene in genes:
                    issues.append(['duplicate_gene', gene])
                kind, method = decide_type(gene, attrs, mapping)
                resolution_counts[method] += 1
                genes[gene] = (fields[6], start, end, kind)
                ledger.append((gene, attrs.get('locus_type', ''), ','.join(sorted(mapping.get(gene, set()))), kind or '', method))
                if not attrs.get('locus_type'):
                    missing_types[kind or 'UNRESOLVED'] += 1
                if kind is None:
                    issues.append([method, gene])
                else:
                    attrs['gene_biotype'] = kind
                    attrs['locus_type'] = kind
            elif feature == 'mRNA':
                tx, parent = attrs['ID'], attrs.get('Parent', '')
                if tx in transcripts or not parent or ',' in parent:
                    issues.append(['transcript_identity', tx])
                transcripts[tx] = (parent, fields[6], start, end)
            elif feature == 'CDS':
                cds.append((attrs.get('Parent', ''), fields[6], start, end, fields[7], number))
            fields[0] = SEQID
            fields[8] = ';'.join(f'{k}={v}' for k, v in attrs.items())
            records.append(fields)
    for tx, (gene, strand, start, end) in transcripts.items():
        if gene not in genes:
            issues.append(['transcript_missing_gene', tx])
        elif genes[gene][0] != strand or not genes[gene][1] <= start < end <= genes[gene][2]:
            issues.append(['transcript_outside_parent', tx])
    for parent, strand, start, end, phase, number in cds:
        if parent not in transcripts:
            issues.append(['CDS_missing_single_mRNA_parent', number])
        elif transcripts[parent][1] != strand or not transcripts[parent][2] <= start < end <= transcripts[parent][3]:
            issues.append(['CDS_outside_parent', number])
        if phase not in {'0', '1', '2'}:
            issues.append(['CDS_phase', number])
    result = {'stage': 'reference_only_no_prediction_access',
              'source_record': 'https://zenodo.org/records/15889110',
              'seqid': SEQID, 'gene_count': len(genes), 'feature_rows': dict(feature_counts),
              'functional_table_chr3_models': len(model_types),
              'table_subfeatures_not_gene_types': subfeatures,
              'table_unknown_objects': unknown,
              'functional_table_chr3_loci': len(mapping),
              'type_resolution': dict(resolution_counts),
              'previously_missing_type_resolution': dict(missing_types),
              'resolved_gene_types': dict(Counter(g[3] for g in genes.values())),
              'issues': issues, 'qualified': not issues,
              'normalization': 'all Chr3 rows, unchanged coordinates/strand/phase/feature; structural attributes only; gene type from explicit GFF or unambiguous same-release table'}
    with (out/'gene_type_resolution.tsv').open('x') as handle:
        handle.write('gene\tGFF_locus_type\ttable_type\tresolved_type\tresolution\n')
        for row in ledger:
            handle.write('\t'.join(row)+'\n')
    with (out/'input_qualification.json').open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    if not issues:
        with (out/'araport_chr3.qualified.gff3').open('x', encoding='utf-8') as handle:
            handle.write('##gff-version 3\n')
            for fields in records:
                handle.write('\t'.join(fields)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
