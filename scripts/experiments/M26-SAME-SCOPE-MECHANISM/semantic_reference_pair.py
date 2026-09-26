#!/usr/bin/env python3
"""Fixed-prediction source pairing after explicit gene-like object qualification."""
import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import reference_pair as R
import qualify_araport as Q
P, S = R.P, R.S
SEQID, LENGTH, SPECIES = R.SEQID, R.LENGTH, R.SPECIES
source_views, background_mask, fasta_target = R.source_views, R.background_mask, R.fasta_target


def classify_type(gene, feature, attrs, mapping):
    kind, evidence = Q.decide_type(gene, attrs, mapping)
    if kind is None and attrs.get('locus_type') == 'mirna' and mapping.get(gene) == {'miRNA_primary_transcript'}:
        return 'mirna', 'explicit_noncoding_vocabulary_disagreement'
    if kind is None and feature == 'pseudogene' and attrs.get('locus_type') == 'pseudogene' and mapping.get(gene) == {'pre_trna'}:
        return 'pseudogene', 'explicit_noncoding_pseudogene_tRNA_disagreement'
    if kind == 'protein_coding' and feature != 'gene':
        return None, 'coding_gene_like_feature_conflict'
    return kind, evidence


def qualify_gene_like(data, out):
    with gzip.open(data/'Araport11_functional_descriptions_20240630.txt.gz', 'rt', encoding='latin-1') as handle:
        mapping, _, subfeatures, unknown = Q.type_table(handle)
    old = json.loads((S.ROOT/'outputs/M26-ARAPORT-INPUT-R6-TYPEVIEW/input_qualification.json').read_text())
    old_parent_ids = {item[1] for item in old['issues'] if item[0] == 'transcript_missing_gene'}
    rows, genes, txs, cds = [], {}, {}, []
    issues = [['unknown_table_object', x['model']] for x in unknown]
    types, resolutions = Counter(), Counter()
    with gzip.open(data/'Araport11.20240701.gff.gz', 'rt', encoding='latin-1') as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if f[0] != 'Chr3':
                continue
            if len(f) != 9:
                raise ValueError(f'GFF columns at {number}')
            attrs, feature = Q.structural_attrs(f[8]), f[2]
            start, end = int(f[3])-1, int(f[4])
            if not (0 <= start < end <= LENGTH) or f[6] not in {'+', '-'}:
                issues.append(['coordinates', number])
            if feature in {'gene', 'pseudogene', 'transposable_element_gene'}:
                gene = attrs['ID']
                if gene in genes:
                    issues.append(['duplicate_gene_like_ID', gene])
                kind, evidence = classify_type(gene, feature, attrs, mapping)
                genes[gene] = {'feature': feature, 'kind': kind, 'strand': f[6], 'start': start, 'end': end}
                types[str(kind)] += 1
                resolutions[evidence] += 1
                if kind is None:
                    issues.append([evidence, gene])
                else:
                    attrs['gene_biotype'] = kind
                    if 'locus_type' not in attrs:
                        attrs['locus_type'] = kind
            elif feature == 'mRNA':
                tx = attrs['ID']
                if tx in txs or not attrs.get('Parent') or ',' in attrs['Parent']:
                    issues.append(['mRNA_identity', tx])
                txs[tx] = {'gene': attrs.get('Parent', ''), 'strand': f[6], 'start': start, 'end': end}
            elif feature == 'CDS':
                cds.append((attrs.get('Parent', ''), f[6], start, end, f[7], number))
            f[0], f[8] = SEQID, ';'.join(f'{k}={v}' for k, v in attrs.items())
            rows.append(f)
    parent_kinds, old_resolution = Counter(), Counter()
    for tx, t in txs.items():
        g = genes.get(t['gene'])
        if g is None:
            issues.append(['mRNA_parent_missing', tx])
            continue
        parent_kinds[g['feature']] += 1
        if tx in old_parent_ids:
            old_resolution[g['feature']] += 1
        if t['strand'] != g['strand'] or not g['start'] <= t['start'] < t['end'] <= g['end']:
            issues.append(['mRNA_parent_coordinates', tx])
    for parent, strand, start, end, phase, number in cds:
        t = txs.get(parent)
        if t is None:
            issues.append(['CDS_missing_single_mRNA', number])
        elif strand != t['strand'] or not t['start'] <= start < end <= t['end']:
            issues.append(['CDS_parent_coordinates', number])
        if phase not in {'0', '1', '2'}:
            issues.append(['CDS_phase', number])
    if sum(old_resolution.values()) != len(old_parent_ids):
        issues.append(['R6_missing_parent_not_all_resolved', len(old_parent_ids)-sum(old_resolution.values())])
    result = {'scope': 'protein_coding_pair_with_other_gene_like_objects_retained',
              'gene_like_features': dict(Counter(g['feature'] for g in genes.values())),
              'gene_like_types': dict(types), 'type_resolution': dict(resolutions),
              'mRNA_parent_feature': dict(parent_kinds),
              'R6_reported_parent_issues_resolved_to': dict(old_resolution),
              'uORF_table_subfeatures': len(subfeatures),
              'CDS_rows': len(cds), 'all_Chr3_feature_rows_retained': len(rows),
              'issues': issues, 'protein_coding_type_uncertainty_loci': sum(g['kind'] is None for g in genes.values()),
              'coding_view_qualified': not issues}
    with (out/'coding_input_qualification.json').open('x') as handle:
        json.dump(result, handle, indent=2)
    with (out/'mRNA_parent_ledger.tsv').open('x') as handle:
        handle.write('mRNA\tgene_like_parent\tparent_feature\tparent_type\tR6_flagged\n')
        for tx, t in txs.items():
            g = genes.get(t['gene'], {})
            handle.write('\t'.join(map(str, [tx,t['gene'],g.get('feature'),g.get('kind'),tx in old_parent_ids]))+'\n')
    if issues:
        raise ValueError('Input not qualified; see coding_input_qualification.json')
    normalized = out/'araport_chr3.coding_qualified.gff3'
    with normalized.open('x') as handle:
        handle.write('##gff-version 3\n')
        for row in rows:
            handle.write('\t'.join(row)+'\n')
    return normalized, result


def locus_primary_class(left_primary, right_primary, left_any, right_any, left_raw=None, right_raw=None):
    if left_primary is None and right_primary is None:
        return 'neither_CDS_assessable'
    if left_primary is None:
        return 'Araport_only_CDS_assessable'
    if right_primary is None:
        return 'RefSeq_only_CDS_assessable'
    if left_primary == right_primary:
        return 'same_primary'
    if left_primary in right_any and right_primary in left_any:
        return 'primary_selection_only'
    if left_raw is not None and left_primary in right_raw and right_primary in left_raw:
        return 'CDS_assessability_difference'
    return 'primary_structure_disagreement'


def structural_partition(source_data, predictions, out):
    views = []
    for source, (path, selected, paired) in source_data.items():
        raw_g, _, _, _ = P.raw_features(open(path), SEQID)
        primary = {raw_g[t['gene_id']].get('locus_tag', t['gene_id']): S.key(t) for t in selected}
        any_keys, raw_keys = defaultdict(set), defaultdict(set)
        for (locus, key), flags in paired.items():
            raw_keys[locus].add(key)
            if flags['CDS_assessable']:
                any_keys[locus].add(key)
        loci = {a.get('locus_tag', g) for g, a in raw_g.items() if a.get('gene_biotype') == 'protein_coding'}
        views.append((source, primary, any_keys, loci, raw_keys))
    (_, lp, la, ll, lraw), (_, rp, ra, rl, rraw) = views
    groups, ledger = defaultdict(list), []
    pred_keys = {m: {S.key(t) for t in ts} for m, ts in predictions.items()}
    for locus in sorted(ll | rl):
        left, right = lp.get(locus), rp.get(locus)
        category = locus_primary_class(left, right, la[locus], ra[locus], lraw[locus], rraw[locus])
        groups[category].append(locus)
        row = {'locus':locus, 'category':category, 'RefSeq_coding_locus':locus in ll,
               'Araport_coding_locus':locus in rl,
               'RefSeq_primary_CDS':repr(left), 'Araport_primary_CDS':repr(right),
               'RefSeq_only_assessable_isoform_chains':len(la[locus]-ra[locus]),
               'Araport_only_assessable_isoform_chains':len(ra[locus]-la[locus]),
               'RefSeq_assessable_chains_absent_from_Araport_raw_CDS':len(la[locus]-rraw[locus]),
               'Araport_assessable_chains_absent_from_RefSeq_raw_CDS':len(ra[locus]-lraw[locus])}
        for method, keys in pred_keys.items():
            row[method+'_RefSeq_primary_exact'] = left in keys if left else False
            row[method+'_Araport_primary_exact'] = right in keys if right else False
            row[method+'_RefSeq_any_exact'] = bool(la[locus] & keys)
            row[method+'_Araport_any_exact'] = bool(ra[locus] & keys)
        ledger.append(row)
    with (out/'semantic_locus_partition.tsv').open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(ledger[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(ledger)
    result = {'coding_loci_RefSeq':len(ll), 'coding_loci_Araport':len(rl),
              'shared_coding_loci':len(ll & rl),
              'loci_with_any_assessable_isoform_set_difference':sum(la[x] != ra[x] for x in ll | rl),
              'RefSeq_only_assessable_isoform_chains':sum(len(la[x]-ra[x]) for x in ll | rl),
              'Araport_only_assessable_isoform_chains':sum(len(ra[x]-la[x]) for x in ll | rl),
              'shared_loci_with_genuinely_source_only_assessable_CDS':sum(bool((la[x]-rraw[x]) | (ra[x]-lraw[x])) for x in ll & rl),
              'residual_structure_method_flips':{},
              'groups':{}}
    residual = {x for x in ll & rl if (la[x]-rraw[x]) or (ra[x]-lraw[x])}
    for method, keys in pred_keys.items():
        result['residual_structure_method_flips'][method] = {
            'residual_loci':len(residual),
            'primary_exact_judgment_flips':sum((lp.get(x) in keys) != (rp.get(x) in keys) for x in residual),
            'any_assessable_exact_judgment_flips':sum(bool(la[x] & keys) != bool(ra[x] & keys) for x in residual)}
    for group, loci in groups.items():
        result['groups'][group] = {'loci':len(loci), 'methods':{}}
        for method, keys in pred_keys.items():
            result['groups'][group]['methods'][method] = {
                'RefSeq_primary_exact_loci':sum(lp.get(x) in keys for x in loci),
                'Araport_primary_exact_loci':sum(rp.get(x) in keys for x in loci),
                'RefSeq_any_assessable_exact_loci':sum(bool(la[x] & keys) for x in loci),
                'Araport_any_assessable_exact_loci':sum(bool(ra[x] & keys) for x in loci)}
    return result


def phase_chain_consistent(tx):
    parts = tx['CDS'] if tx['strand'] == '+' else list(reversed(tx['CDS']))
    if not parts:
        return True
    expected = int(parts[0][2])
    for start, end, phase in parts:
        if int(phase) != expected:
            return False
        expected = (expected-(end-start)) % 3
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', required=True)
    out = Path(parser.parse_args().output_dir)
    out.mkdir(parents=True, exist_ok=True)
    if (out/'reference_pair.json').exists():
        raise FileExistsError('Do not overwrite completed results')
    data = S.ROOT/'data/m26_reference_pair_20260926'
    sequence = S.V.read_fasta(S.ROOT/'data/m1_screen'/SPECIES/'genome.fa')[SEQID].upper()
    assert len(sequence) == LENGTH
    assert fasta_target(data/'TAIR10_chr_all.fas.gz', 'Chr3') == sequence, 'FASTA mismatch: stop coordinate comparison'
    normalized, qualification = qualify_gene_like(data, out)
    paths = {'RefSeq': S.ROOT/'data/m1_screen'/SPECIES/'reference.gff3', 'Araport11_20240701': normalized}
    previous = json.loads((S.ROOT/'outputs/M26-SAME-SCOPE-MECHANISM-R1/summary.json').read_text())
    predictions = {m: S.E.primary_transcripts(S.E.parse_annotation(S.ROOT/t.format(species=SPECIES), {SEQID:LENGTH}), complete_only=False) for m,t in S.METHODS.items()}
    result = {'status':'development_reference_representation_control_shared_annotation_ancestry',
              'species':SPECIES, 'seqid':SEQID, 'length':LENGTH, 'FASTA_exact_uppercase_match':True,
              'setaria_access':False, 'sources':{}}
    paired_sources = {}
    source_data = {}
    for source, path in paths.items():
        views, coding, paired, stats = source_views(path, sequence)
        phase_issues = [t['id'] for t in coding['transcripts'].values() if not phase_chain_consistent(t)]
        stats['phase_chain_inconsistent_transcripts'] = phase_issues
        if phase_issues:
            with (out/(source+'_phase_issues.json')).open('x') as handle:
                json.dump(phase_issues, handle)
            raise ValueError('Coding phase chain inconsistent; do not score silently')
        paired_sources[source] = paired
        source_data[source] = (path, views['CDS_assessable'], paired)
        if source == 'RefSeq':
            assert len(views['parent_filter']) == 4151
        all_coding = np.zeros(LENGTH, bool)
        for g in coding['genes'].values():
            all_coding[g['start']:g['end']] = True
        stats['fixed_all_coding_gene_background_bp'] = int((~all_coding).sum())
        stats['policies'] = {}
        for policy, refs in views.items():
            background = ~background_mask(refs)
            ref_keys = {S.key(t) for t in refs}
            per_method = {}
            for method, txs in predictions.items():
                keys = {S.key(t) for t in txs}
                span, cds = P.coverage_masks(txs, LENGTH)
                metrics = S.prf(len(keys & ref_keys), len(keys), len(ref_keys))
                counts = {'span_fp_bp':int((span & background).sum()), 'background_bp':int(background.sum()),
                          'CDS_coverage_on_background_bp':int((cds & background).sum())}
                fpr = counts['span_fp_bp']/counts['background_bp']
                if source == 'RefSeq' and policy == 'parent_filter':
                    old = previous['methods'][method]['per_species'][SPECIES]
                    assert metrics == old['exact_chain']
                    assert abs(fpr-old['frozen_metrics']['intergenic_FPR']) < 1e-12
                per_method[method] = {'chain':metrics, 'policy_span_background_counts':counts,
                                      'policy_span_FPR':fpr,
                                      'all_coding_background_span_bp':int((span & ~all_coding).sum()),
                                      'all_coding_background_CDS_bp':int((cds & ~all_coding).sum())}
            stats['policies'][policy] = {'reference_unique_chains':len(ref_keys), 'methods':per_method}
        result['sources'][source] = stats
    result['input_qualification'] = qualification
    result['semantic_locus_partition'] = structural_partition(source_data, predictions, out)
    n, a = paired_sources.values()
    shared = n.keys() & a.keys()
    same_cds = {k for k in shared if n[k]['CDS_assessable'] and a[k]['CDS_assessable']}
    discordant = {k for k in same_cds if not n[k]['parent_eligible'] and a[k]['parent_eligible']}
    result['shared_chain_control'] = {
        'same_locus_and_CDS_chain_records':len(shared),
        'CDS_assessable_in_both':len(same_cds),
        'CDS_assessable_parent_eligibility_RefSeq_Araport':dict(Counter(f'{int(n[k]["parent_eligible"])},{int(a[k]["parent_eligible"])}' for k in same_cds)),
        'RefSeq_excluded_Araport_retained_shared_chains':len(discordant),
        'RefSeq_excluded_Araport_retained_shared_loci':len({k[0] for k in discordant}),
        'methods':{}}
    for method, txs in predictions.items():
        keys = {S.key(t) for t in txs}
        matched = {k for k in discordant if k[1] in keys}
        result['shared_chain_control']['methods'][method] = {'exact_discordant_chains':len(matched), 'exact_discordant_loci':len({k[0] for k in matched})}
    with (out/'shared_chain_eligibility.tsv').open('x') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['locus','strand','CDS_0based_halfopen','RefSeq_parent','Araport_parent','CDS_assessable_both'])
        for locus, key in sorted(shared):
            writer.writerow([locus,key[1],';'.join(f'{x}-{y}' for x,y in key[2]),n[(locus,key)]['parent_eligible'],a[(locus,key)]['parent_eligible'],(locus,key) in same_cds])
    with (out/'reference_pair.json').open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
