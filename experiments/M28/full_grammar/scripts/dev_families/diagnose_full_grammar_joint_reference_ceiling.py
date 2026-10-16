#!/usr/bin/env python3
"""Reference-informed structural oracle on saved original TRAIN records only."""
from collections import Counter, defaultdict
import itertools
import json
import os
from pathlib import Path
import time

import diagnose_native_train_reference as D

ROOT = D.ROOT
SOURCE = ROOT / 'outputs/M28-FULL-GRAMMAR-NATIVE-GAP-SUPERVISION'
OUT = ROOT / 'outputs/M28-FULL-GRAMMAR-JOINT-REFERENCE-CEILING'
PLAN = 'reports/M28-METHOD-RESTART/full_grammar_train_energy_diagnostic.md'


def span(row):
    return row['CDS'][0][0], row['CDS'][-1][1]


def maximum_compatible(rows):
    """Unit-weight interval scheduling; opposite strands have separate budgets."""
    groups = defaultdict(list)
    for row in rows:
        groups[row['species'], row['seqid'], row['strand']].append(row)
    selected = set()
    for part in groups.values():
        last_end = -1
        for row in sorted(part, key=lambda r: (span(r)[1], span(r)[0], D.record_id(r))):
            start, end = span(row)
            if start >= last_end:
                selected.add(D.record_id(row))
                last_end = end
    return selected


def brute_count(rows):
    best = 0
    for mask in itertools.product((False, True), repeat=len(rows)):
        chosen = [r for r, take in zip(rows, mask) if take]
        if all((D.record_id(a)[:3] != D.record_id(b)[:3]
                or span(a)[1] <= span(b)[0] or span(b)[1] <= span(a)[0])
               for a, b in itertools.combinations(chosen, 2)):
            best = max(best, len(chosen))
    return best


def fixture():
    def row(name, start, end, strand='+', seqid='q', supported=True):
        return dict(species='s', seqid=seqid, strand=strand, id=name,
                    CDS=[(start, end, 0)], single_chain_supported=supported)
    cases = [[], [row('one', 0, 3)],
        [row('a', 0, 3), row('b', 3, 6)],
        [row('long', 0, 12), row('a', 1, 3), row('b', 3, 6), row('c', 7, 11)],
        [row('same', 0, 6), row('same', 0, 6, '-'), row('same', 0, 6, seqid='r')],
        [row('x', 0, 4), row('y', 2, 6), row('z', 4, 8)],
        [row('a', 0, 3), row('b', 0, 3), row('c', 3, 6)],
        [row('unsupported', 0, 2, supported=False), row('a', 1, 3), row('b', 3, 6)]]
    for rows in cases:
        assert len(maximum_compatible(rows)) == brute_count(rows)
        allowed = [r for r in rows if r['single_chain_supported']]
        assert len(maximum_compatible(allowed)) == brute_count(allowed)
        assert maximum_compatible(rows) == maximum_compatible(list(reversed(rows)))
    assert len(maximum_compatible(cases[2])) == 2
    assert len(maximum_compatible(cases[3])) == 3
    assert len(maximum_compatible(cases[4])) == 3
    return dict(cases=len(cases), exhaustive_small_set_comparison=True,
                half_open_touching_nested_opposite_strand_distinct_contig_and_ties=True)


def describe(rows):
    n = len(rows)
    allowed = [r for r in rows if r['single_chain_supported']]
    geometric = len(maximum_compatible(rows))
    structural = len(maximum_compatible(allowed))
    projections = {}
    for method in (D.MODEL,) + D.METHODS:
        hits = [r for r in rows if (r['exact_primary_hit'] if method == D.MODEL
                                   else r['native_exact_hits'][method])]
        supported_hits = [r for r in hits if r['single_chain_supported']]
        projections[method] = dict(exact_reference_hits=len(hits),
            individually_supported_hits=len(supported_hits),
            maximum_span_compatible_hits=len(maximum_compatible(hits)),
            maximum_supported_span_compatible_hits=len(maximum_compatible(supported_hits)))
    anchor = projections[D.MODEL]
    if not (anchor['exact_reference_hits'] == anchor['individually_supported_hits']
            == anchor['maximum_span_compatible_hits']
            == anchor['maximum_supported_span_compatible_hits'] <= structural):
        raise ValueError('Original caller hits violate saved grammar or same-strand span compatibility')
    return dict(reference_chains=n, individually_supported=len(allowed),
        geometry_only_maximum_hits=geometric, supported_geometry_maximum_hits=structural,
        geometry_only_forced_omissions=n-geometric,
        single_chain_unsupported=n-len(allowed),
        extra_forced_omissions_among_supported=len(allowed)-structural,
        total_structurally_forced_omissions=n-structural,
        additional_loss_of_maximum_after_grammar_filter=geometric-structural,
        oracle_precision=1.0 if structural else None,
        structural_recall_ceiling=structural/n if n else None,
        structural_F1_ceiling=2*structural/(n+structural) if n else None,
        completed_hit_set_projections=projections)


def main():
    began = time.perf_counter()
    checked = fixture()  # Algorithm failures stop before reading research records.
    source, paired = D.completed(SOURCE), D.completed(D.OUT)
    if ((SOURCE/'JOBID').read_text().strip() != '13462169'
            or source['job_id'] != '13462169' or paired['job_id'] != '13461637'
            or (D.OUT/'JOBID').read_text().strip() != '13461637'
            or source['source_jobs']['native'] != paired['job_id']
            or source['reference_policy'] != paired['reference_policy']
            or source['reference_primary_chains'] != 22688
            or not source['one_to_one_complete_primary_join']
            or not source['full_and_species_counts_reproduce_completed_sources']
            or any(s['DEV_used'] or s['sealed_data_used'] or s['optimizer_steps_this_run'] != 0
                   for s in (source, paired))):
        raise ValueError('Completed original TRAIN supervision and paired reference sources required')
    rows = D.read_rows(SOURCE/'references.jsonl')
    if (len(rows) != 22688 or Counter(r['seqid'] for r in rows) != Counter(D.EXPECTED)
            or len({D.record_id(r) for r in rows}) != len(rows)
            or len({D.reference_key(r) for r in rows}) != len(rows)
            or sum(r['single_chain_supported'] for r in rows) != 22479
            or D.summarize(rows) != {k:v for k,v in paired['pooled'].items() if k != 'strata'}):
        raise ValueError('Complete original reference geometries, support or frozen hit memberships differ')
    for row in rows:
        cds = row['CDS']
        if (row['strand'] not in ('+', '-') or len(cds) != row['CDS_segments']
                or any(not 0 <= a < b for a,b,_ in cds)
                or any(cds[i][1] > cds[i+1][0] for i in range(len(cds)-1))
                or span(row)[1]-span(row)[0] != row['span_bp']
                or type(row['single_chain_supported']) is not bool):
            raise ValueError('Saved sorted half-open CDS geometry or support flag differs')
    pooled = describe(rows)
    species = {s:describe([r for r in rows if r['species'] == s])
               for s in source['per_species']}
    chromosomes = {q:describe([r for r in rows if r['seqid'] == q]) for q in D.EXPECTED}
    for field in ('reference_chains', 'individually_supported', 'geometry_only_maximum_hits',
                  'supported_geometry_maximum_hits', 'total_structurally_forced_omissions'):
        if pooled[field] != sum(v[field] for v in species.values()) or pooled[field] != sum(
                v[field] for v in chromosomes.values()):
            raise ValueError('Independent chromosome/species oracle counts do not sum')
    geometrical = maximum_compatible(rows)
    structural = maximum_compatible([r for r in rows if r['single_chain_supported']])
    report = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], plan=PLAN,
        status='reference_informed_original_TRAIN_structural_oracle_not_prediction_performance',
        fixture=checked, source_supervision_job=source['job_id'], source_paired_reference_job=paired['job_id'],
        source_observation_job=source['source_jobs']['original'],
        reference_policy=source['reference_policy'], reference_primary_chains=len(rows),
        reference_filters_added=False, all_original_reference_denominators_retained=True,
        structural_definition='cap128_individually_DNA_supported_CDS_spans_nonoverlapping_per_contig_and_strand',
        interval_definition='outer_CDS_zero_based_half_open_not_UTR_or_reference_gene_feature_span',
        algorithm='unit_weight_earliest_finish_interval_scheduling_per_species_contig_strand',
        tie_policy='CDS_end_then_CDS_start_then_composite_reference_identity',
        pooled=pooled, per_species=species, chromosomes=chromosomes,
        macro_species_structural_F1_ceiling=sum(v['structural_F1_ceiling'] for v in species.values())/len(species),
        worst_species_structural_F1_ceiling=min(v['structural_F1_ceiling'] for v in species.values()),
        DNA_read=False, real_reference_GFF_read=False, prediction_GFF_read=False, checkpoint_read=False,
        saved_reference_profiles_read=True, new_model_inference=False, predictions_modified=False,
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
        scientific_winner_or_DEV_candidate_selected=False,
        limitations=[
            'Reference-informed structural oracle, not a raw-FASTA prediction or model performance result.',
            'Fine/coarse observation availability, learned scores and neural realizability are not oracle constraints.',
            'Selected identities form one maximum-cardinality set; exclusion is not per-gene irrecoverability.',
            'Unsupported and overlapping references remain in all original evaluation denominators.',
            'Native hit projections maximize each fixed hit subset separately; they are not an altered baseline evaluation.',
            'Original TRAIN and upstream exposure cannot establish independent applicability or generalization.'])
    with (OUT/'oracle_membership.jsonl').open('x') as handle:
        for row in rows:
            start,end = span(row)
            handle.write(json.dumps(dict(species=row['species'], seqid=row['seqid'], strand=row['strand'],
                id=row['id'], CDS_span=[start,end], single_chain_supported=row['single_chain_supported'],
                geometry_only_oracle_selected=D.record_id(row) in geometrical,
                supported_geometry_oracle_selected=D.record_id(row) in structural), allow_nan=False)+'\n')
    report['seconds'] = time.perf_counter()-began
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(completed=True, seconds=report['seconds'], pooled=pooled)), flush=True)


if __name__ == '__main__':
    main()
