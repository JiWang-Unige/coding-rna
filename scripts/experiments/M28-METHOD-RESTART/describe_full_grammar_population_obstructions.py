#!/usr/bin/env python3
"""Classify saved original-TRAIN fine supervision gaps, without DNA or model reads."""
from collections import Counter
import json
import os
from pathlib import Path
import time

ROOT = Path('/home/users/j/jwang/coding-rna')
SOURCE = ROOT/'outputs/M28-FULL-GRAMMAR-POPULATION-SUPERVISION'
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-POPULATION-OBSTRUCTIONS'
PLAN = 'reports/M28-METHOD-RESTART/full_grammar_fresh_draws_plan.md'
BLOCKED = 'containing_cores_without_complete_matching_fine'
AVAILABLE = 'complete_matching_fine_core_available'
UNSUPPORTED = 'outside_fixed_single_chain_grammar'
NO_CORE = ('span_exceeds_existing_core_width', 'no_containing_existing_core', UNSUPPORTED)
COUNT_FIELDS = ('matching_fine_rows', 'disagreeing_fine_rows', 'coarse_only_rows', 'fully_unknown_rows')
SCOPE = {
    'arabidopsis_thaliana': {'NC_003070.9':7135, 'NC_003071.7':4293, 'NC_003076.8':6311},
    'oryza_sativa': {'NC_089043.1':1638, 'NC_089044.1':1601, 'NC_089046.1':1710}}


def classify(row):
    category = row['population_supervision_category']
    counts = row['population_containing_core_row_counts']
    span = row['CDS'][-1][1]-row['CDS'][0][0]
    if (span != row['span_bp'] or len(counts) != row['population_containing_core_count']
            or len({v['window_id'] for v in counts}) != len(counts)
            or row['single_chain_supported'] != (category != UNSUPPORTED)):
        raise ValueError('Saved span, containing-core count or grammar status differs')
    for value in counts:
        if (any(type(value[k]) is not int or value[k] < 0 for k in COUNT_FIELDS)
                or sum(value[k] for k in COUNT_FIELDS) != span
                or value['complete_matching_fine'] != (value['matching_fine_rows'] == span)):
            raise ValueError('Every containing-core observation must partition the same reference span')
    complete = [v['window_id'] for v in counts if v['complete_matching_fine']]
    if set(complete) != set(row['population_matching_core_window_ids']):
        raise ValueError('Saved complete fine qualification differs from row counts')
    compatible = [v for v in counts if v['disagreeing_fine_rows'] == 0]
    gaps = set()
    for value in compatible:
        coarse, unknown = value['coarse_only_rows'], value['fully_unknown_rows']
        gaps.add('mixed' if coarse and unknown else 'coarse_only' if coarse else 'unknown_only' if unknown else 'none')
    if category in NO_CORE:
        if counts or complete:
            raise ValueError('Non-evaluated geometry/grammar strata gained a containing-core record')
        obstruction = category
    elif category == AVAILABLE:
        if not complete:
            raise ValueError('Available reference lost every complete fine core')
        obstruction = AVAILABLE
    elif category == BLOCKED:
        if not counts or complete or 'none' in gaps:
            raise ValueError('Blocked reference must retain containing but incomplete fine cores')
        obstruction = 'all_containing_cores_have_fine_conflict' if not compatible else 'fine_compatible_gaps:'+'|'.join(sorted(gaps))
    else:
        raise ValueError('Unexpected population supervision category')
    return dict(**{k:row[k] for k in ('species','seqid','strand','id','single_chain_supported')},
        population_supervision_category=category, row_state_obstruction=obstruction,
        containing_core_windows=len(counts), fine_compatible_core_windows=len(compatible),
        compatible_window_ids=[v['window_id'] for v in compatible], compatible_gap_types=sorted(gaps),
        fine_conflict_rows_range=[min(v['disagreeing_fine_rows'] for v in counts),
                                  max(v['disagreeing_fine_rows'] for v in counts)] if counts else None,
        compatible_coarse_only_rows_range=[min(v['coarse_only_rows'] for v in compatible),
                                           max(v['coarse_only_rows'] for v in compatible)] if compatible else None,
        compatible_fully_unknown_rows_range=[min(v['fully_unknown_rows'] for v in compatible),
                                             max(v['fully_unknown_rows'] for v in compatible)] if compatible else None,
        original_exact_hit=row['exact_primary_hit'], native_exact_hits=row['native_exact_hits'])


def describe(rows):
    blocked = [r for r in rows if r['population_supervision_category'] == BLOCKED]
    hard = sum(r['fine_compatible_core_windows'] == 0 for r in blocked)
    return dict(reference_chains=len(rows), supported=sum(r['single_chain_supported'] for r in rows),
        population_categories=dict(Counter(r['population_supervision_category'] for r in rows)),
        row_state_obstructions=dict(Counter(r['row_state_obstruction'] for r in rows)),
        containing_but_without_complete_fine=len(blocked), all_containing_cores_have_fine_conflict=hard,
        has_fine_compatible_containing_core=len(blocked)-hard,
        compatible_gap_patterns=dict(Counter('|'.join(r['compatible_gap_types']) for r in blocked
                                            if r['fine_compatible_core_windows'] > 0)),
        any_containing_core_has_fine_conflict=sum(r['fine_conflict_rows_range'] is not None
            and r['fine_conflict_rows_range'][1] > 0 for r in rows))


def fixture():
    def make(values, category=BLOCKED):
        counts = [dict(window_id='w'+str(i), **dict(zip(COUNT_FIELDS,v)),
            complete_matching_fine=v[0] == 9) for i,v in enumerate(values)]
        return dict(species='fixture',seqid='q',strand='+',id='t',CDS=[[0,9,0]],span_bp=9,
            single_chain_supported=category != UNSUPPORTED, population_supervision_category=category,
            population_containing_core_row_counts=counts, population_containing_core_count=len(counts),
            population_matching_core_window_ids=[v['window_id'] for v in counts if v['complete_matching_fine']],
            exact_primary_hit=False,native_exact_hits={m:False for m in ('ANNEVO','Tiberius','Helixer')})
    examples = [
        (make([(9,0,0,0)],AVAILABLE), AVAILABLE),
        (make([(6,0,3,0)]), 'fine_compatible_gaps:coarse_only'),
        (make([(6,0,0,3)]), 'fine_compatible_gaps:unknown_only'),
        (make([(6,0,1,2)]), 'fine_compatible_gaps:mixed'),
        (make([(6,0,3,0),(6,0,0,3)]), 'fine_compatible_gaps:coarse_only|unknown_only'),
        (make([(5,1,0,3),(5,2,2,0)]), 'all_containing_cores_have_fine_conflict'),
        (make([(5,1,0,3),(6,0,0,3)]), 'fine_compatible_gaps:unknown_only')]
    examples += [(make([],cat),cat) for cat in NO_CORE]
    for row,expected in examples:
        assert classify(row)['row_state_obstruction'] == expected
    for bad in (make([(6,0,0,2)]), make([(9,0,0,0)])):
        try:
            classify(bad)
        except ValueError:
            pass
        else:
            raise AssertionError('Invalid span conservation or blocked/complete qualification accepted')
    return dict(valid_cases=len(examples), invalid_cases_rejected=2,
        definition='all_cores_vs_any_compatible_core_coarse_unknown_mixed_combinations_unsupported_and_geometry_retained')


def main():
    began = time.perf_counter()
    checked = fixture()
    source = json.loads((SOURCE/'summary.json').read_text())
    if ((SOURCE/'STATUS').read_text().strip() != 'COMPLETED' or not source['completed']
            or (SOURCE/'JOBID').read_text().strip() != '13462343' or source['job_id'] != '13462343'
            or source['reference_primary_chains'] != 22688 or source['full_existing_TRAIN_windows'] != 25512
            or source['reference_filters_added'] or source['labels_changed'] or source['DEV_used']
            or source['sealed_data_used'] or source['optimizer_steps_this_run'] != 0):
        raise ValueError('Completed fixed original-six TRAIN population census required')
    records,ids,scope = [],set(),Counter()
    with (SOURCE/'references.jsonl').open() as handle:
        for line in handle:
            row = json.loads(line)
            key = tuple(row[k] for k in ('species','seqid','strand','id'))
            if key in ids or row['seqid'] not in SCOPE.get(row['species'],{}):
                raise ValueError('Duplicate or out-of-scope primary reference')
            ids.add(key)
            scope[row['species'],row['seqid']] += 1
            records.append(classify(row))
    if scope != Counter({(s,q):n for s,contigs in SCOPE.items() for q,n in contigs.items()}):
        raise ValueError('All six unchanged reference denominators required')
    pooled = describe(records)
    per_species = {s:describe([r for r in records if r['species'] == s]) for s in SCOPE}
    for observed,expected in [(pooled,source['pooled']['all_references'])]+[
            (per_species[s],source['per_species'][s]['all_references']) for s in SCOPE]:
        if any(observed[k] != expected[k] for k in ('reference_chains','supported')) or observed['population_categories'] != expected['categories']:
            raise ValueError('Whole population and species eligibility strata must reproduce the census')
    if pooled['containing_but_without_complete_fine'] != 6220:
        raise ValueError('Every declared containing-but-incomplete reference required')
    native_gaps = [r for r in records if not r['original_exact_hit'] and any(r['native_exact_hits'].values())]
    consensus_gaps = [r for r in native_gaps if all(r['native_exact_hits'].values())]
    if len(native_gaps) != 7831 or len(consensus_gaps) != 4954:
        raise ValueError('Frozen native diagnostic group memberships differ')
    report = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],plan=PLAN,fixture=checked,fixture_only=False,
        status='completed_saved_original_TRAIN_fine_row_state_obstruction_not_model_performance',
        source_census_job=source['job_id'],reference_policy=source['reference_policy'],reference_primary_chains=len(records),
        all_saved_containing_core_row_partitions_verified=True,source_population_and_species_categories_reproduced=True,
        pooled=pooled,per_species=per_species,
        fixed_any_native_anchor_gap=describe(native_gaps),fixed_all_three_native_anchor_gap=describe(consensus_gaps),
        saved_population_reference_counts_read=True,DNA_read=False,real_reference_GFF_read=False,
        R4_labels_rebuilt=False,checkpoint_read=False,feature_tensors_read=False,new_model_inference=False,
        labels_or_predictions_modified=False,optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,
        sampling_or_candidate_selected=False,baseline_or_model_outcomes_used_for_obstruction_classification=False,
        limitations=[
            'Fine compatibility is confined to the saved reference span, not the full outside-span fine event or gene-set support.',
            'Gap types describe missing fine observations, not why R4 masks were generated or whether coarse regions agree.',
            'Multiple overlapping containing-core observations are not independent reference or genomic samples.',
            'Coarse channels are positive at epsilon0.1; their presence is not absent gradient or a hard fine contradiction.',
            'Native-gap groups are descriptive fixed outcome subsets, never training-selection rules.',
            'No label repair, new fit, frozen DEV change, causal mechanism or generalization is justified by this classification alone.'],
        seconds=time.perf_counter()-began)
    with (OUT/'references.jsonl').open('x') as handle:
        for row in records:
            handle.write(json.dumps(row,allow_nan=False)+'\n')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(completed=True,reference_primary_chains=len(records),
        containing_but_without_complete_fine=6220)),flush=True)


if __name__ == '__main__':
    main()
