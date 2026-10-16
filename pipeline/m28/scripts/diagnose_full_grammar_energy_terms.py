#!/usr/bin/env python3
"""Read existing fixed TRAIN chain-score terms; no forward, fitting or reweighting."""
import argparse
from collections import Counter
import json
import math
import os
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT/'outputs/M28-FULL-GRAMMAR-TRAIN-ENERGY-DIAGNOSTIC'
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-TRAIN-ENERGY-TERMS'
GROUPS = dict(coding_body=(1, 2, 3), intron_body=(4, 5, 6),
              splice_coding_boundary=(7, 8, 9, 10, 11, 12), start_stop=(13, 14))
SCOPE = {'NC_003070.9':'arabidopsis_thaliana', 'NC_003071.7':'arabidopsis_thaliana',
         'NC_003076.8':'arabidopsis_thaliana', 'NC_089043.1':'oryza_sativa',
         'NC_089044.1':'oryza_sativa', 'NC_089046.1':'oryza_sativa'}
TOLERANCE = 1e-5


def terms(row):
    gains, counts = row['gain_by_role'], row['rows_by_role']
    if len(gains) != 15 or len(counts) != 15 or counts[0] != 0 or gains[0] != 0:
        raise ValueError('Complete supported chain requires the existing fifteen-role path')
    if not math.isclose(math.fsum(gains), row['reference_gain_vs_background'],
                        rel_tol=1e-9, abs_tol=TOLERANCE):
        raise ValueError('Role sums do not reconstruct the fixed reference score')
    return {name:dict(gain=math.fsum(gains[i] for i in indices),
                      rows=sum(counts[i] for i in indices)) for name, indices in GROUPS.items()}


def summarize(rows):
    eligible = [r for r in rows if r['single_chain_supported']]
    grouped = [terms(r) for r in eligible]
    result = dict(reference=len(rows), individually_supported_reference=len(eligible),
                  unsupported_reference=len(rows)-len(eligible), categories={})
    for category in sorted({r['category'] for r in eligible}):
        pairs = [(r, g) for r, g in zip(eligible, grouped) if r['category'] == category]
        selected = [r for r, _ in pairs]
        group_results = {}
        for name in GROUPS:
            gain = math.fsum(g[name]['gain'] for _, g in pairs)
            count = sum(g[name]['rows'] for _, g in pairs)
            group_results[name] = dict(total_gain_vs_background=gain, base_rows=count,
                row_weighted_mean_gain=gain/count if count else None,
                mean_gain_per_chain=gain/len(pairs),
                negative_gain_chains=sum(g[name]['gain'] < -TOLERANCE for _, g in pairs),
                positive_gain_chains=sum(g[name]['gain'] > TOLERANCE for _, g in pairs))
        patterns = Counter()
        isolated_losses = Counter()
        for r, g in pairs:
            negative = [name for name in GROUPS if g[name]['gain'] < -TOLERANCE]
            patterns['+'.join(negative) or 'no_negative_group'] += 1
            if r['reference_gain_vs_background'] < -TOLERANCE:
                for name in negative:
                    other = math.fsum(value['gain'] for key, value in g.items() if key != name)
                    if other > TOLERANCE:
                        isolated_losses[name] += 1
        result['categories'][category] = dict(reference=len(pairs),
            total_reference_gain=math.fsum(r['reference_gain_vs_background'] for r in selected),
            role_groups=group_results, negative_group_sign_patterns=dict(patterns),
            per_role={str(i):dict(
                total_gain_vs_background=math.fsum(r['gain_by_role'][i] for r in selected),
                base_rows=sum(r['rows_by_role'][i] for r in selected),
                negative_gain_chains=sum(r['gain_by_role'][i] < -TOLERANCE for r in selected),
                positive_gain_chains=sum(r['gain_by_role'][i] > TOLERANCE for r in selected))
                for i in range(1, 15)},
            negative_whole_chain_with_this_group_negative_and_other_groups_positive=dict(isolated_losses))
    return result


def fixture():
    gains, counts = [0.]*15, [0]*15
    for index, gain, count in ((1, 2., 2), (4, -5., 5), (7, 1., 1), (13, 1., 3)):
        gains[index], counts[index] = gain, count
    row = dict(single_chain_supported=True, gain_by_role=gains, rows_by_role=counts,
               reference_gain_vs_background=-1., category='supported_reference_negative_gain_vs_background')
    value = summarize([row, dict(single_chain_supported=False)])
    category = value['categories'][row['category']]
    assert value['reference'] == 2 and value['unsupported_reference'] == 1
    assert category['role_groups']['intron_body']['row_weighted_mean_gain'] == -1.
    assert category['negative_group_sign_patterns'] == {'intron_body':1}
    assert category['negative_whole_chain_with_this_group_negative_and_other_groups_positive'] == {'intron_body':1}
    row['reference_gain_vs_background'] = 1.
    try:
        terms(row)
    except ValueError:
        return 'role_partition_score_reconstruction_signs_and_unsupported_denominator'
    raise AssertionError('Wrong role sum was accepted')


def main(seqid):
    began = time.perf_counter()
    checked = fixture()
    seqids = [seqid] if seqid else list(SCOPE)
    all_rows, source_jobs = [], {}
    for identifier in seqids:
        folder = SOURCE/identifier
        source = json.loads((folder/'summary.json').read_text())
        if ((folder/'STATUS').read_text().strip() != 'COMPLETED' or not source['completed']
                or source['seqid'] != identifier or source['species'] != SCOPE[identifier]
                or source['DEV_used'] or source['sealed_data_used']
                or source['optimizer_steps_this_run'] or source['predictions_modified']):
            raise ValueError('Completed unchanged fixed TRAIN source required')
        rows = [json.loads(line) for line in (folder/'reference_energy.jsonl').read_text().splitlines()]
        if (len(rows) != source['counts']['reference']
                or Counter(r['category'] for r in rows) != Counter(source['counts']['categories'])
                or any((r['seqid'], r['species']) != (identifier, SCOPE[identifier]) for r in rows)
                or sum(r['single_chain_supported'] for r in rows) != source['counts']['individually_supported_reference']):
            raise ValueError('Existing reference denominator, category or identity differs')
        all_rows.extend(rows)
        source_jobs[identifier] = source['job_id']
    target = OUT/(seqid or 'all_original_TRAIN')
    target.mkdir(parents=True, exist_ok=False)
    report = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], fixture=checked,
        scope='one_completed_TRAIN_chromosome' if seqid else 'six_completed_original_TRAIN_chromosomes',
        source_jobs=source_jobs, role_groups={k:list(v) for k, v in GROUPS.items()},
        pooled=summarize(all_rows),
        per_species={species:summarize([r for r in all_rows if r['species'] == species])
                     for species in sorted({r['species'] for r in all_rows})},
        checkpoint_read=False, DNA_or_reference_annotation_read=False, native_forward=False,
        predictions_modified=False, optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
        role_reweighting=False, near_zero_reporting_tolerance=TOLERANCE,
        term_partition='coding_body excludes splice_coding_boundary and start_stop; intron_body roles4-6',
        counts_are_chain_or_base_row_exposures_not_unique_genomic_coverage=True,
        competing_selected_role_terms_unavailable=True,
        interpretation='descriptive_reference_score_arithmetic_not_counterfactual_predictions_or_causal_attribution',
        seconds=time.perf_counter()-began)
    with (target/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write('\n')
    (target/'JOBID').write_text(os.environ['SLURM_JOB_ID']+'\n')
    (target/'STATUS').write_text('COMPLETED\n')
    print(json.dumps(report, allow_nan=False), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seqid', choices=list(SCOPE))
    main(parser.parse_args().seqid)
