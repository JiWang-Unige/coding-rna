#!/usr/bin/env python3
"""Recompute the frozen A13 TRAIN gate from paired actual candidate records."""
import json
import math
from collections import Counter
from itertools import zip_longest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OLD = ROOT/'outputs/M28-A12-QUALITY-TRAIN-DIAGNOSTIC'
FIT = ROOT/'outputs/M28-A13-CHAIN-COMPETITION-FIT'
NEW = ROOT/'outputs/M28-A13-CHAIN-COMPETITION-TRAIN-READOUT'
EXPECTED_JOBS = {OLD:'13360650',FIT:'13360834',NEW:'13360836'}
GROUPS = ('positive','free_positive','known_nonmatch','five_or_more_positive','unknown')


def rows(path):
    with path.open() as handle:
        for line in handle:
            yield json.loads(line)


def compact(candidates):
    groups = dict(positive=[r for r in candidates if r['label'] == 1],
        free_positive=[r for r in candidates if r['label'] == 1 and r['origin'] == 'free'],
        known_nonmatch=[r for r in candidates if r['label'] == 0],
        five_or_more_positive=[r for r in candidates if r['label'] == 1 and len(r['chain']) >= 5],
        unknown=[r for r in candidates if r['label'] == -1])
    out = {}
    for name,values in groups.items():
        conflicts = [r for r in values if r['label'] == 1 and r['max_overlapping_known_nonmatch_gain'] is not None]
        out[name] = dict(candidates=len(values),positive_gain=sum(r['gain'] > 0 for r in values),
            nonpositive_gain=sum(r['gain'] <= 0 for r in values),
            positive_candidates_with_known_overlap=len(conflicts),
            positive_candidates_not_above_max_known_overlap=sum(
                r['gain'] <= r['max_overlapping_known_nonmatch_gain'] for r in conflicts))
    return out


def compare(actual, reported):
    for group in GROUPS:
        if any(reported[group][key] != value for key,value in actual[group].items()):
            raise ValueError('Actual TRAIN profile differs: '+group)


def main():
    for folder,job in EXPECTED_JOBS.items():
        if (folder/'STATUS').read_text().strip() != 'COMPLETED' or (folder/'JOBID').read_text().strip() != job:
            raise ValueError('Matching terminal TRAIN source required')
    fitted = json.loads((FIT/'summary.json').read_text())
    report = json.loads((NEW/'summary.json').read_text())
    if (fitted['completed'] is not True or fitted['optimizer_steps'] != 4608 or
        [e['steps'] for e in fitted['epochs']] != [1536]*3 or
        fitted['final_checkpoint'] != 'outputs/M28-A13-CHAIN-COMPETITION-FIT/step_004608.pt' or
        report['checkpoint'] != fitted['final_checkpoint'] or report['configuration'] != fitted['configuration'] or
        report['protocol'] != fitted['protocol'] or report['completed'] is not True or
        report['source_generator'] != 'outputs/M28-A3-FIRST-FIT/step_004608.pt' or
        report['source_A6'] != str(OLD.relative_to(ROOT)) or report['optimizer_steps_this_run'] != 0 or
        report['threshold'] != 0. or report['candidate_geometry_labels_and_origins_unchanged'] is not True or
        any(report[k] is not False for k in ('new_candidate_search','new_backbone_forward','DEV_forward_used',
            'test_setaria_used','threshold_scan','head_weights_updated','checkpoints_written'))):
        raise ValueError('Fixed zero-update TRAIN readout provenance differs')
    all_rows = {'A6':[],'A13':[]}; count = 0; labels = Counter()
    for old,new in zip_longest(rows(OLD/'windows.jsonl'),rows(NEW/'windows.jsonl')):
        if old is None or new is None or old['window'] != new['window'] or old['counts'] != new['counts']:
            raise ValueError('Actual ordered original TRAIN window/count identity changed')
        count += 1
        for a,b in zip_longest(old['candidates'],new['candidates']):
            if a is None or b is None or any(a[k] != b[k] for k in
                ('window_id','species','chain','origin','label','exon_class')):
                raise ValueError('Actual candidate geometry/labels/origins changed')
            labels[b['label']] += 1
        for arm,record in (('A6',old),('A13',new)):
            local = record['candidates']
            for row in local:
                if (not all(math.isfinite(row[k]) for k in ('gain','additive','nonadditive','null')) or
                    abs(row['gain']-(row['additive']+row['nonadditive']-row['null'])) > 1e-4):
                    raise ValueError('Actual gain/component identity differs')
                a,b = row['chain'][0][0],row['chain'][-1][1]
                scores = [r['gain'] for r in local if r['label'] == 0 and
                    r['chain'][0][0] < b and r['chain'][-1][1] > a] if row['label'] == 1 else []
                if row['max_overlapping_known_nonmatch_gain'] != (max(scores) if scores else None):
                    raise ValueError('Actual known span-conflict maximum differs')
            all_rows[arm].extend(local)
    if (count != 1491 or report['windows'] != count or report['original_draws'] != 1536 or
        len(all_rows['A13']) != 16857 or report['candidate_window_pairs'] != 16857 or
        dict(labels) != {1:2916,0:12566,-1:1375}):
        raise ValueError('Frozen original TRAIN population differs')
    pooled = {arm:compact(values) for arm,values in all_rows.items()}
    for arm in all_rows:
        compare(pooled[arm],report['pooled'][arm])
    per_species = {}; checks = {}
    for species in ('arabidopsis_thaliana','oryza_sativa'):
        values = {arm:compact([r for r in rs if r['species'] == species]) for arm,rs in all_rows.items()}
        for arm in all_rows:
            compare(values[arm],report['per_species'][species][arm])
        per_species[species] = values; a,b = values['A6'],values['A13']
        checks[species+'_free_correct_positive_gain_increased'] = b['free_positive']['positive_gain'] > a['free_positive']['positive_gain']
        checks[species+'_all_correct_positive_gain_increased'] = b['positive']['positive_gain'] > a['positive']['positive_gain']
        checks[species+'_known_nonmatch_positive_gain_not_increased'] = b['known_nonmatch']['positive_gain'] <= a['known_nonmatch']['positive_gain']
    a,b = pooled['A6'],pooled['A13']
    checks['pooled_five_or_more_correct_positive_gain_increased'] = b['five_or_more_positive']['positive_gain'] > a['five_or_more_positive']['positive_gain']
    checks['pooled_correct_not_above_known_conflict_not_increased'] = b['positive']['positive_candidates_not_above_max_known_overlap'] <= a['positive']['positive_candidates_not_above_max_known_overlap']
    if checks != report['predeclared_TRAIN_preparation_gate'] or all(checks.values()) != report['TRAIN_preparation_supported']:
        raise ValueError('Frozen TRAIN gate differs from actual records')
    print(json.dumps(dict(actual_records_verified=True,windows=count,candidate_window_pairs=len(all_rows['A13']),
        labels=dict(labels),per_species=per_species,pooled=pooled,gate=checks,
        TRAIN_preparation_supported=all(checks.values()),optimizer_steps_this_run=0,DEV_used=False,
        limitation='In-sample TRAIN gate, not final gene F1/background or independent generalization.')))


if __name__ == '__main__':
    main()
