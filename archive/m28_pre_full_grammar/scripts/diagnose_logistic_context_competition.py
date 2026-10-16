#!/usr/bin/env python3
"""Actual selected blockers of closed A38 positive-gain missed TRAIN chains."""
from collections import Counter
from itertools import zip_longest
import json
import diagnose_logistic_context_attrition as A
import audit_chain_displacement as G

D, R, I = A.D, A.R, A.I
OUT = I.ROOT/'outputs/M28-A38-LOGISTIC-CONTEXT-COMPETITION'
ATTRITION_JOB = '13414584'


def diagnose_window(new, old, lab):
    baseline, labels = old['variants']['A31_A3'], lab['labels']
    if (new['window'] != old['window'] or new['window'] != lab['window'] or
            new['chains'] != baseline['chains'] or new['reference_used_in_scoring'] or
            old['reference_used_in_scoring'] or new['window']['split'] != 'train' or
            any(y not in (-1, 0, 1) for y in labels)):
        raise ValueError('Ordered reference-free same-pool TRAIN identity and protected labels required')
    counts = {arm: R.selected_counts(value, labels) for arm, value in
              (('A31_A3', baseline), ('A38', new))}
    if any(counts[arm] != lab['counts'][arm] for arm in counts):
        raise ValueError('Actual per-window baseline and A38 selections must reproduce closed readout')
    local = [dict(chain=c, gain=g, label=y, origin='free') for c, g, y in
             zip(new['chains'], new['scores'], labels)]
    chosen, before = set(new['selected']), set(baseline['selected'])
    missed = [i for i, r in enumerate(local) if r['label'] == 1 and r['gain'] > 0 and i not in chosen]
    cases = []
    for i in missed:
        value = G.exchange(local, new['selected'], [i])
        blockers = value['selected_conflict_indices']
        if not blockers:
            raise ValueError('Positive missed chain has no actual selected conflict')
        cases.append(dict(candidate_index=i, chain=local[i]['chain'], gain=local[i]['gain'],
            baseline_selected=i in before, five_or_more=len(local[i]['chain']) >= 5,
            **value, all_individual_blockers_strictly_weaker=all(local[b]['gain'] < local[i]['gain'] for b in blockers)))
    joint = G.exchange(local, new['selected'], missed)
    return dict(window=new['window'], counts=counts,
        eligible_references=lab['eligible_references'], eligible_five_or_more=lab['eligible_five_or_more'],
        cases=cases, joint_positive_unselected_exchange=joint)


def profile(items):
    return dict(pairs=len(items), multi_blocker_pairs=sum(len(v['selected_conflict_indices']) >= 2 for v in items),
        all_individual_blockers_strictly_weaker_pairs=sum(v['all_individual_blockers_strictly_weaker'] for v in items),
        blocker_label_combinations=dict(Counter('+'.join(str(y) for y in (1, 0, -1)
            if v['conflict_counts_by_label'][str(y)]) for v in items)),
        blocker_members_by_label={str(y): sum(v['conflict_counts_by_label'][str(y)] for v in items)
                                  for y in (1, 0, -1)})


def profiles(items):
    return dict(all_positive_unselected=profile(items),
        baseline_lost_positive=profile([v for v in items if v['baseline_selected']]),
        all_positive_unselected_five_or_more=profile([v for v in items if v['five_or_more']]),
        baseline_lost_positive_five_or_more=profile([v for v in items if v['baseline_selected'] and v['five_or_more']]))


def census(saved, previous, labelled):
    records, items = [], {s: [] for s in R.C.R.ALLOW}
    totals = {s: {a: Counter() for a in ('A31_A3', 'A38')} for s in items}
    refs = {s: Counter(eligible=0, five_or_more=0) for s in items}
    candidates = 0
    for new, old, lab in zip_longest(saved, previous, labelled):
        if any(v is None for v in (new, old, lab)):
            raise ValueError('Complete three-way original windows required')
        record = diagnose_window(new, old, lab); records.append(record)
        species = new['window']['species']; items[species].extend(record['cases'])
        for arm in totals[species]:
            totals[species][arm].update(record['counts'][arm])
        refs[species].update(eligible=record['eligible_references'], five_or_more=record['eligible_five_or_more'])
        candidates += len(new['chains'])
    return records, dict(windows=len(records), natural_candidates=candidates,
        selected_counts_per_species={s: {a: dict(v) for a, v in t.items()} for s, t in totals.items()},
        eligible_reference_pairs={s: dict(v) for s, v in refs.items()},
        per_species={s: profiles(v) for s, v in items.items()},
        pooled=profiles([v for local in items.values() for v in local]))


def reconcile(result, report, attrition):
    if (result['windows'] != 1491 or result['natural_candidates'] != 16321 or
            result['eligible_reference_pairs'] != report['eligible_reference_pairs'] or
            any(result['selected_counts_per_species'][s][a] != report['selected_counts_per_species'][s][a]
                for s in result['selected_counts_per_species'] for a in ('A31_A3', 'A38'))):
        raise ValueError('Full population, closed selections and original denominators must reproduce')
    for species, current in {**result['per_species'], 'pooled': result['pooled']}.items():
        old = (attrition['pooled'] if species == 'pooled' else attrition['per_species'][species])['A31_A3_to_A38']
        expected = dict(all_positive_unselected=old['1']['lost_positive_not_selected']+old['1']['neither_positive_not_selected'],
            baseline_lost_positive=old['1']['lost_positive_not_selected'],
            all_positive_unselected_five_or_more=old['correct_at_least_5_exons']['lost_positive_not_selected']+
                old['correct_at_least_5_exons']['neither_positive_not_selected'],
            baseline_lost_positive_five_or_more=old['correct_at_least_5_exons']['lost_positive_not_selected'])
        if any(current[k]['pairs'] != v for k, v in expected.items()):
            raise ValueError('All-positive and baseline-lost strata must reproduce completed attrition')


def main():
    D.F.A.require_jobs([I.FIT_JOB, I.READOUT_JOB, I.CHECK_JOB, ATTRITION_JOB])
    fitted = R.completed(D.F.OUT, I.FIT_JOB)
    report = R.completed(D.OUT, I.READOUT_JOB); checked = R.completed(D.CHECK, I.CHECK_JOB)
    gates = A.validate_completed(fitted, report, checked)
    attrition = R.completed(A.OUT, ATTRITION_JOB)
    if (all(gates.values()) or attrition['scientific_status'] != 'A38_TRAIN_FAIL_closed_no_rescue' or
            attrition['configuration'] != D.CONFIG or attrition['protocol'] != D.PROTOCOL or
            attrition['source_jobs'] != dict(fit=I.FIT_JOB, readout=I.READOUT_JOB, CPU=I.CHECK_JOB,
                baseline_readout='13398521', context_readout='13410739', context_CPU='13410742') or
            attrition['selected_counts_per_species'] != report['selected_counts_per_species'] or
            attrition['eligible_reference_pairs'] != report['eligible_reference_pairs']):
        raise ValueError('Actual closed A38 and complete matching attrition required before score access')
    R.completed(R.READOUT, '13398521')
    records, result = census(R.C.R.rows(D.OUT/'natural_scores.jsonl'),
        R.C.R.rows(R.READOUT/'natural_scores.jsonl'), R.C.R.rows(D.OUT/'windows.jsonl'))
    reconcile(result, report, attrition)
    with (OUT/'windows.jsonl').open('x') as handle:
        for record in records:
            handle.write(json.dumps(record, allow_nan=False)+'\n')
    R.C.R.write_summary(OUT, dict(completed=True, experiment='M28-A38-LOGISTIC-CONTEXT-COMPETITION',
        protocol=D.PROTOCOL, configuration=D.CONFIG,
        source_jobs=dict(fit=I.FIT_JOB, readout=I.READOUT_JOB, CPU=I.CHECK_JOB, attrition=ATTRITION_JOB,
            baseline_readout='13398521'), scientific_status='A38_TRAIN_FAIL_closed_no_rescue', **result,
        all_individual_and_joint_actual_selected_exchanges_verified=True,
        exchange_logic='unchanged_audit_chain_displacement.exchange', exchange_relative_roundoff_allowance=1e-6,
        gain_threshold=0., threshold_scan=False, model_or_backbone_forward=False, model_or_feature_tensor_opened=False,
        optimizer_steps_this_run=0, new_candidate_search=False, new_R4_labelling=False,
        TRAIN_catalog_or_reference_opened=False, DEV_used=False, test_setaria_used=False, A33_outputs_opened=False,
        limitation='Actual selected blocker composition and feasible exchange, not a best forced-set optimum, '
            'single-cause proof, biological false-positive determination, threshold rescue or gene F1. '
            'Unknown blockers stay protected; repeated blockers across missed pairs are not unique genes.'))


if __name__ == '__main__':
    main()
