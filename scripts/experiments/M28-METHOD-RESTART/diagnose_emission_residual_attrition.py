#!/usr/bin/env python3
"""Closed A41 saved TRAIN score transitions; no fit, tensor access or rescue."""
from collections import Counter, defaultdict
from itertools import zip_longest
import math

import infer_emission_residual as I
from diagnose_structured_emission_attrition import transitions

D, R = I.D, I.D.R
OUT = I.ROOT/'outputs/M28-A41-EMISSION-RESIDUAL-ATTRITION'


def validate_closed(fitted, report, checked):
    D.validate_fit(fitted)
    for value in (report, checked):
        if (value['completed'] is not True or value['configuration'] != D.CONFIG or
                value['protocol'] != D.PROTOCOL or value['windows'] != 1491 or
                value['natural_candidates'] != 16321 or value['potential_receipts'] != 1491 or
                value['optimizer_steps_this_run'] != 0 or
                any(value[k] is not False for k in ('new_candidate_search', 'new_R4_labelling',
                    'DEV_used', 'test_setaria_used', 'A33_outputs_opened'))):
            raise ValueError('Complete original TRAIN-only A41 result required before score access')
    if (report['source_fit_job'] != I.FIT_JOB or report['source_checks_job'] != '13420720' or
            report['source_preflight_job'] != D.F.PREFLIGHT_JOB or
            report['checkpoint'] != fitted['final_checkpoint'] or report['original_draws'] != 1536 or
            report['reference_injected_predictions'] != 0 or report['label_catalog_opened'] is not False or
            checked['source_jobs'] != dict(fit=I.FIT_JOB, readout=I.READOUT_JOB) or
            any(checked[k] is not True for k in I.CPU_PROOFS)):
        raise ValueError('Actual A41 source jobs and six CPU reconciliations required before score access')
    counts, refs = report['selected_counts_per_species'], report['eligible_reference_pairs']
    expected_refs = dict(arabidopsis_thaliana=dict(eligible=2172, five_or_more=854),
                         oryza_sativa=dict(eligible=744, five_or_more=232))
    if (set(counts) != set(expected_refs) or refs != expected_refs or
            checked['selected_counts_per_species'] != counts or checked['eligible_reference_pairs'] != refs):
        raise ValueError('Both species and original complete paired counts and denominators required')
    gates = D.gate(counts)
    if (report['predeclared_TRAIN_gate'] != gates or checked['predeclared_TRAIN_gate'] != gates or
            report['TRAIN_preparation_supported'] is not False or checked['TRAIN_preparation_supported'] is not False or
            all(gates.values())):
        raise ValueError('Actual closed A41 scientific FAIL required, not raw execution readiness')
    return gates


def lost_positive_blockers(chains, labels, old, new):
    """Label composition of actual selected overlaps, not a counterfactual decode."""
    groups = {k: Counter() for k in ('correct', 'correct_at_least_5_exons')}
    before, after = set(old['selected']), set(new['selected'])
    for i in sorted(before - after):
        if labels[i] != 1 or new['scores'][i] <= 0:
            continue
        a, b = chains[i][0][0], chains[i][-1][1]
        blockers = [j for j in new['selected'] if chains[j][0][0] < b and chains[j][-1][1] > a]
        gain = math.fsum(new['scores'][j] for j in blockers)
        if new['scores'][i] - gain > 1e-8 + 1e-12*abs(gain):
            raise ValueError('Positive lost chain would improve the saved original global DP objective')
        combination = ','.join(map(str, sorted({labels[j] for j in blockers}))) if blockers else 'no_selected_overlap'
        keys = ['correct'] + (['correct_at_least_5_exons'] if len(chains[i]) >= 5 else [])
        for key in keys:
            groups[key].update(chains=1, **{combination: 1})
    return groups


def census(saved, previous, labelled):
    groups = {s: defaultdict(Counter) for s in R.C.R.ALLOW}
    blockers = {s: {k: Counter() for k in ('correct', 'correct_at_least_5_exons')} for s in groups}
    totals = {s: {a: Counter() for a in ('A31_A3', 'A41')} for s in groups}
    refs = {s: Counter(eligible=0, five_or_more=0) for s in groups}
    windows = candidates = 0; maximum = 0.
    for new, old, lab in zip_longest(saved, previous, labelled):
        if any(v is None for v in (new, old, lab)):
            raise ValueError('Complete ordered original three-way windows required')
        baseline = old['variants']['A31_A3']; labels = lab['labels']
        if (new['window'] != old['window'] or new['window'] != lab['window'] or
                new['chains'] != baseline['chains'] or new['window']['split'] != 'train' or
                new['reference_used_in_scoring'] is not False or old['reference_used_in_scoring'] is not False):
            raise ValueError('Original reference-free same-natural-pool TRAIN identity differs')
        I.verify_score_components(dict(score=new['scores'], frozen_path_gain=new['frozen_path_gain'],
                                       residual_path_gain=new['residual_path_gain'], null=new['null']))
        for a, b in zip(new['frozen_path_gain'], baseline['scores']):
            if not math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-8):
                raise ValueError('Actual frozen component differs from original A31 score')
            maximum = max(maximum, abs(a-b))
        counts = dict(A31_A3=R.selected_counts(baseline, labels), A41=R.selected_counts(new, labels))
        if counts != lab['counts']:
            raise ValueError('Actual original two-arm per-window counts differ')
        species = new['window']['species']
        for arm, value in counts.items():
            totals[species][arm].update(value)
        refs[species].update(eligible=lab['eligible_references'], five_or_more=lab['eligible_five_or_more'])
        for key, value in transitions(new['chains'], labels, baseline, new).items():
            groups[species][key].update(value)
        for key, value in lost_positive_blockers(new['chains'], labels, baseline, new).items():
            blockers[species][key].update(value)
        windows += 1; candidates += len(labels)
    pooled = defaultdict(Counter); pooled_blockers = {k: Counter() for k in ('correct', 'correct_at_least_5_exons')}
    for species in groups:
        for key, value in groups[species].items():
            pooled[key].update(value)
        for key, value in blockers[species].items():
            pooled_blockers[key].update(value)
    return dict(windows=windows, natural_candidates=candidates,
        selected_counts_per_species={s: {a: dict(c) for a, c in v.items()} for s, v in totals.items()},
        eligible_reference_pairs={s: dict(c) for s, c in refs.items()},
        per_species={s: {k: dict(c) for k, c in v.items()} for s, v in groups.items()},
        pooled={k: dict(c) for k, c in pooled.items()},
        lost_positive_correct_blocker_label_sets={s: {k: dict(c) for k, c in v.items()} for s, v in blockers.items()},
        pooled_lost_positive_correct_blocker_label_sets={k: dict(c) for k, c in pooled_blockers.items()},
        max_original_A31_frozen_component_error=maximum)


def main():
    D.F.A.A.require_jobs([I.FIT_JOB, I.READOUT_JOB, I.CHECK_JOB, '13398521', '13398524'])
    fitted = R.completed(D.F.OUT, I.FIT_JOB)
    report = R.completed(D.OUT, I.READOUT_JOB); checked = R.completed(D.CHECK, I.CHECK_JOB)
    gates = validate_closed(fitted, report, checked)
    baseline = R.completed(R.READOUT, '13398521')
    result = census(R.C.R.rows(D.OUT/'natural_scores.jsonl'), R.C.R.rows(R.READOUT/'natural_scores.jsonl'),
                    R.C.R.rows(D.OUT/'windows.jsonl'))
    counts = result['selected_counts_per_species']
    if (result['windows'] != 1491 or result['natural_candidates'] != 16321 or
            counts != report['selected_counts_per_species'] or
            result['eligible_reference_pairs'] != report['eligible_reference_pairs'] or
            any(counts[s]['A31_A3'] != baseline['selected_counts_per_species'][s]['A31_A3'] for s in counts)):
        raise ValueError('Complete transition census must reproduce two closed arms and full denominators')
    R.C.R.write_summary(OUT, dict(completed=True, experiment='M28-A41-EMISSION-RESIDUAL-ATTRITION',
        protocol=D.PROTOCOL, configuration=D.CONFIG,
        source_jobs=dict(fit=I.FIT_JOB, readout=I.READOUT_JOB, CPU=I.CHECK_JOB, baseline_readout='13398521'),
        scientific_status='A41_TRAIN_FAIL_closed_no_rescue',
        grain='original_unique_TRAIN_window_natural_chain_pair_not_unique_gene',
        **result, predeclared_TRAIN_gate=gates,
        transition_logic='unchanged_diagnose_structured_emission_attrition.transitions',
        exact_saved_separate_FP64_sum_reconciled=True,
        actual_selected_overlap_label_sets_and_replacement_objective_reconciled=True,
        actual_saved_two_arm_counts_and_full_denominators_recomputed=True,
        model_or_backbone_forward=False, model_or_feature_tensor_opened=False, optimizer_steps_this_run=0,
        new_candidate_search=False, new_R4_labelling=False, TRAIN_catalog_or_reference_opened=False,
        DEV_used=False, test_setaria_used=False, A33_outputs_opened=False,
        limitation='Descriptive closed score transitions and actual selected overlaps, not training-loss causality, '
            'role ablation, new gene/biological false-positive claim, threshold rescue, raw annotation or generalization.'))


if __name__ == '__main__':
    main()
