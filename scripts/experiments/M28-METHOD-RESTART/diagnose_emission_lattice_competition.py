#!/usr/bin/env python3
"""Closed A42 saved chain competition, without search, tensors or new labels."""
from collections import Counter
from itertools import zip_longest
import math

import run_emission_lattice_train as D
import verify_emission_lattice_train as V

READOUT = D.ROOT/'outputs/M28-A42-EMISSION-LATTICE-TRAIN-READOUT-R1'
CHECK = D.ROOT/'outputs/M28-A42-EMISSION-LATTICE-RESULT-CHECK-R1'
OUT = D.ROOT/'outputs/M28-A42-EMISSION-LATTICE-COMPETITION'
JOBS = dict(potentials='13422268', decode='13422269', readout='13423713', CPU='13423721')
CPU_PROOFS = (
    'actual_full_ordered_1491_TRAIN_windows_original_1536_draws_identity_verified',
    'actual_frozen_A3_saved_potential_identity_and_no_update_receipts_verified',
    'actual_original_A31_FP64_scores_DP_and_stronger_baseline_recomputed',
    'actual_new_chains_DNA_ORF_motif_cap_span_and_FP64_gain_recomputed',
    'actual_complete_32_shards_natural_search_closed_before_labels_verified',
    'actual_frozen_R4_labels_unknown_protection_and_denominators_recomputed',
    'actual_per_species_counts_and_five_conditions_recomputed')


def validate_closed(report, checked):
    for value in (report, checked):
        if (value['completed'] is not True or value['configuration'] != D.CONFIG or
                value['protocol'] != D.PROTOCOL or value['windows'] != 1491 or
                value['original_draws'] != 1536 or value['original_candidates'] != 16321 or
                value['optimizer_steps_this_run'] != 0 or value['new_geometry_labels_computed'] is not True or
                any(value[k] is not False for k in ('model_or_backbone_forward', 'new_candidate_search',
                    'R4_policy_or_catalog_changed', 'DEV_used', 'test_setaria_used', 'A33_outputs_opened'))):
            raise ValueError('Complete original TRAIN-only A42 result required before saved rows')
    if (report['source_jobs'] != {k: JOBS[k] for k in ('potentials', 'decode')} or
            checked['source_jobs'] != {k: JOBS[k] for k in ('potentials', 'decode', 'readout')} or
            any(checked[k] is not True for k in CPU_PROOFS) or
            report['all_natural_search_and_score_files_closed_before_labels'] is not True or
            report['original_R4_labels_exact'] is not True or report['reference_injected_predictions'] != 0 or
            checked['independent_neural_forward_replay'] is not False):
        raise ValueError('All seven CPU proofs and actual source identities required before saved rows')
    refs = dict(arabidopsis_thaliana=dict(eligible=2172, five_or_more=854),
                oryza_sativa=dict(eligible=744, five_or_more=232))
    counts = report['selected_counts_per_species']; gates = D.gate(counts)
    if (set(counts) != set(refs) or report['eligible_reference_pairs'] != refs or
            checked['eligible_reference_pairs'] != refs or checked['selected_counts_per_species'] != counts or
            any(v['predeclared_TRAIN_gate'] != gates or v['TRAIN_preparation_supported'] != all(gates.values())
                for v in (report, checked))):
        raise ValueError('Closed counts, original denominators or five-condition arithmetic differ')
    return gates


def key(chain):
    return tuple(map(tuple, chain))


def close(a, b):
    return math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-8)


def census(saved, previous, labelled):
    groups = {s: {k: Counter(retained=0, lost=0, gained=0) for k in ('-1', '0', '1', 'correct_at_least_5_exons')}
              for s in D.R.C.R.ALLOW}
    overlaps = {s: {k: Counter() for k in ('correct', 'correct_at_least_5_exons')} for s in groups}
    totals = {s: {a: Counter() for a in D.ARMS} for s in groups}
    refs = {s: Counter(eligible=0, five_or_more=0) for s in groups}
    patterns = {s: Counter() for s in groups}
    windows = original_candidates = new_selected = lost_replacements = 0
    minimum_replacement_margin = None
    for new, prior, lab in zip_longest(saved, previous, labelled):
        if any(v is None for v in (new, prior, lab)):
            raise ValueError('Complete ordered original three-way windows required')
        old = prior['variants']['A31_A3']; labels = lab['labels']; w = new['window']
        if (w != prior['window'] or w != lab['window'] or w['split'] != 'train' or
                new['original_ordinal'] != windows or new['reference_used_in_scoring'] is not False or
                prior['reference_used_in_scoring'] is not False or new['selected'] != list(range(len(new['chains'])))):
            raise ValueError('Actual original ordered reference-free window identity differs')
        counts = {a: V.counts(v, labels[a]) for a, v in (('A31_A3', old), ('A42', new))}
        if counts != lab['counts']:
            raise ValueError('Saved two-arm selected counts or labels differ')
        before_pool = {key(c): (old['scores'][i], labels['A31_A3'][i]) for i, c in enumerate(old['chains'])}
        before = {key(old['chains'][i]): before_pool[key(old['chains'][i])] for i in old['selected']}
        after = {key(c): (new['scores'][i], labels['A42'][i]) for i, c in enumerate(new['chains'])}
        if len(before_pool) != len(old['chains']) or len(after) != len(new['chains']):
            raise ValueError('Saved natural geometry must be unique before identity joins')
        for geometry in before_pool.keys() & after.keys():
            if before_pool[geometry][1] != after[geometry][1] or not close(before_pool[geometry][0], after[geometry][0]):
                raise ValueError('Common geometry frozen labels or FP64 score differs')
        old_gain = math.fsum(v[0] for v in before.values()); new_gain = math.fsum(v[0] for v in after.values())
        tolerance = 1e-8 + 1e-12*max(abs(old_gain), abs(new_gain))
        if (not close(old_gain, new['original_baseline_gain']) or not close(new_gain, new['gain']) or
                new_gain + tolerance < old_gain or any(v[0] <= 0 for v in (*before.values(), *after.values()))):
            raise ValueError('Actual unchanged positive gains or stronger objective differ')
        s = w['species']
        for arm, value in counts.items(): totals[s][arm].update(value)
        refs[s].update(eligible=lab['eligible_references'], five_or_more=lab['eligible_five_or_more'])
        delta = counts['A42']['selected_correct'] - counts['A31_A3']['selected_correct']
        patterns[s].update({('objective_increased' if new_gain-old_gain > tolerance else 'objective_equal') +
            ('_correct_increased' if delta > 0 else '_correct_decreased' if delta < 0 else '_correct_equal'): 1})
        for transition, geometries, source in (
                ('retained', before.keys() & after.keys(), after),
                ('lost', before.keys() - after.keys(), before),
                ('gained', after.keys() - before.keys(), after)):
            for geometry in geometries:
                label = source[geometry][1]; groups[s][str(label)][transition] += 1
                if label == 1 and len(geometry) >= 5:
                    groups[s]['correct_at_least_5_exons'][transition] += 1
                if transition != 'lost': continue
                a, b = geometry[0][0], geometry[-1][1]
                blockers = [v for c, v in after.items() if c[0][0] < b and c[-1][1] > a]
                blocker_gain = math.fsum(v[0] for v in blockers)
                margin = blocker_gain-source[geometry][0]
                if margin < -(1e-8 + 1e-12*max(abs(blocker_gain), abs(source[geometry][0]))):
                    raise ValueError('Lost legal positive chain would improve full-grammar selected objective')
                lost_replacements += 1
                minimum_replacement_margin = margin if minimum_replacement_margin is None else min(minimum_replacement_margin, margin)
                if label == 1:
                    combination = ','.join(map(str, sorted({v[1] for v in blockers})))
                    for category in ['correct'] + (['correct_at_least_5_exons'] if len(geometry) >= 5 else []):
                        overlaps[s][category].update(chains=1, **{combination: 1})
        windows += 1; original_candidates += len(old['chains']); new_selected += len(after)
    pooled = {k: Counter() for k in ('-1', '0', '1', 'correct_at_least_5_exons')}
    pooled_overlaps = {k: Counter() for k in ('correct', 'correct_at_least_5_exons')}
    for s in groups:
        for category, value in groups[s].items():
            old_key = 'selected_five_or_more_correct' if category == 'correct_at_least_5_exons' else {
                '-1': 'selected_unknown', '0': 'selected_known_nonmatch', '1': 'selected_correct'}[category]
            if (value['retained']+value['lost'] != totals[s]['A31_A3'][old_key] or
                    value['retained']+value['gained'] != totals[s]['A42'][old_key]):
                raise ValueError('Exact geometry transition decomposition differs from selected counts')
            pooled[category].update(value)
        for category, value in overlaps[s].items(): pooled_overlaps[category].update(value)
    return dict(windows=windows, original_candidates=original_candidates, new_selected_chains=new_selected,
        selected_counts_per_species={s: {a: dict(v) for a, v in arms.items()} for s, arms in totals.items()},
        eligible_reference_pairs={s: dict(v) for s, v in refs.items()},
        per_species={s: {k: dict(v) for k, v in categories.items()} for s, categories in groups.items()},
        pooled={k: dict(v) for k, v in pooled.items()},
        paired_window_pattern={s: dict(v) for s, v in patterns.items()},
        lost_correct_selected_overlap_label_sets={s: {k: dict(v) for k, v in categories.items()} for s, categories in overlaps.items()},
        pooled_lost_correct_selected_overlap_label_sets={k: dict(v) for k, v in pooled_overlaps.items()},
        lost_selected_positive_replacement_checks=lost_replacements,
        minimum_selected_overlap_gain_minus_lost_chain_gain=minimum_replacement_margin)


def main():
    D.B.F.A.require_jobs([JOBS['readout'], JOBS['CPU'], '13398521', '13398524'])
    report = D.R.completed(READOUT, JOBS['readout']); checked = D.R.completed(CHECK, JOBS['CPU'])
    gates = validate_closed(report, checked)
    result = census(D.R.C.R.rows(READOUT/'natural_scores.jsonl'),
                    D.R.C.R.rows(D.R.READOUT/'natural_scores.jsonl'), D.R.C.R.rows(READOUT/'windows.jsonl'))
    if (result['windows'] != 1491 or result['original_candidates'] != 16321 or
            result['selected_counts_per_species'] != report['selected_counts_per_species'] or
            result['eligible_reference_pairs'] != report['eligible_reference_pairs']):
        raise ValueError('All original windows, two-arm counts and full denominators must reconcile')
    D.R.C.R.write_summary(OUT, dict(completed=True, experiment='M28-A42-EMISSION-LATTICE-COMPETITION',
        protocol=D.PROTOCOL, configuration=D.CONFIG, source_jobs=JOBS, **result,
        predeclared_TRAIN_gate=gates, TRAIN_preparation_supported=all(gates.values()),
        actual_saved_chain_identity_transitions_and_full_denominators_reconciled=True,
        actual_common_geometry_frozen_scores_and_labels_reconciled=True,
        actual_lost_selected_chain_local_replacement_objective_reconciled=True,
        grain='original_unique_TRAIN_window_chain_pair_not_unique_gene',
        model_or_backbone_forward=False, model_or_feature_tensor_opened=False, optimizer_steps_this_run=0,
        new_candidate_search=False, new_geometry_labels_computed=False, TRAIN_catalog_or_reference_opened=False,
        DEV_used=False, independent_species_used=False, test_setaria_used=False, A33_outputs_opened=False,
        limitation='Saved competition evidence, not training-loss causality, biological false-positive/new-gene proof, '
            'new scoring, search rescue, raw annotation, external SOTA or generalization.'))


if __name__ == '__main__':
    main()
