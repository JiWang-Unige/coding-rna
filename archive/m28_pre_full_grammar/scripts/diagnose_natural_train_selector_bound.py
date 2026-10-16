#!/usr/bin/env python3
"""Label-assisted fixed A3 TRAIN pool bounds, not model scores or a rescue."""
from collections import Counter
import json
import subprocess
from itertools import zip_longest
import run_local_path_caller as R

OUT = R.ROOT/'outputs/M28-A31-NATURAL-TRAIN-SELECTOR-BOUND'
PROTOCOL = 'reports/M28-METHOD-RESTART/a31_score_pool_train_diagnostic_protocol.md'
JOBS = dict(readout='13398521', CPU='13398524')
MODES = ('maximum_correct_then_long', 'maximum_long_only')


def source_jobs():
    result = subprocess.run(['sacct', '-X', '-j', ','.join(JOBS.values()), '-n', '-P',
        '--format=JobID%30,State%30,ExitCode'], check=True, capture_output=True, text=True)
    actual = {p[0].strip(): [v.strip() for v in p[1:3]]
              for line in result.stdout.splitlines() if len(p := line.split('|')) >= 3}
    if any(actual.get(job) != ['COMPLETED', '0:0'] for job in JOBS.values()):
        raise ValueError('Successful actual closed A31 source roots required')


def validate_sources(prepared, checked):
    for value in (prepared, checked):
        if (value['completed'] is not True or value['configuration'] != R.CONFIG or
                value['protocol'] != R.PROTOCOL or value['windows'] != 1491):
            raise ValueError('Original complete unchanged TRAIN source required')
    if (checked['source_jobs'] != dict(readout=JOBS['readout']) or
            checked['actual_saved_natural_pool_and_greedy_subset_identity_verified'] is not True or
            checked['actual_R4_labels_selection_reference_denominators_and_proxy_metrics_recomputed'] is not True or
            checked['selected_counts_per_species'] != prepared['selected_counts_per_species'] or
            checked['eligible_reference_pairs'] != prepared['eligible_reference_pairs']):
        raise ValueError('Actual original CPU reconciliation required before label access')


def bounds(chains, labels):
    if len(chains) != len(labels) or any(y not in (-1, 0, 1) for y in labels):
        raise ValueError('Original ordered complete protected labels required')
    # One extra correct outweighs all possible long-chain tie bonuses in this window.
    base = len(chains)+1
    weights = dict(maximum_correct_then_long=[base+int(len(c) >= 5) if y == 1 else 0.
                   for c, y in zip(chains, labels)],
                   maximum_long_only=[float(y == 1 and len(c) >= 5) for c, y in zip(chains, labels)])
    result = {}
    for mode, scores in weights.items():
        selected = R.select_nonoverlapping(chains, scores)
        result[mode] = dict(selected=selected, counts=R.selected_counts(
            dict(chains=chains, scores=scores, selected=selected), labels))
    return result


def interpretation(totals):
    individual = {}; witness = {}
    for species, arms in totals.items():
        old, new = arms['A31_A3'], arms[MODES[0]]
        key = species+'_selected_correct_strictly_increased'
        individual[key] = new['selected_correct'] > old['selected_correct']
        witness[key] = individual[key]
        witness[species+'_known_nonmatch_not_increased'] = new['selected_known_nonmatch'] <= old['selected_known_nonmatch']
    old_long = sum(v['A31_A3']['selected_five_or_more_correct'] for v in totals.values())
    individual['pooled_long_correct_strictly_increased'] = sum(
        v[MODES[1]]['selected_five_or_more_correct'] for v in totals.values()) > old_long
    witness['pooled_long_correct_strictly_increased'] = sum(
        v[MODES[0]]['selected_five_or_more_correct'] for v in totals.values()) > old_long
    status = ('label_assisted_witness_meets_all_five_not_model_PASS' if all(witness.values()) else
              'individual_upper_bound_rules_out_all_five_on_this_fixed_pool' if not all(individual.values()) else
              'individual_bounds_sufficient_but_joint_feasibility_not_established')
    return dict(individual_upper_bound_checks=individual, maximum_correct_witness_checks=witness, status=status)


def main():
    source_jobs()
    prepared = R.completed(R.READOUT, JOBS['readout']); checked = R.completed(R.CHECK, JOBS['CPU'])
    validate_sources(prepared, checked)
    ids = list(dict.fromkeys(d['window_id'] for d in R.C.R.rows(R.Q.DRAWS)))
    totals = {s: {a: Counter() for a in ('A31_A3', *MODES)} for s in R.C.R.ALLOW}
    refs = {s: Counter(eligible=0, five_or_more=0) for s in totals}
    seen = []; candidates = 0
    with (OUT/'windows.jsonl').open('x') as handle:
        for free, lab in zip_longest(R.C.R.rows(R.READOUT/'natural_scores.jsonl'),
                                    R.C.R.rows(R.READOUT/'windows.jsonl')):
            if (free is None or lab is None or free['window'] != lab['window'] or
                    free['reference_used_in_scoring'] is not False or free['window']['split'] != 'train' or
                    free['window']['seqid'] not in R.C.R.ALLOW[free['window']['species']]):
                raise ValueError('Complete ordered allowed original TRAIN identity required')
            old = free['variants']['A31_A3']; labels = lab['labels']['A31_A3']
            counts = R.selected_counts(old, labels)
            if counts != lab['counts']['A31_A3']:
                raise ValueError('Actual original saved score/DP selection counts differ')
            oracle = bounds(old['chains'], labels); species = free['window']['species']
            totals[species]['A31_A3'].update(counts)
            for mode, value in oracle.items(): totals[species][mode].update(value['counts'])
            refs[species].update(eligible=lab['eligible_references'], five_or_more=lab['eligible_five_or_more'])
            seen.append(free['window']['id']); candidates += len(labels)
            handle.write(json.dumps(dict(window=free['window'], label_assisted=True,
                baseline_counts=counts, oracle=oracle), allow_nan=False)+'\n')
    totals = {s: {a: dict(v) for a, v in arms.items()} for s, arms in totals.items()}
    refs = {s: dict(v) for s, v in refs.items()}
    if (len(ids) != 1491 or seen != ids or candidates != 16321 or
            refs != prepared['eligible_reference_pairs'] or sum(v['eligible'] for v in refs.values()) != 2916 or
            sum(v['five_or_more'] for v in refs.values()) != 1086 or
            any(v['A31_A3'] != prepared['selected_counts_per_species'][s]['A31_A3'] for s, v in totals.items()) or
            [sum(v['A31_A3'][k] for v in totals.values()) for k in
             ('pool_correct', 'pool_known_nonmatch', 'pool_unknown')] != [2380, 12566, 1375]):
        raise ValueError('Original full counts, labels and reference denominators must reproduce')
    R.C.R.write_summary(OUT, dict(completed=True, experiment='M28-A31-NATURAL-TRAIN-SELECTOR-BOUND',
        protocol=PROTOCOL, source_jobs=JOBS, source_configuration=R.CONFIG,
        windows=1491, natural_candidates=candidates, selected_counts_per_species=totals,
        eligible_reference_pairs=refs, label_assisted_interpretation=interpretation(totals),
        grain='original_TRAIN_window_chain_pair_not_unique_gene',
        label_assisted_scores=True, original_scores_or_labels_changed=False,
        original_TRAIN_gates_unchanged=True, optimizer_steps_this_run=0,
        model_or_backbone_forward=False, feature_or_checkpoint_opened=False,
        new_candidate_search=False, new_R4_labelling=False, TRAIN_catalog_or_reference_opened=False,
        A40_outputs_opened=False, DEV_used=False, test_setaria_used=False, A33_outputs_opened=False,
        limitation='Exploratory label-assisted arbitrary-score bounds on the unchanged natural TRAIN pool. '
            'Separate maxima need not jointly occur; a missing witness is not infeasibility. '
            'No attainable neural-score, actual model PASS, raw-FASTA gene F1, biological truth or generalization claim.'))


if __name__ == '__main__':
    main()
