#!/usr/bin/env python3
"""Same-score scalar objective diagnostics on closed A37/A38 natural TRAIN."""
import argparse
from collections import Counter, defaultdict
from itertools import zip_longest
import json
import statistics
from pathlib import Path

import torch
from src.m28.chain_set import chain_set_loss
import diagnose_logistic_context_competition as C
import infer_context_chain as I37

D, R, I = C.D, C.R, C.I
OUT = I.ROOT/'outputs/M28-A38-SAVED-LOSS-GRADIENT'
COMPETITION_JOB = '13414782'


def source_reports():
    D.F.A.require_jobs([I.FIT_JOB, I.READOUT_JOB, I.CHECK_JOB, C.ATTRITION_JOB,
                       COMPETITION_JOB, '13410061', I37.READOUT_JOB, I37.CHECK_JOB, '13398521'])
    fitted = R.completed(D.F.OUT, I.FIT_JOB)
    report = R.completed(D.OUT, I.READOUT_JOB)
    checked = R.completed(D.CHECK, I.CHECK_JOB)
    gates = C.A.validate_completed(fitted, report, checked)
    competition = R.completed(C.OUT, COMPETITION_JOB)
    if (all(gates.values()) or competition['scientific_status'] != 'A38_TRAIN_FAIL_closed_no_rescue' or
            competition['source_jobs'] != dict(fit=I.FIT_JOB, readout=I.READOUT_JOB, CPU=I.CHECK_JOB,
                attrition=C.ATTRITION_JOB, baseline_readout='13398521') or
            competition['configuration'] != D.CONFIG or competition['protocol'] != D.PROTOCOL or
            competition['all_individual_and_joint_actual_selected_exchanges_verified'] is not True or
            competition['eligible_reference_pairs'] != report['eligible_reference_pairs']):
        raise ValueError('Closed, actually reconciled A38 competition required')
    fitted37 = R.completed(I37.D.F.OUT, '13410061'); I37.D.validate_fit(fitted37)
    report37 = R.completed(I37.D.OUT, I37.READOUT_JOB)
    checked37 = R.completed(I37.D.CHECK, I37.CHECK_JOB)
    if (report37['configuration'] != I37.D.CONFIG or checked37['configuration'] != I37.D.CONFIG or
            report37['protocol'] != I37.D.PROTOCOL or checked37['protocol'] != I37.D.PROTOCOL or
            checked37['source_jobs'] != dict(fit='13410061', readout=I37.READOUT_JOB) or
            any(checked37[k] is not True for k in I37.CPU_PROOFS) or
            checked37['selected_counts_per_species'] != report37['selected_counts_per_species'] or
            checked37['eligible_reference_pairs'] != report37['eligible_reference_pairs'] or
            report37['TRAIN_preparation_supported'] is not False or
            checked37['TRAIN_preparation_supported'] is not False or
            I37.D.gate(report37['selected_counts_per_species']) != report37['predeclared_TRAIN_gate'] or
            checked37['predeclared_TRAIN_gate'] != report37['predeclared_TRAIN_gate']):
        raise ValueError('Actual closed numerically reconciled A37 required')
    R.completed(R.READOUT, '13398521')
    return report, competition, report37


def scalar_values(row, labels):
    gain = torch.tensor(row['scores'], dtype=torch.float64, requires_grad=True)
    result = chain_set_loss(dict(logits=gain, null=gain.new_tensor(0.)), row['chains'], labels)
    q = torch.autograd.grad(result['log_partition'], gain)[0].detach()
    p = gain.detach().sigmoid()
    if (not bool(torch.isfinite(q).all()) or bool((q < -1e-10).any()) or
            any(float(q[i]) > float(p[i])+1e-10 for i, y in enumerate(labels) if y >= 0) or
            any(float(q[i]) != 0. for i, y in enumerate(labels) if y == -1)):
        raise ValueError('Diagnostic float64 occupancy bound or protected direct gradient differs')
    n = result['known_candidates']
    values = [None if y == -1 else dict(p=float(p[i]), q=float(q[i]),
                BCE=(float(p[i])-y)/n, set_NLL=(float(q[i])-y)/n,
                composite=((float(p[i])+float(q[i]))/2-y)/n, label=y)
              for i, y in enumerate(labels)]
    shift = sum(float(p[i]-q[i]) for i, y in enumerate(labels) if y >= 0)/(2*n) if n else 0.
    return values, shift


def summarize(values):
    result = dict(pairs=len(values))
    for key in ('p', 'q', 'BCE', 'set_NLL', 'composite'):
        result[key+'_mean'] = statistics.fmean(v[key] for v in values)
    if all(v['label'] == 0 for v in values):
        ratios = [v['q']/v['p'] for v in values if v['p'] > 0.]
        result.update(set_over_BCE_direct_gradient_ratio_median=statistics.median(ratios) if ratios else None,
                      ratio_pairs=len(ratios), zero_BCE_probability_pairs=len(values)-len(ratios))
    return result


def main(out):
    report, competition, report37 = source_reports()
    species_names = tuple(report['selected_counts_per_species'])
    counts = {s: {a: Counter() for a in ('A31_A3', 'A37', 'A38')} for s in species_names}
    refs = {s: Counter(eligible=0, five_or_more=0) for s in species_names}
    groups = defaultdict(list); shifts = defaultdict(list); cases = {s: [] for s in species_names}
    unique_blockers = {k: set() for k in ('all_known0', 'baseline_lost_known0')}
    windows = candidates = 0
    sources = (R.C.R.rows(D.OUT/'natural_scores.jsonl'), R.C.R.rows(R.READOUT/'natural_scores.jsonl'),
               R.C.R.rows(D.OUT/'windows.jsonl'), R.C.R.rows(I37.D.OUT/'natural_scores.jsonl'),
               R.C.R.rows(C.OUT/'windows.jsonl'))
    with (out/'windows.jsonl').open('x') as handle:
        for new, old, lab, ctx, blockers in zip_longest(*sources):
            if any(v is None for v in (new, old, lab, ctx, blockers)):
                raise ValueError('Complete five-way original TRAIN records required')
            baseline, labels = old['variants']['A31_A3'], lab['labels']
            if (any(v['window'] != new['window'] for v in (old, lab, ctx, blockers)) or
                    new['chains'] != baseline['chains'] or new['chains'] != ctx['chains'] or
                    any(v['reference_used_in_scoring'] is not False for v in (new, old, ctx)) or
                    new['window']['split'] != 'train'):
                raise ValueError('Original ordered same-pool TRAIN identity differs')
            species = new['window']['species']
            arms = dict(A31_A3=baseline, A37=ctx, A38=new)
            actual = {a: R.selected_counts(row, labels) for a, row in arms.items()}
            if actual != lab['counts'] or any(blockers['counts'][a] != actual[a] for a in ('A31_A3', 'A38')):
                raise ValueError('Actual original per-window selections differ')
            for a, values in actual.items():
                counts[species][a].update(values)
            refs[species].update(eligible=lab['eligible_references'], five_or_more=lab['eligible_five_or_more'])
            probabilities = {}; local = {}
            for arm, row in (('A37', ctx), ('A38', new)):
                values, shift = scalar_values(row, labels); local[arm] = values
                probabilities[arm] = dict(known_occupancy=[None if v is None else v['q'] for v in values],
                    composite_minus_BCE_shared_null_scalar_derivative=shift)
                shifts[(species, arm)].append(shift)
                for i, value in enumerate(values):
                    if value is None:
                        continue
                    key = 'known_positive' if labels[i] == 1 else 'known_nonmatch'
                    groups[(species, arm, key)].append(value)
                    if labels[i] == 0 and i in row['selected']:
                        groups[(species, arm, 'selected_known_nonmatch')].append(value)
            for case in blockers['cases']:
                i = case['candidate_index']; chain = new['chains'][i]
                expected = {b for b in new['selected'] if chain[0][0] < new['chains'][b][-1][1]
                            and new['chains'][b][0][0] < chain[-1][1]}
                if (labels[i] != 1 or new['scores'][i] <= 0. or i in new['selected'] or
                        chain != case['chain'] or new['scores'][i] != case['gain'] or
                        case['baseline_selected'] != (i in baseline['selected']) or
                        set(case['selected_conflict_indices']) != expected or
                        len(case['selected_conflict_indices']) != len(expected)):
                    raise ValueError('Actual positive-unselected blocker identity differs')
                for b in case['selected_conflict_indices']:
                    if labels[b] != 0:
                        continue
                    unique_blockers['all_known0'].add((new['window']['id'], b))
                    if case['baseline_selected']:
                        unique_blockers['baseline_lost_known0'].add((new['window']['id'], b))
                    for arm in ('A37', 'A38'):
                        groups[(species, arm, 'A38_actual_known0_blocker_memberships')].append(local[arm][b])
                        if case['baseline_selected']:
                            groups[(species, arm, 'A38_baseline_lost_known0_blocker_memberships')].append(local[arm][b])
            cases[species].extend(blockers['cases'])
            handle.write(json.dumps(dict(window=new['window'], known=labels.count(0)+labels.count(1),
                protected_unknown=labels.count(-1), diagnostic=probabilities), allow_nan=False)+'\n')
            windows += 1; candidates += len(labels)
    actual = {s: {a: dict(v) for a, v in values.items()} for s, values in counts.items()}
    references = {s: dict(v) for s, v in refs.items()}
    profiles = {s: C.profiles(v) for s, v in cases.items()}
    pooled_profiles = C.profiles([case for values in cases.values() for case in values])
    if (windows != 1491 or candidates != 16321 or actual != report['selected_counts_per_species'] or
            references != report['eligible_reference_pairs'] or profiles != competition['per_species'] or
            pooled_profiles != competition['pooled'] or
            any(actual[s]['A37'] != report37['selected_counts_per_species'][s]['A37'] for s in species_names)):
        raise ValueError('Complete original population, actual selections, blocker strata and denominators differ')
    pooled = defaultdict(list)
    for (species, arm, group), values in groups.items():
        pooled[(arm, group)].extend(values)
    summary = dict(completed=True, experiment='M28-A38-SAVED-LOSS-GRADIENT',
        source_jobs=dict(A38_fit=I.FIT_JOB, A38_readout=I.READOUT_JOB, A38_CPU=I.CHECK_JOB,
            A37_fit='13410061', A37_readout=I37.READOUT_JOB, A37_CPU=I37.CHECK_JOB,
            competition=COMPETITION_JOB, baseline_readout='13398521'),
        scientific_status='closed_A37_A38_descriptive_same_score_counterfactual_not_A39_result',
        windows=windows, natural_candidates=candidates, selected_counts_per_species=actual,
        eligible_reference_pairs=references, actual_blocker_profiles=pooled_profiles,
        unique_known0_blockers={k: len(v) for k, v in unique_blockers.items()},
        per_species={s: {a: {g: summarize(v) for (sp, ar, g), v in groups.items() if sp == s and ar == a}
                        for a in ('A37', 'A38')} for s in species_names},
        pooled={a: {g: summarize(v) for (ar, g), v in pooled.items() if ar == a} for a in ('A37', 'A38')},
        per_species_window_mean_null_derivative_shift={s: {a: statistics.fmean(shifts[(s, a)])
            for a in ('A37', 'A38')} for s in species_names},
        diagnostic_arithmetic='float64_scalar_autograd_original_known_partition_at_saved_FP32_gains',
        diagnostic_occupancy_bound_absolute_tolerance=1e-10, clipping_or_probability_floor=False,
        draw_importance_or_repetition_applied=False, original_loss_implementations_unchanged=True,
        full_population_labels_DP_references_and_actual_blockers_reconciled=True,
        model_or_backbone_forward=False, model_or_feature_tensor_opened=False,
        optimizer_steps_this_run=0, new_candidate_search=False, new_R4_labelling=False,
        TRAIN_catalog_or_reference_opened=False, DEV_used=False, test_setaria_used=False,
        A33_outputs_opened=False, A39_weights_or_scores_opened=False,
        limitation='Unique TRAIN window/candidate census; blocker memberships may repeat candidates. '
        'Counterfactual scalar objectives at fixed old scores, not neural gradients, AdamW steps, '
        'causal failure explanation, biological negatives, full gene F1, generalization or threshold rescue.')
    R.C.R.write_summary(out, summary)
    print(json.dumps(dict(windows=windows, natural_candidates=candidates,
                         unique_known0_blockers=summary['unique_known0_blockers']), allow_nan=False), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(); torch.set_num_threads(2)
    main(args.output_dir)
