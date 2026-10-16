#!/usr/bin/env python3
"""Complete A39 TRAIN transitions at fixed saved scores, whether PASS or FAIL."""
from collections import Counter, defaultdict
from itertools import zip_longest
import infer_null_anchored_set as I
from diagnose_structured_emission_attrition import transitions

D, R = I.D, I.D.R
OUT = I.ROOT/'outputs/M28-A39-NULL-ANCHORED-SET-ATTRITION'
COMPARISONS = {'A31_A3_to_A39': 'A31_A3', 'A37_to_A39': 'A37', 'A38_to_A39': 'A38'}


def validate_completed(fitted, report, checked):
    D.validate_fit(fitted)
    for value in (report, checked):
        if (value['completed'] is not True or value['configuration'] != D.CONFIG or
                value['protocol'] != D.PROTOCOL or value['windows'] != 1491 or
                value['natural_candidates'] != 16321 or value['reference_injected_predictions'] != 0 or
                value['optimizer_steps_this_run'] != 0 or any(value[k] is not False for k in
                    ('new_candidate_search', 'new_R4_labelling', 'DEV_used', 'test_setaria_used', 'A33_outputs_opened'))):
            raise ValueError('Complete unchanged TRAIN-only A39 result required before score access')
    if (report['source_fit_job'] != I.FIT_JOB or report['source_checks_job'] != I.CHECKS_JOB or
            report['source_preflight_job'] != D.F.PREFLIGHT_JOB or
            report['checkpoint'] != fitted['final_checkpoint'] or
            report['source_context_comparator_jobs'] != dict(fit='13410061', readout='13410739', CPU='13410742') or
            report['source_logistic_comparator_jobs'] != dict(fit='13412630', readout='13412701', CPU='13412702') or
            report['source_baseline_jobs'] != dict(readout='13398521', CPU='13398524') or
            checked['source_jobs'] != dict(fit=I.FIT_JOB, readout=I.READOUT_JOB) or
            any(checked[k] is not True for k in I.CPU_PROOFS)):
        raise ValueError('Actual A39 source jobs and five CPU reconciliations required before score access')
    if (any(report[k] is not True for k in ('natural_score_file_closed_before_saved_label_arrays_open',
            'original_natural_geometry_order_sources_unchanged', 'original_saved_labels_reused',
            'model_weights_unchanged', 'all_parameters_frozen_and_gradients_absent')) or
            report['label_catalog_opened'] is not False or report['new_backbone_forward'] is not False or
            checked['model_or_backbone_forward'] is not False or checked['independent_neural_forward_replay'] is not False):
        raise ValueError('Actual unchanged closed score provenance required')
    counts, refs = report['selected_counts_per_species'], report['eligible_reference_pairs']
    if (set(counts) != set(R.C.R.ALLOW) or set(refs) != set(R.C.R.ALLOW) or
            checked['selected_counts_per_species'] != counts or checked['eligible_reference_pairs'] != refs or
            sum(v['eligible'] for v in refs.values()) != 2916 or sum(v['five_or_more'] for v in refs.values()) != 1086):
        raise ValueError('Both species, paired counts and original full reference denominators required')
    gates = D.gate(counts)
    if (report['predeclared_TRAIN_gate'] != gates or checked['predeclared_TRAIN_gate'] != gates or
            report['TRAIN_preparation_supported'] is not all(gates.values()) or
            checked['TRAIN_preparation_supported'] is not all(gates.values())):
        raise ValueError('Actual original five-condition arithmetic required, whether PASS or FAIL')
    return gates


def census(saved, previous, labelled, context, logistic):
    groups = {s: {k: defaultdict(Counter) for k in COMPARISONS} for s in R.C.R.ALLOW}
    counts = {s: {a: Counter() for a in ('A31_A3', 'A37', 'A38', 'A39')} for s in groups}
    refs = {s: Counter(eligible=0, five_or_more=0) for s in groups}
    windows = candidates = 0
    for new, old, lab, prior_context, prior_logistic in zip_longest(saved, previous, labelled, context, logistic):
        if any(v is None for v in (new, old, lab, prior_context, prior_logistic)):
            raise ValueError('Complete five-way paired original windows required')
        baseline, labels = old['variants']['A31_A3'], lab['labels']
        if (any(v['window'] != new['window'] for v in (old, lab, prior_context, prior_logistic)) or
                any(v['chains'] != new['chains'] for v in (baseline, prior_context, prior_logistic)) or
                any(v['reference_used_in_scoring'] is not False for v in (new, old, prior_context, prior_logistic)) or
                new['window']['split'] != 'train'):
            raise ValueError('Original ordered reference-free same-natural-pool TRAIN identity differs')
        species = new['window']['species']
        arms = dict(A31_A3=baseline, A37=prior_context, A38=prior_logistic, A39=new)
        actual = {a: R.selected_counts(row, labels) for a, row in arms.items()}
        if actual != lab['counts']:
            raise ValueError('Actual original four-arm per-window selections differ')
        for arm, value in actual.items():
            counts[species][arm].update(value)
        refs[species].update(eligible=lab['eligible_references'], five_or_more=lab['eligible_five_or_more'])
        for comparison, source in COMPARISONS.items():
            for key, value in transitions(new['chains'], labels, arms[source], new).items():
                groups[species][comparison][key].update(value)
        windows += 1; candidates += len(labels)
    pooled = {k: defaultdict(Counter) for k in COMPARISONS}
    for species in groups.values():
        for comparison, values in species.items():
            for key, value in values.items():
                pooled[comparison][key].update(value)
    return dict(windows=windows, natural_candidates=candidates,
        selected_counts_per_species={s: {a: dict(c) for a, c in v.items()} for s, v in counts.items()},
        eligible_reference_pairs={s: dict(c) for s, c in refs.items()},
        per_species={s: {a: {k: dict(c) for k, c in v.items()} for a, v in species.items()}
            for s, species in groups.items()},
        pooled={a: {k: dict(c) for k, c in v.items()} for a, v in pooled.items()})


def main():
    D.F.A.require_jobs([I.FIT_JOB, I.READOUT_JOB, I.CHECK_JOB, '13398521', '13410742', '13412702'])
    fitted = R.completed(D.F.OUT, I.FIT_JOB)
    report = R.completed(D.OUT, I.READOUT_JOB); checked = R.completed(D.CHECK, I.CHECK_JOB)
    gates = validate_completed(fitted, report, checked)
    sources = dict(A31_A3=R.completed(R.READOUT, '13398521'),
                   A37=R.completed(D.L.F.A.D.CHECK, '13410742'), A38=R.completed(D.L.CHECK, '13412702'))
    result = census(R.C.R.rows(D.OUT/'natural_scores.jsonl'), R.C.R.rows(R.READOUT/'natural_scores.jsonl'),
        R.C.R.rows(D.OUT/'windows.jsonl'), R.C.R.rows(D.L.F.A.D.OUT/'natural_scores.jsonl'),
        R.C.R.rows(D.L.OUT/'natural_scores.jsonl'))
    counts = result['selected_counts_per_species']
    if (result['windows'] != 1491 or result['natural_candidates'] != 16321 or
            counts != report['selected_counts_per_species'] or
            result['eligible_reference_pairs'] != report['eligible_reference_pairs'] or
            any(counts[s][a] != source['selected_counts_per_species'][s][a]
                for s in counts for a, source in sources.items())):
        raise ValueError('Full transition census must reproduce four completed arms and original denominators')
    R.C.R.write_summary(OUT, dict(completed=True, experiment='M28-A39-NULL-ANCHORED-SET-ATTRITION',
        protocol=D.PROTOCOL, configuration=D.CONFIG,
        source_jobs=dict(fit=I.FIT_JOB, readout=I.READOUT_JOB, CPU=I.CHECK_JOB,
            baseline_readout='13398521', context_readout='13410739', context_CPU='13410742',
            logistic_readout='13412701', logistic_CPU='13412702'),
        scientific_status=('A39_TRAIN_gate_PASS_not_raw_FASTA_or_generalization_evidence' if all(gates.values())
            else 'A39_TRAIN_FAIL_closed_no_rescue'),
        grain='original_unique_TRAIN_window_natural_chain_pair_not_unique_gene',
        **result, predeclared_TRAIN_gate=gates,
        transition_logic='unchanged_diagnose_structured_emission_attrition.transitions',
        model_or_backbone_forward=False, model_or_feature_tensor_opened=False, optimizer_steps_this_run=0,
        new_candidate_search=False, new_R4_labelling=False, TRAIN_catalog_or_reference_opened=False,
        DEV_used=False, test_setaria_used=False, A33_outputs_opened=False,
        limitation='Descriptive complete saved-score transitions against the stronger same-pool baseline '
            'and matched-recipe A37/A38; not neural causality, threshold rescue, gene F1 or generalization. '
            'Protected unknown stays separate; reference-relative nonmatch is not biological false positive.'))


if __name__ == '__main__':
    main()
