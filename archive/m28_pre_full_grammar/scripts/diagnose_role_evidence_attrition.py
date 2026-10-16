#!/usr/bin/env python3
"""Closed A40 TRAIN transitions at unchanged gains, whether its gate PASSes or FAILs."""
from collections import Counter, defaultdict
from itertools import zip_longest
import infer_role_evidence_context as I
from diagnose_structured_emission_attrition import transitions

D, R = I.D, I.D.R
OUT = I.ROOT/'outputs/M28-A40-ROLE-EVIDENCE-CONTEXT-ATTRITION'
COMPARISONS = {'A31_A3_to_A40': 'A31_A3', 'A39_to_A40': 'A39'}


def validate_completed(fitted, report, checked):
    D.validate_fit(fitted)
    for value in (report, checked):
        if (value['completed'] is not True or value['configuration'] != D.CONFIG or
                value['protocol'] != D.PROTOCOL or value['windows'] != 1491 or
                value['natural_candidates'] != 16321 or value['local_logit_receipts'] != 1491 or
                value['reference_injected_predictions'] != 0 or value['optimizer_steps_this_run'] != 0 or
                any(value[k] is not False for k in ('new_candidate_search', 'new_R4_labelling',
                    'DEV_used', 'test_setaria_used', 'A33_outputs_opened'))):
            raise ValueError('Complete unchanged TRAIN-only A40 result required before score access')
    if (report['source_fit_job'] != I.FIT_JOB or report['checkpoint'] != fitted['final_checkpoint'] or
            report['source_baseline_jobs'] != dict(readout='13398521', CPU='13398524') or
            report['source_previous_arm_jobs'] != dict(fit='13415276', readout='13415309', CPU='13415319') or
            checked['source_jobs'] != dict(fit=I.FIT_JOB, readout=I.READOUT_JOB) or
            any(checked[k] is not True for k in I.CPU_PROOFS)):
        raise ValueError('Actual A40 source jobs and seven CPU reconciliations required before score access')
    counts, refs = report['selected_counts_per_species'], report['eligible_reference_pairs']
    if (set(counts) != set(I.C.R.ALLOW) or set(refs) != set(I.C.R.ALLOW) or
            checked['selected_counts_per_species'] != counts or checked['eligible_reference_pairs'] != refs or
            sum(v['eligible'] for v in refs.values()) != 2916 or
            sum(v['five_or_more'] for v in refs.values()) != 1086):
        raise ValueError('Both species and original complete paired counts and denominators required')
    gates = D.gate(counts)
    if (report['predeclared_TRAIN_gate'] != gates or checked['predeclared_TRAIN_gate'] != gates or
            report['TRAIN_preparation_supported'] is not all(gates.values()) or
            checked['TRAIN_preparation_supported'] is not all(gates.values())):
        raise ValueError('Actual original five-condition arithmetic required, whether PASS or FAIL')
    return gates


def census(saved, previous, labelled, prior):
    groups = {s: {k: defaultdict(Counter) for k in COMPARISONS} for s in I.C.R.ALLOW}
    counts = {s: {a: Counter() for a in ('A31_A3', 'A39', 'A40')} for s in groups}
    refs = {s: Counter(eligible=0, five_or_more=0) for s in groups}
    windows = candidates = 0
    for new, old, lab, last in zip_longest(saved, previous, labelled, prior):
        if any(v is None for v in (new, old, lab, last)):
            raise ValueError('Complete four-way ordered original windows required')
        baseline, labels = old['variants']['A31_A3'], lab['labels']
        if (any(v['window'] != new['window'] for v in (old, lab, last)) or
                any(v['chains'] != new['chains'] for v in (baseline, last)) or
                any(v['reference_used_in_scoring'] is not False for v in (new, old, last)) or
                new['window']['split'] != 'train'):
            raise ValueError('Original ordered reference-free same-natural TRAIN identity differs')
        species = new['window']['species']
        arms = dict(A31_A3=baseline, A39=last, A40=new)
        actual = {a: R.selected_counts(row, labels) for a, row in arms.items()}
        if actual != lab['counts']:
            raise ValueError('Actual original three-arm per-window selections differ')
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
    D.F.A.A.require_jobs([I.FIT_JOB, I.READOUT_JOB, I.CHECK_JOB, '13398521', '13415319'])
    fitted = R.completed(D.F.OUT, I.FIT_JOB)
    report = R.completed(D.OUT, I.READOUT_JOB); checked = R.completed(D.CHECK, I.CHECK_JOB)
    gates = validate_completed(fitted, report, checked)
    baseline = R.completed(R.READOUT, '13398521')
    previous = R.completed(D.F.A.M.D.CHECK, '13415319')
    result = census(R.C.R.rows(D.OUT/'natural_scores.jsonl'), R.C.R.rows(R.READOUT/'natural_scores.jsonl'),
        R.C.R.rows(D.OUT/'windows.jsonl'), R.C.R.rows(D.F.A.M.D.OUT/'natural_scores.jsonl'))
    counts = result['selected_counts_per_species']
    if (result['windows'] != 1491 or result['natural_candidates'] != 16321 or
            counts != report['selected_counts_per_species'] or
            result['eligible_reference_pairs'] != report['eligible_reference_pairs'] or
            any(counts[s]['A31_A3'] != baseline['selected_counts_per_species'][s]['A31_A3'] or
                counts[s]['A39'] != previous['selected_counts_per_species'][s]['A39'] for s in counts)):
        raise ValueError('Complete transition census must reproduce three closed arms and original denominators')
    R.C.R.write_summary(OUT, dict(completed=True, experiment='M28-A40-ROLE-EVIDENCE-CONTEXT-ATTRITION',
        protocol=D.PROTOCOL, configuration=D.CONFIG,
        source_jobs=dict(fit=I.FIT_JOB, readout=I.READOUT_JOB, CPU=I.CHECK_JOB,
            baseline_readout='13398521', previous_readout='13415309', previous_CPU='13415319'),
        scientific_status=('A40_TRAIN_gate_PASS_not_raw_FASTA_or_generalization_evidence' if all(gates.values())
            else 'A40_TRAIN_FAIL_closed_no_rescue'),
        grain='original_unique_TRAIN_window_natural_chain_pair_not_unique_gene',
        **result, predeclared_TRAIN_gate=gates,
        transition_logic='unchanged_diagnose_structured_emission_attrition.transitions',
        model_or_backbone_forward=False, model_or_feature_tensor_opened=False, optimizer_steps_this_run=0,
        new_candidate_search=False, new_R4_labelling=False, TRAIN_catalog_or_reference_opened=False,
        DEV_used=False, test_setaria_used=False, A33_outputs_opened=False,
        limitation='Descriptive full saved-score transitions against stronger A31_A3 and closed A39, '
            'not role-feature causality, matched capacity/exposure, threshold rescue, gene F1 or generalization. '
            'Protected unknown remains separate; reference-relative nonmatch is not biological false positive.'))


if __name__ == '__main__':
    main()
