#!/usr/bin/env python3
"""Count new A42 selected geometry relative to the closed original natural pool.

Read only the same three saved JSON files as the closed A42 competition
diagnosis, not A12 injected gold, references, features, DNA or models.
"""
from collections import Counter
from itertools import zip_longest

import diagnose_emission_lattice_competition as A42

OUT = A42.D.ROOT/'outputs/M28-A43-NATURAL-COMPETITOR-SUPPORT-CENSUS'
PROTOCOL = 'reports/M28-METHOD-RESTART/a43_natural_competitor_preparation.md'


def census(original, natural, labelled):
    names = {-1: 'unknown', 0: 'known_nonmatch', 1: 'correct'}
    groups = {s: {name: Counter(present=0, new=0) for name in names.values()}
              for s in A42.D.R.C.R.ALLOW}
    old_counts = {s: Counter() for s in groups}
    refs = {s: Counter(eligible=0, five_or_more=0) for s in groups}
    windows = original_candidates = selected = 0
    for prior, new, lab in zip_longest(original, natural, labelled):
        if any(value is None for value in (prior, new, lab)):
            raise ValueError('Complete ordered three-way saved population required')
        w = new['window']; old = prior['variants']['A31_A3']
        a, b = lab['labels']['A31_A3'], lab['labels']['A42']
        if (w != prior['window'] or w != lab['window'] or w['split'] != 'train' or
                new['original_ordinal'] != windows or new['reference_used_in_scoring'] is not False or
                prior['reference_used_in_scoring'] is not False or
                new['selected'] != list(range(len(new['chains']))) or
                len(old['chains']) != len(a) or len(new['chains']) != len(b) or
                any(y not in names for y in (*a, *b))):
            raise ValueError('Original ordered reference-free geometry or labels differ')
        old_keys = [A42.key(c) for c in old['chains']]
        new_keys = [A42.key(c) for c in new['chains']]
        if len(set(old_keys)) != len(old_keys) or len(set(new_keys)) != len(new_keys):
            raise ValueError('Closed natural geometry must be unique')
        old_labels = dict(zip(old_keys, a)); s = w['species']
        old_counts[s].update(names[y] for y in a)
        refs[s].update(eligible=lab['eligible_references'], five_or_more=lab['eligible_five_or_more'])
        for geometry, label in zip(new_keys, b):
            present = geometry in old_labels
            if present and old_labels[geometry] != label:
                raise ValueError('Common geometry frozen labels differ')
            groups[s][names[label]]['present' if present else 'new'] += 1
        windows += 1; original_candidates += len(old_keys); selected += len(new_keys)
    pooled = {name: Counter() for name in names.values()}
    original_pooled = Counter()
    for s in groups:
        original_pooled.update(old_counts[s])
        for name in pooled: pooled[name].update(groups[s][name])
    return dict(windows=windows, original_natural_candidates=original_candidates,
        A42_natural_selected_candidates=selected,
        original_natural_labels_per_species={s: dict(v) for s, v in old_counts.items()},
        original_natural_labels_pooled=dict(original_pooled),
        eligible_reference_pairs={s: dict(v) for s, v in refs.items()},
        per_species={s: {name: dict(v) for name, v in categories.items()} for s, categories in groups.items()},
        pooled={name: dict(v) for name, v in pooled.items()})


def main():
    D = A42.D
    D.B.F.A.require_jobs([A42.JOBS['readout'], A42.JOBS['CPU'], '13398521', '13398524'])
    report = D.R.completed(A42.READOUT, A42.JOBS['readout'])
    checked = D.R.completed(A42.CHECK, A42.JOBS['CPU'])
    gates = A42.validate_closed(report, checked)
    result = census(D.R.C.R.rows(D.R.READOUT/'natural_scores.jsonl'),
        D.R.C.R.rows(A42.READOUT/'natural_scores.jsonl'), D.R.C.R.rows(A42.READOUT/'windows.jsonl'))
    if (result['windows'] != 1491 or result['original_natural_candidates'] != 16321 or
            result['A42_natural_selected_candidates'] != 3630 or
            result['original_natural_labels_pooled'] != dict(correct=2380, known_nonmatch=12566, unknown=1375) or
            result['eligible_reference_pairs'] != report['eligible_reference_pairs']):
        raise ValueError('Full original natural support and reference denominators differ')
    for s, values in result['per_species'].items():
        for name, count_key in (('correct', 'selected_correct'), ('known_nonmatch', 'selected_known_nonmatch'),
                                ('unknown', 'selected_unknown')):
            if sum(values[name].values()) != report['selected_counts_per_species'][s]['A42'][count_key]:
                raise ValueError('Actual A42 selected label counts differ')
    D.R.C.R.write_summary(OUT, dict(completed=True, experiment='M28-A43-NATURAL-COMPETITOR-SUPPORT-CENSUS',
        protocol=PROTOCOL, source_A42_protocol=D.PROTOCOL, source_A42_configuration=D.CONFIG,
        source_jobs={**A42.JOBS, 'original_readout': '13398521', 'original_CPU': '13398524'}, **result,
        actual_original_ordered_support_and_saved_labels_reconciled=True,
        actual_A42_selected_labels_and_original_denominators_reconciled=True,
        new_known_competitor_geometry_present=result['pooled']['known_nonmatch']['new'] > 0,
        original_A12_injected_gold_opened=False, model_or_feature_tensor_opened=False,
        model_or_backbone_forward=False, optimizer_steps_this_run=0, new_candidate_search=False,
        new_geometry_labels_computed=False, TRAIN_catalog_or_reference_opened=False,
        DEV_used=False, independent_species_used=False, test_setaria_used=False, A33_outputs_opened=False,
        source_predeclared_TRAIN_gate=gates, source_TRAIN_preparation_supported=all(gates.values()),
        grain='original_unique_TRAIN_window_chain_pair_not_unique_gene',
        limitation='Finite source-pool membership only; not A12 gold/union support verification, '
            'training authorization, model improvement, training causality, raw annotation or generalization.'))


if __name__ == '__main__':
    main()
