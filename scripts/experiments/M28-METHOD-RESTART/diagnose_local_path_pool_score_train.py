#!/usr/bin/env python3
"""Three saved A3 TRAIN cells; diagnostic only, never a new selection gate."""
from collections import Counter
import json
import torch
import run_local_path_caller as R

OUT = R.ROOT/'outputs/M28-A31-SCORE-POOL-TRAIN-DIAGNOSTIC'
PROTOCOL = 'reports/M28-METHOD-RESTART/a31_score_pool_train_diagnostic_protocol.md'
ARMS = ('greedy_query', 'greedy_local', 'union_local')


def three_cells(natural, labelled):
    greedy = natural['variants']['A3_greedy']; union = natural['variants']['A31_A3']
    labels = labelled['labels']
    index = {tuple(map(tuple, c)): i for i, c in enumerate(union['chains'])}
    try:
        positions = [index[tuple(map(tuple, c))] for c in greedy['chains']]
    except KeyError as exc:
        raise ValueError('Greedy geometry is not a subset of the verified natural union') from exc
    subset_labels = [labels['A31_A3'][i] for i in positions]
    if subset_labels != labels['A3_greedy']:
        raise ValueError('Same greedy geometry must retain the same saved R4 labels')
    scores = [union['scores'][i] for i in positions]
    local = dict(chains=greedy['chains'], scores=scores,
                 selected=R.select_nonoverlapping(greedy['chains'], scores))
    cells = dict(greedy_query=greedy, greedy_local=local, union_local=union)
    targets = dict(greedy_query=labels['A3_greedy'], greedy_local=subset_labels, union_local=labels['A31_A3'])
    counts = {a: R.selected_counts(cells[a], targets[a]) for a in ARMS}
    return counts, dict(chains=local['chains'], scores=scores, selected=local['selected'], labels=subset_labels)


def contrast(old, new, old_metric, new_metric):
    # Rate differences telescope along this ordered path, not causal shares.
    keys = ('selected_correct', 'selected_known_nonmatch', 'selected_unknown', 'selected_five_or_more_correct')
    return dict(count_changes={k: new[k]-old[k] for k in keys},
                proxy_change_percentage_points={k: 100*(new_metric[k]-old_metric[k]) for k in ('precision', 'recall', 'F1')})


def main():
    prepared = R.completed(R.READOUT, 13398521); checked = R.completed(R.CHECK, 13398524)
    if (prepared['configuration'] != R.CONFIG or checked['configuration'] != R.CONFIG or
            prepared['protocol'] != R.PROTOCOL or checked['protocol'] != R.PROTOCOL or
            prepared['windows'] != 1491 or checked['windows'] != 1491 or
            checked['source_jobs'] != {'readout': '13398521'} or
            not checked['actual_saved_natural_pool_and_greedy_subset_identity_verified'] or
            not checked['actual_R4_labels_selection_reference_denominators_and_proxy_metrics_recomputed'] or
            prepared['selected_counts_per_species'] != checked['selected_counts_per_species'] or
            prepared['eligible_reference_pairs'] != checked['eligible_reference_pairs']):
        raise ValueError('Original closed TRAIN readout and actual CPU reconciliation required')
    # Neither variant's PASS/FAIL is an execution condition for this diagnostic.
    natural = R.C.R.rows(R.READOUT/'natural_scores.jsonl')
    labelled = R.C.R.rows(R.READOUT/'windows.jsonl')
    ids = list(dict.fromkeys(d['window_id'] for d in R.C.R.rows(R.Q.DRAWS)))
    if (len(ids) != 1491 or [r['window']['id'] for r in natural] != ids or
            [r['window']['id'] for r in labelled] != ids):
        raise ValueError('Complete ordered original TRAIN window population required')
    totals = {s: {a: Counter() for a in ARMS} for s in R.C.R.ALLOW}
    references = {s: Counter(eligible=0, five_or_more=0) for s in R.C.R.ALLOW}
    with (OUT/'greedy_local.jsonl').open('x') as handle:
        for free, saved in zip(natural, labelled):
            if any(free['window'][k] != saved['window'][k] for k in R.C.R.GEOMETRY):
                raise ValueError('Saved paired window geometry differs')
            counts, local = three_cells(free, saved); species = free['window']['species']
            for a in ARMS:
                totals[species][a].update(counts[a])
            references[species].update(eligible=saved['eligible_references'], five_or_more=saved['eligible_five_or_more'])
            handle.write(json.dumps(dict(window=free['window'], greedy_local=local, counts=counts), allow_nan=False)+'\n')
    totals = {s: {a: dict(c) for a, c in v.items()} for s, v in totals.items()}
    references = {s: dict(c) for s, c in references.items()}
    if references != prepared['eligible_reference_pairs']:
        raise ValueError('Original complete 2916 window/reference denominator differs')
    for s in totals:
        if (totals[s]['greedy_query'] != prepared['selected_counts_per_species'][s]['A3_greedy'] or
                totals[s]['union_local'] != prepared['selected_counts_per_species'][s]['A31_A3']):
            raise ValueError('The two original cells must reproduce the source counts exactly')
    metrics = R.proxies(totals, references)
    effects = {s: dict(
        score_change_within_greedy=contrast(totals[s]['greedy_query'], totals[s]['greedy_local'],
                                           metrics[s]['greedy_query'], metrics[s]['greedy_local']),
        pool_expansion_at_local_score=contrast(totals[s]['greedy_local'], totals[s]['union_local'],
                                              metrics[s]['greedy_local'], metrics[s]['union_local'])) for s in totals}
    R.C.R.write_summary(OUT, dict(completed=True, experiment='A31_score_pool_TRAIN_diagnostic', protocol=PROTOCOL,
        source_jobs=dict(readout='13398521', CPU='13398524'), windows=1491,
        generator_checkpoint=R.CONFIG['generators']['A3'], threshold=0.,
        selected_counts_per_species=totals, eligible_reference_pairs=references,
        known_label_window_pair_proxies=metrics, ordered_conditional_changes=effects,
        source_two_cells_and_reference_denominators_reproduced=True,
        original_TRAIN_gates_unchanged=True, diagnostic_used_for_model_selection=False,
        original_TRAIN_preparation_supported=prepared['TRAIN_preparation_supported'],
        optimizer_steps_this_run=0, new_parameters=0, model_or_backbone_forward=False,
        new_candidate_search=False, new_R4_labelling=False, feature_cache_used=False,
        DEV_used=False, test_setaria_used=False,
        limitation='Same frozen A3, in-sample window/chain diagnostic. Greedy score contrast and local-score pool '
            'contrast are conditional; missing union-query cell prevents interaction estimation. '
            'No biological causality, independent generalization, gene F1 or new continuation gate.'))


if __name__ == '__main__':
    torch.set_num_threads(2); main()
