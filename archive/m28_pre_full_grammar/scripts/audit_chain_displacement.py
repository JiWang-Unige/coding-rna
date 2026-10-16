#!/usr/bin/env python3
"""A22 read-only displacement audit of frozen natural TRAIN candidate sets."""
from collections import Counter
import json
import math
from pathlib import Path
import time

from src.m28.core import select_nonoverlapping

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT/'outputs/M28-A22-CHAIN-DISPLACEMENT-AUDIT'
PROTOCOL = 'reports/M28-METHOD-RESTART/a22_chain_displacement_protocol.md'
SOURCES = {
    'A12': ('M28-A12-QUALITY-TRAIN-DIAGNOSTIC', '13360650'),
    'A15': ('M28-A15-NATURAL-CHAIN-SET-TRAIN-READOUT-RECOVERY-R1', '13387947'),
    'A15_check': ('M28-A15-NATURAL-CHAIN-SET-RESULT-CHECK-RECOVERY-R1', '13387977'),
    'A21': ('M28-A21-JUNCTION-CHAIN-TRAIN-READOUT', '13390188'),
    'A21_check': ('M28-A21-JUNCTION-CHAIN-RESULT-CHECK', '13390189'),
}
ARMS = ('A6', 'A15', 'A21')
SPECIES = ('arabidopsis_thaliana', 'oryza_sativa')
IDENTITY = ('window_id', 'species', 'chain', 'origin', 'label', 'exon_class')
CATEGORIES = ('no_original_positive_overlap', 'no_natural_positive_overlap',
              'natural_positive_all_nonpositive', 'positive_gain_set_competition')
COUNT_KEYS = ('selected', 'selected_correct', 'selected_known_nonmatch',
              'selected_unknown', 'selected_five_or_more_correct')


def rows(path):
    with path.open() as handle:
        return [json.loads(line) for line in handle]


def overlap(a, b):
    return a[0][0] < b[-1][1] and b[0][0] < a[-1][1]


def require_compatible(chains):
    ordered = sorted(chains, key=lambda c: c[0][0])
    if any(overlap(a, b) for a, b in zip(ordered, ordered[1:])):
        raise ValueError('Positive or selected set is not span-compatible')


def exchange(local, selected, positives):
    """Actual removed set when inserting a compatible subset of missed gold."""
    require_compatible([local[i]['chain'] for i in positives])
    blockers = [i for i in selected if any(
        overlap(local[i]['chain'], local[p]['chain']) for p in positives)]
    inserted = sum(local[i]['gain'] for i in positives)
    removed = sum(local[i]['gain'] for i in blockers)
    if inserted > removed + 1e-6*max(1., abs(inserted), abs(removed)):
        raise ValueError('Saved natural selection violates a feasible gold-set exchange')
    return dict(positive_indices=positives, selected_conflict_indices=blockers,
        inserted_gain=inserted, removed_gain=removed,
        conflict_counts_by_label={str(y): sum(local[i]['label'] == y for i in blockers)
                                  for y in (1, 0, -1)},
        conflict_gain_by_label={str(y): sum(local[i]['gain'] for i in blockers if local[i]['label'] == y)
                                for y in (1, 0, -1)})


def audit_window(local, original):
    if any(r['origin'] != 'free' or not math.isfinite(r['gain']) for r in local):
        raise ValueError('Only finite unchanged natural candidates may enter selection')
    gold = [r for r in original if r['label'] == 1]
    natural_gold = [i for i, r in enumerate(local) if r['label'] == 1]
    require_compatible([r['chain'] for r in gold])
    selected = select_nonoverlapping([r['chain'] for r in local], [r['gain'] for r in local])
    require_compatible([local[i]['chain'] for i in selected])
    chosen = set(selected)
    if any(local[i]['gain'] <= 0 for i in selected):
        raise ValueError('Nonpositive chain selected at the frozen null operating point')
    bad = []
    for i in selected:
        r = local[i]
        if r['label'] != 0:
            continue
        available = [p for p in natural_gold if overlap(r['chain'], local[p]['chain'])]
        if not any(overlap(r['chain'], g['chain']) for g in gold):
            category = CATEGORIES[0]
        elif not available:
            category = CATEGORIES[1]
        elif all(local[p]['gain'] <= 0 for p in available):
            category = CATEGORIES[2]
        else:
            category = CATEGORIES[3]
        bad.append(dict(candidate_index=i, category=category,
                        overlapping_natural_positive_indices=available))
    missed = [i for i in natural_gold if i not in chosen]
    positive_missed = [i for i in missed if local[i]['gain'] > 0]
    individual = [exchange(local, selected, [i]) for i in positive_missed]
    joint = exchange(local, selected, positive_missed)
    counts = dict(selected=len(selected),
        selected_correct=sum(local[i]['label'] == 1 for i in selected),
        selected_known_nonmatch=len(bad),
        selected_unknown=sum(local[i]['label'] == -1 for i in selected),
        selected_five_or_more_correct=sum(local[i]['label'] == 1 and len(local[i]['chain']) >= 5 for i in selected))
    return dict(selected_indices=selected, counts=counts, selected_known_nonmatch=bad,
                missed_natural_positive_indices=missed,
                missed_nonpositive_natural_positive_indices=[i for i in missed if local[i]['gain'] <= 0],
                positive_gain_missed_exchanges=individual, joint_positive_gain_missed_exchange=joint)


def category_counts(bad):
    return {key: sum(r['category'] == key for r in bad) for key in CATEGORIES}


def aggregate(records):
    result = {s: {} for s in SPECIES}
    transitions = {s: {} for s in SPECIES}
    for species in SPECIES:
        subset = [r for r in records if r['window']['species'] == species]
        for arm in ARMS:
            counts, categories, missed, blockers = Counter(), Counter(), Counter(), Counter()
            for rec in subset:
                a = rec['arms'][arm]
                counts.update(a['counts'])
                categories.update(category_counts(a['selected_known_nonmatch']))
                missed.update(dict(nonpositive=len(a['missed_nonpositive_natural_positive_indices']),
                                   positive_gain_competition=len(a['positive_gain_missed_exchanges'])))
                for item in a['positive_gain_missed_exchanges']:
                    labels = [y for y in (1, 0, -1) if item['conflict_counts_by_label'][str(y)]]
                    blockers.update(['+'.join(map(str, labels))])
            result[species][arm] = dict(counts=dict(counts), selected_known_nonmatch_categories=dict(categories),
                missed_natural_positives=dict(missed), positive_gain_missed_blocker_label_combinations=dict(blockers))
        for baseline in ('A6', 'A15'):
            gained, dropped, common = [], [], 0
            for rec in subset:
                a = {r['candidate_index']: r for r in rec['arms']['A21']['selected_known_nonmatch']}
                b = {r['candidate_index']: r for r in rec['arms'][baseline]['selected_known_nonmatch']}
                gained.extend(a[i] for i in a.keys()-b.keys())
                dropped.extend(b[i] for i in b.keys()-a.keys())
                common += len(a.keys() & b.keys())
            transitions[species]['A21_vs_'+baseline] = dict(common=common, gained=len(gained), dropped=len(dropped),
                net=len(gained)-len(dropped), gained_categories=category_counts(gained), dropped_categories=category_counts(dropped))
    return result, transitions


def main():
    began = time.perf_counter()
    summaries, data = {}, {}
    for key, (name, job) in SOURCES.items():
        folder = ROOT/'outputs'/name
        if (folder/'STATUS').read_text().strip() != 'COMPLETED' or (folder/'JOBID').read_text().strip() != job:
            raise ValueError('Wrong or incomplete source '+key)
        summaries[key] = json.loads((folder/'summary.json').read_text())
        if not summaries[key]['completed']:
            raise ValueError('Source not completed '+key)
        if not key.endswith('_check'):
            data[key] = rows(folder/'windows.jsonl')
    if (not summaries['A21_check']['actual_last_head_shapes_finite_tensors_components_selection_profiles_and_gate_verified']
        or summaries['A21_check']['TRAIN_preparation_supported']
        or not summaries['A15_check']['actual_paired_selection_profiles_and_frozen_A6_gate_recomputed']):
        raise ValueError('Closed, verified failed A21 and paired A15 required')
    ids = [r['window']['id'] for r in data['A12']]
    if len(ids) != 1491 or len(set(ids)) != 1491 or any(
        [r['window']['id'] for r in data[key]] != ids for key in ('A15', 'A21')):
        raise ValueError('Original ordered unique TRAIN population changed')
    diagnostics = []
    for original, a15, a21 in zip(data['A12'], data['A15'], data['A21']):
        w = original['window']
        if w['split'] != 'train' or w['species'] not in SPECIES or a15['window'] != w or a21['window'] != w:
            raise ValueError('TRAIN window identity changed')
        free = [r for r in original['candidates'] if r['origin'] == 'free']
        if any(r['origin'] != 'free' and (r['origin'] != 'reference_injected' or r['label'] != 1)
               for r in original['candidates']):
            raise ValueError('Original injected positive semantics changed')
        identity = [{k: r[k] for k in IDENTITY} for r in free]
        for rec in (a15, a21):
            if [{k: r[k] for k in IDENTITY} for r in rec['candidates']] != identity:
                raise ValueError('Paired natural candidate identity changed')
        arms = {arm: audit_window(local, original['candidates']) for arm, local in
                zip(ARMS, (free, a15['candidates'], a21['candidates']))}
        for arm, audit in arms.items():
            if any(audit['counts'][k] != a21['natural_selection'][arm][k] for k in COUNT_KEYS):
                raise ValueError('Recomputed selection differs from closed A21 paired readout')
        diagnostics.append(dict(window=w, arms=arms))
    per_species, transitions = aggregate(diagnostics)
    for species in SPECIES:
        for arm in ARMS:
            if any(per_species[species][arm]['counts'][k] !=
                   summaries['A21_check']['natural_selection_per_species'][species][arm][k] for k in COUNT_KEYS):
                raise ValueError('Actual paired selection totals changed')
    pairs = sum(len(r['candidates']) for r in data['A21'])
    if pairs != 16321:
        raise ValueError('Natural candidate count changed')
    with (OUT/'windows.jsonl').open('x') as handle:
        for rec in diagnostics:
            handle.write(json.dumps(rec, allow_nan=False)+'\n')
    summary = dict(experiment='M28-A22-CHAIN-DISPLACEMENT-AUDIT', completed=True, protocol=PROTOCOL,
        source_jobs={k: v[1] for k, v in SOURCES.items()}, windows=len(diagnostics), candidate_window_pairs=pairs,
        original_injected_positive_metadata=sum(r['origin'] == 'reference_injected'
            for rec in data['A12'] for r in rec['candidates']),
        source_identity_and_actual_selection_reconciled=True, all_individual_and_joint_exchanges_passed=True,
        per_species=per_species, transitions=transitions, threshold=0., threshold_scan=False,
        reference_injected_into_selection=False, optimizer_steps_this_run=0, new_candidate_search=False,
        new_head_or_backbone_inference=False, DEV_forward_used=False, test_setaria_used=False,
        A21_scientific_gate_unchanged=True, seconds=time.perf_counter()-began,
        limitation='Frozen in-sample TRAIN diagnostic, not independent gene F1, biological negative truth or causal isolation.')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(summary, handle, indent=2, allow_nan=False)
    print(json.dumps(dict(completed=True, windows=len(diagnostics), seconds=summary['seconds'])), flush=True)


if __name__ == '__main__':
    main()
