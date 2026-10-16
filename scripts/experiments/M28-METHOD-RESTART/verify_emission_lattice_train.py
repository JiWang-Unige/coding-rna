#!/usr/bin/env python3
"""A42 CPU saved-potential, geometry, frozen R4 and five-condition reconciliation."""
import argparse
from collections import Counter
import json
import math
import time

import torch
import run_emission_lattice_train as D
from src.m28.labels import candidate_targets
from src.m28.phase_evidence_chain import candidate_roles


def verify_scores(local, chains, expected):
    if (local.dtype != torch.float32 or local.ndim != 2 or local.shape[1] != 15 or
            not bool(torch.isfinite(local).all())):
        raise ValueError('Finite saved FP32 local emission matrix required')
    # Full-window role assignment, separately from the decoder prefix arithmetic.
    logits = local.double(); positions = torch.arange(len(logits))
    actual = [float((logits[positions, candidate_roles(c, len(logits))]-logits[:, 0]).sum()) for c in chains]
    torch.testing.assert_close(torch.tensor(actual, dtype=torch.float64),
                               torch.tensor(expected, dtype=torch.float64), rtol=1e-12, atol=1e-8)
    return actual, max((abs(a-b) for a, b in zip(actual, expected)), default=0.)


def counts(value, labels):
    chains, scores, selected = value['chains'], value['scores'], value['selected']
    if (len(chains) != len(scores) or len(chains) != len(labels) or
            any(y not in (-1, 0, 1) for y in labels) or not all(math.isfinite(s) for s in scores) or
            selected != D.P.select_nonoverlapping(chains, scores)):
        raise ValueError('Actual labels, finite scores and original DP selection must reconcile')
    return dict(selected=len(selected), selected_correct=sum(labels[i] == 1 for i in selected),
        selected_known_nonmatch=sum(labels[i] == 0 for i in selected),
        selected_unknown=sum(labels[i] == -1 for i in selected),
        selected_five_or_more_correct=sum(labels[i] == 1 and len(chains[i]) >= 5 for i in selected))


@torch.inference_mode()
def main(potential_job, decode_job, readout_job):
    began = time.perf_counter(); D.B.F.A.require_jobs([readout_job])
    report = D.R.completed(D.OUT, readout_job); potentials = D.require_potentials(potential_job)
    if (report['configuration'] != D.CONFIG or report['protocol'] != D.PROTOCOL or
            report['source_jobs'] != dict(potentials=potential_job, decode=decode_job) or
            report['windows'] != 1491 or report['original_draws'] != 1536 or report['original_candidates'] != 16321 or
            not report['all_natural_search_and_score_files_closed_before_labels'] or
            not report['new_geometry_labels_computed'] or not report['original_R4_labels_exact'] or
            report['R4_policy_or_catalog_changed'] or report['optimizer_steps_this_run'] != 0 or
            report['reference_injected_predictions'] != 0 or report['DEV_used'] or
            report['test_setaria_used'] or report['A33_outputs_opened']):
        raise ValueError('Actual complete fixed A42 TRAIN readout required')
    draws, free, prior, _, old = D.B.population()
    shards = D.closed_shards(decode_job, potential_job)
    natural = D.R.C.R.rows(D.OUT/'natural_scores.jsonl')
    labelled = D.R.C.R.rows(D.OUT/'windows.jsonl')
    if (len(natural) != 1491 or len(labelled) != 1491 or len(shards) != 1491 or
            [r['original_ordinal'] for r in natural] != list(range(1491)) or
            sorted(r['original_ordinal'] for r in shards) != list(range(1491)) or
            sorted(shards, key=lambda r: r['original_ordinal']) != natural or
            len(list(D.POTENTIALS.glob('window_*.pt'))) != 1491):
        raise ValueError('All actual ordered 1491 records, 32 shards and potentials required')
    loaded_draws, _, windows, catalogs = D.R.Q.load_train_inputs()
    if loaded_draws != draws or set(windows) != {p['window']['id'] for p in free}:
        raise ValueError('Original TRAIN target population differs')
    source_labels = D.R.C.R.rows(D.R.READOUT/'windows.jsonl')
    if len(source_labels) != 1491:
        raise ValueError('All original label rows required')
    totals = {s: {a: Counter() for a in D.ARMS} for s in D.R.C.R.ALLOW}
    refs = {s: Counter(eligible=0, five_or_more=0) for s in D.R.C.R.ALLOW}
    maximum_old_error = 0.; maximum_new_error = 0.
    for ordinal, row in enumerate(natural):
        pool = free[ordinal]; baseline = prior[ordinal]['variants']['A31_A3']
        value = D.load_potential(ordinal, potentials['records'][ordinal], potentials)
        if (row['window'] != pool['window'] or row['window'] != value['window'] or
                row['reference_used_in_scoring'] or row['potential_job'] != potential_job or
                row['shard'] != ordinal % 32 or row['feature_revision'] != potentials['feature_revision'] or
                row['potential_artifact'] != str(D.potential_path(ordinal).relative_to(D.ROOT)) or
                value['original_chains'] != pool['chains'] or value['original_chains'] != baseline['chains'] or
                value['original_scores'] != baseline['scores'] or value['original_selected'] != baseline['selected'] or
                row['selected'] != list(range(len(row['chains'])))):
            raise ValueError('Actual source, original pool or new ordered window identity differs')
        sequence = value['sequence'].upper(); local = value['local_logits']; previous_end = 0
        for chain in value['original_chains']:
            D.P.check_chain(sequence, chain)
        for chain in row['chains']:
            D.P.check_chain(sequence, chain)
            if chain[0][0] < previous_end:
                raise ValueError('Actual selected complete-gene spans overlap')
            previous_end = chain[-1][1]
        old_scores, error = verify_scores(local, baseline['chains'], baseline['scores'])
        maximum_old_error = max(maximum_old_error, error)
        new_scores, error = verify_scores(local, row['chains'], row['scores'])
        maximum_new_error = max(maximum_new_error, error)
        if D.P.select_nonoverlapping(baseline['chains'], old_scores) != baseline['selected']:
            raise ValueError('Original baseline DP selection differs')
        old_gain = float(torch.tensor([old_scores[i] for i in baseline['selected']], dtype=torch.float64).sum())
        new_gain = torch.tensor(new_scores, dtype=torch.float64).sum()
        torch.testing.assert_close(new_gain, torch.tensor(row['gain'], dtype=torch.float64), rtol=1e-12, atol=1e-8)
        if (any(s <= 0 for s in row['scores']) or row['gain'] + 1e-8 + 1e-12*abs(old_gain) < old_gain or
                abs(row['original_baseline_gain']-old_gain) > 1e-8+1e-12*abs(old_gain)):
            raise ValueError('Full grammar objective or positive selected gains differ')
        w = windows[row['window']['id']]
        if {k: w[k] for k in D.R.C.R.GEOMETRY} != row['window']:
            raise ValueError('Actual frozen target window differs')
        targets = D.R.Q.build_targets(w, catalogs[w['species']], sequence)
        old_labels = candidate_targets(baseline['chains'], targets)
        new_labels = candidate_targets(row['chains'], targets)
        original = source_labels[ordinal]
        old_counts = counts(baseline, old_labels); new_counts = counts(row, new_labels)
        eligible = len(targets['positive_chains']); long = sum(len(c) >= 5 for c in targets['positive_chains'])
        actual = dict(window=row['window'], labels={'A31_A3': old_labels, 'A42': new_labels},
                      counts={'A31_A3': old_counts, 'A42': new_counts},
                      eligible_references=eligible, eligible_five_or_more=long)
        if (actual != labelled[ordinal] or original['window'] != row['window'] or
                old_labels != original['labels']['A31_A3'] or
                old_counts != D.selected_only(original['counts']['A31_A3']) or
                eligible != original['eligible_references'] or long != original['eligible_five_or_more']):
            raise ValueError('Actual frozen R4 labels, unknowns, counts or denominators differ')
        for arm, count in actual['counts'].items():
            totals[w['species']][arm].update(count)
        refs[w['species']].update(eligible=eligible, five_or_more=long)
        del value, local, targets
        if (ordinal+1) % 128 == 0 or ordinal+1 == 1491:
            print(json.dumps(dict(verified_windows=ordinal+1, total=1491)), flush=True)
    selected = {s: {a: dict(c) for a, c in v.items()} for s, v in totals.items()}
    references = {s: dict(c) for s, c in refs.items()}
    if (selected != report['selected_counts_per_species'] or references != report['eligible_reference_pairs'] or
            sum(v['eligible'] for v in references.values()) != 2916 or
            sum(v['five_or_more'] for v in references.values()) != 1086 or
            any(selected[s]['A31_A3'] != D.selected_only(old['selected_counts_per_species'][s]['A31_A3']) for s in selected)):
        raise ValueError('Actual per-species aggregates, baseline or full denominators differ')
    checks = {}
    for s, values in selected.items():
        checks[s+'_natural_selected_correct_strictly_increased'] = values['A42']['selected_correct'] > values['A31_A3']['selected_correct']
        checks[s+'_natural_selected_known_nonmatch_not_increased'] = values['A42']['selected_known_nonmatch'] <= values['A31_A3']['selected_known_nonmatch']
    checks['pooled_natural_selected_five_or_more_correct_strictly_increased'] = (
        sum(v['A42']['selected_five_or_more_correct'] for v in selected.values()) >
        sum(v['A31_A3']['selected_five_or_more_correct'] for v in selected.values()))
    if checks != report['predeclared_TRAIN_gate'] or all(checks.values()) != report['TRAIN_preparation_supported']:
        raise ValueError('Actual original five-condition arithmetic differs')
    D.R.C.R.write_summary(D.CHECK, dict(completed=True, configuration=D.CONFIG, protocol=D.PROTOCOL,
        source_jobs=dict(potentials=potential_job, decode=decode_job, readout=readout_job),
        windows=1491, original_draws=1536, original_candidates=16321,
        actual_full_ordered_1491_TRAIN_windows_original_1536_draws_identity_verified=True,
        actual_frozen_A3_saved_potential_identity_and_no_update_receipts_verified=True,
        actual_original_A31_FP64_scores_DP_and_stronger_baseline_recomputed=True,
        actual_new_chains_DNA_ORF_motif_cap_span_and_FP64_gain_recomputed=True,
        actual_complete_32_shards_natural_search_closed_before_labels_verified=True,
        actual_frozen_R4_labels_unknown_protection_and_denominators_recomputed=True,
        actual_per_species_counts_and_five_conditions_recomputed=True,
        max_original_A31_score_error=maximum_old_error, max_new_gain_error=maximum_new_error,
        selected_counts_per_species=selected, eligible_reference_pairs=references,
        predeclared_TRAIN_gate=checks, TRAIN_preparation_supported=all(checks.values()),
        model_or_backbone_forward=False, independent_neural_forward_replay=False,
        optimizer_steps_this_run=0, new_candidate_search=False, new_geometry_labels_computed=True,
        R4_policy_or_catalog_changed=False, DEV_used=False, test_setaria_used=False, A33_outputs_opened=False,
        seconds=time.perf_counter()-began,
        limitation='CPU checks saved potentials, geometry and frozen window labels, not neural replay, raw annotation or generalization.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--potential-job', required=True); parser.add_argument('--decode-job', required=True)
    parser.add_argument('--readout-job', required=True)
    args = parser.parse_args(); torch.set_num_threads(2)
    main(args.potential_job, args.decode_job, args.readout_job)
