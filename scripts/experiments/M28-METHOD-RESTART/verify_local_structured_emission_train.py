#!/usr/bin/env python3
"""Actual paired A44 exposure, saved predictions and frozen comparison gates."""
import argparse
from collections import Counter
import json
import math
import time
import torch

import readout_local_structured_emission_train as D
import verify_emission_lattice_train as V42
import verify_structured_emission_train as V36


def verify_training_rows(arm, rows, draws, records, windows, local_positions):
    if len(rows) != 3*len(draws):
        raise ValueError('Complete fixed paired three-pass receipts required')
    coefficients = D.F.CONFIG['loss_coefficients_by_arm'][arm]; expected = Counter()
    for draw in draws: expected.update(records[draw['window_id']]['counts'])
    finite = ('weighted_loss', 'local_loss', 'structured_loss', 'structured_candidate_mean_loss',
              'log_partition', 'gold_gain', 'gradient_norm_before_clip')
    for i, row in enumerate(rows):
        epoch, ordinal = divmod(i, len(draws)); draw = draws[ordinal]
        ident = draw['window_id']; rec = records[ident]; w = windows[ident]
        known = rec['counts']['label_1']+rec['counts']['label_0'] if arm == 'A44_JOINT' else 0
        unknown = rec['counts']['label_-1'] if arm == 'A44_JOINT' else 0
        weight = draw['importance_weight_uniform_within_stratum']
        if (row['arm'] != arm or row['step'] != i+1 or row['epoch'] != epoch+1
                or row['draw_ordinal'] != ordinal or row['window_id'] != ident or row['species'] != draw['species']
                or row['valid_bases'] != w['valid_bases'] or row['importance_weight'] != weight
                or row['loss_coefficients'] != coefficients
                or row['structured_known_candidates_used'] != known
                or row['structured_unknown_chains_excluded'] != unknown
                or row['local_supervised_positions'] != local_positions[ident]
                or any(row[k] != v for k, v in rec['counts'].items())
                or any(not math.isfinite(row[k]) for k in finite)):
            raise ValueError('Actual paired draw, loss source, masks or finite update differs')
        numerator = row['log_partition']-row['gold_gain']
        if not known and any(row[k] != 0 for k in ('log_partition', 'gold_gain', 'structured_loss', 'structured_candidate_mean_loss')):
            raise ValueError('Empty paired structural input must retain exact zero')
        structural = numerator/w['valid_bases']
        old_mean = numerator/known if known else 0.
        total = (row['local_loss']+coefficients[1]*structural)*weight
        if (not math.isclose(structural, row['structured_loss'], rel_tol=1e-12, abs_tol=1e-8)
                or not math.isclose(old_mean, row['structured_candidate_mean_loss'], rel_tol=1e-12, abs_tol=1e-8)
                or not math.isclose(total, row['weighted_loss'], rel_tol=1e-12, abs_tol=1e-8)
                or (not local_positions[ident] and row['local_loss'] != 0)):
            raise ValueError('Actual per-valid-base NLL, local zero or once-only importance differs')
    return dict(expected)


@torch.inference_mode()
def verify_checkpoints(arm, fitted, training, replay):
    if len(fitted['epochs']) != 3 or len(training) != 4608:
        raise ValueError('Actual fixed paired exposure incomplete')
    model = D.F.A.A36.new_model('cpu'); params = list(model.parameters())
    for epoch, step in enumerate((1536, 3072, 4608)):
        path = D.F.output_root(arm)/('step_'+str(step).zfill(6)+'.pt'); report = fitted['epochs'][epoch]
        if (report['epoch'] != epoch+1 or report['steps'] != 1536 or report['counts'] != replay
                or report['checkpoint'] != str(path.relative_to(D.ROOT)) or not path.stat().st_size
                or not math.isclose(report['mean_weighted_loss'],
                    sum(r['weighted_loss'] for r in training[epoch*1536:(epoch+1)*1536])/1536,
                    rel_tol=1e-12, abs_tol=1e-8)):
            raise ValueError('Actual paired checkpoint exposure or recorded epoch mean differs')
        state = torch.load(path, map_location='cpu', weights_only=True)
        if (state['arm'] != arm or state['configuration'] != D.F.CONFIG or state['protocol'] != D.F.PROTOCOL
                or state['step'] != step or state['revision'] != D.R.C.I.F.REVISION
                or sum(v.numel() for v in state['model'].values()) != 366191
                or any(v.dtype != torch.float32 or not bool(torch.isfinite(v).all()) for v in state['model'].values())):
            raise ValueError('Actual finite same-capacity paired checkpoint differs')
        model.load_state_dict(state['model'], strict=True)
        optimizer = state['optimizer']; groups = optimizer['param_groups']
        if (len(groups) != 1 or len(groups[0]['params']) != len(params) or len(optimizer['state']) != len(params)
                or groups[0]['lr'] != D.F.CONFIG['lr'] or groups[0]['weight_decay'] != D.F.CONFIG['weight_decay']
                or list(groups[0]['betas']) != D.F.CONFIG['betas'] or groups[0]['eps'] != D.F.CONFIG['eps']):
            raise ValueError('Actual paired AdamW population or hyperparameters differ')
        for ident, parameter in zip(groups[0]['params'], params):
            item = optimizer['state'][ident]
            if (float(item['step']) != step or item['exp_avg'].shape != parameter.shape
                    or item['exp_avg_sq'].shape != parameter.shape
                    or any(v.dtype != torch.float32 or not bool(torch.isfinite(v).all())
                           for v in (item['exp_avg'], item['exp_avg_sq']))):
                raise ValueError('Actual finite paired AdamW parameter exposure differs')
        del state, optimizer


@torch.inference_mode()
def main(potential_job, decode_job, readout_job):
    began = time.perf_counter(); D.F.A.B.D.B.F.A.require_jobs([readout_job])
    report = D.R.completed(D.OUT, readout_job); potentials = D.require_potentials(potential_job)
    if (report['configuration'] != D.CONFIG or report['protocol'] != D.PROTOCOL
            or report['source_jobs'] != dict(fits=D.FIT_JOBS, potentials=potential_job, decode=decode_job)
            or report['windows_per_arm'] != 1491 or report['original_draws'] != 1536
            or report['original_candidates'] != 16321 or report['reference_injected_predictions'] != 0
            or not report['all_both_arm_natural_search_and_score_files_closed_before_labels']
            or not report['original_R4_labels_exact'] or not report['new_geometry_labels_computed']
            or report['R4_policy_or_catalog_changed'] or report['optimizer_steps_this_run'] != 0
            or any(report[k] for k in ('DEV_used', 'test_setaria_used', 'independent_species_used', 'A33_outputs_opened'))):
        raise ValueError('Actual complete fixed paired natural TRAIN readout required')
    fitted = D.require_fits(); draws, windows, _ = D.D.population()
    support_draws, _, records, support, targets_windows, catalogs, _ = D.F.P.population('13424089')
    old_draws, free, prior, _, old = D.D.A42.B.population()
    shards = D.closed_shards(decode_job, potential_job)
    natural = {a: D.R.C.R.rows(D.OUT/(a+'_natural_scores.jsonl')) for a in D.ARMS}
    labelled = D.R.C.R.rows(D.OUT/'windows.jsonl'); source_labels = D.R.C.R.rows(D.R.READOUT/'windows.jsonl')
    if (draws != support_draws or draws != old_draws or len(labelled) != 1491 or len(source_labels) != 1491
            or any(len(natural[a]) != 1491 or sorted(shards[a], key=lambda r: r['original_ordinal']) != natural[a]
                   or [r['original_ordinal'] for r in natural[a]] != list(range(1491))
                   or len(list((D.POTENTIALS/a).glob('window_*.pt'))) != 1491 for a in D.ARMS)):
        raise ValueError('Complete matching paired predictions, sources and potentials required')
    rows = []; local_positions = {}; max_error = {a: 0. for a in D.ARMS}; min_margin = {a: math.inf for a in D.ARMS}
    for ordinal in range(1491):
        values = {a: D.load_potential(a, ordinal, potentials['records'][a][ordinal], potentials) for a in D.ARMS}
        sequence = values[D.ARMS[0]]['sequence'].upper(); w = windows[ordinal]['window']
        if (values[D.ARMS[1]]['sequence'].upper() != sequence or any(values[a]['window'] != w for a in D.ARMS)
                or w != free[ordinal]['window'] or w != prior[ordinal]['window']):
            raise ValueError('Actual paired geometry or original oriented DNA differs')
        baseline = prior[ordinal]['variants']['A31_A3']
        if free[ordinal]['chains'] != baseline['chains']:
            raise ValueError('Original stronger natural pool geometry differs')
        for chain in baseline['chains']: D.D.A42.P.check_chain(sequence, chain)
        paired = {a: natural[a][ordinal] for a in D.ARMS}
        for arm in D.ARMS:
            row = paired[arm]; local = values[arm]['local_logits']
            if (row['arm'] != arm or row['window'] != w or row['potential_job'] != potential_job
                    or row['source_fit_job'] != D.FIT_JOBS[arm]
                    or row['potential_artifact'] != str(D.potential_path(arm, ordinal).relative_to(D.ROOT))
                    or row['feature_revision'] != potentials['feature_revision']
                    or row['task'] != D.ARMS.index(arm)*32+ordinal % 32 or row['shard'] != ordinal % 32
                    or row['reference_used_in_scoring'] or row['selected'] != list(range(len(row['chains'])))):
                raise ValueError('Actual paired prediction source, arm or selection differs')
            D.D.validate_prediction(values[arm], dict(chains=row['chains'], chain_gains=row['scores'],
                                                     gain=row['gain'], reference_used=False))
            old_scores = D.D.A42.P.local_path_scores(local, baseline['chains']).tolist()
            max_error[arm] = max(max_error[arm], V36.verify_direct_scores(local,
                baseline['chains']+row['chains'], old_scores+row['scores']))
            selected = D.D.A42.P.select_nonoverlapping(baseline['chains'], old_scores)
            finite_gain = float(torch.tensor([old_scores[i] for i in selected], dtype=torch.float64).sum())
            if row['gain']+1e-8+1e-12*abs(finite_gain) < finite_gain:
                raise ValueError('Full grammar objective below same paired model original-natural-pool objective')
            min_margin[arm] = min(min_margin[arm], row['gain']-finite_gain)
        target = targets_windows[w['id']]; D.D.A42.validate_target_geometry(target, w)
        targets = D.R.Q.build_targets(target, catalogs[w['species']], sequence)
        fine, region = targets['seg'], targets['region']
        local_positions[w['id']] = int(((fine >= 0) | ((fine < 0) & (region >= 0) & (region < 3))).sum())
        actual = D.label_paired_record(paired, prior[ordinal], source_labels[ordinal], targets)
        for arm, prediction in (('A31_A3', baseline), *((a, paired[a]) for a in D.ARMS)):
            if V42.counts(prediction, actual['labels'][arm]) != actual['counts'][arm]:
                raise ValueError('Actual paired selected counts or null0 selection differs')
        if actual != labelled[ordinal]:
            raise ValueError('Actual paired R4 labels, unknowns or complete denominators differ')
        rows.append(actual); del values, local, targets
        if (ordinal+1) % 128 == 0 or ordinal+1 == 1491: print(json.dumps(dict(paired_verified_windows=ordinal+1, total=1491)), flush=True)
    replay = {}
    for arm in D.ARMS:
        training = D.R.C.R.rows(D.F.output_root(arm)/'training.jsonl')
        replay[arm] = verify_training_rows(arm, training, draws, records, targets_windows, local_positions)
        if replay[arm] != support['merged_replayed_draw_counts']:
            raise ValueError('Actual paired original support exposure differs')
        verify_checkpoints(arm, fitted[arm], training, replay[arm])
    counts, refs = D.aggregate(rows); D.D.A42.check_baseline(counts, refs, old); checks, supported = D.gates(counts)
    if (counts != report['selected_counts_per_species'] or refs != report['eligible_reference_pairs']
            or checks != report['predeclared_TRAIN_gates'] or supported != report['TRAIN_preparation_supported']):
        raise ValueError('Actual paired aggregates or fixed three five-condition comparisons differ')
    D.R.C.R.write_summary(D.CHECK, dict(completed=True, configuration=D.CONFIG, protocol=D.PROTOCOL,
        source_jobs=dict(fits=D.FIT_JOBS, potentials=potential_job, decode=decode_job, readout=readout_job),
        windows_per_arm=1491, original_draws=1536, original_candidates=16321, candidate_window_pairs=19189,
        actual_both_original_draw_weights_loss_units_masks_and_4608_updates_verified=True,
        actual_six_finite_paired_checkpoints_and_AdamW_exposure_verified=True,
        actual_paired_last_model_potential_identity_and_complete_64_tasks_verified=True,
        actual_new_DNA_ORF_motif_cap_span_and_full_window_FP64_log_probability_scores_verified=True,
        actual_both_full_grammar_objectives_not_below_same_model_original_natural_pool=True,
        actual_frozen_R4_unknowns_stronger_baseline_and_complete_denominators_recomputed=True,
        actual_per_species_counts_and_fixed_three_five_condition_comparisons_recomputed=True,
        replayed_original_draw_counts=replay, max_direct_log_probability_error=max_error,
        minimum_full_grammar_minus_rescored_original_pool_objective=min_margin,
        selected_counts_per_species=counts, eligible_reference_pairs=refs,
        predeclared_TRAIN_gates=checks, TRAIN_preparation_supported=supported,
        model_or_backbone_forward=False, independent_neural_forward_replay=False,
        independent_full_grammar_DP_replay=False, new_candidate_search=False,
        new_geometry_labels_computed=True, reference_injected_predictions=0, seconds=time.perf_counter()-began,
        limitation='Actual paired TRAIN fit and saved natural decisions, not independent neural/DP replay, '
            'raw assembly, biological false positives, matched total supervision or generalization.', **D.SCOPE))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--potential-job', required=True); parser.add_argument('--decode-job', required=True)
    parser.add_argument('--readout-job', required=True)
    args = parser.parse_args(); torch.set_num_threads(2); main(args.potential_job, args.decode_job, args.readout_job)
