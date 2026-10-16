#!/usr/bin/env python3
"""Actual A43 training exposure and saved full-grammar TRAIN decisions, CPU."""
import argparse
from collections import Counter
import json
import math
import time
import torch

import readout_natural_competitor_train as D
import verify_emission_lattice_train as V42
import verify_structured_emission_train as V36


def verify_training_rows(rows, draws, records):
    if len(rows) != 3*len(draws):
        raise ValueError('Complete fixed three-pass training receipts required')
    expected = Counter()
    for draw in draws: expected.update(records[draw['window_id']]['counts'])
    for i, row in enumerate(rows):
        epoch, ordinal = divmod(i, len(draws)); draw = draws[ordinal]; rec = records[draw['window_id']]
        weight = draw['importance_weight_uniform_within_stratum']
        if (row['step'] != i+1 or row['epoch'] != epoch+1 or row['draw_ordinal'] != ordinal or
                row['window_id'] != draw['window_id'] or row['species'] != draw['species'] or
                row['importance_weight'] != weight or any(row[k] != v for k, v in rec['counts'].items()) or
                any(not math.isfinite(row[k]) for k in ('weighted_loss', 'log_partition', 'gold_gain', 'gradient_norm_before_clip'))):
            raise ValueError('Actual original draw/weight/union exposure or finite update differs')
        known = rec['counts']['label_1']+rec['counts']['label_0']
        loss = (row['log_partition']-row['gold_gain'])/known*weight if known else 0.
        if not math.isclose(loss, row['weighted_loss'], rel_tol=1e-12, abs_tol=1e-8):
            raise ValueError('Actual known-only union NLL and original importance differ')
    return dict(expected)


@torch.inference_mode()
def verify_fit():
    fitted = D.require_fit()
    draws, _, records, support = D.F.A.population('13424089')
    training = D.R.C.R.rows(D.F.OUT/'training.jsonl')
    if len(training) != 4608 or len(fitted['epochs']) != 3:
        raise ValueError('Actual fixed A43 exposure incomplete')
    replay = verify_training_rows(training, draws, records)
    if replay != support['merged_replayed_draw_counts']:
        raise ValueError('Actual three-pass support replay differs')
    model = D.F.A.A36.new_model('cpu')
    for epoch, step in enumerate((1536, 3072, 4608)):
        path = D.F.OUT/('step_'+str(step).zfill(6)+'.pt'); report = fitted['epochs'][epoch]
        if (report['epoch'] != epoch+1 or report['steps'] != 1536 or report['counts'] != replay or
                report['checkpoint'] != str(path.relative_to(D.ROOT)) or not path.stat().st_size or
                not math.isclose(report['mean_weighted_loss'],
                    sum(r['weighted_loss'] for r in training[epoch*1536:(epoch+1)*1536])/1536,
                    rel_tol=1e-12, abs_tol=1e-8)):
            raise ValueError('Actual saved epoch/checkpoint exposure or recorded mean differs')
        state = torch.load(path, map_location='cpu', weights_only=True)
        if (state['arm'] != 'A43' or state['configuration'] != D.F.CONFIG or state['protocol'] != D.F.PROTOCOL or
                state['step'] != step or state['revision'] != D.R.C.I.F.REVISION or
                sum(v.numel() for v in state['model'].values()) != 366191 or
                any(v.dtype != torch.float32 or not bool(torch.isfinite(v).all()) for v in state['model'].values())):
            raise ValueError('Matching finite fixed A43 checkpoint required')
        model.load_state_dict(state['model'], strict=True)
        opt = state['optimizer']; groups = opt['param_groups']; params = list(model.parameters())
        if (len(groups) != 1 or len(groups[0]['params']) != len(opt['state']) or len(opt['state']) != len(params) or
                groups[0]['lr'] != D.F.CONFIG['lr'] or groups[0]['weight_decay'] != D.F.CONFIG['weight_decay'] or
                list(groups[0]['betas']) != D.F.CONFIG['betas'] or groups[0]['eps'] != D.F.CONFIG['eps']):
            raise ValueError('Actual fixed AdamW population or hyperparameters differ')
        for ident, parameter in zip(groups[0]['params'], params):
            item = opt['state'][ident]
            if (float(item['step']) != step or item['exp_avg'].shape != parameter.shape or
                    item['exp_avg_sq'].shape != parameter.shape or
                    any(v.dtype != torch.float32 or not bool(torch.isfinite(v).all())
                        for v in (item['exp_avg'], item['exp_avg_sq']))):
                raise ValueError('Actual finite AdamW parameter exposure differs')
        del state, opt
    return replay


@torch.inference_mode()
def main(potential_job, decode_job, readout_job):
    began = time.perf_counter(); D.F.A.B.D.B.F.A.require_jobs([readout_job])
    report = D.R.completed(D.OUT, readout_job); potentials = D.require_potentials(potential_job)
    if (report['configuration'] != D.CONFIG or report['protocol'] != D.PROTOCOL or
            report['source_jobs'] != dict(fit=D.FIT_JOB, potentials=potential_job, decode=decode_job) or
            report['windows'] != 1491 or report['original_draws'] != 1536 or report['original_candidates'] != 16321 or
            not report['all_natural_search_and_score_files_closed_before_labels'] or
            not report['original_R4_labels_exact'] or not report['new_geometry_labels_computed'] or
            report['R4_policy_or_catalog_changed'] or report['optimizer_steps_this_run'] != 0 or
            report['reference_injected_predictions'] != 0 or report['DEV_used'] or report['test_setaria_used'] or
            report['independent_species_used'] or report['A33_outputs_opened']):
        raise ValueError('Actual complete fixed full-grammar A43 readout required')
    replay = verify_fit()
    draws, windows, _ = D.population(); old_draws, free, prior, _, old = D.A42.B.population()
    shards = D.closed_shards(decode_job, potential_job)
    natural = D.R.C.R.rows(D.OUT/'natural_scores.jsonl'); labelled = D.R.C.R.rows(D.OUT/'windows.jsonl')
    if (draws != old_draws or len(natural) != 1491 or len(labelled) != 1491 or
            sorted(shards, key=lambda r: r['original_ordinal']) != natural or
            [r['original_ordinal'] for r in natural] != list(range(1491)) or
            len(list(D.POTENTIALS.glob('window_*.pt'))) != 1491):
        raise ValueError('Actual complete original ordered searches and potentials required')
    loaded_draws, _, target_windows, catalogs = D.R.Q.load_train_inputs()
    if loaded_draws != draws or set(target_windows) != {r['window']['id'] for r in windows}:
        raise ValueError('Original TRAIN-only target identity differs')
    source_labels = D.R.C.R.rows(D.R.READOUT/'windows.jsonl')
    if len(source_labels) != 1491:
        raise ValueError('Complete original saved label population required')
    actual_rows = []; max_error = 0.; minimum_margin = math.inf
    for ordinal, row in enumerate(natural):
        value = D.load_potential(ordinal, potentials['records'][ordinal], potentials)
        w = row['window']; local = value['local_logits']; sequence = value['sequence'].upper()
        baseline = prior[ordinal]['variants']['A31_A3']
        if (w != windows[ordinal]['window'] or w != free[ordinal]['window'] or w != value['window'] or
                w != prior[ordinal]['window'] or row['reference_used_in_scoring'] or
                row['source_fit_job'] != D.FIT_JOB or row['potential_job'] != potential_job or
                row['potential_artifact'] != str(D.potential_path(ordinal).relative_to(D.ROOT)) or
                row['feature_revision'] != potentials['feature_revision'] or row['shard'] != ordinal % 32 or
                row['selected'] != list(range(len(row['chains']))) or free[ordinal]['chains'] != baseline['chains']):
            raise ValueError('Actual source, original baseline or natural geometry differs')
        # The original finite pool is first opened after all searches close;
        # rescoring is a numerical lower bound, never a prediction candidate set.
        for chain in baseline['chains']: D.A42.P.check_chain(sequence, chain)
        D.validate_prediction(value, dict(chains=row['chains'], chain_gains=row['scores'],
                                         gain=row['gain'], reference_used=row['reference_used_in_scoring']))
        baseline_scores = D.A42.P.local_path_scores(local, baseline['chains']).tolist()
        max_error = max(max_error, V36.verify_direct_scores(local, baseline['chains']+row['chains'], baseline_scores+row['scores']))
        selected = D.A42.P.select_nonoverlapping(baseline['chains'], baseline_scores)
        finite_gain = float(torch.tensor([baseline_scores[i] for i in selected], dtype=torch.float64).sum())
        if row['gain']+1e-8+1e-12*abs(finite_gain) < finite_gain:
            raise ValueError('Full grammar objective below same A43 original-natural-pool objective')
        minimum_margin = min(minimum_margin, row['gain']-finite_gain)
        target = target_windows[w['id']]; D.A42.validate_target_geometry(target, w)
        targets = D.R.Q.build_targets(target, catalogs[w['species']], sequence)
        actual = D.label_record(row, prior[ordinal], source_labels[ordinal], targets)
        for arm, candidate in (('A31_A3', baseline), ('A43', row)):
            if V42.counts(candidate, actual['labels'][arm]) != actual['counts'][arm]:
                raise ValueError('Actual selected counts or null0 DP differ')
        if actual != labelled[ordinal]:
            raise ValueError('Actual frozen R4 labels, unknowns or complete denominators differ')
        actual_rows.append(actual)
        del value, local, targets
        if (ordinal+1) % 128 == 0 or ordinal+1 == 1491:
            print(json.dumps(dict(verified_windows=ordinal+1, total=1491)), flush=True)
    counts, refs = D.aggregate(actual_rows); D.A42.check_baseline(counts, refs, old); checks = D.gate(counts)
    if (counts != report['selected_counts_per_species'] or refs != report['eligible_reference_pairs'] or
            checks != report['predeclared_TRAIN_gate'] or all(checks.values()) != report['TRAIN_preparation_supported']):
        raise ValueError('Actual aggregates or original five-condition arithmetic differ')
    D.R.C.R.write_summary(D.CHECK, dict(completed=True, configuration=D.CONFIG, protocol=D.PROTOCOL,
        source_jobs=dict(fit=D.FIT_JOB, potentials=potential_job, decode=decode_job, readout=readout_job),
        windows=1491, original_draws=1536, original_candidates=16321, candidate_window_pairs=19189,
        actual_original_draw_weights_union_NLL_and_4608_updates_verified=True,
        actual_three_finite_A43_checkpoints_and_AdamW_exposure_verified=True,
        actual_frozen_last_model_potential_identity_and_complete_32_shards_verified=True,
        actual_new_DNA_ORF_motif_cap_span_and_full_window_FP64_log_probability_scores_verified=True,
        actual_full_grammar_objective_not_below_same_model_original_natural_pool=True,
        actual_frozen_R4_unknowns_stronger_baseline_and_complete_denominators_recomputed=True,
        actual_per_species_counts_and_original_five_conditions_recomputed=True,
        replayed_original_draw_counts=replay, max_direct_log_probability_error=max_error,
        minimum_full_grammar_minus_rescored_original_pool_objective=minimum_margin,
        selected_counts_per_species=counts, eligible_reference_pairs=refs,
        predeclared_TRAIN_gate=checks, TRAIN_preparation_supported=all(checks.values()),
        model_or_backbone_forward=False, independent_neural_forward_replay=False,
        independent_full_grammar_DP_replay=False, optimizer_steps_this_run=0, new_candidate_search=False,
        new_geometry_labels_computed=True, R4_policy_or_catalog_changed=False, reference_injected_predictions=0,
        DEV_used=False, test_setaria_used=False, independent_species_used=False, A33_outputs_opened=False,
        seconds=time.perf_counter()-began,
        limitation='Actual finite-support fit and saved natural TRAIN decisions, not independent neural/DP replay, '
            'raw assembly, biological false positives, causal ablation or generalization.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--potential-job', required=True); parser.add_argument('--decode-job', required=True)
    parser.add_argument('--readout-job', required=True); args = parser.parse_args(); torch.set_num_threads(2)
    main(args.potential_job, args.decode_job, args.readout_job)
