#!/usr/bin/env python3
"""CPU verification of actual A39 branch exposure, saved components and gates."""
import argparse
from collections import Counter
import math
import time
import torch
import readout_null_anchored_set_train as D
from verify_context_chain_train import verify_components, expected_optimizer_steps


@torch.inference_mode()
def verify_fit():
    fitted = D.R.completed(D.F.OUT, D.FIT_JOB); D.validate_fit(fitted)
    draws, _, _, records, _, _, _ = D.F.A.population()
    training = D.R.C.R.rows(D.F.OUT/'training.jsonl')
    if len(training) != 4608 or len(fitted['epochs']) != 3:
        raise ValueError('Complete actual original three-pass fit required')
    expected_counts = Counter()
    for draw in draws:
        expected_counts.update(records[draw['window_id']]['counts'])
    for i, row in enumerate(training):
        epoch, ordinal = divmod(i, 1536); draw = draws[ordinal]; rec = records[draw['window_id']]
        if (row['step'] != i+1 or row['epoch'] != epoch+1 or row['draw_ordinal'] != ordinal or
                row['window_id'] != draw['window_id'] or row['species'] != draw['species'] or
                row['importance_weight'] != draw['importance_weight_uniform_within_stratum'] or
                any(row[k] != v for k, v in rec['counts'].items()) or
                any(not math.isfinite(row[k]) for k in ('weighted_loss','unweighted_mean_known_BCE','unweighted_mean_known_set_NLL',
                    'unweighted_fixed_composite_loss','log_partition','gold_gain','gradient_norm_before_clip')) or
                (rec['counts']['zero_known_windows'] and any(row[k] != 0 for k in
                    ('weighted_loss','unweighted_mean_known_BCE','unweighted_mean_known_set_NLL',
                     'unweighted_fixed_composite_loss','log_partition','gold_gain','gradient_norm_before_clip')))):
            raise ValueError('Actual draw/order/weight/counts or finite protected update record differs')
        combined = (row['unweighted_mean_known_BCE']+row['unweighted_mean_known_set_NLL'])/2
        if (not math.isclose(row['unweighted_fixed_composite_loss'], combined, rel_tol=1e-6, abs_tol=1e-7) or
                not math.isclose(row['weighted_loss'], row['unweighted_fixed_composite_loss']*row['importance_weight'],
                                 rel_tol=1e-6, abs_tol=1e-7)):
            raise ValueError('Actual fixed equal composite and original importance arithmetic differs')
    model = D.F.A.new_model('cpu'); named = list(model.named_parameters()); final_exposure = None
    for epoch, step in enumerate((1536,3072,4608)):
        path = D.F.OUT/('step_'+str(step).zfill(6)+'.pt'); report = fitted['epochs'][epoch]
        if (report['epoch'] != epoch+1 or report['steps'] != 1536 or
                report['counts'] != dict(expected_counts) or report['checkpoint'] != str(path.relative_to(D.ROOT)) or
                not path.stat().st_size):
            raise ValueError('Actual saved epoch/checkpoint exposure differs')
        state = torch.load(path, map_location='cpu', weights_only=True)
        if (state['arm'] != 'A39' or state['configuration'] != D.CONFIG or state['protocol'] != D.F.PROTOCOL or
                state['step'] != step or state['revision'] != D.F.A.R.C.I.F.REVISION or
                sum(v.numel() for v in state['model'].values()) != 612324 or
                any(v.dtype != torch.float32 or not bool(torch.isfinite(v).all()) for v in state['model'].values())):
            raise ValueError('Finite complete fixed A39 checkpoint structure required')
        model.load_state_dict(state['model'], strict=True)
        opt = state['optimizer']; groups = opt['param_groups']
        if (len(groups) != 1 or len(groups[0]['params']) != len(named) or
                len(set(groups[0]['params'])) != len(named) or groups[0]['lr'] != D.CONFIG['lr'] or
                groups[0]['weight_decay'] != D.CONFIG['weight_decay'] or list(groups[0]['betas']) != D.CONFIG['betas'] or
                groups[0]['eps'] != D.CONFIG['eps']):
            raise ValueError('Actual fixed AdamW configuration differs')
        exposure = expected_optimizer_steps([k for k, _ in named], draws, records, epoch+1)
        expected_ids = {ident for ident, (name, _) in zip(groups[0]['params'], named) if exposure[name]}
        if set(opt['state']) != expected_ids:
            raise ValueError('Actual called/uncalled branch optimizer population differs')
        for ident, (name, param) in zip(groups[0]['params'], named):
            if not exposure[name]:
                continue
            value = opt['state'][ident]
            if (float(value['step']) != exposure[name] or
                    any(v.shape != param.shape or v.dtype != torch.float32 or not bool(torch.isfinite(v).all())
                        for v in (value['exp_avg'], value['exp_avg_sq']))):
                raise ValueError('Actual finite per-parameter optimizer exposure differs: '+name)
        final_exposure = exposure
        del state, opt
    if fitted['optimizer_steps_per_parameter'] != final_exposure:
        raise ValueError('Final actual branch-step summary differs from checkpoint and original pool')
    return dict(expected_counts), final_exposure


@torch.inference_mode()
def main(readout_job):
    began = time.perf_counter(); D.F.A.require_jobs([D.FIT_JOB, str(readout_job)])
    report = D.R.completed(D.OUT, readout_job)
    if (report['protocol'] != D.PROTOCOL or report['configuration'] != D.CONFIG or
            report['source_fit_job'] != D.FIT_JOB or report['windows'] != 1491 or
            report['natural_candidates'] != 16321 or report['reference_injected_predictions'] != 0 or
            any(report[k] is not True for k in ('natural_score_file_closed_before_saved_label_arrays_open',
                'original_natural_geometry_order_sources_unchanged','original_saved_labels_reused',
                'model_weights_unchanged','all_parameters_frozen_and_gradients_absent')) or
            report['optimizer_steps_this_run'] != 0 or any(report[k] is not False for k in
                ('label_catalog_opened','new_R4_labelling','new_candidate_search','new_backbone_forward',
                 'DEV_used','test_setaria_used','A33_outputs_opened'))):
        raise ValueError('Actual complete unchanged natural A39 readout required')
    replay, exposure = verify_fit()
    _, free, previous, _, baseline = D.B.population()
    saved = D.R.C.R.rows(D.OUT/'natural_scores.jsonl')
    labelled = D.R.C.R.rows(D.OUT/'windows.jsonl')
    if len(saved) != 1491 or len(labelled) != 1491:
        raise ValueError('All1491 saved natural scores and labelled records required')
    for row, pool in zip(saved, free):
        if (row['window'] != pool['window'] or row['chains'] != pool['chains'] or
                row['decoder_sources'] != pool['decoder_sources']):
            raise ValueError('Original saved natural identity/order/sources differ')
        verify_components(row)
    source_labels = D.R.C.R.rows(D.R.READOUT/'windows.jsonl')
    context = D.R.C.R.rows(D.L.F.A.D.OUT/'natural_scores.jsonl')
    logistic = D.R.C.R.rows(D.L.OUT/'natural_scores.jsonl')
    actual, counts, references = D.paired_readout(saved, previous, source_labels, context, logistic)
    if (actual != labelled or counts != report['selected_counts_per_species'] or
            references != report['eligible_reference_pairs'] or sum(v['eligible'] for v in references.values()) != 2916 or
            sum(v['five_or_more'] for v in references.values()) != 1086 or
            any(counts[s]['A31_A3'] != baseline['selected_counts_per_species'][s]['A31_A3'] for s in counts)):
        raise ValueError('Actual labels/DP/counts/denominators/stronger baseline differ')
    checks = D.gate(counts)
    if checks != report['predeclared_TRAIN_gate'] or all(checks.values()) != report['TRAIN_preparation_supported']:
        raise ValueError('Actual original five-condition arithmetic differs')
    D.R.C.R.write_summary(D.CHECK, dict(completed=True, experiment='M28-A39-NULL-ANCHORED-SET-RESULT-CHECK',
        protocol=D.PROTOCOL, configuration=D.CONFIG, source_jobs=dict(fit=D.FIT_JOB, readout=str(readout_job)),
        windows=1491, natural_candidates=16321, reference_injected_predictions=0,
        actual_original_draw_weights_counts_and_4608_updates_verified=True,
        actual_saved_fixed_composite_loss_and_importance_arithmetic_verified=True,
        actual_three_model_checkpoint_tensors_and_branch_optimizer_exposure_verified=True,
        actual_saved_natural_identity_and_FP32_score_components_verified=True,
        actual_saved_labels_DP_baseline_denominators_and_five_conditions_recomputed=True,
        original_replayed_per_pass_counts=replay, actual_optimizer_steps_per_parameter=exposure,
        selected_counts_per_species=counts, eligible_reference_pairs=references,
        predeclared_TRAIN_gate=checks, TRAIN_preparation_supported=all(checks.values()),
        model_or_backbone_forward=False, independent_neural_forward_replay=False, optimizer_steps_this_run=0,
        new_candidate_search=False, new_R4_labelling=False, DEV_used=False, test_setaria_used=False,
        A33_outputs_opened=False, seconds=time.perf_counter()-began,
        limitation='Exact saved FP32 component arithmetic and original TRAIN decisions, not independent '
            'neural replay, biological causality, genome geneF1 or independent generalization.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--readout-job', required=True)
    parser.add_argument('--fit-job', required=True)
    args = parser.parse_args(); D.FIT_JOB = args.fit_job
    torch.set_num_threads(2); main(args.readout_job)
