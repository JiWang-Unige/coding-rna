#!/usr/bin/env python3
"""CPU reconciliation of actual A36 exposure, saved logits and fixed five gates."""
import argparse
from collections import Counter
import json
import math
import time
import torch
import readout_structured_emission_train as D


@torch.inference_mode()
def verify_direct_scores(logits, chains, expected):
    if logits.ndim != 2 or logits.shape[1] != 15 or logits.dtype != torch.float32 or not bool(torch.isfinite(logits).all()):
        raise ValueError('Finite saved valid-base FP32 emission matrix required')
    logp = logits.double().log_softmax(-1)
    positions = torch.arange(len(logp))
    values = [float((logp[positions, D.R.candidate_roles(c, len(logp))]-logp[:, 0]).sum()) for c in chains]
    torch.testing.assert_close(torch.tensor(values, dtype=torch.float64),
        torch.tensor(expected, dtype=torch.float64), rtol=1e-12, atol=1e-8)
    return max((abs(a-b) for a, b in zip(values, expected)), default=0.)


def verify_fit():
    fitted = D.R.completed(D.F.OUT, D.FIT_JOB); D.validate_fit(fitted)
    draws, _, _, records, old_rows = D.F.A.A14.population()
    if len(draws) != 1536 or len(records) != 1491 or len(old_rows) != 16857:
        raise ValueError('Unchanged complete original training population required')
    train = D.R.C.R.rows(D.F.OUT/'training.jsonl')
    if len(train) != 4608 or len(fitted['epochs']) != 3:
        raise ValueError('Actual complete fixed three-pass training required')
    expected_counts = Counter()
    for draw in draws:
        expected_counts.update(records[draw['window_id']]['counts'])
    for i, row in enumerate(train):
        epoch, ordinal = divmod(i, 1536); draw = draws[ordinal]; rec = records[draw['window_id']]
        if (row['step'] != i+1 or row['epoch'] != epoch+1 or row['draw_ordinal'] != ordinal or
                row['window_id'] != draw['window_id'] or row['species'] != draw['species'] or
                row['importance_weight'] != draw['importance_weight_uniform_within_stratum'] or
                any(row[k] != v for k, v in rec['counts'].items()) or
                any(not math.isfinite(row[k]) for k in ('weighted_loss', 'log_partition', 'gold_gain', 'gradient_norm_before_clip'))):
            raise ValueError('Actual draw/order/weight/exposure or finite update record differs')
    model = D.F.A.new_model('cpu')
    for epoch, step in enumerate((1536, 3072, 4608)):
        path = D.F.OUT/('step_'+str(step).zfill(6)+'.pt')
        report = fitted['epochs'][epoch]
        if (report['epoch'] != epoch+1 or report['steps'] != 1536 or report['counts'] != dict(expected_counts) or
                report['checkpoint'] != str(path.relative_to(D.ROOT)) or not path.stat().st_size):
            raise ValueError('Actual saved epoch/checkpoint exposure differs')
        state = torch.load(path, map_location='cpu', weights_only=True)
        if (state['arm'] != 'A36' or state['configuration'] != D.CONFIG or state['protocol'] != D.F.PROTOCOL or
                state['step'] != step or state['revision'] != D.F.A.C.I.F.REVISION or
                sum(v.numel() for v in state['model'].values()) != 366191 or
                any(v.dtype != torch.float32 or not bool(torch.isfinite(v).all()) for v in state['model'].values())):
            raise ValueError('Finite matching actual fixed checkpoint structure required')
        model.load_state_dict(state['model'], strict=True)
        opt = state['optimizer']; group = opt['param_groups']
        if (len(group) != 1 or len(group[0]['params']) != len(opt['state']) or
                len(opt['state']) != len(list(model.parameters())) or
                group[0]['lr'] != D.CONFIG['lr'] or group[0]['weight_decay'] != D.CONFIG['weight_decay'] or
                list(group[0]['betas']) != D.CONFIG['betas'] or group[0]['eps'] != D.CONFIG['eps'] or
                any(float(v['step']) != step or not bool(torch.isfinite(v['exp_avg']).all()) or
                    not bool(torch.isfinite(v['exp_avg_sq']).all()) for v in opt['state'].values())):
            raise ValueError('Actual complete finite fixed AdamW exposure required')
        del state, opt
    return dict(expected_counts)


@torch.inference_mode()
def main(readout_job):
    began = time.perf_counter(); D.F.A.require_jobs([D.FIT_JOB, str(readout_job)])
    report = D.R.completed(D.OUT, readout_job)
    if (report['protocol'] != D.PROTOCOL or report['configuration'] != D.CONFIG or
            report['source_fit_job'] != D.FIT_JOB or report['windows'] != 1491 or
            report['natural_candidates'] != 16321 or report['reference_injected_predictions'] != 0 or
            not report['natural_score_file_closed_before_saved_label_arrays_open'] or
            not report['original_natural_geometry_order_sources_unchanged'] or not report['original_saved_labels_reused'] or
            not report['model_weights_unchanged'] or not report['all_parameters_frozen_and_gradients_absent'] or
            report['optimizer_steps_this_run'] != 0 or report['label_catalog_opened'] or report['new_R4_labelling'] or
            report['new_candidate_search'] or report['new_backbone_forward'] or report['DEV_used'] or
            report['test_setaria_used'] or report['A33_outputs_opened']):
        raise ValueError('Actual complete unchanged natural-only A36 readout required')
    replay = verify_fit()
    draws, free, prior, _, old = D.population()
    saved = D.R.C.R.rows(D.OUT/'natural_scores.jsonl')
    labelled = D.R.C.R.rows(D.OUT/'windows.jsonl')
    source_labels = D.R.C.R.rows(D.R.READOUT/'windows.jsonl')
    if len(saved) != 1491 or len(labelled) != 1491 or len(list((D.OUT/'logits').glob('*.pt'))) != 1491:
        raise ValueError('All actual 1491 natural records and saved logit matrices required')
    error = 0.
    for ordinal, (row, pool) in enumerate(zip(saved, free)):
        path = D.OUT/'logits'/('%04d.pt' % ordinal)
        if (row['window'] != pool['window'] or row['chains'] != pool['chains'] or
                row['decoder_sources'] != pool['decoder_sources'] or
                row['logits_artifact'] != str(path.relative_to(D.ROOT)) or not path.stat().st_size):
            raise ValueError('Actual original natural identity/order/sources or saved artifact differs')
        tensor = torch.load(path, map_location='cpu', weights_only=True)
        if tensor['window'] != row['window'] or tuple(tensor['logits'].shape) != (row['window']['valid_bases'], 15):
            raise ValueError('Actual saved valid-base emission geometry differs')
        error = max(error, verify_direct_scores(tensor['logits'], row['chains'], row['scores']))
        del tensor
        if (ordinal+1) % 128 == 0 or ordinal+1 == len(saved):
            print(json.dumps(dict(verified_windows=ordinal+1, total=1491)), flush=True)
    actual, counts, references = D.paired_readout(saved, prior, source_labels)
    if (actual != labelled or counts != report['selected_counts_per_species'] or
            references != report['eligible_reference_pairs'] or sum(v['eligible'] for v in references.values()) != 2916 or
            any(counts[s]['A31_A3'] != old['selected_counts_per_species'][s]['A31_A3'] for s in counts)):
        raise ValueError('Saved labels/selection/denominators or stronger baseline do not reconcile')
    checks = D.gate(counts)
    if checks != report['predeclared_TRAIN_gate'] or all(checks.values()) != report['TRAIN_preparation_supported']:
        raise ValueError('Actual original five-condition arithmetic differs')
    D.R.C.R.write_summary(D.CHECK, dict(completed=True, experiment='M28-A36-STRUCTURED-EMISSION-RESULT-CHECK',
        protocol=D.PROTOCOL, configuration=D.CONFIG, source_jobs=dict(fit=D.FIT_JOB, readout=str(readout_job)),
        windows=1491, natural_candidates=16321, reference_injected_predictions=0,
        actual_original_training_draw_weights_and_4608_updates_verified=True,
        actual_three_checkpoint_tensors_and_optimizer_exposure_verified=True,
        actual_saved_natural_identity_and_complete_logits_verified=True,
        actual_saved_full_window_log_probability_scores_verified=True,
        actual_saved_labels_selections_stronger_baseline_and_five_conditions_recomputed=True,
        original_replayed_per_pass_counts=replay, max_direct_log_probability_error=error,
        selected_counts_per_species=counts, eligible_reference_pairs=references,
        predeclared_TRAIN_gate=checks, TRAIN_preparation_supported=all(checks.values()),
        model_or_backbone_forward=False, optimizer_steps_this_run=0, new_candidate_search=False,
        new_R4_labelling=False, DEV_used=False, test_setaria_used=False, A33_outputs_opened=False,
        seconds=time.perf_counter()-began,
        limitation='CPU verifies saved new-model emissions and the original TRAIN window/chain decision; '
            'no independent neural forward, gene F1, biological causality or generalization.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--readout-job', required=True)
    args = parser.parse_args(); torch.set_num_threads(2); main(args.readout_job)
