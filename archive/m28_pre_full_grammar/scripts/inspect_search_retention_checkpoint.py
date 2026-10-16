#!/usr/bin/env python3
"""Read-only A30 checkpoint preparation, never a training/recovery command."""
from collections import Counter
import json
import math
from pathlib import Path
import subprocess
import torch
import run_query_search_retention as R

SOURCE = R.FIT
SOURCE_JOB = '13398266'
OUT = R.ROOT/'outputs/M28-A30-CHECKPOINT-PREPARATION'
INSPECTION = SOURCE/'step_001024.pt'


def root_receipt(text):
    roots = [line.strip().split('|') for line in text.splitlines()
             if line.strip().split('|')[0] == SOURCE_JOB]
    if len(roots) != 1 or len(roots[0]) < 3:
        raise ValueError('Authoritative original A30 root receipt required')
    _, state, exit_code = roots[0][:3]
    return dict(job=SOURCE_JOB, state=state.strip(), exit_code=exit_code.strip())


def source_receipt():
    if (SOURCE/'JOBID').read_text().strip() != SOURCE_JOB:
        raise ValueError('Original source job differs')
    result = subprocess.run(['sacct', '-j', SOURCE_JOB, '-n', '-P',
        '--format=JobID,State,ExitCode'], check=True, capture_output=True, text=True)
    return root_receipt(result.stdout)


def require_source_timeout():
    receipt = source_receipt()
    if receipt['state'] != 'TIMEOUT':
        raise ValueError('Only an actual source TIMEOUT permits recovery consideration')
    return receipt


def latest_timeout_checkpoint():
    """Inspect highest saved step first; completeness is a separate check, not a metric."""
    require_source_timeout()  # Must precede examining or selecting source checkpoints.
    paths = list(SOURCE.glob('step_*.pt'))
    if not paths:
        raise ValueError('No persisted A30 model/optimizer checkpoint')
    return max(paths, key=lambda p: int(p.stem.split('_')[1]))


def load_optimizer(model, saved, logical_step):
    optimizer = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=.01,
                                 betas=(.9, .999), eps=1e-8)
    optimizer.load_state_dict(saved)
    if len(optimizer.param_groups) != 1:
        raise ValueError('Original single AdamW parameter group required')
    group = optimizer.param_groups[0]
    expected = dict(lr=3e-4, weight_decay=.01, betas=(.9, .999), eps=1e-8,
                    amsgrad=False, maximize=False)
    if any(group[k] != v for k, v in expected.items()):
        raise ValueError('Frozen AdamW hyperparameters differ')
    if logical_step > 0 and not optimizer.state:
        raise ValueError('Persisted AdamW moments required, not fresh optimizer')
    for parameter, state in optimizer.state.items():
        if set(state) != {'step', 'exp_avg', 'exp_avg_sq'}:
            raise ValueError('Original AdamW moment state required')
        step = float(state['step'])
        if not math.isfinite(step) or not step.is_integer() or not 0 <= step <= logical_step:
            raise ValueError('AdamW applied step is inconsistent with checkpoint')
        for key in ('exp_avg', 'exp_avg_sq'):
            value = state[key]
            if value.shape != parameter.shape or value.device != parameter.device or not bool(torch.isfinite(value).all()):
                raise ValueError('AdamW moment shape, device or finiteness differs')
        if not bool((state['exp_avg_sq'] >= 0).all()):
            raise ValueError('Invalid AdamW second moment')
    return optimizer


def restore_checkpoint(path):
    state = torch.load(path, map_location='cpu', weights_only=True)
    step = int(path.stem.split('_')[1])
    if (state['arm'] != 'A30' or state['configuration'] != R.CONFIG or
            state['protocol'] != R.PROTOCOL or state['revision'] != R.A.F.REVISION or
            state['step'] != step or not 0 < step <= 4608 or step % 256):
        raise ValueError('Matching original A30 engineering checkpoint required')
    model = R.new_model('cpu')
    model.load_state_dict(state['model'], strict=True)
    if sum(p.numel() for p in model.parameters()) != 649457 or any(
            not bool(torch.isfinite(v).all()) for v in model.state_dict().values()):
        raise ValueError('Original finite complete A30 model required')
    # Source saves no RNG state. Current fixed path has no stochastic train module.
    for module in model.modules():
        if ((isinstance(module, torch.nn.Dropout) and module.p != 0.) or
                (isinstance(module, (torch.nn.MultiheadAttention, torch.nn.RNNBase)) and module.dropout != 0.)):
            raise ValueError('Unsaved stochastic training state cannot be assumed restored')
    optimizer = load_optimizer(model, state['optimizer'], step)
    return model, optimizer, step


def validated_prefix(path, step, draws, windows):
    """Use only checkpoint-matching rows; a discarded trailing partial row is irrelevant."""
    with path.open() as handle:
        for i in range(step):
            line = handle.readline()
            if not line:
                raise ValueError('Checkpoint prefix is missing logged updates')
            row = json.loads(line)
            epoch, ordinal = divmod(i, 1536)
            draw = draws[ordinal]; window = windows[draw['window_id']]
            if ((row['step'], row['epoch'], row['draw_ordinal'], row['window_id'], row['species']) !=
                    (i+1, epoch+1, ordinal, window['id'], window['species']) or
                    row['importance_weight'] != draw['importance_weight_uniform_within_stratum'] or
                    window['split'] != 'train' or window['seqid'] not in R.A.ALLOW[window['species']] or
                    row['reference_injected_candidates'] != 0):
                raise ValueError('Original fixed TRAIN prefix identity or weight differs')
            if set(row['losses']) != {'local', 'existence_bce', 'boundary_ce', 'termination_bce'} or any(
                    not math.isfinite(v) for v in (*row['losses'].values(), row['weighted_total'],
                        row['weighted_original_total'], row['search_retention_margin'], row['gradient_norm_before_clip'])):
                raise ValueError('Finite original four-loss prefix required')
            yield row


def inspect():
    if (SOURCE/'JOBID').read_text().strip() != SOURCE_JOB:
        raise ValueError('Original source job differs')
    model, optimizer, step = restore_checkpoint(INSPECTION)
    draws = R.Q.read_rows(R.Q.DRAWS)
    ids = {d['window_id'] for d in draws}; windows = {}
    for species in R.A.ALLOW:
        windows.update({w['id']: w for w in R.Q.train_records(R.Q.LABELS/species/'windows_train.jsonl')
                        if w['id'] in ids})
    if len(draws) != 1536 or len(ids) != 1491 or set(windows) != ids:
        raise ValueError('Original TRAIN draw/window metadata required')
    counts = Counter(); next_draw = draws[step % 1536]
    for row in validated_prefix(SOURCE/'training.jsonl', step, draws, windows):
        counts[row['species']] += 1
    applied = [int(s['step']) for s in optimizer.state.values()]
    receipt = source_receipt()
    R.C.R.write_summary(OUT, dict(completed=True, arm='A30_checkpoint_preparation',
        configuration=R.CONFIG, protocol=R.PROTOCOL, source_job=SOURCE_JOB, source_receipt=receipt,
        inspected_checkpoint=str(INSPECTION.relative_to(R.ROOT)), inspected_logical_step=step,
        checkpoint_matching_TRAIN_prefix_records=sum(counts.values()), prefix_draws_by_species=dict(counts),
        next_original_draw=dict(global_step=step+1, epoch=step//1536+1,
            draw_ordinal=step % 1536, window_id=next_draw['window_id'],
            importance_weight=next_draw['importance_weight_uniform_within_stratum']),
        original_model_structure_and_tensors_finite=True, trainable_parameters=649457,
        original_AdamW_loaded_shapes_moments_and_hyperparameters_verified=True,
        AdamW_parameters_with_state=len(applied), AdamW_applied_step_range=[min(applied), max(applied)],
        original_draw_weight_and_prefix_identity_verified=True,
        source_RNG_state_not_saved=True, current_training_modules_have_no_enabled_dropout=True,
        actual_TRAIN_optimizer_steps_this_run=0, actual_model_forward=False,
        actual_candidate_search=False, new_backbone_forward=False, DEV_forward_used=False, test_setaria_used=False,
        source_artifacts_modified=False, source_JOBID_STATUS_or_resources_modified=False,
        recovery_training_implemented=False, recovery_training_submitted=False,
        inspection_checkpoint_is_future_recovery_selection=False,
        limitation='Read-only preparation, not a completed recovery, bitwise uninterrupted-GPU proof or science gate. '
            'If the source actually times out, inspect the latest persisted own checkpoint and matching prefix anew; '
            'never use this older inspection checkpoint as a metric-based fallback.'))


if __name__ == '__main__':
    torch.set_num_threads(2)
    inspect()
