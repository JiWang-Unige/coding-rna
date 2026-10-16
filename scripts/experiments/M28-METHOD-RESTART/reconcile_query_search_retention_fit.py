#!/usr/bin/env python3
"""A30 full final artifacts under root TIMEOUT: reconcile, then read out only."""
import argparse
import json
import subprocess
import time
import torch
import resume_query_search_retention as S

K, R, V = S.K, S.R, S.V
OUT = R.ROOT/'outputs/M28-A30-TIMEOUT-FIT-RECONCILIATION'
READOUT = R.ROOT/'outputs/M28-A30-SEARCH-RETENTION-TRAIN-READOUT-RECOVERY-R1'
CHECK = R.ROOT/'outputs/M28-A30-SEARCH-RETENTION-RESULT-CHECK-RECOVERY-R1'


def reconciliation_receipt(job):
    result = subprocess.run(['sacct', '-j', job, '-n', '-P',
        '--format=JobID,State,ExitCode'], check=True, capture_output=True, text=True)
    roots = [line.strip().split('|') for line in result.stdout.splitlines()
             if line.strip().split('|')[0] == job]
    if len(roots) != 1 or roots[0][1:3] != ['COMPLETED', '0:0']:
        raise ValueError('Successful actual reconciliation root required, not a batch marker')


def reconcile():
    began = time.perf_counter(); receipt = K.require_source_timeout()
    fitted = json.loads((K.SOURCE/'summary.json').read_text()); R.validate_fit(fitted)
    if ((K.SOURCE/'STATUS').read_text().strip() != 'COMPLETED' or
            fitted['supervised_checkpoint_inheritance'] or fitted['new_backbone_forward'] or
            fitted['reference_injected_candidates'] != 0 or fitted['trainable_parameters'] != 649457 or
            fitted['main_evaluation'] != 'last_step_004608_only' or len(fitted['epochs']) != 3):
        raise ValueError('Complete original fixed final fit artifacts required')
    checkpoint = K.SOURCE/'step_004608.pt'
    model, optimizer, step = K.restore_checkpoint(checkpoint)
    if step != 4608:
        raise ValueError('Only original final step 004608, not a partial recovery')
    for saved_step in range(256, 4609, 256):
        if (K.SOURCE/('step_'+str(saved_step).zfill(6)+'.pt')).stat().st_size == 0:
            raise ValueError('Original saved checkpoint prefix missing')
    draws = R.Q.read_rows(R.Q.DRAWS); ids = {d['window_id'] for d in draws}; windows = {}
    for species in R.A.ALLOW:
        windows.update({w['id']: w for w in R.Q.train_records(R.Q.LABELS/species/'windows_train.jsonl')
                        if w['id'] in ids})
    if len(draws) != 1536 or len(ids) != 1491 or set(windows) != ids:
        raise ValueError('Original TRAIN draw/window metadata required')
    with (K.SOURCE/'training.jsonl').open() as handle:
        if sum(1 for _ in handle) != 4608:
            raise ValueError('Exactly 4608 logged original updates required')
    saved = S.accumulator()
    for row in K.validated_prefix(K.SOURCE/'training.jsonl', step, draws, windows):
        S.accumulate(saved, row)
    original = json.loads((R.A.STAGE_DIRS['fit']/'summary.json').read_text())
    R.A.check_last_identity(original)
    if saved['total'] != 4608 or saved['steps'] != [1536]*3 or len(original['epochs']) != 3:
        raise ValueError('Original three-epoch exposure required')
    for i, epoch in enumerate(fitted['epochs']):
        if (epoch['epoch'] != i+1 or epoch['steps'] != 1536 or
                epoch['checkpoint'] != str((K.SOURCE/('step_'+str((i+1)*1536).zfill(6)+'.pt')).relative_to(R.ROOT)) or
                epoch['supervision_by_species'] != {s: dict(c) for s, c in saved['counts'][i].items()} or
                epoch['supervision_by_species'] != original['epochs'][i]['supervision_by_species'] or
                epoch['mean_losses'] != {k: v/1536 for k, v in saved['losses'][i].items()}):
            raise ValueError('Original final summary and fixed ordered log disagree')
    if (dict(saved['margin_steps']) != fitted['positive_search_margin_steps_by_species'] or
            any(saved['margin_steps'][s] == 0 for s in R.A.ALLOW)):
        raise ValueError('Original both-species search-signal accounting differs')
    applied = [int(s['step']) for s in optimizer.state.values()]
    R.write(OUT, dict(completed=True, arm='A30_timeout_fit_reconciliation',
        source_job=K.SOURCE_JOB, source_root_receipt=receipt, source_root_state_not_reclassified=True,
        original_final_fit_artifacts_reconciled=True,
        checkpoint=str(checkpoint.relative_to(R.ROOT)), optimizer_steps_in_original_fit=4608,
        actual_4608_draw_weight_and_epoch_identity_verified=True,
        original_model_structure_and_all_tensors_finite=True, trainable_parameters=649457,
        original_AdamW_shapes_moments_devices_hyperparameters_verified=True,
        AdamW_parameters_with_state=len(applied), AdamW_applied_step_range=[min(applied), max(applied)],
        source_STATUS_and_JOBID_modified=False, source_artifacts_modified=False,
        actual_TRAIN_optimizer_updates_this_run=0, actual_model_forward=False, new_candidate_search=False,
        new_backbone_forward=False, DEV_forward_used=False, test_setaria_used=False,
        seconds=time.perf_counter()-began,
        limitation='Full final artifact reconciliation under actual root TIMEOUT; not a successful '
            'Slurm fit root, training retry, natural-pool gate, selected gene F1 or generalization.'))


def require_reconciled(job):
    K.require_source_timeout(); reconciliation_receipt(job)
    if ((OUT/'JOBID').read_text().strip() != job or
            (OUT/'STATUS').read_text().strip() != 'COMPLETED'):
        raise ValueError('Matching completed reconciliation artifacts required')
    report = json.loads((OUT/'summary.json').read_text())
    if (not report['completed'] or not report['original_final_fit_artifacts_reconciled'] or
            report['source_job'] != K.SOURCE_JOB or
            report['source_root_receipt']['state'] != 'TIMEOUT' or
            report['configuration'] != R.CONFIG or report['protocol'] != R.PROTOCOL or
            report['optimizer_steps_in_original_fit'] != 4608 or
            report['checkpoint'] != str((K.SOURCE/'step_004608.pt').relative_to(R.ROOT)) or
            report['actual_TRAIN_optimizer_updates_this_run'] != 0):
        raise ValueError('Original final checkpoint reconciliation differs')


def readout(job):
    require_reconciled(job)
    R.READOUT = READOUT  # Keep R.FIT bound to the original completed artifacts.
    R.readout()


def check(job, readout_job):
    require_reconciled(job); R.READOUT, V.OUT = READOUT, CHECK
    V.main(K.SOURCE_JOB, readout_job)
    path = CHECK/'summary.json'; result = json.loads(path.read_text())
    result.update(source_fit_root_state='TIMEOUT', source_fit_artifact_reconciliation_job=job,
        actual_TRAIN_optimizer_updates_this_recovery=0, source_fit_root_state_not_reclassified=True,
        source_artifacts_modified=False)
    with path.open('w') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('reconcile', 'readout', 'check'), required=True)
    parser.add_argument('--reconciliation-job'); parser.add_argument('--readout-job')
    args = parser.parse_args(); torch.set_num_threads(2)
    if args.mode == 'reconcile':
        reconcile()
    elif not args.reconciliation_job:
        raise ValueError('Actual CPU reconciliation job ID required')
    elif args.mode == 'readout':
        readout(args.reconciliation_job)
    elif not args.readout_job:
        raise ValueError('Actual recovery readout job ID required')
    else:
        check(args.reconciliation_job, args.readout_job)
