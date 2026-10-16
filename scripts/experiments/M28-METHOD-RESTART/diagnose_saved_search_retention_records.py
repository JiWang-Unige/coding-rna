#!/usr/bin/env python3
"""Read-only A30 supervision-location counts from already reconciled fit logs."""
from collections import Counter, defaultdict
import json
import math
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[3]
FIT = ROOT/'outputs/M28-A30-SEARCH-RETENTION-FIT'
CHECK = ROOT/'outputs/M28-A30-SEARCH-RETENTION-RESULT-CHECK-RECOVERY-R1'
OUT = ROOT/'outputs/M28-A30-SAVED-SEARCH-LOSS-DIAGNOSTIC'
SOURCE_JOBS = ('13399750', '13399752', '13399755')


def main():
    receipts = subprocess.run(['sacct', '-X', '-j', ','.join(SOURCE_JOBS), '-n', '-P',
        '--format=JobID%30,State%30,ExitCode'], check=True, capture_output=True, text=True)
    actual = {}
    for line in receipts.stdout.splitlines():
        parts = [p.strip() for p in line.split('|')]
        if len(parts) >= 3 and parts[0] in SOURCE_JOBS:
            actual[parts[0]] = parts[1:3]
    if any(actual.get(job) != ['COMPLETED', '0:0'] for job in SOURCE_JOBS):
        raise ValueError('Actual closed reconciliation/readout/CPU roots required')
    if ((CHECK/'STATUS').read_text().strip() != 'COMPLETED' or
            (CHECK/'JOBID').read_text().strip() != '13399755' or OUT.exists()):
        raise ValueError('Closed original source identity and unique diagnostic output required')
    checked = json.loads((CHECK/'summary.json').read_text())
    if (not checked['completed'] or checked['arm'] != 'A30' or checked['TRAIN_generator_supported'] or
            not checked['actual_current_search_saved_negative_qualification_and_margin_recomputed'] or
            not checked['actual_4608_draw_weight_and_A3_supervision_identity_verified'] or
            checked['source_fit_root_state'] != 'TIMEOUT' or
            checked['source_fit_artifact_reconciliation_job'] != '13399750' or
            checked['actual_TRAIN_optimizer_updates_this_recovery'] != 0 or checked['test_setaria_used']):
        raise ValueError('Accepted closed failed A30 source-only diagnostic required')
    fitted = json.loads((FIT/'summary.json').read_text())
    if not fitted['completed'] or fitted['optimizer_steps'] != 4608:
        raise ValueError('All original logical training updates required')
    cells = defaultdict(lambda: dict(steps=0, positive_query_records=0, margin_sum_over_steps=0.,
        dispositions=Counter(), first_loss_stages_all=Counter(), safe_first_loss_stages=Counter(),
        safe_target_closed_exon_counts=Counter(), safe_full_reference_exon_counts=Counter(),
        all_positive_full_reference_exon_counts=Counter()))
    epoch_sums = Counter(); dispositions = Counter(); records = steps = safe = 0
    with (FIT/'training.jsonl').open() as handle:
        for ordinal, line in enumerate(handle):
            row = json.loads(line); epoch = ordinal//1536+1
            if row['step'] != ordinal+1 or row['epoch'] != epoch:
                raise ValueError('Complete ordered original step/epoch population required')
            cell = cells[str(epoch)+'|'+row['species']]
            cell['steps'] += 1; steps += 1
            cell['margin_sum_over_steps'] += row['search_retention_margin']
            epoch_sums[epoch] += row['search_retention_margin']
            collection = row['collection_records']
            if len(collection) != row['positive_queries']:
                raise ValueError('Original assigned-positive-query population differs')
            for record in collection:
                disposition = record['supervision_disposition']
                cell['positive_query_records'] += 1; records += 1
                cell['dispositions'][disposition] += 1; dispositions[disposition] += 1
                cell['all_positive_full_reference_exon_counts'][str(len(record['chain']))] += 1
                loss = record.get('first_loss')
                if loss is not None:
                    cell['first_loss_stages_all'][loss['stage']] += 1
                if disposition == 'known_nonmatch_margin':
                    if loss is None or record['cutoff_annotation']['label'] != 0:
                        raise ValueError('Accepted saved known-nonmatch record required')
                    safe += 1; cell['safe_first_loss_stages'][loss['stage']] += 1
                    cell['safe_target_closed_exon_counts'][str(len(record['target_node']['chain']))] += 1
                    cell['safe_full_reference_exon_counts'][str(len(record['chain']))] += 1
    if (steps != 4608 or records != sum(checked['positive_query_dispositions'].values()) or
            dict(dispositions) != checked['positive_query_dispositions'] or
            safe != checked['safe_known_nonmatch_pairs']):
        raise ValueError('Closed CPU population and qualification totals must reproduce')
    for epoch in fitted['epochs']:
        mean = epoch_sums[epoch['epoch']]/epoch['steps']
        if not math.isclose(mean, epoch['mean_losses']['search_retention_margin'], rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError('Original epoch mean margin must reproduce')
    for cell in cells.values():
        cell['mean_margin_over_steps'] = cell.pop('margin_sum_over_steps')/cell['steps']
    report = dict(completed=True, experiment='M28-A30-SAVED-SEARCH-LOSS-DIAGNOSTIC',
        source_jobs=SOURCE_JOBS, original_fit_job='13398266', original_fit_root_state='TIMEOUT',
        original_optimizer_steps=4608, positive_query_records=records, safe_known_nonmatch_records=safe,
        by_epoch_species=dict(cells), original_epoch_mean_margins={str(e['epoch']):
            e['mean_losses']['search_retention_margin'] for e in fitted['epochs']},
        actual_closed_query_population_qualification_totals_and_epoch_means_reproduced=True,
        model_forward=False, new_search=False, new_R4_labelling=False, optimizer_updates_this_run=0,
        DEV_outputs_opened=False, test_setaria_used=False, original_gate_or_stop_changed=False,
        interpretation='Saved repeated TRAIN query observations, not independent genes or endpoint performance. '
            'Target closed-exon count is the recorded node history, not a claim that its chain completes. '
            'Stage coverage or falling loss does not identify one causal defect or rescue the original FAIL.')
    OUT.mkdir(exist_ok=False)
    (OUT/'JOBID').write_text(os.environ['SLURM_JOB_ID']+'\n')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False); handle.write('\n')
    (OUT/'STATUS').write_text('COMPLETED\n')
    print(json.dumps(report, allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
