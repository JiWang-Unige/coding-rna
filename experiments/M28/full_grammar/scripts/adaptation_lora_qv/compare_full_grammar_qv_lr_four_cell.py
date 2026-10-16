#!/usr/bin/env python3
"""Fixed four-cell TRAIN LR contrasts; completed JSON receipts only, no model/data reads."""
import argparse
import json
import math
import os
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[3]
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_qv_lr_controls_plan.md'
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-QV-LR-FOUR-CELL'
CHECKS = ROOT/'outputs/M28-FULL-GRAMMAR-QV-LR-FOUR-CELL-CHECKS'
SPECIES = ('arabidopsis_thaliana', 'oryza_sativa')
SCOPE = (
    (SPECIES[0], 'NC_003070.9', 30427671),
    (SPECIES[0], 'NC_003071.7', 19698289),
    (SPECIES[0], 'NC_003076.8', 26975502),
    (SPECIES[1], 'NC_089043.1', 27474823),
    (SPECIES[1], 'NC_089044.1', 23931887),
    (SPECIES[1], 'NC_089046.1', 28271460))
ENDPOINTS = ('exact_chain', 'CDS_base_unstranded_union')
PRF = ('precision', 'recall', 'f1')
BG_COUNTS = ('bp', 'predicted_span_bp', 'predicted_CDS_bp', 'wholly_background_unique_chains')
RECORD_COUNTS = ('prediction_records', 'duplicate_chain_records', 'exact_other_assessable_isoform_not_primary_TP')
CELLS = {
    'lora_low': ('JOINT-LORA-LR1E5', '13545578', 1e-5, PROTOCOL),
    'lora_high': ('JOINT-LORA', '13499817', 3e-4, 'reports/M28-METHOD-RESTART/full_grammar_joint_lora_fit.md'),
    'dense_low': ('DENSE-QV', '13510324', 1e-5, 'reports/M28-METHOD-RESTART/full_grammar_dense_qv_plan.md'),
    'dense_high': ('DENSE-QV-LR3E4', '13545584', 3e-4, PROTOCOL)}
FIXED_CONFIG = ('encoder_model', 'encoder_revision', 'encoder_blocks_bp', 'encoder_token_bp',
    'seed', 'draws_per_pass', 'unique_windows', 'passes', 'fine', 'coarse', 'unknown', 'epsilon',
    'max_exons', 'score_accumulation', 'lr', 'weight_decay', 'betas', 'eps', 'clip_norm',
    'starting_fit_job', 'starting_optimizer_steps', 'additional_optimizer_steps',
    'sampling_seed', 'preparation_job', 'support_job')
LIMITATION = ('Single fixed fit/seed per cell, original six TRAIN only; nonlinear F1 and distinct '
    'parameterization/rounding/joint clipping/optimizer trajectories prevent a pure rank causal '
    'claim. No DEV selection, independent generalization, or publication-success claim.')


def close(a, b):
    return math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12)


def check_row(row):
    for endpoint in ENDPOINTS:
        r = row[endpoint]
        tp, predicted, reference = (r[k] for k in ('tp', 'predicted', 'reference'))
        if not 0 <= tp <= min(predicted, reference):
            raise ValueError('Invalid endpoint counts')
        expected = (tp/predicted if predicted else 0., tp/reference if reference else 0.,
                    2*tp/(predicted+reference) if predicted+reference else 0.)
        if any(not close(r[k], v) for k, v in zip(PRF, expected)):
            raise ValueError('PRF does not reconstruct counts')
    bg = row['background']
    if (bg['bp'] <= 0 or not 0 <= bg['predicted_CDS_bp'] <= bg['predicted_span_bp'] <= bg['bp']
            or bg['wholly_background_unique_chains'] < 0
            or not close(bg['span_FPR'], bg['predicted_span_bp']/bg['bp'])
            or not close(bg['wholly_background_chains_per_Mb'], bg['wholly_background_unique_chains']*1e6/bg['bp'])):
        raise ValueError('Background counts, denominator, or rates differ')
    chain = row['exact_chain']
    if (row['prediction_records'] != chain['predicted']+row['duplicate_chain_records']
            or row['duplicate_chain_records'] < 0
            or not 0 <= row['exact_other_assessable_isoform_not_primary_TP'] <= chain['predicted']-chain['tp']):
        raise ValueError('Primary/nonprimary/duplicate accounting differs')


def aggregate_matches(row, parts):
    for endpoint in ENDPOINTS:
        if any(row[endpoint][k] != sum(p[endpoint][k] for p in parts) for k in ('tp', 'predicted', 'reference')):
            raise ValueError('Aggregate endpoint counts differ from chromosomes')
    if (any(row['background'][k] != sum(p['background'][k] for p in parts) for k in BG_COUNTS)
            or any(row[k] != sum(p[k] for p in parts) for k in RECORD_COUNTS)):
        raise ValueError('Aggregate background or record counts differ from chromosomes')


def measures(row, phase_count):
    values = {f'{e}.{k}': row[e][k] for e in ENDPOINTS for k in ('tp', 'predicted', 'reference')}
    values.update({f'{e}.{k}_percent': 100*row[e][k] for e in ENDPOINTS for k in PRF})
    values.update({f'background.{k}': row['background'][k] for k in BG_COUNTS})
    values['background.span_FPR_percent'] = 100*row['background']['span_FPR']
    values['background.wholly_background_chains_per_Mb'] = row['background']['wholly_background_chains_per_Mb']
    values.update({k: row[k] for k in RECORD_COUNTS})
    values['phase_checked_prediction_records'] = phase_count
    values['nonprimary_unique_chains'] = row['exact_chain']['predicted']-row['exact_chain']['tp']
    values['neither_primary_nor_assessable_other_isoform_unique_chains'] = (
        values['nonprimary_unique_chains']-row['exact_other_assessable_isoform_not_primary_TP'])
    return values


def contrasts(arms):
    lora = {k: arms['lora_high'][k]-arms['lora_low'][k] for k in arms['lora_low']}
    dense = {k: arms['dense_high'][k]-arms['dense_low'][k] for k in arms['dense_low']}
    return dict(lora_high_minus_low=lora, dense_high_minus_low=dense,
                difference_in_differences={k: dense[k]-lora[k] for k in lora})


def summarize(reports):
    if set(reports) != set(CELLS):
        raise ValueError('All four declared cells required, including negative results')
    first = reports['lora_low']
    seqids = {seqid for _, seqid, _ in SCOPE}
    for report in reports.values():
        if (report['reference_policy'] != first['reference_policy']
                or set(report['chromosomes']) != seqids or set(report['per_species']) != set(SPECIES)):
            raise ValueError('Common reference policy and original six TRAIN geometry required')
        for species, seqid, length in SCOPE:
            row, common = report['chromosomes'][seqid], first['chromosomes'][seqid]
            if (row['species'], row['seqid'], row['length']) != (species, seqid, length):
                raise ValueError('Chromosome geometry changed')
            if (any(row[k] != common[k] for k in ('reference_primary_records',
                    'reference_assessable_isoform_records', 'retained_unsupported_reference_strata'))
                    or any(row[e]['reference'] != common[e]['reference'] for e in ENDPOINTS)
                    or row['background']['bp'] != common['background']['bp']):
                raise ValueError('Retained reference population or background denominator changed')
            check_row(row)
            if (row['phase_checked_prediction_records'] != row['prediction_records']
                    or row['prediction_nonexact_primary'] != row['exact_chain']['predicted']-row['exact_chain']['tp']
                    or row['reference_primary_records'] != row['exact_chain']['reference']):
                raise ValueError('Phase or unique primary/nonprimary accounting differs')
        for scope in (*SPECIES, 'pooled'):
            parts = [r for r in report['chromosomes'].values() if scope == 'pooled' or r['species'] == scope]
            row = report['pooled'] if scope == 'pooled' else report['per_species'][scope]
            check_row(row)
            aggregate_matches(row, parts)
            expected_reference = 22688 if scope == 'pooled' else 17739 if scope == SPECIES[0] else 4949
            if row['exact_chain']['reference'] != expected_reference:
                raise ValueError('Keep original full TRAIN primary denominator')
        species_f1 = [report['per_species'][s]['exact_chain']['f1'] for s in SPECIES]
        if (not close(report['macro_species_exact_chain_F1'], sum(species_f1)/2)
                or not close(report['worst_species_exact_chain_F1'], min(species_f1))):
            raise ValueError('Saved macro/worst F1 differs')
    scopes = {}
    for scope in (*(seqid for _, seqid, _ in SCOPE), *SPECIES, 'pooled'):
        arms = {}
        for name, report in reports.items():
            row = report['pooled'] if scope == 'pooled' else report['per_species'][scope] if scope in SPECIES else report['chromosomes'][scope]
            parts = [r for r in report['chromosomes'].values() if scope == 'pooled' or r['species'] == scope or r['seqid'] == scope]
            arms[name] = measures(row, sum(r['phase_checked_prediction_records'] for r in parts))
        scopes[scope] = dict(arms=arms, contrasts=contrasts(arms))
    rate_keys = [f'{e}.{k}_percent' for e in ENDPOINTS for k in PRF]
    bg_keys = ('background.span_FPR_percent', 'background.wholly_background_chains_per_Mb')
    for scope in ('macro_species', 'worst_species'):
        arms = {}
        for name in CELLS:
            arms[name] = {k: sum(scopes[s]['arms'][name][k] for s in SPECIES)/2 if scope == 'macro_species'
                else min(scopes[s]['arms'][name][k] for s in SPECIES) for k in rate_keys}
            arms[name].update({k: sum(scopes[s]['arms'][name][k] for s in SPECIES)/2 if scope == 'macro_species'
                else max(scopes[s]['arms'][name][k] for s in SPECIES) for k in bg_keys})
        scopes[scope] = dict(arms=arms, contrasts=contrasts(arms))
    return scopes


def completed(root, job):
    value = json.loads((root/'summary.json').read_text())
    if ((root/'STATUS').read_text().strip() != 'COMPLETED' or (root/'JOBID').read_text().strip() != job
            or not value['completed'] or value['job_id'] != job or value['DEV_used']
            or value['sealed_data_used']):
        raise ValueError('Completed matching original TRAIN source receipt required')
    return value


def load_cell(name, evaluation_job, paired_job):
    stem, fit_job, native_lr, protocol = CELLS[name]
    base = ROOT/'outputs'/('M28-FULL-GRAMMAR-'+stem)
    fit = completed(Path(str(base)+'-FIT'), fit_job)
    evaluation = completed(Path(str(base)+'-TRAIN-DIAGNOSTIC'), evaluation_job)
    paired = completed(Path(str(base)+'-PAIRED'), paired_job)
    inputs = json.loads((Path(str(base)+'-TRAIN-DIAGNOSTIC')/'paired_input_match.json').read_text())
    config = fit['configuration']
    age_key = 'adapter_optimizer_steps' if name.startswith('lora') else 'native_qv_optimizer_steps'
    if (fit['protocol'] != protocol or fit['source_fit_job'] != '13434603'
            or fit['cumulative_head_optimizer_steps'] != 9216 or fit[age_key] != 4608
            or fit['optimizer_steps_this_run'] != 4608 or fit['fresh_draws'] != 4608
            or fit['fresh_unique_windows'] != 4219 or config['lr'] != 3e-4
            or config['sampling_seed'] != 1 or fit['cached_features_used']
            or not fit['saved_source_weights_and_optimizer_unchanged']
            or not fit['fixed_fresh_draw_order_and_unit_importance_preserved']
            or not fit['all_fixed_draws_retained'] or fit['original_targets_changed']
            or fit['context_or_structural_grammar_changed']):
        raise ValueError('Fixed original4608-to9216 scientific source differs')
    if name == 'lora_high':
        if config['optimizer_state'] != 'original_head_AdamW4608_restored_new_adapter_AdamW0_same_hyperparameters':
            raise ValueError('Historical E native group must use verified original shared hyperparameters')
    elif config['native_adapter_lr' if name.startswith('lora') else 'native_qv_lr'] != native_lr:
        raise ValueError('Native LR cell differs')
    checkpoint = str(Path(str(base)+'-FIT')/'step_009216.pt')
    if (evaluation['checkpoint'] != checkpoint or fit['final_checkpoint'] != str(Path(checkpoint).relative_to(ROOT))
            or evaluation['optimizer_steps_this_run'] or evaluation['reference_used_in_inference']
            or not evaluation['unsupported_reference_strata_retained']
            or not evaluation['training_eligibility_does_not_filter_reference']):
        raise ValueError('Strict raw final9216 evaluation source differs')
    if (inputs['job_id'] != evaluation_job or inputs['protocol'] != protocol
            or inputs['adapted_fit_job'] != fit_job or inputs['baseline_fit_job'] != '13494895'
            or not inputs['original_six_TRAIN_inputs_directly_equal']
            or not inputs['comparison_before_reference_access'] or inputs['DEV_used'] or inputs['sealed_data_used']
            or paired['source_jobs'] != dict(adapted=evaluation_job, baseline='13494898')
            or paired['source_fit_jobs'] != dict(adapted=fit_job, baseline='13494895')
            or paired['protocol'] != protocol or paired['optimizer_steps_this_run'] or paired['native_forward']
            or not paired['all_paired_caller_inputs_and_normalized_original_DNA_directly_equal']
            or not paired['DNA_comparison_before_adapted_reference_access'] or not paired['counts_and_PRF_reconciled']):
        raise ValueError('All four cells must already directly match same B and original DNA before reference access')
    return evaluation, config


def checks():
    import unittest
    import test_m28_qv_lr_four_cell as T
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(T))
    # Exercise old completed schemas only; no native forward or DNA/GFF/checkpoint read.
    if result.wasSuccessful():
        load_cell('lora_high', '13499823', '13499824')
        load_cell('dense_low', '13510351', '13510352')
    report = dict(completed=result.wasSuccessful(), job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL,
        tests_run=result.testsRun, errors=len(result.errors), failures=len(result.failures),
        existing_completed_E_M_JSON_receipts_read=result.wasSuccessful(), old_sources_structural_schema_verified=result.wasSuccessful(),
        DNA_read=False, real_reference_GFF_read=False, checkpoint_read=False, native_forward=False,
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False)
    (CHECKS/os.environ['SLURM_JOB_ID']/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    if not result.wasSuccessful():
        raise RuntimeError('Four-cell accounting checks failed')


def main(args):
    began = time.perf_counter()
    qualification = completed(CHECKS/args.checks_job, args.checks_job)
    if qualification['tests_run'] != 3 or qualification['errors'] or qualification['failures'] or not qualification['old_sources_structural_schema_verified']:
        raise ValueError('Three qualified accounting groups and old source schema match required')
    jobs = dict(lora_low=(args.low_evaluation_job, args.low_paired_job), lora_high=('13499823', '13499824'),
                dense_low=('13510351', '13510352'), dense_high=(args.high_evaluation_job, args.high_paired_job))
    reports, configurations = {}, {}
    for name, (evaluation_job, paired_job) in jobs.items():
        reports[name], configurations[name] = load_cell(name, evaluation_job, paired_job)
    if any(any(c[k] != configurations['lora_low'][k] for k in FIXED_CONFIG) for c in configurations.values()):
        raise ValueError('Common scientific recipe fields differ across cells')
    scopes = summarize(reports)
    report = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL, checks_job=args.checks_job,
        source_evaluation_jobs={name: pair[0] for name, pair in jobs.items()},
        source_paired_jobs={name: pair[1] for name, pair in jobs.items()},
        source_fit_jobs={name: cell[1] for name, cell in CELLS.items()},
        native_learning_rates={name: cell[2] for name, cell in CELLS.items()}, head_learning_rate=3e-4,
        reference_policy=reports['lora_low']['reference_policy'], all_four_declared_cells_included=True,
        source_DNA_equality_before_reference_receipts_verified=True, common_scientific_recipe_verified=True,
        counts_PRF_background_phase_nonprimary_macro_worst_reconciled=True,
        retained_full_primary_chains=22688, scopes=scopes,
        units='PRF and span FPR in percent; their contrasts in percentage points; all counts and chains_per_Mb retain native units',
        worst_species_policy='min PRF separately per cell and metric; max harmful background rate separately per cell; then contrast',
        DNA_read=False, real_reference_GFF_read=False, checkpoint_read=False, native_forward=False,
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False, frozen_DEV_family_changed=False,
        predictions_modified=False, seconds=time.perf_counter()-began, limitation=LIMITATION)
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False); handle.write('\n')
    lines = ['# 原 TRAIN 的 q/v 参数化与学习率四格', '',
        '| 范围 | LoRA 低lr F1 % | LoRA 高lr F1 % | Dense 低lr F1 % | Dense 高lr F1 % | LoRA 高减低 pp | Dense 高减低 pp | 差中差 pp |',
        '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for scope in (*SPECIES, 'pooled', 'macro_species', 'worst_species'):
        row = scopes[scope]
        key = 'exact_chain.f1_percent'
        values = [row['arms'][name][key] for name in CELLS]
        effects = [row['contrasts'][name][key] for name in ('lora_high_minus_low', 'dense_high_minus_low', 'difference_in_differences')]
        lines.append('| '+scope+' | '+' | '.join(f'{v:.6f}' for v in values+effects)+' |')
    lines += ['', '全部六染色体、base、背景、phase、其他isoform及duplicate计数与对照保存在summary.json；负值不省略。',
        '', '非主链／未匹配annotated assessable isoform不等于新基因或已证实的假阳性。', '', LIMITATION]
    with (OUT/'comparison.md').open('x') as handle:
        handle.write('\n'.join(lines)+'\n')
    print(json.dumps(dict(completed=True, job_id=report['job_id'], source_evaluation_jobs=report['source_evaluation_jobs'])), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checks', action='store_true')
    for name in ('low-evaluation-job', 'high-evaluation-job', 'low-paired-job', 'high-paired-job', 'checks-job'):
        parser.add_argument('--'+name)
    args = parser.parse_args()
    if args.checks:
        checks()
    elif any(getattr(args, k.replace('-', '_')) is None for k in ('low-evaluation-job', 'high-evaluation-job', 'low-paired-job', 'high-paired-job', 'checks-job')):
        parser.error('All new source and checks job IDs required')
    else:
        main(args)
