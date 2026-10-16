#!/usr/bin/env python3
"""Prespecified four-cell TRAIN context x native-CDS-prior descriptive comparison."""
import copy
import json
import math
import os
import time
import aggregate_full_grammar_continuation_train as C

E = C.E
ROOT = C.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-CONTEXT-PRIOR-TRAIN-COMPARISON'
CHECKS = ROOT/'outputs/M28-FULL-GRAMMAR-CONTEXT-PRIOR-TRAIN-CHECKS'
PLAN = 'reports/M28-METHOD-RESTART/full_grammar_joint24_native_prior_plan.md'
STATUS = 'complete_original_TRAIN_context_prior_four_cell_diagnostic_not_DEV_or_generalization'
SPECIES = C.SPECIES
REFERENCE_COUNTS = C.REFERENCE_COUNTS
SPECS = {
    'F00': C.SPECS['M28_STRUCTURED_CONTINUE_final9216'],
    'F01': C.SPECS['M28_NATIVE_CDS_PRIOR_final9216'],
    'F10': C.SPECS['M28_JOINT24576_final9216'],
    'F11': (ROOT/'outputs/M28-FULL-GRAMMAR-JOINT24-NATIVE-PRIOR-TRAIN-DIAGNOSTIC',
        ROOT/'outputs/M28-FULL-GRAMMAR-JOINT24-NATIVE-PRIOR-FIT/step_009216.pt', '13457852', 9216)}
CELL_DEFINITION = {
    'F00': dict(encoder_context_bp=6144, native_CDS_prior=False, added_parameters=0),
    'F01': dict(encoder_context_bp=6144, native_CDS_prior=True, added_parameters=28),
    'F10': dict(encoder_context_bp=24576, native_CDS_prior=False, added_parameters=0),
    'F11': dict(encoder_context_bp=24576, native_CDS_prior=True, added_parameters=28)}
CONTRASTS = {
    'context_without_prior': {'F10': 1, 'F00': -1},
    'prior_without_joint_context': {'F01': 1, 'F00': -1},
    'context_with_prior': {'F11': 1, 'F01': -1},
    'prior_with_joint_context': {'F11': 1, 'F10': -1},
    'combined_minus_control': {'F11': 1, 'F00': -1},
    'descriptive_interaction': {'F11': 1, 'F10': -1, 'F01': -1, 'F00': 1}}


def validate_report(cell, report, control):
    if (not report['completed'] or report['status'] != 'original_TRAIN_end_to_end_diagnostic_not_DEV_or_generalization'
            or report['checkpoint'] != str(SPECS[cell][1]) or report['job_id'] != SPECS[cell][2]
            or report['reference_policy'] != 'M26_CDS_assessable_longest_eligible_CDS_tie_transcript_ID'
            or report['reference_used_in_inference'] or not report['training_eligibility_does_not_filter_reference']
            or not report['unsupported_reference_strata_retained'] or report['DEV_used']
            or report['sealed_data_used'] or report['optimizer_steps_this_run'] != 0
            or set(report['chromosomes']) != set(REFERENCE_COUNTS) or set(report['per_species']) != set(SPECIES)):
        raise ValueError('Fixed report identity or complete TRAIN scope differs: '+cell)
    for species, seqid, length in E.I.SCOPE:
        row = report['chromosomes'][seqid]
        if ((row['species'], row['seqid'], row['length'], row['reference_primary_records'], row['exact_chain']['reference'])
                != (species, seqid, length, REFERENCE_COUNTS[seqid], REFERENCE_COUNTS[seqid])
                or C.reference_metadata(row) != C.reference_metadata(control['chromosomes'][seqid])):
            raise ValueError('Frozen reference denominator, background or retained strata differ: '+cell)
        for field in ('exact_chain', 'CDS_base_unstranded_union'):
            counts = row[field]
            if counts != E.C.S.prf(*(counts[k] for k in ('tp', 'predicted', 'reference'))):
                raise ValueError('Per-chromosome count arithmetic differs: '+cell)
        if (row['phase_checked_prediction_records'] != row['prediction_records'] or
                row['prediction_records'] != row['exact_chain']['predicted']+row['duplicate_chain_records']):
            raise ValueError('Prediction, duplicate or phase-check counts differ: '+cell)
    for species in SPECIES:
        rows = [report['chromosomes'][seqid] for s, seqid, _ in E.I.SCOPE if s == species]
        if report['per_species'][species] != E.aggregate(rows):
            raise ValueError('Per-species aggregation differs: '+cell)
    f1s = [report['per_species'][s]['exact_chain']['f1'] for s in SPECIES]
    if (report['pooled'] != E.aggregate([report['chromosomes'][seqid] for _, seqid, _ in E.I.SCOPE]) or
            not math.isclose(report['macro_species_exact_chain_F1'], sum(f1s)/len(f1s), rel_tol=1e-12, abs_tol=1e-12) or
            report['worst_species_exact_chain_F1'] != min(f1s)):
        raise ValueError('Micro counts, macro or worst-species aggregation differs: '+cell)


def f1_vector(report):
    return dict(per_species={s: report['per_species'][s]['exact_chain']['f1'] for s in SPECIES},
        pooled_micro=report['pooled']['exact_chain']['f1'],
        macro_species=report['macro_species_exact_chain_F1'],
        worst_species=report['worst_species_exact_chain_F1'])


def build_comparison(reports):
    if set(reports) != set(SPECS):
        raise ValueError('All four prespecified cells, not a selected subset, required')
    for cell in SPECS:
        validate_report(cell, reports[cell], reports['F00'])
    measurements = {cell: dict(checkpoint=reports[cell]['checkpoint'], cumulative_optimizer_steps=SPECS[cell][3],
        **{k: copy.deepcopy(reports[cell][k]) for k in ('chromosomes', 'per_species', 'pooled',
           'macro_species_exact_chain_F1', 'worst_species_exact_chain_F1')}) for cell in SPECS}
    vectors = {cell: f1_vector(reports[cell]) for cell in SPECS}
    contrasts = {}
    for name, weights in CONTRASTS.items():
        contrasts[name] = dict(
            per_species={s: 100*sum(w*vectors[c]['per_species'][s] for c, w in weights.items()) for s in SPECIES},
            **{k: 100*sum(w*vectors[c][k] for c, w in weights.items())
               for k in ('pooled_micro', 'macro_species', 'worst_species')})
    return dict(cells=list(SPECS), cell_definition=copy.deepcopy(CELL_DEFINITION), measurements=measurements,
        matched_primary_control='F00', contrast_weights=copy.deepcopy(CONTRASTS),
        exact_chain_F1_contrasts_percentage_points=contrasts,
        all_declared_cells_included=True, cumulative_optimizer_steps_matched=True,
        reference_metadata_directly_equal=True, reference_ID_or_geometry_equality_retested=False,
        reference_primary_chains={s: sum(REFERENCE_COUNTS[q] for t, q, _ in E.I.SCOPE if t == s) for s in SPECIES},
        pooled_reference_primary_chains=sum(REFERENCE_COUNTS.values()),
        worst_species_contrast_definition='contrast_of_each_cell_species_minimum_not_minimum_species_contrast')


def load_completed_reports():
    reports = {}
    for cell, (folder, _, job, _) in SPECS.items():
        if ((folder/'STATUS').read_text().strip() != 'COMPLETED' or
                (folder/'JOBID').read_text().strip() != job):
            raise ValueError('All four exact evaluation jobs must complete before comparison')
        reports[cell] = json.loads((folder/'summary.json').read_text())
    return reports


def main():
    began = time.perf_counter()
    reports = load_completed_reports()
    value = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], plan=PLAN, status=STATUS,
        source_evaluation_jobs={cell: report['job_id'] for cell, report in reports.items()}, **build_comparison(reports),
        research_DNA_annotations_predictions_or_checkpoints_read=False, completed_measurement_summaries_only=True,
        new_predictions_or_baseline_runs_created=False, optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
        uncertainty='single_fit_seed_original_TRAIN_shared_upstream_exposure_descriptive_nonlinear_F1_scale',
        statistical_significance_or_mechanism_causality_claimed=False,
        scientific_winner_or_DEV_candidate_selected=False, seconds=time.perf_counter()-began)
    with (OUT/'summary.json').open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')
    lines = ['# 原 TRAIN 联合上下文 × 双头 CDS 先验的四格比较', '',
        '全部六条原 TRAIN 染色体及22688条primary CDS参考链；只复用四格末9216的已完成同尺评价。',
        '单seed TRAIN描述性结果，不是显著性、机制因果、干净泛化或基线竞争力证据。', '',
        '| 格 | 编码 bp | 先验参数 | 范围 | 链 P % | 链 R % | 链 F1 % | CDS base F1 % | 纯背景链数 |',
        '| --- | --- | --- | --- | --- | --- | --- | --- | --- |']
    for cell in SPECS:
        for scope in (*SPECIES, 'pooled'):
            row = reports[cell]['pooled'] if scope == 'pooled' else reports[cell]['per_species'][scope]
            rates = [100*row['exact_chain'][k] for k in ('precision', 'recall', 'f1')]
            lines.append('| %s | %d | %d | %s | %.4f | %.4f | %.4f | %.4f | %d |' %
                (cell, CELL_DEFINITION[cell]['encoder_context_bp'], CELL_DEFINITION[cell]['added_parameters'], scope,
                 *rates, 100*row['CDS_base_unstranded_union']['f1'], row['background']['wholly_background_unique_chains']))
    lines += ['', '| F1差值，百分点 | AT | Rice | pooled micro | macro | worst |',
        '| --- | --- | --- | --- | --- | --- |']
    for name, contrast in value['exact_chain_F1_contrasts_percentage_points'].items():
        rates = [contrast['per_species'][s] for s in SPECIES]
        rates += [contrast[k] for k in ('pooled_micro', 'macro_species', 'worst_species')]
        lines.append('| '+name+' | '+' | '.join('%+.4f'%v for v in rates)+' |')
    lines += ['', '描述性交互量为F11−F10−F01+F00；pooled F1先汇总TP／预测／参考数后计算，不能平均染色体F1。',
        'worst先取每格的物种最小值，再作对比；最差物种可能切换，不是各物种差值的最小值。',
        '逐染色体、完整P/R/F1计数、unsupported、isoform、重复预测、phase检查数量及背景见summary.json。',
        'CDS base为unstranded union；phase检查数量不是phase准确率；运行计时不是公平速度基准。',
        '增加28参数和共同clip改变优化过程；匹配draws／更新次数不等于有效优化幅度完全匹配。',
        '不选择获胜方法、不设置置信区间，不改变七权重TRAIN汇总、八／十方法DEV家族或封存数据边界。']
    with (OUT/'comparison.md').open('x') as handle:
        handle.write('\n'.join(lines)+'\n')


if __name__ == '__main__':
    main()
