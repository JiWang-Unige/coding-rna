#!/usr/bin/env python3
"""Fixed seven-checkpoint TRAIN comparison using completed unchanged-ruler reports."""
import copy
import json
import math
import os
import time
import evaluate_full_grammar_original_train as E

ROOT = E.I.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-CONTINUATION-TRAIN-FAMILY'
CHECKS = ROOT/'outputs/M28-FULL-GRAMMAR-CONTINUATION-TRAIN-CHECKS'
PLAN = 'reports/M28-METHOD-RESTART/full_grammar_objective_contrast_contract.md'
ANCHOR = 'M28_final4608'
CONTROL = 'M28_STRUCTURED_CONTINUE_final9216'
SPECIES = tuple(dict.fromkeys(s for s,_,_ in E.I.SCOPE))
REFERENCE_COUNTS = dict(zip((seqid for _,seqid,_ in E.I.SCOPE), (7135,4293,6311,1638,1601,1710)))
SPECS = {
    ANCHOR:(E.OUT, E.I.CHECKPOINT, '13443329', 4608),
    **{'M28_'+arm+'_final9216':(
        ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-TRAIN-DIAGNOSTIC'/arm,
        ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST'/arm/'step_009216.pt', job, 9216)
        for arm,job in (('STRUCTURED_CONTINUE','13451625'),('LOCAL_CHANNEL','13451623'),
                        ('HALF_STRUCTURED_HALF_LOCAL','13451624'))},
    'M28_NATIVE_CDS_PRIOR_final9216':(
        ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CDS-PRIOR-TRAIN-DIAGNOSTIC',
        ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CDS-PRIOR-FIT/step_009216.pt','13453184',9216),
    'M28_FRESH_DRAWS_final9216':(
        ROOT/'outputs/M28-FULL-GRAMMAR-FRESH-DRAWS-TRAIN-DIAGNOSTIC',
        ROOT/'outputs/M28-FULL-GRAMMAR-FRESH-DRAWS-FIT/step_009216.pt','13453583',9216),
    'M28_JOINT24576_final9216':(
        ROOT/'outputs/M28-FULL-GRAMMAR-JOINT24-TRAIN-DIAGNOSTIC',
        ROOT/'outputs/M28-FULL-GRAMMAR-JOINT24-FIT/step_009216.pt','13453780',9216)}


def reference_metadata(row):
    return dict(species=row['species'],seqid=row['seqid'],length=row['length'],
        primary_records=row['reference_primary_records'],
        primary_chains=row['exact_chain']['reference'],
        assessable_isoform_records=row['reference_assessable_isoform_records'],
        CDS_base_union=row['CDS_base_unstranded_union']['reference'],
        retained_strata=row['retained_unsupported_reference_strata'],
        background_definition=row['background']['definition'],background_bp=row['background']['bp'])


def validate_report(label, report, anchor):
    if (not report['completed'] or report['status'] != 'original_TRAIN_end_to_end_diagnostic_not_DEV_or_generalization'
            or report['checkpoint'] != str(SPECS[label][1]) or report['job_id'] != SPECS[label][2]
            or report['reference_policy'] != 'M26_CDS_assessable_longest_eligible_CDS_tie_transcript_ID'
            or report['reference_used_in_inference'] or not report['training_eligibility_does_not_filter_reference']
            or not report['unsupported_reference_strata_retained'] or report['DEV_used']
            or report['sealed_data_used'] or report['optimizer_steps_this_run'] != 0
            or set(report['chromosomes']) != set(REFERENCE_COUNTS) or set(report['per_species']) != set(SPECIES)):
        raise ValueError('Fixed report identity or complete TRAIN scope differs: '+label)
    for species,seqid,length in E.I.SCOPE:
        row = report['chromosomes'][seqid]
        if ((row['species'],row['seqid'],row['length'],row['reference_primary_records'],row['exact_chain']['reference'])
                != (species,seqid,length,REFERENCE_COUNTS[seqid],REFERENCE_COUNTS[seqid])
                or reference_metadata(row) != reference_metadata(anchor['chromosomes'][seqid])):
            raise ValueError('Frozen reference denominator, background or retained strata differ: '+label)
        for field in ('exact_chain','CDS_base_unstranded_union'):
            counts = row[field]
            if counts != E.C.S.prf(*(counts[k] for k in ('tp','predicted','reference'))):
                raise ValueError('Per-chromosome count arithmetic differs: '+label)
        if (row['phase_checked_prediction_records'] != row['prediction_records'] or
                row['prediction_records'] != row['exact_chain']['predicted']+row['duplicate_chain_records']):
            raise ValueError('Prediction, duplicate or phase-check counts differ: '+label)
    for species in SPECIES:
        rows = [report['chromosomes'][seqid] for s,seqid,_ in E.I.SCOPE if s == species]
        if report['per_species'][species] != E.aggregate(rows):
            raise ValueError('Per-species aggregation differs: '+label)
    f1s = [report['per_species'][s]['exact_chain']['f1'] for s in SPECIES]
    if (report['pooled'] != E.aggregate([report['chromosomes'][seqid] for _,seqid,_ in E.I.SCOPE]) or
            not math.isclose(report['macro_species_exact_chain_F1'],sum(f1s)/len(f1s),rel_tol=1e-12,abs_tol=1e-12) or
            report['worst_species_exact_chain_F1'] != min(f1s)):
        raise ValueError('Micro counts, macro or worst-species aggregation differs: '+label)


def build_comparison(reports):
    if set(reports) != set(SPECS):
        raise ValueError('All seven fixed checkpoints, not a selected subset, required')
    anchor = reports[ANCHOR]
    for label in SPECS:
        validate_report(label,reports[label],anchor)
    measurements = {label:dict(checkpoint=report['checkpoint'],cumulative_optimizer_steps=SPECS[label][3],
        **{k:copy.deepcopy(report[k]) for k in ('chromosomes','per_species','pooled',
           'macro_species_exact_chain_F1','worst_species_exact_chain_F1')}) for label,report in reports.items()}
    contrasts = {}
    for label,report in reports.items():
        if label == CONTROL:
            continue
        other = reports[CONTROL]
        contrasts[label] = dict(cumulative_optimizer_steps_matched=SPECS[label][3] == SPECS[CONTROL][3],
            per_species={s:100*(report['per_species'][s]['exact_chain']['f1']-
                other['per_species'][s]['exact_chain']['f1']) for s in SPECIES},
            pooled_micro=100*(report['pooled']['exact_chain']['f1']-other['pooled']['exact_chain']['f1']),
            macro_species=100*(report['macro_species_exact_chain_F1']-other['macro_species_exact_chain_F1']),
            worst_species=100*(report['worst_species_exact_chain_F1']-other['worst_species_exact_chain_F1']))
    return dict(methods=list(SPECS),measurements=measurements,matched_primary_control=CONTROL,
        differences_from_control_percentage_points=contrasts,all_declared_candidates_included=True,
        reference_metadata_directly_equal=True,reference_ID_or_geometry_equality_retested=False,
        reference_primary_chains={s:sum(REFERENCE_COUNTS[q] for t,q,_ in E.I.SCOPE if t == s) for s in SPECIES},
        pooled_reference_primary_chains=sum(REFERENCE_COUNTS.values()))


def main():
    began = time.perf_counter()
    reports = {}
    for label,(folder,_,_,_) in SPECS.items():
        if (folder/'STATUS').read_text().strip() != 'COMPLETED':
            raise ValueError('Every fixed TRAIN evaluation must complete before comparison')
        reports[label] = json.loads((folder/'summary.json').read_text())
    value = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],plan=PLAN,
        status='complete_original_TRAIN_seven_checkpoint_diagnostic_not_DEV_or_generalization',
        source_evaluation_jobs={label:r['job_id'] for label,r in reports.items()},**build_comparison(reports),
        research_DNA_annotations_predictions_or_checkpoints_read=False,completed_measurement_summaries_only=True,
        new_predictions_or_baseline_runs_created=False,optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,
        uncertainty='one_fit_seed_original_TRAIN_shared_upstream_exposure_no_independent_replicates_or_generalization',
        seconds=time.perf_counter()-began)
    with (OUT/'summary.json').open('x') as handle:
        json.dump(value,handle,indent=2,allow_nan=False); handle.write('\n')
    lines = ['# 原 TRAIN 的七个固定末权重比较','',
        '全部原六染色体及22688条primary CDS参考链；复用已经完成的同尺评价。',
        'TRAIN描述性结果，不是独立泛化、基线竞争力或论文成功证据。','',
        '| 末权重 | AT F1 % | Rice F1 % | pooled F1 % | macro F1 % | worst F1 % |',
        '| --- | --- | --- | --- | --- | --- |']
    for label in SPECS:
        report = reports[label]
        values = [report['per_species'][s]['exact_chain']['f1'] for s in SPECIES]
        values += [report['pooled']['exact_chain']['f1'],report['macro_species_exact_chain_F1'],
                   report['worst_species_exact_chain_F1']]
        lines.append('| '+label+' | '+' | '.join('%.4f'%(100*v) for v in values)+' |')
    lines += ['', '| 末权重 | AT 纯背景链数 | Rice 纯背景链数 | pooled 背景 span FPR % |',
        '| --- | --- | --- | --- |']
    for label in SPECS:
        report = reports[label]
        counts = [report['per_species'][s]['background']['wholly_background_unique_chains'] for s in SPECIES]
        lines.append('| %s | %d | %d | %.4f |'%(label,*counts,100*report['pooled']['background']['span_FPR']))
    lines += ['',
        '主要控制为STRUCTURED_CONTINUE末9216；final4608少4608次更新，不能视为纯目标因果对照。',
        '七方法按预定顺序保留，不排名选模；逐染色体、CDS base、isoform、重复记录、phase检查数量和计时均见summary.json。',
        'CDS base是unstranded union，phase检查数量不是phase准确率；运行计时不是公平速度基准。',
        '双头先验增加28个参数；新抽样改变提议、有限曝光及梯度方差；联合编码改变窗口可见输入。',
        '匹配累计更新次数不能自动升级为所有因素均匹配的因果消融，不据TRAIN汇总宣称SOTA或干净泛化。',
        '不改变预定末权重、阈值、DEV候选或冻结评价；不采用gene级伪独立置信区间。']
    with (OUT/'comparison.md').open('x') as handle:
        handle.write('\n'.join(lines)+'\n')


if __name__ == '__main__':
    main()
