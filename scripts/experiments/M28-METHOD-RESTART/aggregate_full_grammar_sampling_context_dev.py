#!/usr/bin/env python3
"""Separate ten-method DEV extension; never overwrite the original eight-method family."""
import copy
import json
import math
import os
import time
import infer_full_grammar_sampling_context_dev as I
import aggregate_full_grammar_continuation_dev as B

E = B.E
OUT = I.ROOT/'outputs/M28-FULL-GRAMMAR-SAMPLING-CONTEXT-DEV-FAMILY'
NEW_SPECS = {method:(I.COMPARISON/branch,train.CHECKPOINT)
    for branch in I.BRANCHES for train,method in (I.candidate(branch),)}
SPECS = dict(B.SPECS,**NEW_SPECS)
METHODS = (*B.METHODS,*NEW_SPECS)


def new_measurement(label,report,anchor):
    if (not report['completed'] or report['status'] != 'fixed_two_DEV_development_comparison_not_independent_generalization'
            or report['checkpoint'] != str(NEW_SPECS[label][1]) or not report['DEV_used']
            or report['sealed_data_used'] or report['optimizer_steps_this_run'] != 0
            or not report['input_records_directly_equal'] or not report['unsupported_reference_strata_retained']
            or report['reference_policy'] != 'M26_CDS_assessable_longest_eligible_CDS_tie_transcript_ID'
            or set(report['species']) != set(E.EXPECTED) or set(report['pooled']) != {label,*E.N.METHODS}):
        raise ValueError('Fixed new-candidate, input or comparison identity differs: '+label)
    for species,seqid,length in E.I.SCOPE:
        item = report['species'][species]
        if ((item['seqid'],item['length'],item['reference_primary_chains']) != (seqid,length,E.EXPECTED[species])
                or item['retained_reference_strata'] != anchor['species'][species]['retained_reference_strata']
                or set(item['methods']) != {label,*E.N.METHODS}):
            raise ValueError('Frozen DEV denominator, contig or retained strata differ')
        for method in E.N.METHODS:
            if item['methods'][method] != anchor['species'][species]['methods'][method]:
                raise ValueError('Fixed baseline per-species measurements differ: '+method)
        for method,value in item['methods'].items():
            row = value['exact_chain']
            if row != E.C.S.prf(row['tp'],row['predicted'],E.EXPECTED[species]):
                raise ValueError('Per-species exact-chain count arithmetic differs')
    for method,pooled in report['pooled'].items():
        rows = [report['species'][s]['methods'][method]['exact_chain'] for s in E.EXPECTED]
        desired = E.C.S.prf(*(sum(r[k] for r in rows) for k in ('tp','predicted','reference')))
        if (pooled['exact_chain'] != desired or not math.isclose(pooled['macro_species_exact_chain_F1'],
                    sum(r['f1'] for r in rows)/2,rel_tol=1e-12,abs_tol=1e-12)
                or pooled['worst_species_exact_chain_F1'] != min(r['f1'] for r in rows)):
            raise ValueError('New-candidate micro, macro or worst-species aggregation differs')
        if method in E.N.METHODS and pooled != anchor['pooled'][method]:
            raise ValueError('Fixed baseline pooled measurements differ: '+method)
    return dict(per_species={s:copy.deepcopy(report['species'][s]['methods'][label]) for s in E.EXPECTED},
        pooled=copy.deepcopy(report['pooled'][label]))


def build_comparison(reports):
    if set(reports) != set(SPECS):
        raise ValueError('All seven fixed M28 reports, not a selected subset, required')
    original = B.build_comparison({label:reports[label] for label in B.SPECS})
    measured = original['measurements']
    for label in NEW_SPECS:
        measured[label] = new_measurement(label,reports[label],reports['M28_final4608'])
    table = []
    for method in METHODS:
        value = measured[method]
        table.append(dict(method=method,pooled_exact_chain=value['pooled']['exact_chain'],
            per_species_exact_chain_F1={s:value['per_species'][s]['exact_chain']['f1'] for s in E.EXPECTED},
            macro_species_exact_chain_F1=value['pooled']['macro_species_exact_chain_F1'],
            worst_species_exact_chain_F1=value['pooled']['worst_species_exact_chain_F1']))
    contrasts = {}
    for label in NEW_SPECS:
        candidate = measured[label]
        contrasts[label] = {}
        for method in METHODS:
            if method == label:
                continue
            other = measured[method]
            contrasts[label][method] = dict(per_species={s:100*(candidate['per_species'][s]['exact_chain']['f1']
                -other['per_species'][s]['exact_chain']['f1']) for s in E.EXPECTED},
                pooled_micro=100*(candidate['pooled']['exact_chain']['f1']-other['pooled']['exact_chain']['f1']),
                macro_species=100*(candidate['pooled']['macro_species_exact_chain_F1']-other['pooled']['macro_species_exact_chain_F1']),
                worst_species=100*(candidate['pooled']['worst_species_exact_chain_F1']-other['pooled']['worst_species_exact_chain_F1']))
    return dict(methods=list(METHODS),measurements=measured,exact_chain_table=table,
        new_candidate_differences_percentage_points=contrasts,matched_primary_control=B.CONTROL,
        original_eight_method_family_preserved=True,original_eight_methods=list(B.METHODS),
        reference_primary_chains=dict(E.EXPECTED),pooled_reference_primary_chains=sum(E.EXPECTED.values()),
        all_declared_candidates_included=True,reused_native_baselines_directly_equal_across_reports=True,
        original_micro_macro_worst_arithmetic_verified=True,
        prospective_extension_not_part_of_original_eight_method_freeze=True)


def main():
    began = time.perf_counter()
    I.require_train_closure()
    reports = {}
    for label,(folder,_) in SPECS.items():
        if (folder/'STATUS').read_text().strip() != 'COMPLETED':
            raise ValueError('Every fixed comparison must complete, not only a preferred candidate')
        reports[label] = json.loads((folder/'summary.json').read_text())
    value = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],plan=I.PLAN,
        status='separate_prospective_ten_method_DEV_extension_not_independent_generalization',
        source_comparison_jobs={label:r['job_id'] for label,r in reports.items()},**build_comparison(reports),
        DEV_used=True,sealed_data_used=False,optimizer_steps_this_run=0,new_predictions_or_baseline_runs_created=False,
        uncertainty='one_fit_seed_two_historically_reused_DEV_contigs_no_independent_replicates_or_clean_zero_shot',
        seconds=time.perf_counter()-began)
    with (OUT/'summary.json').open('x') as handle:
        json.dump(value,handle,indent=2,allow_nan=False); handle.write('\n')
    species = list(E.EXPECTED)
    lines = ['# 新抽样和联合编码的独立扩展 DEV 比较','',
        '原八方法不变，另加两个预定末9216候选；7728条primary参考链、同一共同尺。',
        '这是历史复用开发染色体上的描述性比较，不是独立泛化或论文成功。','',
        '| 模型 | AT F1 % | Rice F1 % | pooled micro F1 % | macro F1 % | worst F1 % |',
        '| --- | --- | --- | --- | --- | --- |']
    for row in value['exact_chain_table']:
        values = [row['per_species_exact_chain_F1'][s]*100 for s in species]
        values += [row['pooled_exact_chain']['f1']*100,row['macro_species_exact_chain_F1']*100,
                   row['worst_species_exact_chain_F1']*100]
        lines.append('| '+row['method']+' | '+' | '.join('%.4f'%v for v in values)+' |')
    lines += ['',
        '两个新候选均以STRUCTURED_CONTINUE末9216为主要匹配控制；差值及完整测量见summary.json。',
        '新抽样改变提议、有限曝光和梯度方差，不能单独归因于覆盖；联合编码只改变每窗口可见输入。',
        'CDS base、背景、phase、isoform及结构层保留，不能由汇总F1代替。',
        '不按DEV选择checkpoint、阈值、上下文或后续训练，不使用gene级伪独立置信区间。']
    with (OUT/'comparison.md').open('x') as handle:
        handle.write('\n'.join(lines)+'\n')


if __name__ == '__main__':
    main()
