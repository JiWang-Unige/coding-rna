#!/usr/bin/env python3
"""Separate six-method table: all three fixed candidates and all native baselines."""
import copy
import json
import os
import time
import infer_full_grammar_joint_controls_dev as I
import evaluate_full_grammar_dev_baselines as E

OUT = I.ROOT/'outputs/M28-FULL-GRAMMAR-JOINT-CONTROLS-DEV-FAMILY-R1'
METHODS = (*I.METHODS, *E.N.METHODS)
CONTROL = I.METHODS[0]


def build_comparison(reports):
    if set(reports) != set(I.METHODS):
        raise ValueError('All three fixed candidates required, not a selected subset')
    anchor = reports[CONTROL]
    measured = {}
    for branch, label in zip(I.BRANCHES,I.METHODS):
        train, _ = I.candidate(branch)
        report = reports[label]
        if (not report['completed']
                or report['status'] != 'fixed_two_DEV_development_comparison_not_independent_generalization'
                or report['checkpoint'] != str(train.CHECKPOINT) or not report['DEV_used']
                or report['sealed_data_used'] or report['optimizer_steps_this_run'] != 0
                or not report['input_records_directly_equal'] or not report['unsupported_reference_strata_retained']
                or report['reference_policy'] != 'M26_CDS_assessable_longest_eligible_CDS_tie_transcript_ID'
                or set(report['species']) != set(E.EXPECTED) or set(report['pooled']) != {label,*E.N.METHODS}):
            raise ValueError('Fixed candidate, source inputs or common reference policy differs')
        for species, seqid, length in E.I.SCOPE:
            item = report['species'][species]
            if ((item['seqid'],item['length'],item['reference_primary_chains']) != (seqid,length,E.EXPECTED[species])
                    or set(item['methods']) != {label,*E.N.METHODS}
                    or item['retained_reference_strata'] != anchor['species'][species]['retained_reference_strata']):
                raise ValueError('Frozen reference denominator, contig or retained strata differs')
            for method, row in item['methods'].items():
                for endpoint in ('exact_chain','CDS_base_unstranded_union'):
                    counts = row[endpoint]
                    if counts != E.C.S.prf(*(counts[k] for k in ('tp','predicted','reference'))):
                        raise ValueError('Per-species count-derived PRF differs')
                if row['exact_chain']['reference'] != E.EXPECTED[species]:
                    raise ValueError('Primary chain denominator differs')
                for counts in row['descriptive_structure_metrics_not_reference_filters'].values():
                    if counts != E.C.S.prf(*(counts[k] for k in ('tp','predicted','reference'))):
                        raise ValueError('Descriptive structure PRF differs')
                if method in E.N.METHODS and row != anchor['species'][species]['methods'][method]:
                    raise ValueError('Fixed baseline per-species measurements differ')
        for method, pooled in report['pooled'].items():
            rows = [report['species'][s]['methods'][method] for s in E.EXPECTED]
            desired = E.aggregate(rows)
            f1s = [row['exact_chain']['f1'] for row in rows]
            desired.update(macro_species_exact_chain_F1=sum(f1s)/2,worst_species_exact_chain_F1=min(f1s))
            if pooled != desired:
                raise ValueError('Counts, base, background, macro or worst aggregation differs')
            if method in E.N.METHODS and pooled != anchor['pooled'][method]:
                raise ValueError('Fixed baseline pooled measurements differ')
        measured[label] = dict(per_species={s:copy.deepcopy(report['species'][s]['methods'][label]) for s in E.EXPECTED},
            pooled=copy.deepcopy(report['pooled'][label]))
    for method in E.N.METHODS:
        measured[method] = dict(per_species={s:copy.deepcopy(anchor['species'][s]['methods'][method]) for s in E.EXPECTED},
            pooled=copy.deepcopy(anchor['pooled'][method]))
    table = []
    contrasts = {}
    for label in I.METHODS:
        candidate = measured[label]
        contrasts[label] = {}
        for other_label in METHODS:
            if other_label == label:
                continue
            other = measured[other_label]
            contrasts[label][other_label] = dict(per_species={s:100*(candidate['per_species'][s]['exact_chain']['f1']
                -other['per_species'][s]['exact_chain']['f1']) for s in E.EXPECTED},
                pooled_micro=100*(candidate['pooled']['exact_chain']['f1']-other['pooled']['exact_chain']['f1']),
                macro_species=100*(candidate['pooled']['macro_species_exact_chain_F1']-other['pooled']['macro_species_exact_chain_F1']),
                worst_species=100*(candidate['pooled']['worst_species_exact_chain_F1']-other['pooled']['worst_species_exact_chain_F1']))
    for method in METHODS:
        value = measured[method]
        table.append(dict(method=method,pooled_exact_chain=value['pooled']['exact_chain'],
            per_species_exact_chain_F1={s:value['per_species'][s]['exact_chain']['f1'] for s in E.EXPECTED},
            macro_species_exact_chain_F1=value['pooled']['macro_species_exact_chain_F1'],
            worst_species_exact_chain_F1=value['pooled']['worst_species_exact_chain_F1']))
    return dict(methods=list(METHODS),measurements=measured,exact_chain_table=table,
        candidate_differences_percentage_points=contrasts,matched_primary_control=CONTROL,
        reference_primary_chains=dict(E.EXPECTED),pooled_reference_primary_chains=sum(E.EXPECTED.values()),
        all_declared_candidates_included=True,reused_native_baselines_directly_equal_across_reports=True,
        counts_base_background_and_micro_macro_worst_reconciled=True)


def main():
    began = time.perf_counter()
    I.require_train_closure()
    reports = {}
    for branch,label in zip(I.BRANCHES,I.METHODS):
        folder = I.COMPARISON/branch
        report = json.loads((folder/'summary.json').read_text())
        if ((folder/'STATUS').read_text().strip() != 'COMPLETED'
                or report['job_id'] != (folder/'JOBID').read_text().strip()):
            raise ValueError('Every fixed comparison must complete before the full table')
        reports[label] = report
    value = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],plan=I.PLAN,
        status='separate_three_fixed_candidates_DEV_family_not_independent_generalization',
        source_comparison_jobs={label:row['job_id'] for label,row in reports.items()},
        source_fit_jobs={label:I.SPECS[branch][1] for branch,label in zip(I.BRANCHES,I.METHODS)},
        **build_comparison(reports),DEV_used=True,sealed_data_used=False,optimizer_steps_this_run=0,
        new_predictions_or_baseline_runs_created=False,
        uncertainty='single_fit_historically_reused_two_species_DEV_not_clean_zero_shot_or_model_selection',
        seconds=time.perf_counter()-began)
    with (OUT/'summary.json').open('x') as handle:
        json.dump(value,handle,indent=2,allow_nan=False)
        handle.write('\n')
    lines = ['# 联合编码三个固定候选的 DEV 与基线比较','',
        '三候选均报告，原7728条primary CDS链和三固定原生基线；旧八／十方法文件不变。',
        'uniform与codon以fresh×joint为匹配控制；这是历史复用开发集描述性比较，不是独立泛化。','',
        '| 模型 | AT F1 % | Rice F1 % | pooled micro F1 % | macro F1 % | worst F1 % |',
        '| --- | --- | --- | --- | --- | --- |']
    for row in value['exact_chain_table']:
        numbers = [100*row['per_species_exact_chain_F1'][s] for s in E.EXPECTED]
        numbers += [100*row['pooled_exact_chain']['f1'],100*row['macro_species_exact_chain_F1'],
            100*row['worst_species_exact_chain_F1']]
        lines.append('| '+row['method']+' | '+' | '.join('%.4f'%n for n in numbers)+' |')
    lines += ['', '完整计数、CDS base、背景、phase、isoform、结构层和差值均保留在summary.json。',
        '不按DEV选择分支、checkpoint、参数或后续训练；新增codon参数年龄和联合裁剪属于该干预。']
    with (OUT/'comparison.md').open('x') as handle:
        handle.write('\n'.join(lines)+'\n')


if __name__ == '__main__':
    main()
