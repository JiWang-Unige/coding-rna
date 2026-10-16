#!/usr/bin/env python3
"""Complete predeclared DEV family, with direct common-baseline reconciliation."""
import copy
import json
import math
import os
import time
import evaluate_full_grammar_dev_baselines as E
import infer_full_grammar_objective_contrast_dev as C
import evaluate_full_grammar_native_cds_prior_dev as N

ROOT = C.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-CONTINUATION-DEV-FAMILY'
PLAN = N.I.PLAN
SPECS = {
    'M28_final4608':(E.OUT, E.I.CHECKPOINT),
    **{'M28_'+arm+'_final9216':(ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-DEV-COMPARISON'/arm,
                             C.final_checkpoint(arm)) for arm in C.A.ARMS},
    N.I.METHOD:(N.OUT, N.I.CHECKPOINT)}
CONTROL = 'M28_STRUCTURED_CONTINUE_final9216'
METHODS = (*SPECS, *E.N.METHODS)


def build_comparison(reports):
    if set(reports) != set(SPECS): raise ValueError('All five fixed M28 reports, not a selected subset, required')
    anchor = reports['M28_final4608']; measured = {}
    for label, report in reports.items():
        if (not report['completed'] or report['status'] != 'fixed_two_DEV_development_comparison_not_independent_generalization'
                or report['checkpoint'] != str(SPECS[label][1]) or not report['DEV_used']
                or report['sealed_data_used'] or report['optimizer_steps_this_run'] != 0
                or not report['input_records_directly_equal'] or not report['unsupported_reference_strata_retained']
                or report['reference_policy'] != 'M26_CDS_assessable_longest_eligible_CDS_tie_transcript_ID'
                or set(report['species']) != set(E.EXPECTED) or set(report['pooled']) != {label, *E.N.METHODS}):
            raise ValueError('Fixed candidate, reference, input or comparison identity differs: '+label)
        for species, seqid, length in E.I.SCOPE:
            item = report['species'][species]
            if (item['seqid'], item['length'], item['reference_primary_chains']) != (seqid, length, E.EXPECTED[species]):
                raise ValueError('Frozen DEV reference denominator or contig differs')
            if (set(item['methods']) != {label, *E.N.METHODS}
                    or item['retained_reference_strata'] != anchor['species'][species]['retained_reference_strata']):
                raise ValueError('Reference strata or complete comparator family differs')
            for method in E.N.METHODS:
                if item['methods'][method] != anchor['species'][species]['methods'][method]:
                    raise ValueError('Fixed baseline per-species measurements differ: '+method)
        for method, pooled in report['pooled'].items():
            rows = [report['species'][species]['methods'][method]['exact_chain'] for species in E.EXPECTED]
            counts = [sum(row[key] for row in rows) for key in ('tp', 'predicted', 'reference')]
            expected = E.C.S.prf(*counts)
            if (any(row['reference'] != E.EXPECTED[s] for row, s in zip(rows, E.EXPECTED))
                    or pooled['exact_chain'] != expected
                    or not math.isclose(pooled['macro_species_exact_chain_F1'], sum(r['f1'] for r in rows)/2,
                                        rel_tol=1e-12, abs_tol=1e-12)
                    or pooled['worst_species_exact_chain_F1'] != min(r['f1'] for r in rows)):
                raise ValueError('Micro counts, macro or worst-species aggregation differs')
            if method in E.N.METHODS and pooled != anchor['pooled'][method]:
                raise ValueError('Fixed baseline pooled measurements differ: '+method)
        measured[label] = dict(per_species={s:copy.deepcopy(report['species'][s]['methods'][label]) for s in E.EXPECTED},
                              pooled=copy.deepcopy(report['pooled'][label]))
    for method in E.N.METHODS:
        measured[method] = dict(per_species={s:copy.deepcopy(anchor['species'][s]['methods'][method]) for s in E.EXPECTED},
                               pooled=copy.deepcopy(anchor['pooled'][method]))
    table = []
    for method in METHODS:
        value = measured[method]; pooled = value['pooled']
        table.append(dict(method=method, pooled_exact_chain=pooled['exact_chain'],
            per_species_exact_chain_F1={s:value['per_species'][s]['exact_chain']['f1'] for s in E.EXPECTED},
            macro_species_exact_chain_F1=pooled['macro_species_exact_chain_F1'],
            worst_species_exact_chain_F1=pooled['worst_species_exact_chain_F1']))
    candidate = measured[N.I.METHOD]
    contrasts = {}
    for method in METHODS:
        if method == N.I.METHOD: continue
        other = measured[method]
        contrasts[method] = dict(per_species={s:100*(candidate['per_species'][s]['exact_chain']['f1']
            -other['per_species'][s]['exact_chain']['f1']) for s in E.EXPECTED},
            pooled_micro=100*(candidate['pooled']['exact_chain']['f1']-other['pooled']['exact_chain']['f1']),
            macro_species=100*(candidate['pooled']['macro_species_exact_chain_F1']-other['pooled']['macro_species_exact_chain_F1']),
            worst_species=100*(candidate['pooled']['worst_species_exact_chain_F1']-other['pooled']['worst_species_exact_chain_F1']))
    return dict(methods=list(METHODS), measurements=measured, exact_chain_table=table,
        native_prior_differences_percentage_points=contrasts, matched_primary_control=CONTROL,
        reference_primary_chains=dict(E.EXPECTED), pooled_reference_primary_chains=sum(E.EXPECTED.values()),
        all_declared_candidates_included=True, reused_native_baselines_directly_equal_across_reports=True,
        original_micro_macro_worst_arithmetic_verified=True)


def main():
    began = time.perf_counter(); N.I.require_train_closure(); reports = {}
    for label, (folder, _) in SPECS.items():
        if (folder/'STATUS').read_text().strip() != 'COMPLETED':
            raise ValueError('Every fixed DEV comparison must complete, not only a preferred candidate')
        reports[label] = json.loads((folder/'summary.json').read_text())
    value = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], plan=PLAN,
        status='complete_frozen_DEV_development_family_not_independent_generalization',
        source_comparison_jobs={label:r['job_id'] for label, r in reports.items()},
        **build_comparison(reports), DEV_used=True, sealed_data_used=False, optimizer_steps_this_run=0,
        new_predictions_or_baseline_runs_created=False,
        uncertainty='one_fit_seed_two_historically_reused_DEV_contigs_no_independent_replicates_or_clean_zero_shot',
        seconds=time.perf_counter()-began)
    with (OUT/'summary.json').open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')
    species = list(E.EXPECTED)
    lines = ['# 冻结 DEV 的完整续训比较', '',
        '原两条开发染色体、7728个primary CDS参考链；同一共同评价尺。不是独立泛化或论文成功证据。', '',
        '| 模型 | AT F1 % | Rice F1 % | pooled micro F1 % | macro F1 % | worst F1 % |',
        '| --- | --- | --- | --- | --- | --- |']
    for row in value['exact_chain_table']:
        values = [row['per_species_exact_chain_F1'][s]*100 for s in species]
        values += [row['pooled_exact_chain']['f1']*100, row['macro_species_exact_chain_F1']*100,
                   row['worst_species_exact_chain_F1']*100]
        lines.append('| '+row['method']+' | '+' | '.join('%.4f'%v for v in values)+' |')
    lines += ['', '双头支路的主要匹配控制为STRUCTURED_CONTINUE末9216。各物种和汇总差值均见summary.json；',
        '完整CDS base、背景假阳性、phase、isoform和结构层保留在measurements，不能由汇总F1代替。', '',
        '单seed、历史复用DEV、方法相关的上游注释暴露；不采用gene级伪独立置信区间，也不以DEV选模。']
    with (OUT/'comparison.md').open('x') as handle: handle.write('\n'.join(lines)+'\n')


if __name__ == '__main__': main()
