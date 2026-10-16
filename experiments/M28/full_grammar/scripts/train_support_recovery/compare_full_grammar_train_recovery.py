#!/usr/bin/env python3
"""Fixed TRAIN reference recovery across seven checkpoints and the separate combination."""
import copy
import json
import os
from collections import Counter
from pathlib import Path
import time
import diagnose_native_train_reference as D

ROOT = D.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-TRAIN-RECOVERY'
CHECKS = ROOT/'outputs/M28-FULL-GRAMMAR-TRAIN-RECOVERY-CHECKS'
PLAN = 'reports/M28-METHOD-RESTART/full_grammar_train_recovery_plan.md'
ANCHOR = 'M28_final4608'
CONTROL = 'M28_STRUCTURED_CONTINUE_final9216'
LOCAL = 'M28_LOCAL_CHANNEL_final9216'
COMBO = 'M28_JOINT24576_NATIVE_CDS_PRIOR_final9216'
METHODS = (ANCHOR,CONTROL,LOCAL,'M28_HALF_STRUCTURED_HALF_LOCAL_final9216',
    'M28_NATIVE_CDS_PRIOR_final9216','M28_FRESH_DRAWS_final9216','M28_JOINT24576_final9216',COMBO)
RAW_ROOTS = {
    **{'M28_'+arm+'_final9216':ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-TRAIN-RAW'/arm
       for arm in ('STRUCTURED_CONTINUE','LOCAL_CHANNEL','HALF_STRUCTURED_HALF_LOCAL')},
    'M28_NATIVE_CDS_PRIOR_final9216':ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CDS-PRIOR-TRAIN-RAW',
    'M28_FRESH_DRAWS_final9216':ROOT/'outputs/M28-FULL-GRAMMAR-FRESH-DRAWS-TRAIN-RAW',
    'M28_JOINT24576_final9216':ROOT/'outputs/M28-FULL-GRAMMAR-JOINT24-TRAIN-RAW',
    COMBO:ROOT/'outputs/M28-FULL-GRAMMAR-JOINT24-NATIVE-PRIOR-TRAIN-RAW'}
GROUPS = ('anchor_hit_any_native','anchor_hit_no_native',
          'anchor_missed_any_native','anchor_missed_no_native')


def fixed_group(row):
    native = any(row['native_exact_hits'].values())
    return ('anchor_hit_' if row['exact_primary_hit'] else 'anchor_missed_')+('any_native' if native else 'no_native')


def pair(rows, label, control=ANCHOR):
    counts = Counter()
    for row in rows:
        a,b = row['model_hits'][control],row['model_hits'][label]
        counts['retained' if a and b else 'gained' if b else 'lost' if a else 'still_missed'] += 1
    n = len(rows)
    hits = counts['retained']+counts['gained']
    return dict(reference_chains=n,control_hits=counts['retained']+counts['lost'],candidate_hits=hits,
        **{k:counts[k] for k in ('retained','gained','lost','still_missed')},
        net_reference_hits=counts['gained']-counts['lost'],reference_recall=hits/n if n else None)


def summarize(rows, label):
    value = pair(rows,label)
    value['fixed_anchor_native_groups'] = {k:pair([r for r in rows if fixed_group(r)==k],label) for k in GROUPS}
    value['all_three_native_exact_anchor_gap'] = pair(
        [r for r in rows if not r['exact_primary_hit'] and all(r['native_exact_hits'].values())],label)
    return value


def describe(rows, label):
    value = summarize(rows,label)
    value['strata'] = {
        'CDS_segments':{k:summarize([r for r in rows if (r['CDS_segments']==1)==single],label)
                       for k,single in (('single',True),('multiple',False))},
        'span_bp':{k:summarize([r for r in rows if (r['span_bp']>6144)==long],label)
                   for k,long in (('at_most_6144',False),('over_6144',True))},
        'original_draw_R4_positive_chain':{
            k:summarize([r for r in rows if r['original_draw_R4_positive_chain']==present],label)
            for k,present in (('present',True),('absent',False))},
        'original_model_error_category':{
            k:summarize([r for r in rows if r['original_model_error_category']==k],label) for k in D.CATEGORIES}}
    return value


def source_specs():
    import aggregate_full_grammar_context_prior_train as X
    C = X.C
    specs = {**C.SPECS,COMBO:X.SPECS['F11']}
    if tuple(specs) != METHODS or C.ANCHOR != ANCHOR or C.CONTROL != CONTROL:
        raise ValueError('Fixed TRAIN candidates differ; no automatic family changes')
    return C,X,specs


def load_reports(labels):
    C,X,specs = source_specs()
    # All declared sources close before any saved reference/GFF read in the full run.
    for label in labels:
        folder,_,job,_ = specs[label]
        if ((folder/'STATUS').read_text().strip() != 'COMPLETED'
                or (folder/'JOBID').read_text().strip() != job):
            raise ValueError('Every declared evaluation must complete: '+label)
    reports = {label:json.loads((specs[label][0]/'summary.json').read_text()) for label in labels}
    for label,report in reports.items():
        if label == COMBO: X.validate_report('F11',report,reports[ANCHOR])
        else: C.validate_report(label,report,reports[ANCHOR])
    return C,specs,reports


def source_records(control, C):
    source = D.completed(D.OUT)
    if (source['job_id'] != '13461637' or source['source_evaluation_job'] != '13461307'
            or source['source_original_evaluation_job'] != control['job_id']
            or source['checkpoint'] != control['checkpoint']
            or source['reference_policy'] != control['reference_policy']
            or source['reference_primary_chains'] != sum(D.EXPECTED.values())
            or source['reference_filters_added'] or not source['exact_hit_counts_reproduce_completed_evaluations']
            or source['DEV_used'] or source['sealed_data_used'] or source['optimizer_steps_this_run'] != 0):
        raise ValueError('Completed paired reference source differs')
    rows = D.read_rows(D.OUT/'references.jsonl')
    scope = {(s,q):D.EXPECTED[q] for s,q,_ in C.E.I.SCOPE}
    if (Counter((r['species'],r['seqid']) for r in rows) != Counter(scope)
            or len({D.reference_key(r) for r in rows}) != len(rows)
            or any(set(r['native_exact_hits']) != set(D.METHODS) for r in rows)
            or D.summarize(rows) != {k:v for k,v in source['pooled'].items() if k!='strata'}):
        raise ValueError('Complete paired reference identities or frozen memberships differ')
    for row in rows: row['model_hits'] = {ANCHOR:row['exact_primary_hit']}
    for species,seqid,_ in C.E.I.SCOPE:
        part = [r for r in rows if (r['species'],r['seqid'])==(species,seqid)]
        if pair(part,ANCHOR)['candidate_hits'] != control['chromosomes'][seqid]['exact_chain']['tp']:
            raise ValueError('Original4608 per-contig hit membership differs')
    return rows,source


def assign_hits(rows, label, predicted, prediction_records, official):
    keys = {D.reference_key(r) for r in rows}
    observed = D.exact_counts(predicted,keys)
    if (observed != {k:official['exact_chain'][k] for k in ('tp','predicted','reference')}
            or prediction_records != official['prediction_records']
            or prediction_records-len(predicted) != official['duplicate_chain_records']):
        raise ValueError('Prediction replay differs from completed unchanged ruler: '+label)
    for row in rows: row['model_hits'][label] = D.reference_key(row) in predicted


def replay(rows, label, report, C):
    for species,seqid,length in C.E.I.SCOPE:
        part = [r for r in rows if (r['species'],r['seqid'])==(species,seqid)]
        official = report['chromosomes'][seqid]
        prediction = ROOT/official['raw_prediction_path']
        if prediction.resolve() != (RAW_ROOTS[label]/seqid/'prediction/full_grammar.gff3').resolve():
            raise ValueError('Completed candidate prediction path differs: '+label)
        txs = C.E.C.S.E.primary_transcripts(C.E.C.S.E.parse_annotation(prediction,{seqid:length}),complete_only=False)
        predicted = {(species,)+C.E.C.S.key(t) for t in txs}
        assign_hits(part,label,predicted,len(txs),official)
    if pair(rows,label)['candidate_hits'] != report['pooled']['exact_chain']['tp']:
        raise ValueError('Pooled candidate reference hits differ: '+label)
    for species in C.SPECIES:
        if pair([r for r in rows if r['species']==species],label)['candidate_hits'] != report['per_species'][species]['exact_chain']['tp']:
            raise ValueError('Species candidate reference hits differ: '+label)


def complete_family(rows, labels):
    if tuple(labels) != METHODS or any(set(r['model_hits']) != set(METHODS) for r in rows):
        raise ValueError('All eight predetermined TRAIN models, not a selected subset, required')
    return {label:describe(rows,label) for label in labels}


def main():
    began = time.perf_counter()
    checked = D.completed(CHECKS)
    if not checked['completed_LOCAL_replay_verified'] or checked['source_LOCAL_evaluation_job'] != '13451623':
        raise ValueError('Synthetic tests and completed LOCAL replay required')
    C,specs,reports = load_reports(METHODS)
    rows,source = source_records(reports[ANCHOR],C)
    for label in METHODS[1:]: replay(rows,label,reports[label],C)
    pooled = complete_family(rows,METHODS)
    per_species = {s:{m:describe([r for r in rows if r['species']==s],m) for m in METHODS} for s in C.SPECIES}
    matched = {scope:{m:pair(part,m,CONTROL) for m in METHODS}
               for scope,part in [('pooled',rows)]+[(s,[r for r in rows if r['species']==s]) for s in C.SPECIES]}
    report = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],plan=PLAN,
        status='complete_original_TRAIN_fixed_model_reference_recovery_not_DEV_or_generalization',
        methods=list(METHODS),source_evaluation_jobs={m:r['job_id'] for m,r in reports.items()},
        source_paired_reference_job=source['job_id'],reference_primary_chains=len(rows),
        reference_policy=reports[ANCHOR]['reference_policy'],reference_filters_added=False,
        all_declared_candidates_included=True,exact_hits_and_unique_predictions_reproduce_completed_evaluations=True,
        key_definition=source['key_definition'],fixed_groups_defined_by_original4608_and_frozen_native_hits=True,
        pooled=pooled,per_species=per_species,matched_structured9216_reference_pairs=matched,
        measurements={m:dict(cumulative_optimizer_steps=specs[m][3],
            **{k:copy.deepcopy(r[k]) for k in ('checkpoint','chromosomes','per_species','pooled',
                'macro_species_exact_chain_F1','worst_species_exact_chain_F1')}) for m,r in reports.items()},
        DNA_read=False,real_reference_GFF_read=False,checkpoint_read=False,prediction_GFF_read=True,
        saved_paired_reference_profile_read=True,new_model_inference=False,predictions_modified=False,
        optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,
        scientific_winner_or_DEV_candidate_selected=False,statistical_significance_or_mechanism_causality_claimed=False,
        limitations=[
            'Reference recovery/loss partitions are descriptive recall accounting, not subset prediction precision or F1.',
            'Original4608 has fewer updates; matched structured9216 counts do not establish all-factor causal equivalence.',
            'Fixed anchor/native groups are outcome-defined diagnostic subsets, not pre-treatment causal strata.',
            'Old R4 positive-chain flags do not describe all training exposure, especially after fresh sampling.',
            'Native-only and model-only hits are annotated primary references, not newly discovered genes.',
            'TRAIN, one fit seed, shared upstream exposure and no independent biological replicates cannot establish generalization.',
            'The separate combination is not added to the existing seven-checkpoint or frozen eight/ten-method DEV families.'])
    report['seconds'] = time.perf_counter()-began
    with (OUT/'references.jsonl').open('x') as handle:
        for row in rows: handle.write(json.dumps(row,allow_nan=False)+'\n')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    lines = ['# 原 TRAIN 固定模型的参考恢复与损失','',
        '全部22688条原primary参考；七权重和独立组合完整保留，不排序选模、不改冻结DEV家族。',
        '恢复和损失以原4608命中集合为锚点；基线缺口为原4608遗漏而至少一个固定原生工具命中的7831条参考。','',
        '| 模型 | 范围 | 完整链F1 % | 恢复旧遗漏 | 丢失旧命中 | 净TP变化 | 恢复基线缺口 | 基线缺口分母 | 恢复三基线一致缺口 |',
        '| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for m in METHODS:
        for scope in (*C.SPECIES,'pooled'):
            stats = pooled[m] if scope=='pooled' else per_species[scope][m]
            metric = reports[m]['pooled'] if scope=='pooled' else reports[m]['per_species'][scope]
            gap = stats['fixed_anchor_native_groups']['anchor_missed_any_native']
            lines.append('| %s | %s | %.4f | %d | %d | %+d | %d | %d | %d |' %
                (m,scope,100*metric['exact_chain']['f1'],stats['gained'],stats['lost'],stats['net_reference_hits'],
                 gap['candidate_hits'],gap['reference_chains'],stats['all_three_native_exact_anchor_gap']['candidate_hits']))
    lines += ['', '完整P/R/F1、CDS-base、背景、isoform、重复、phase检查数量及全部结构分层见summary.json。',
        '恢复若伴随旧命中流失或背景增加，不能单独用恢复数声称净性能改善。',
        '6144跨度层是原单块编码长度，不是硬grammar限制；原R4标记缺失不是独立未见。',
        '这是单seed原TRAIN描述性计数，不是显著性、模块因果、新基因或泛化证据。']
    with (OUT/'comparison.md').open('x') as handle: handle.write('\n'.join(lines)+'\n')


if __name__ == '__main__': main()
