#!/usr/bin/env python3
"""Certify all frozen complete-core-proposal TRAIN events without model fitting."""
import argparse
from collections import Counter
import json
import os
import re
import time

import torch
import prepare_full_grammar_complete_core_draws as P
import check_full_grammar_fresh_draws_support as S

OUT = P.ROOT/'outputs/M28-FULL-GRAMMAR-COMPLETE-CORE-DRAWS-SUPPORT'
AGGREGATE = P.ROOT/'outputs/M28-FULL-GRAMMAR-COMPLETE-CORE-DRAWS-SUPPORT-AGGREGATE'
SHARDS = 16
CONFIG = dict(preparation_configuration=P.CONFIG, preparation_job='13462399',
    scope='all4107_fixed_complete_core_proposal_unique_TRAIN_events_no_resampling',
    shards=SHARDS, partition='first_draw_occurrence_ordinal_modulo16',
    max_exons=128, checkpoint_bp=256, epsilon=.1, potentials='zero_CPU_FP32',
    support='exact_fine_only_coarse_channel_likelihood_strictly_positive_for_every_role',
    reuse='matching_original_and_fixed_FRESH_finite_events_same_R4_builder_and_inputs',
    optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False)
FIELDS = ('fine_positions','coarse_positions','fully_unknown_positions',
          'genomic_context','core_oriented_interval')


def population():
    ready = json.loads((P.OUT/'summary.json').read_text())
    if ((P.OUT/'STATUS').read_text().strip() != 'COMPLETED' or not ready['completed']
            or ready['job_id'] != CONFIG['preparation_job'] or ready['configuration'] != P.CONFIG
            or ready['prospective_draws'] != 4608 or ready['selected_unique_windows'] != 4107
            or ready['optimizer_steps_this_run'] != 0 or ready['DEV_used'] or ready['sealed_data_used']):
        raise ValueError('Completed fixed complete-core-proposal TRAIN preparation required')
    plan = list(P.P.rows(P.OUT/'training_plan.jsonl'))
    unique = {}
    for row in plan:
        unique.setdefault(row['window_id'],row)
    if len(plan) != 4608 or len(unique) != 4107:
        raise ValueError('All 4608 fixed draws and 4107 unique events must remain')
    return ready,list(unique.values())


def certificates():
    result = {ident:(row,job,'REUSED_ORIGINAL_FINITE_EVENT')
              for ident,(row,job) in S.original_supports().items()}
    fresh = json.loads((S.AGGREGATE/'summary.json').read_text())
    if ((S.AGGREGATE/'STATUS').read_text().strip() != 'COMPLETED' or not fresh['completed']
            or fresh['job_id'] != '13453408' or fresh['array_job_id'] != '13453388'
            or fresh['configuration'] != S.CONFIG or fresh['unique_events'] != 4219
            or not fresh['all_fixed_draws_retained'] or not fresh['all_fresh_exact_fine_events_feasible']
            or fresh['infeasible_events'] != 0 or fresh['optimizer_steps_this_run'] != 0
            or fresh['DEV_used'] or fresh['sealed_data_used']):
        raise ValueError('Completed same-input fixed FRESH event support required')
    _,unique = S.population()
    seen = []
    for shard in range(S.SHARDS):
        root = S.OUT/f'shard_{shard}'
        report = json.loads((root/'summary.json').read_text())
        expected = [r['window_id'] for i,r in enumerate(unique) if i % S.SHARDS == shard]
        if ((root/'STATUS').read_text().strip() != 'COMPLETED' or not report['completed']
                or report['configuration'] != S.CONFIG or report['array_job_id'] != fresh['array_job_id']
                or report['shard'] != shard or report['job_id'] != fresh['source_jobs'][shard]
                or report['optimizer_steps_this_run'] != 0 or report['DEV_used'] or report['sealed_data_used']
                or [r['window_id'] for r in report['records']] != expected
                or any(not r['feasible'] for r in report['records'])):
            raise ValueError('Complete same-input FRESH certificate shard required')
        for row in report['records']:
            seen.append(row['window_id'])
            result.setdefault(row['window_id'],(row,report['job_id'],'REUSED_FRESH_FINITE_EVENT'))
    if len(seen) != 4219 or len(set(seen)) != 4219:
        raise ValueError('All fixed FRESH events must be present exactly once')
    return result


def run(shard):
    began = time.perf_counter()
    _,unique = population()
    prior = certificates()
    windows,_,catalogs,_ = P.P.load_population()
    sequences,descriptors = P.P.F.P.read_train_DNA()
    out = OUT/f'shard_{shard}'
    P.P.write(out/'started.json',dict(configuration=CONFIG,shard=shard,
        job_id=os.environ['SLURM_JOB_ID'],array_job_id=os.environ['SLURM_ARRAY_JOB_ID']))
    records = []
    with torch.inference_mode(),(out/'windows.jsonl').open('x') as handle:
        for ordinal,expected in enumerate(unique):
            if ordinal % SHARDS != shard:
                continue
            started = time.perf_counter()
            w = windows[expected['window_id']]
            chromosome = sequences[w['species']][w['seqid']]
            geometry = P.P.F.G.context_geometry(w['start'],w['end'],len(chromosome))
            a,b = geometry['genomic_context']
            raw = chromosome[w['start']:w['end']]
            core = raw if w['strand'] == '+' else P.P.F.Q.M.rc(raw)
            targets = P.P.O.build_targets(w,catalogs[w['species']],core)
            n = w['valid_bases']
            fine,region = (torch.as_tensor(targets[k])[:n] for k in ('seg','region'))
            offset = w['start']-a if w['strand'] == '+' else b-w['end']
            row = dict(ordinal=ordinal,window_id=w['id'],species=w['species'],
                fine_positions=int((fine >= 0).sum()),
                coarse_positions=int(((fine < 0) & (region >= 0)).sum()),
                fully_unknown_positions=int(((fine < 0) & (region < 0)).sum()),
                genomic_context=[a,b],core_oriented_interval=[offset,offset+n],
                feasible=True,source='NEW_EXACT_FINE_EVENT_FORWARD',source_job=os.environ['SLURM_JOB_ID'])
            if any(row[k] != expected[k] for k in FIELDS):
                raise ValueError('Frozen preparation observations or original view changed')
            if w['id'] in prior:
                certificate,job,source = prior[w['id']]
                if any(row[k] != certificate[k] for k in FIELDS):
                    raise ValueError('Reused same-R4-input event certificate metadata differs')
                row.update(source=source,source_job=job)
            else:
                context = re.sub('[^ACGT]','N',chromosome[a:b])
                if w['strand'] == '-':
                    context = P.P.F.Q.M.rc(context)
                try:
                    S.fine_event(context,offset,n,targets)
                except S.InfeasibleObservation as error:
                    row.update(feasible=False,reason=str(error))
            row['seconds'] = time.perf_counter()-started
            records.append(row)
            handle.write(json.dumps(row,allow_nan=False)+'\n')
            handle.flush()
            if len(records) % 32 == 0:
                print(json.dumps(dict(shard=shard,checked=len(records),elapsed=time.perf_counter()-began)),flush=True)
    P.P.write(out/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],
        array_job_id=os.environ['SLURM_ARRAY_JOB_ID'],configuration=CONFIG,shard=shard,records=records,
        unique_events=len(records),infeasible_events=sum(not r['feasible'] for r in records),
        original_TRAIN_DNA_sources=descriptors,original_labels_changed=False,all_fixed_draws_retained=True,
        model_forward=False,checkpoint_read=False,optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,
        seconds=time.perf_counter()-began))


def aggregate(array_job):
    began = time.perf_counter()
    _,unique = population()
    records,source_jobs = [],[]
    for shard in range(SHARDS):
        root = OUT/f'shard_{shard}'
        report = json.loads((root/'summary.json').read_text())
        expected = [r['window_id'] for i,r in enumerate(unique) if i % SHARDS == shard]
        if ((root/'STATUS').read_text().strip() != 'COMPLETED' or not report['completed']
                or report['array_job_id'] != array_job or report['shard'] != shard or report['configuration'] != CONFIG
                or report['optimizer_steps_this_run'] != 0 or report['DEV_used'] or report['sealed_data_used']
                or [r['window_id'] for r in report['records']] != expected):
            raise ValueError('All 16 completed fixed-population support shards required')
        records.extend(report['records'])
        source_jobs.append(report['job_id'])
    records.sort(key=lambda r:r['ordinal'])
    if len(records) != 4107 or [r['window_id'] for r in records] != [r['window_id'] for r in unique]:
        raise ValueError('All 4107 frozen events must remain exactly once')
    failed = [r['window_id'] for r in records if not r['feasible']]
    P.P.write(AGGREGATE/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],
        array_job_id=array_job,source_jobs=source_jobs,configuration=CONFIG,protocol=P.PROTOCOL,
        unique_events=len(records),prospective_draws=4608,all_fixed_draws_retained=True,
        source_counts=dict(Counter(r['source'] for r in records)),infeasible_events=len(failed),
        infeasible_window_ids=failed,all_selected_exact_fine_events_feasible=not failed,
        formal_fit_eligible_after_native_preflight=not failed,original_labels_changed=False,
        model_forward=False,checkpoint_read=False,optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,
        seconds=time.perf_counter()-began,interpretation='TRAIN_event_support_not_model_performance'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--shard',type=int)
    parser.add_argument('--aggregate')
    args = parser.parse_args()
    torch.set_num_threads(2)
    if args.aggregate is not None:
        aggregate(args.aggregate)
    elif args.shard is not None and 0 <= args.shard < SHARDS:
        run(args.shard)
    else:
        parser.error('One fixed shard or complete-array aggregate required')
