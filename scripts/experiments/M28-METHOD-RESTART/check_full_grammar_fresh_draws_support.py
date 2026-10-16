#!/usr/bin/env python3
"""All fixed fresh TRAIN fine events in the unchanged complete-gene universe."""
import argparse
from collections import Counter
import json
import math
import os
import re
import time

import torch
import prepare_full_grammar_fresh_draws as P
import probe_full_grammar_train_core_support as S
import probe_full_grammar_coarse_channel_train as N
from src.m28.full_grammar_boundary_compiled import compiled_boundary_full_grammar_log_partition
from src.m28.full_grammar_coarse_channel import coarse_channel_nll
from src.m28.full_grammar_observed_partition import InfeasibleObservation

ROOT = P.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-FRESH-DRAWS-SUPPORT'
AGGREGATE = ROOT/'outputs/M28-FULL-GRAMMAR-FRESH-DRAWS-SUPPORT-AGGREGATE'
CHECKS = ROOT/'outputs/M28-FULL-GRAMMAR-FRESH-DRAWS-SUPPORT-CHECKS'
SHARDS = 16
CONFIG = dict(preparation_configuration=P.CONFIG, preparation_job='13453348',
    scope='all4219_fixed_fresh_unique_TRAIN_core_events_no_resampling', shards=SHARDS,
    partition='first_draw_occurrence_ordinal_modulo16', max_exons=128,
    checkpoint_bp=256, epsilon=.1, potentials='zero_CPU_FP32',
    support='exact_fine_only_coarse_channel_likelihood_strictly_positive_for_every_role',
    reuse='matching_original_finite_hard_event_or_finite_coarse_channel_event',
    optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False)


def population():
    ready = json.loads((P.OUT/'summary.json').read_text())
    if ((P.OUT/'STATUS').read_text().strip() != 'COMPLETED' or not ready['completed']
            or ready['job_id'] != '13453348' or ready['configuration'] != P.CONFIG
            or ready['prospective_draws'] != 4608 or ready['fresh_unique_windows'] != 4219
            or ready['optimizer_steps_this_run'] != 0 or ready['DEV_used'] or ready['sealed_data_used']):
        raise ValueError('Completed unchanged fixed fresh TRAIN schedule required')
    plan = list(P.rows(P.OUT/'training_plan.jsonl'))
    unique = {}
    for row in plan:
        unique.setdefault(row['window_id'], row)
    if len(plan) != 4608 or len(unique) != 4219:
        raise ValueError('All fixed fresh draws and unique events must remain')
    return ready, list(unique.values())


def original_supports():
    hard = json.loads((S.OUT/'summary.json').read_text())
    channel = json.loads((N.OUT/'summary.json').read_text())
    for root, value, job, config in ((S.OUT, hard, '13429444', S.CONFIG),
                                   (N.OUT, channel, '13429625', N.CONFIG)):
        if ((root/'STATUS').read_text().strip() != 'COMPLETED' or not value['completed']
                or value['job_id'] != job or value['configuration'] != config
                or value['optimizer_steps_this_run'] != 0 or value['DEV_used']):
            raise ValueError('Matching completed original support certificates required')
    negative = {r['window_id'] for r in hard['records'] if not r['support_feasible']}
    corrected = {r['window_id']:r for r in channel['records']}
    if len(hard['records']) != 1491 or len(negative) != 24 or set(corrected) != negative:
        raise ValueError('Original full pool and all24 coarse-channel certificates required')
    result = {}
    for row in hard['records']:
        ident = row['window_id']
        if row['support_feasible']:
            result[ident] = (row, '13429444')
        elif corrected[ident]['candidate_event_feasible']:
            result[ident] = (corrected[ident], '13429625')
        else:
            raise ValueError('Original exact-fine event is not certified feasible')
    return result


def fine_event(sequence, offset, n, targets, cap=128):
    fine = torch.as_tensor(targets['seg'])[:n]
    exact = fine >= 0
    allowed = torch.ones((len(sequence), 15), dtype=torch.bool)
    indices = torch.where(exact)[0]+offset
    allowed[indices] = False
    allowed[indices, fine[exact].long()] = True
    cuts = tuple(x for x in (offset, offset+n) if 0 < x < len(sequence))
    result = compiled_boundary_full_grammar_log_partition(sequence,
        torch.zeros((len(sequence), 15), dtype=torch.float32), cap,
        allowed=allowed, cuts=cuts, checkpoint_bp=256)
    mass = float(result['log_partition'])
    if not math.isfinite(mass) or mass < -1e-8:
        raise ValueError('Finite nonnegative zero-potential fine-event log count required')


def fixture():
    sequence = 'ATGAAATAA'
    # Wrong coarse background cannot make the positive-epsilon event empty.
    targets = dict(seg=torch.tensor([13,13,13,-100,-100,-100,14,14,14]),
                   region=torch.zeros(9, dtype=torch.int64))
    fine_event(sequence, 0, 9, targets, 4)
    value = coarse_channel_nll(sequence, torch.zeros((9,15)), 4, 0, 9, targets, .1)
    assert bool(torch.isfinite(value['loss']))
    targets['seg'][0] = 14
    for operation in (lambda:fine_event(sequence, 0, 9, targets, 4),
                      lambda:coarse_channel_nll(sequence, torch.zeros((9,15)), 4, 0, 9, targets, .1)):
        try:
            operation()
        except InfeasibleObservation:
            pass
        else:
            raise AssertionError('Contradictory exact fine observation must remain infeasible')
    return 'coarse_positive_support_and_contradictory_exact_fine_retention'


def run(shard):
    began = time.perf_counter()
    checked = json.loads((CHECKS/'summary.json').read_text())
    if ((CHECKS/'STATUS').read_text().strip() != 'COMPLETED' or not checked['completed']
            or checked['configuration'] != CONFIG or checked['research_data_read']):
        raise ValueError('Successful same-configuration synthetic support checks required')
    _, unique = population()
    old = original_supports()
    windows, _, catalogs, _ = P.load_population()
    sequences, descriptors = P.F.P.read_train_DNA()
    out = OUT/f'shard_{shard}'
    P.write(out/'started.json', dict(configuration=CONFIG, shard=shard,
        job_id=os.environ['SLURM_JOB_ID'], array_job_id=os.environ['SLURM_ARRAY_JOB_ID']))
    records = []
    with torch.inference_mode(), (out/'windows.jsonl').open('x') as handle:
        for ordinal, expected in enumerate(unique):
            if ordinal % SHARDS != shard:
                continue
            started = time.perf_counter()
            w = windows[expected['window_id']]
            chromosome = sequences[w['species']][w['seqid']]
            geometry = P.F.G.context_geometry(w['start'], w['end'], len(chromosome))
            a, b = geometry['genomic_context']
            raw = chromosome[w['start']:w['end']]
            core = raw if w['strand'] == '+' else P.F.Q.M.rc(raw)
            targets = P.O.build_targets(w, catalogs[w['species']], core)
            n = w['valid_bases']
            fine, region = (torch.as_tensor(targets[k])[:n] for k in ('seg', 'region'))
            counts = (int((fine >= 0).sum()), int(((fine < 0) & (region >= 0)).sum()),
                      int(((fine < 0) & (region < 0)).sum()))
            offset = w['start']-a if w['strand'] == '+' else b-w['end']
            if (counts != tuple(expected[k] for k in ('fine_positions', 'coarse_positions', 'fully_unknown_positions'))
                    or [a,b] != expected['genomic_context'] or [offset,offset+n] != expected['core_oriented_interval']):
                raise ValueError('Fixed fresh observations or original view differ')
            context = re.sub('[^ACGT]', 'N', chromosome[a:b])
            if w['strand'] == '-':
                context = P.F.Q.M.rc(context)
            row = dict(ordinal=ordinal, window_id=w['id'], species=w['species'],
                genomic_context=[a,b], core_oriented_interval=[offset,offset+n],
                fine_positions=counts[0], coarse_positions=counts[1], fully_unknown_positions=counts[2],
                feasible=True, source='NEW_EXACT_FINE_EVENT_FORWARD', source_job=os.environ['SLURM_JOB_ID'])
            if w['id'] in old:
                prior, job = old[w['id']]
                if ([a,b] != prior['genomic_context'] or [offset,offset+n] != prior['core_oriented_interval']
                        or counts != tuple(prior[k] for k in ('fine_positions','coarse_positions','fully_unknown_positions'))):
                    raise ValueError('Reused original certificate differs from actual fresh event')
                row.update(source='REUSED_ORIGINAL_FINITE_EVENT', source_job=job)
            else:
                try:
                    fine_event(context, offset, n, targets)
                except InfeasibleObservation as error:
                    row.update(feasible=False, reason=str(error))
            row['seconds'] = time.perf_counter()-started
            records.append(row)
            handle.write(json.dumps(row, allow_nan=False)+'\n')
            handle.flush()
            if len(records) % 32 == 0:
                print(json.dumps(dict(shard=shard, checked=len(records), elapsed=time.perf_counter()-began)), flush=True)
    P.write(out/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'],
        array_job_id=os.environ['SLURM_ARRAY_JOB_ID'], shard=shard, configuration=CONFIG,
        checks_job=checked['job_id'], fixture=checked['fixture'], records=records, unique_events=len(records),
        infeasible_events=sum(not r['feasible'] for r in records), original_TRAIN_DNA_sources=descriptors,
        original_labels_changed=False, all_fixed_draws_retained=True,
        model_forward=False, checkpoints_read=False, optimizer_steps_this_run=0,
        DEV_used=False, sealed_data_used=False, seconds=time.perf_counter()-began))


def aggregate(array_job):
    began = time.perf_counter()
    _, unique = population()
    records = []
    source_jobs = []
    for shard in range(SHARDS):
        out = OUT/f'shard_{shard}'
        report = json.loads((out/'summary.json').read_text())
        expected = [r['window_id'] for i,r in enumerate(unique) if i % SHARDS == shard]
        if ((out/'STATUS').read_text().strip() != 'COMPLETED' or not report['completed']
                or report['array_job_id'] != array_job or report['shard'] != shard
                or report['configuration'] != CONFIG or report['optimizer_steps_this_run'] != 0
                or report['DEV_used'] or report['sealed_data_used']
                or [r['window_id'] for r in report['records']] != expected):
            raise ValueError('All16 completed unchanged fixed-population shards required')
        records.extend(report['records'])
        source_jobs.append(report['job_id'])
    records.sort(key=lambda r:r['ordinal'])
    if len(records) != 4219 or [r['window_id'] for r in records] != [r['window_id'] for r in unique]:
        raise ValueError('All4219 fresh events retained exactly once')
    failed = [r['window_id'] for r in records if not r['feasible']]
    P.write(AGGREGATE/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'],
        array_job_id=array_job, source_jobs=source_jobs, configuration=CONFIG, protocol=P.PROTOCOL,
        unique_events=len(records), prospective_draws=4608, all_fixed_draws_retained=True,
        source_counts=dict(Counter(r['source'] for r in records)),
        infeasible_events=len(failed), infeasible_window_ids=failed,
        all_fresh_exact_fine_events_feasible=not failed, formal_fit_eligible_after_native_preflight=not failed,
        original_labels_changed=False, model_forward=False, checkpoints_read=False,
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
        seconds=time.perf_counter()-began,
        interpretation='actual_fixed_TRAIN_event_support_not_fitting_or_performance'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--shard', type=int)
    parser.add_argument('--aggregate')
    parser.add_argument('--checks', action='store_true')
    args = parser.parse_args()
    torch.set_num_threads(2)
    if args.checks:
        checked = fixture()
        P.write(CHECKS/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'],
            configuration=CONFIG, fixture=checked, cases=2, research_data_read=False,
            native_forward=False, optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False))
    elif args.aggregate is not None:
        aggregate(args.aggregate)
    elif args.shard is not None and 0 <= args.shard < SHARDS:
        run(args.shard)
    else:
        parser.error('One fixed shard or whole-array aggregate required')
