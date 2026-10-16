#!/usr/bin/env python3
"""Freeze three fresh TRAIN draw schedules and measure their actual supervision."""
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
import json
import os
import time

import numpy as np
import torch
import diagnose_full_grammar_training_observations as O
import query_training_census as Q

F, E = O.F, O.E
ROOT = F.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-FRESH-DRAWS-PREPARATION'
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_fresh_draws_plan.md'
CONFIG = dict(scope='original_two_species_six_TRAIN_chromosomes_existing25512_R4_windows',
    sampling='uniform_with_replacement_within_original_stratum_fresh_each_pass',
    sampling_seed=1, passes=3, draws_per_pass=1536, additional_optimizer_steps=4608,
    stratum_order_and_counts='same_as_original1536_draw_schedule_each_pass',
    importance_weight_uniform_within_stratum=1., labels='unchanged_R4_build_targets',
    starting_fit_job='13434603', starting_optimizer_step=4608,
    prospective_final_optimizer_step=9216,
    model_and_loss='unchanged_C0_and_structured_coarse_channel_epsilon0.1',
    source_fit_configuration=F.CONFIG,
    context='unchanged_core24576_plus_one_stride12288_neighbor_each_side',
    formal_checkpoint='prospective_last9216_only_no_DEV_selection',
    train_scope_expanded=False, DEV_used=False, sealed_data_used=False)


def rows(path):
    with path.open() as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def write(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def stratum(row):
    return row['species'], row['seqid'], row['strand']


def fresh_schedule(original, strata, passes=3, seed=1):
    """Only proposal draws change; retain original per-update stratum order."""
    rng = np.random.default_rng(seed)
    schedule = []
    for epoch in range(passes):
        for ordinal, old in enumerate(original):
            key = old['species'], *old['stratum']
            pool = strata[key]
            window = pool[int(rng.integers(len(pool)))]
            schedule.append(dict(additional_step=len(schedule)+1, pass_index=epoch+1,
                draw_ordinal=ordinal, window_id=window['id'], species=key[0],
                stratum=list(key[1:]), stratum_size=len(pool),
                stratum_probability=1./len(pool), importance_weight=1.))
    return schedule


def sampler_check():
    # Two unequal pools: unweighted uniform means must equal the
    # importance-corrected gene-enriched means, not the enriched raw means.
    values = np.array([2., 11., -4.])
    q = np.array([.1, .2, .7])
    corrected = np.sum(q * (1./(len(q)*q)) * values)
    assert np.isclose(corrected, np.mean(values))
    assert not np.isclose(np.sum(q*values), np.mean(values))
    old = [dict(species='a', stratum=['q', '+']),
           dict(species='b', stratum=['r', '-'])]*4
    pools = {('a', 'q', '+'): [dict(id='a'+str(i)) for i in range(3)],
             ('b', 'r', '-'): [dict(id='b'+str(i)) for i in range(5)]}
    result = fresh_schedule(old, pools)
    assert len(result) == 24
    assert [r['additional_step'] for r in result] == list(range(1, 25))
    for epoch in range(3):
        part = result[epoch*8:(epoch+1)*8]
        assert [(r['species'], r['stratum']) for r in part] == [
            (r['species'], r['stratum']) for r in old]
        assert all(r['importance_weight'] == 1. for r in part)
    assert result == fresh_schedule(old, pools)
    assert len({tuple(r['window_id'] for r in result[e*8:(e+1)*8]) for e in range(3)}) > 1
    return 'unequal_pool_expectation_and_fresh_pass_stratum_schedule'


def load_population():
    windows = {}
    strata = defaultdict(list)
    catalogs = {}
    lengths = {(s, q):n for s, q, n in E.I.SCOPE}
    for species in Q.ALLOW:
        local = list(Q.train_records(Q.LABELS/species/'windows_train.jsonl'))
        for w in local:
            n = w['valid_bases']
            if (w['species'] != species or w['seqid'] not in Q.ALLOW[species]
                    or w['strand'] not in ('+', '-') or w['id'] in windows
                    or w['label_policy'] != 'CDS_assessable_primary_R4'
                    or w['width'] != 24576 or w['end']-w['start'] != n
                    or not 0 <= w['start'] < w['end'] <= lengths[species, w['seqid']]
                    or not 0 < n <= 24576):
                raise ValueError('Existing original-six TRAIN R4 window identity differs')
            windows[w['id']] = w
            strata[stratum(w)].append(w)
        needed = {tid for w in local for tid in w['overlapping_transcript_ids']}
        catalog = {}
        for tx in Q.train_records(Q.LABELS/species/'transcripts.jsonl'):
            if tx['id'] in needed:
                if tx['id'] in catalog or tx['seqid'] not in Q.ALLOW[species]:
                    raise ValueError('Original TRAIN transcript identity differs')
                catalog[tx['id']] = tx
        if set(catalog) != needed:
            raise ValueError('All original TRAIN transcript joins required')
        for w in local:
            if not set(w['positive_ids']) <= set(w['overlapping_transcript_ids']):
                raise ValueError('Positive metadata must join overlapping TRAIN transcripts')
            if any((catalog[tid]['seqid'], catalog[tid]['strand']) != (w['seqid'], w['strand'])
                   for tid in w['overlapping_transcript_ids']):
                raise ValueError('TRAIN transcript/window strand or chromosome differs')
        catalogs[species] = catalog
    original = list(rows(Q.DRAWS))
    expected = {(s, q, strand) for s, q, _ in E.I.SCOPE for strand in ('+', '-')}
    if (len(windows) != 25512 or set(strata) != expected or len(original) != 1536
            or len({d['window_id'] for d in original}) != 1491
            or Counter(d['species'] for d in original) != Counter({s:768 for s in Q.ALLOW})):
        raise ValueError('Full natural TRAIN pool and original draw schedule differ')
    for d in original:
        if stratum(windows[d['window_id']]) != (d['species'], *d['stratum']):
            raise ValueError('Original draw stratum differs')
    for pool in strata.values():
        pool.sort(key=lambda w:w['start'])
    return windows, strata, catalogs, original


def main():
    began = time.perf_counter()
    fixture = sampler_check()
    original_summary = json.loads((O.OUT/'summary.json').read_text())
    if ((O.OUT/'STATUS').read_text().strip() != 'COMPLETED'
            or original_summary['job_id'] != '13450304' or not original_summary['completed']
            or original_summary['optimizer_steps_this_run'] != 0
            or original_summary['DEV_used'] or original_summary['sealed_data_used']
            or original_summary['pooled']['reference'] != 22688
            or (E.PROFILE/'STATUS').read_text().strip() != 'COMPLETED'):
        raise ValueError('Completed original observation and fixed reference sources required')
    windows, strata, catalogs, original = load_population()
    schedule = fresh_schedule(original, strata)
    counts = Counter(r['window_id'] for r in schedule)
    by_stratum = Counter((r['species'], *r['stratum']) for r in schedule)
    if by_stratum != Counter({k:3*v for k, v in Counter(
            (r['species'], *r['stratum']) for r in original).items()}):
        raise ValueError('Fixed per-stratum update budget differs')
    write(OUT/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        protocol=PROTOCOL, selected_unique_windows=len(counts)))
    sequences, descriptors = F.P.read_train_DNA()
    reference_key = lambda r:(r['seqid'], r['id'], r['strand'])
    profiles = {reference_key(r):r for r in rows(E.PROFILE/'reference_chain_profile.jsonl')}
    old_refs = {reference_key(r):{k:r[k] for k in (
        'species', 'seqid', 'id', 'strand', 'coverage_class', 'single_chain_supported')}
        for r in rows(O.OUT/'reference_observations.jsonl')}
    if set(profiles) != set(old_refs) or len(profiles) != 22688:
        raise ValueError('Preserve all fixed original primary references')
    selected = defaultdict(list)
    for ident in counts:
        selected[stratum(windows[ident])].append(windows[ident])
    for group in selected.values():
        group.sort(key=lambda w:w['start'])
    exposure = {s:Counter() for s in Q.ALLOW}
    fine_hist = {s:np.zeros(15, dtype=np.int64) for s in Q.ALLOW}
    coarse_hist = {s:np.zeros(3, dtype=np.int64) for s in Q.ALLOW}
    plan_by_id = {}
    reference_rows = []
    with (OUT/'reference_coverage.jsonl').open('x') as handle:
        for species, seqid, length in E.I.SCOPE:
            chromosome = sequences[species][seqid]
            if len(chromosome) != length:
                raise ValueError('Original TRAIN chromosome length differs')
            for strand in ('+', '-'):
                observation = np.zeros(length, dtype=np.uint32)
                fine_by_id = {}
                group = selected[species, seqid, strand]
                starts = [w['start'] for w in group]
                for w in group:
                    raw = chromosome[w['start']:w['end']]
                    core = raw if strand == '+' else F.Q.M.rc(raw)
                    target = O.build_targets(w, catalogs[species], core)
                    n = w['valid_bases']
                    fine, region = target['seg'][:n], target['region'][:n]
                    exact = fine >= 0
                    coarse = (~exact) & (region >= 0)
                    fine_by_id[w['id']] = fine.astype(np.int16)
                    bits = O.encode(fine, region)
                    observation[w['start']:w['end']] |= bits if strand == '+' else bits[::-1]
                    copies = counts[w['id']]
                    fine_hist[species] += np.bincount(fine[exact], minlength=15)*copies
                    coarse_hist[species] += np.bincount(region[coarse], minlength=3)*copies
                    stats = exposure[species]
                    stats['unknown_exposures'] += int((~exact & ~coarse).sum())*copies
                    stats['unique_windows'] += 1
                    stats['valid_base_exposures'] += n*copies
                    geometry = F.G.context_geometry(w['start'], w['end'], length)
                    a, b = geometry['genomic_context']
                    offset = w['start']-a if strand == '+' else b-w['end']
                    denominator = max(1, int(exact.sum()+coarse.sum()))
                    plan_by_id[w['id']] = dict(fine_positions=int(exact.sum()),
                        coarse_positions=int(coarse.sum()), fully_unknown_positions=int((~exact & ~coarse).sum()),
                        known_observation_denominator=denominator, objective_multiplier=1./denominator,
                        genomic_context=[a, b], core_oriented_interval=[offset, offset+n])
                stats = exposure[species]
                stats['genomic_union_any_observation_rows'] += int(np.count_nonzero(observation))
                stats['genomic_union_any_fine_rows'] += int(np.count_nonzero(observation & O.FINE))
                stats['genomic_union_any_coarse_rows'] += int(np.count_nonzero(observation & O.COARSE))
                for key, old in old_refs.items():
                    if (old['species'], old['seqid'], old['strand']) != (species, seqid, strand):
                        continue
                    row = dict(old, original_coverage_class=old['coverage_class'], span=None,
                        full_fine_single_core_unique_windows=0, full_fine_single_core_draw_exposures=0)
                    if old['single_chain_supported']:
                        ref = profiles[key]
                        parts = E.oriented_parts(ref, length)
                        a, b = parts[0][0], parts[-1][1]
                        roles = E.role_path(parts)
                        bits = observation[a:b] if strand == '+' else observation[length-b:length-a][::-1]
                        span = O.observation_counts(bits, roles)
                        row.update(span=span, coverage_class=O.coverage_class(span))
                        lo, hi = ref['CDS'][0][0], ref['CDS'][-1][1]
                        for w in group[bisect_left(starts, hi-24576):bisect_right(starts, lo)]:
                            if w['end'] < hi:
                                continue
                            offset = lo-w['start'] if strand == '+' else w['end']-hi
                            if np.array_equal(fine_by_id[w['id']][offset:offset+len(roles)], roles):
                                row['full_fine_single_core_unique_windows'] += 1
                                row['full_fine_single_core_draw_exposures'] += counts[w['id']]
                    reference_rows.append(row)
                    handle.write(json.dumps(row, allow_nan=False)+'\n')
                del observation, fine_by_id
    with (OUT/'training_plan.jsonl').open('x') as handle:
        for draw in schedule:
            handle.write(json.dumps(dict(draw, **plan_by_id[draw['window_id']]), allow_nan=False)+'\n')
    summaries = {}
    for label, group in [('pooled', reference_rows)]+[
            (s, [r for r in reference_rows if r['species'] == s]) for s in Q.ALLOW]:
        summaries[label] = dict(reference=len(group), supported=sum(r['single_chain_supported'] for r in group),
            original_coverage=dict(Counter(r['original_coverage_class'] for r in group)),
            fresh_coverage=dict(Counter(r['coverage_class'] for r in group)),
            coverage_transitions=[dict(original=a, fresh=b, reference=n) for (a, b), n in sorted(Counter(
                (r['original_coverage_class'], r['coverage_class']) for r in group).items())],
            any_full_fine_single_core_reference=sum(r['full_fine_single_core_unique_windows'] > 0 for r in group),
            fine_disagreement_reference=sum(bool(r['span'] and r['span']['disagreeing_fine_rows']) for r in group),
            coarse_disagreement_reference=sum(bool(r['span'] and r['span']['disagreeing_coarse_rows']) for r in group))
    if (len(reference_rows) != 22688 or len({reference_key(r) for r in reference_rows}) != 22688
            or summaries['pooled']['supported'] != 22479
            or len(schedule) != 4608 or set(plan_by_id) != set(counts)):
        raise ValueError('Reference denominator or prospective draw population differs')
    for species in Q.ALLOW:
        stats = exposure[species]
        if int(fine_hist[species].sum()+coarse_hist[species].sum())+stats['unknown_exposures'] != stats['valid_base_exposures']:
            raise ValueError('Every prospective core row must remain fine, coarse or unknown')
    report = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        protocol=PROTOCOL, fixture=fixture, source_observation_job='13450304',
        original_unique_windows=1491, fresh_unique_windows=len(counts),
        newly_sampled_unique_windows=len(set(counts)-{d['window_id'] for d in original}),
        existing_natural_TRAIN_windows=len(windows), prospective_draws=len(schedule),
        by_pass=[dict(pass_index=e+1, unique_windows=len({r['window_id'] for r in schedule if r['pass_index'] == e+1})) for e in range(3)],
        exposure={s:dict(exposure[s], fine_exposures_by_role=fine_hist[s].tolist(),
            coarse_exposures_by_region=coarse_hist[s].tolist()) for s in Q.ALLOW},
        original_exposure=original_summary['exposure'], reference_coverage=summaries,
        original_TRAIN_DNA_sources=descriptors, original_R4_labels_changed=False,
        primary_reference_changed=False, prospective_fit_submitted=False,
        all_new_exact_fine_events_feasible='NOT_YET_VERIFIED', optimizer_steps_this_run=0,
        feature_tensors_read=False, checkpoint_read=False, native_forward=False,
        DEV_used=False, sealed_data_used=False, predictions_modified=False,
        coverage_is_not_partition_gradient_or_encoder_exposure=True,
        coverage_comparison_is_not_model_performance_or_causal_attribution=True,
        seconds=time.perf_counter()-began)
    write(OUT/'summary.json', report)
    print(json.dumps(dict(completed=True, job_id=report['job_id'], fresh_unique_windows=len(counts),
        newly_sampled_unique_windows=report['newly_sampled_unique_windows'],
        reference_coverage=summaries, seconds=report['seconds']), allow_nan=False), flush=True)


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
