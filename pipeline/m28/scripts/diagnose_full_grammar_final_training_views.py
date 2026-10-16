#!/usr/bin/env python3
"""Final4608 fixed-weight original TRAIN view NLL/MAP diagnosis; zero updates."""
import argparse
from collections import Counter
import json
import math
import os
import re
import time

import numpy as np
import torch
import fit_full_grammar_coarse_channel as F
import diagnose_full_grammar_train_energy as E
from src.m28.full_grammar_coarse_channel import coarse_channel_nll
from src.m28.full_grammar_skyline import decode_full_grammar_skyline
from src.m28.labels import REGION_COLUMNS, build_targets

I = E.I
OUT = F.ROOT/'outputs/M28-FULL-GRAMMAR-FINAL-TRAIN-VIEWS'
CHECKS = F.ROOT/'outputs/M28-FULL-GRAMMAR-FINAL-TRAIN-VIEW-CHECKS'
OBSERVATIONS = F.ROOT/'outputs/M28-FULL-GRAMMAR-TRAIN-OBSERVATION-DIAGNOSTIC'
REGION = np.empty(15, dtype=np.int64)
for region, columns in enumerate(REGION_COLUMNS):
    REGION[list(columns)] = region


def write(path, data):
    with path.open('x') as handle:
        json.dump(data, handle, indent=2, allow_nan=False)
        handle.write('\n')


def map_roles(chains, length):
    roles = np.zeros(length, dtype=np.int64)
    previous = 0
    for chain in chains:
        a, b = chain[0][0], chain[-1][1]
        if not previous <= a < b <= length:
            raise ValueError('Finite-view selected chains must be nonoverlapping and in bounds')
        roles[a:b] = E.role_path(chain)
        previous = b
    return roles


def core_agreement(roles, targets):
    fine = np.asarray(targets['seg'])[:len(roles)]
    region = np.asarray(targets['region'])[:len(roles)]
    exact = fine >= 0
    coarse = (~exact) & (region >= 0)
    fine_table = np.bincount(15*fine[exact]+roles[exact], minlength=225).reshape(15, 15)
    coarse_table = np.bincount(3*region[coarse]+REGION[roles[coarse]], minlength=9).reshape(3, 3)
    return dict(fine_confusion=fine_table.tolist(), coarse_confusion=coarse_table.tolist(),
        fine_positions=int(exact.sum()), fine_matches=int(np.trace(fine_table)),
        coarse_positions=int(coarse.sum()), coarse_matches=int(np.trace(coarse_table)),
        fully_unknown_positions=int((~exact & ~coarse).sum()))


def checks():
    from check_full_grammar_coarse_channel import reference, targets
    began = time.perf_counter()
    sequence = 'CCCATGAAATAACCCATGCCCTAGCCC'
    expected = np.zeros(len(sequence), dtype=np.int64)
    expected[3:6] = 13
    expected[6:9] = [1, 2, 3]
    expected[9:12] = 14
    logits = torch.full((len(sequence), 15), -4., dtype=torch.float64)
    logits[:, 0] = 0.
    for position in range(3, 12):
        logits[position, expected[position]] = 2.
    with torch.no_grad():
        result = decode_full_grammar_skyline(sequence, logits, 2, 3)
        actual = map_roles(result['chains'], len(sequence))
        assert np.array_equal(actual, expected)
        assert result['gain'] == 18.
        t = targets(21)
        t['seg'][0], t['seg'][1], t['seg'][9] = 13, 14, 0
        t['region'][1] = 1  # Fine takes precedence over this coarse label.
        t['region'][4], t['region'][16], t['region'][20] = 1, 1, 0
        observed = core_agreement(actual[3:24], t)
        assert (observed['fine_positions'], observed['fine_matches'],
                observed['coarse_positions'], observed['coarse_matches'],
                observed['fully_unknown_positions']) == (3, 2, 3, 2, 15)
        # Use feasible fine observations for independently enumerated NLL.
        t['seg'][1] = 13
        nll = coarse_channel_nll(sequence, logits, 2, 3, 24, t, .1, cuts=(3, 24), checkpoint_bp=3)
        enumerated = reference(logits, [(3, 12), (15, 24)], 3, 24, t, torch.tensor(.1, dtype=torch.float64))
        assert math.isclose(float(nll['loss']), float(enumerated), rel_tol=0, abs_tol=1e-9)
    assert I.R.genomic_chain(((3, 12),), 100, 127, '-') == [(115, 124)]
    write(CHECKS/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'],
        checks='independent_single_gene_MAP_roles_fine_first_coarse_confusion_enumerated_NLL_minus_coordinates',
        research_data_read=False, checkpoint_read=False, native_forward=False,
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
        seconds=time.perf_counter()-began))


def source_inputs():
    prep = json.loads((F.READY/'summary.json').read_text())
    fit = json.loads((F.FIT/'summary.json').read_text())
    obs = json.loads((OBSERVATIONS/'summary.json').read_text())
    if (any((root/'STATUS').read_text().strip() != 'COMPLETED' for root in (F.READY, F.FIT, OBSERVATIONS))
            or not all(s['completed'] for s in (prep, fit, obs))
            or prep['configuration'] != F.CONFIG or fit['configuration'] != F.CONFIG
            or fit['optimizer_steps'] != 4608 or not obs['original_targets_reconstructed']
            or obs['original_targets_changed'] or obs['DEV_used'] or obs['sealed_data_used']):
        raise ValueError('Completed unchanged original preparation, final fit and observations required')
    plan = [json.loads(line) for line in (F.READY/'training_plan.jsonl').read_text().splitlines()]
    draws, _, windows, catalogs = F.Q.load_train_inputs()
    if len(plan) != 1536 or len(windows) != 1491 or [d['window_id'] for d in draws] != [r['window_id'] for r in plan]:
        raise ValueError('Original TRAIN order and population required')
    return prep, fit, obs, plan, windows, catalogs


@torch.no_grad()
def task(number):
    species, seqid, length = I.SCOPE[number]
    out = OUT/seqid
    began = time.perf_counter()
    checked = json.loads((CHECKS/'summary.json').read_text())
    if not checked['completed'] or (CHECKS/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Completed new view diagnosis checks required')
    prep, fit, obs, plan, windows, catalogs = source_inputs()
    order = [w for w in windows.values() if (w['species'], w['seqid']) == (species, seqid)]
    expected = {r['window_id']:r for r in plan}
    multiplicity = Counter(r['window_id'] for r in plan)
    importance = Counter()
    for r in plan:
        importance[r['window_id']] += r['importance_weight']
    sequences, descriptors = F.P.read_train_DNA()
    chromosome = sequences[species][seqid]
    if len(chromosome) != length:
        raise ValueError('Original TRAIN chromosome length changed')
    del sequences
    references = [r for line in (E.PROFILE/'reference_chain_profile.jsonl').read_text().splitlines()
                  for r in [json.loads(line)] if r['seqid'] == seqid]
    score_rows = [json.loads(line) for line in (E.OUT/seqid/'reference_energy.jsonl').read_text().splitlines()]
    supported = {(r['id'], r['strand']) for r in score_rows if r['single_chain_supported']}
    if len(references) != len(score_rows):
        raise ValueError('Frozen reference population changed')
    write(out/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], task=number,
        species=species, seqid=seqid, length=length, checkpoint=str(I.CHECKPOINT),
        unique_windows=len(order), scope='original_TRAIN_fixed_final_finite_views',
        DEV_used=False, sealed_data_used=False, optimizer_steps_this_run=0))
    saved = torch.load(I.CHECKPOINT, map_location='cpu', weights_only=True)
    configuration = I.R.validate_checkpoint(saved)
    if configuration != F.CONFIG:
        raise ValueError('Final checkpoint differs from actual training configuration')
    model = I.R.M28Core('C0', 2048, 128)
    model.load_state_dict(saved['model'], strict=True)
    del saved
    model.eval().requires_grad_(False).to('cuda')
    backbone, tokenizer = F.G.load_native()
    torch.cuda.reset_peak_memory_stats()
    with (out/'windows.jsonl').open('x') as handle, (out/'reference_views.jsonl').open('x') as joins:
        for ordinal, w in enumerate(order):
            started = time.perf_counter()
            original = expected[w['id']]
            n = w['valid_bases']
            geometry = F.G.context_geometry(w['start'], w['end'], length)
            a, b = geometry['genomic_context']
            raw = chromosome[w['start']:w['end']]
            core = raw if w['strand'] == '+' else F.Q.M.rc(raw)
            targets = build_targets(w, catalogs[species], core)
            offset = w['start']-a if w['strand'] == '+' else b-w['end']
            if [a, b] != original['genomic_context'] or [offset, offset+n] != original['core_oriented_interval']:
                raise ValueError('Original finite view geometry changed')
            context = re.sub('[^ACGT]', 'N', chromosome[a:b])
            if w['strand'] == '-':
                context = F.Q.M.rc(context)
            tensors = []
            for left, right in geometry['genomic_windows']:
                dna_string = re.sub('[^ACGT]', 'N', chromosome[left:right])
                if w['strand'] == '-':
                    dna_string = F.Q.M.rc(dna_string)
                features = F.G.encode_window(backbone, tokenizer, dna_string, 'cuda').clone()
                dna = F.G.one_hot(dna_string, F.Q.M.WIDTH)[None].to('cuda')
                tensors.append(model(features, dna, [right-left])['segmentation_logits'][0, :right-left])
            joined = F.G.owned_context_logits(tensors, geometry, w['strand']).cpu()
            cuts = tuple(x for x in (offset, offset+n) if 0 < x < len(context))
            likelihood = coarse_channel_nll(context, joined, 128, offset, offset+n, targets, .1,
                cuts=cuts, checkpoint_bp=256)
            masses = {name:float(likelihood[name]) for name in
                      ('loss', 'log_partition_all', 'log_partition_observed', 'background_correction')}
            if not all(math.isfinite(v) for v in masses.values()) or masses['loss'] < -1e-8:
                raise ValueError('Finite nonnegative unchanged observation NLL required')
            # MAP uses only raw view emissions and DNA, never observations or references.
            decoded = decode_full_grammar_skyline(context, joined, 128)
            roles = map_roles(decoded['chains'], len(context))
            agreement = core_agreement(roles[offset:offset+n], targets)
            if any(agreement[name] != original[name] for name in
                   ('fine_positions', 'coarse_positions', 'fully_unknown_positions')):
                raise ValueError('Actual observation census differs from the original executed plan')
            selected = {tuple(I.R.genomic_chain(chain, a, b, w['strand'])) for chain in decoded['chains']}
            fine = np.asarray(targets['seg'])[:n]
            for ref in references:
                if ref['strand'] != w['strand'] or (ref['id'], ref['strand']) not in supported:
                    continue
                start, end = ref['CDS'][0][0], ref['CDS'][-1][1]
                if not a <= start < end <= b:
                    continue
                parts = E.oriented_parts(ref, length)
                global_offset = a if w['strand'] == '+' else length-b
                parts = [(x-global_offset, y-global_offset) for x, y in parts]
                x, y = parts[0][0], parts[-1][1]
                lo, hi = max(x, offset), min(y, offset+n)
                expected_roles = E.role_path(parts)
                observed_fine = fine[lo-offset:hi-offset] if lo < hi else np.empty(0, dtype=np.int64)
                desired = expected_roles[lo-x:hi-x] if lo < hi else np.empty(0, dtype=np.int64)
                matches = int(np.count_nonzero((observed_fine >= 0) & (observed_fine == desired)))
                conflicts = int(np.count_nonzero((observed_fine >= 0) & (observed_fine != desired)))
                record = dict(species=species, seqid=seqid, id=ref['id'], strand=w['strand'],
                    window_id=w['id'], reference_span_rows=y-x, matching_core_fine_rows=matches,
                    disagreeing_core_fine_rows=conflicts,
                    whole_span_matching_fine_in_this_core=matches == y-x and conflicts == 0,
                    exact_chain_in_unrestricted_view_MAP=tuple((x, y) for x, y, _ in ref['CDS']) in selected)
                joins.write(json.dumps(record, allow_nan=False)+'\n')
            row = dict(window_id=w['id'], species=species, seqid=seqid, strand=w['strand'],
                original_draw_multiplicity=multiplicity[w['id']], original_importance_sum=importance[w['id']],
                known_observation_denominator=original['known_observation_denominator'],
                genomic_context=[a, b],
                core_oriented_interval=[offset, offset+n], **agreement,
                raw_NLL=masses.pop('loss'), **masses, seconds=time.perf_counter()-started)
            row['normalized_NLL'] = row['raw_NLL']/original['known_observation_denominator']
            row['weighted_total'] = row['normalized_NLL']*importance[w['id']]/multiplicity[w['id']]
            handle.write(json.dumps(row, allow_nan=False)+'\n')
            if (ordinal+1) % 16 == 0 or ordinal+1 == len(order):
                handle.flush(); joins.flush()
                print(json.dumps(dict(stage='original_views', completed=ordinal+1, total=len(order),
                    recent_window_id=w['id'], recent_seconds=row['seconds'])), flush=True)
            del tensors, features, dna, joined, likelihood, decoded, roles, targets
    report = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], task=number, species=species,
        seqid=seqid, length=length, unique_windows=len(order), checkpoint=str(I.CHECKPOINT),
        configuration=configuration, original_TRAIN_DNA_sources=descriptors,
        source_preparation_job=prep['job_id'], source_fit_job=fit['job_id'], source_observation_job=obs['job_id'],
        scope='original_TRAIN_fixed_final_finite_views', native_forward=True,
        original_targets_changed=False, observation_model_changed=False,
        labels_used_for_NLL_and_posthoc_agreement=True, labels_or_references_used_by_encoder_head_or_MAP=False,
        cached_features_used=False, potentials_saved=False, predictions_modified=False,
        optimizer_constructed=False, optimizer_steps_this_run=0, checkpoints_written=False,
        DEV_used=False, sealed_data_used=False,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(), seconds=time.perf_counter()-began,
        interpretation='fixed_final_original_TRAIN_finite_view_fit_diagnosis_not_new_training_or_primary_evaluation')
    write(out/'summary.json', report)
    print(json.dumps(report, allow_nan=False), flush=True)


def totals(rows, repeated):
    count = sum(r['original_draw_multiplicity'] if repeated else 1 for r in rows)
    weights = [r['original_draw_multiplicity'] if repeated else 1 for r in rows]
    fine = sum((np.asarray(r['fine_confusion'], dtype=np.int64)*w for r, w in zip(rows, weights)), np.zeros((15, 15), dtype=np.int64))
    coarse = sum((np.asarray(r['coarse_confusion'], dtype=np.int64)*w for r, w in zip(rows, weights)), np.zeros((3, 3), dtype=np.int64))
    return dict(views_or_draw_exposures=count,
        mean_loss={name:math.fsum(r[name]*w for r, w in zip(rows, weights))/count
                   for name in ('raw_NLL', 'normalized_NLL', 'weighted_total')},
        fine_confusion=fine.tolist(), fine_positions=int(fine.sum()), fine_matches=int(np.trace(fine)),
        coarse_confusion=coarse.tolist(), coarse_positions=int(coarse.sum()), coarse_matches=int(np.trace(coarse)))


def complete_core_map_pattern(reference):
    if not reference['single_chain_supported']:
        return 'unsupported_reference_not_role_joined'
    observed = reference['complete_matching_core_views']
    exact = reference['exact_MAP_complete_matching_core_views']
    if not observed:
        return 'no_single_complete_matching_core_view'
    if not exact:
        return 'no_complete_matching_core_view_MAP_exact'
    if exact == observed:
        return 'all_complete_matching_core_views_MAP_exact'
    return 'some_but_not_all_complete_matching_core_views_MAP_exact'


def aggregate():
    began = time.perf_counter()
    _, _, _, plan, windows, _ = source_inputs()
    rows, joins, sources = [], [], {}
    for number, (species, seqid, length) in enumerate(I.SCOPE):
        root = OUT/seqid
        report = json.loads((root/'summary.json').read_text())
        if ((root/'STATUS').read_text().strip() != 'COMPLETED' or not report['completed']
                or report['task'] != number or (report['species'], report['seqid'], report['length']) != (species, seqid, length)
                or report['checkpoint'] != str(I.CHECKPOINT) or report['DEV_used'] or report['sealed_data_used']
                or report['optimizer_steps_this_run'] != 0):
            raise ValueError('All six completed unchanged final-view sources required')
        rows.extend(json.loads(line) for line in (root/'windows.jsonl').read_text().splitlines())
        joins.extend(json.loads(line) for line in (root/'reference_views.jsonl').read_text().splitlines())
        sources[seqid] = report['job_id']
    by_id = {r['window_id']:r for r in rows}
    if len(rows) != 1491 or len(by_id) != 1491 or set(by_id) != set(windows):
        raise ValueError('Every original unique core exactly once required')
    exposure = Counter(r['window_id'] for r in plan)
    importance = Counter()
    for r in plan:
        importance[r['window_id']] += r['importance_weight']
    if any(r['original_draw_multiplicity'] != exposure[r['window_id']]
           or r['original_importance_sum'] != importance[r['window_id']] for r in rows):
        raise ValueError('Original repeated draw exposure required')
    refrows = [json.loads(line) for line in (OBSERVATIONS/'reference_observations.jsonl').read_text().splitlines()]
    key = lambda r: (r['seqid'], r['id'], r['strand'])
    stats = {key(r):dict(species=r['species'], seqid=r['seqid'], id=r['id'], strand=r['strand'],
        whole_contig_category=r['category'], original_union_coverage=r['coverage_class'],
        single_chain_supported=r['single_chain_supported'], contained_original_views=0,
        exact_MAP_contained_views=0, complete_matching_core_views=0, exact_MAP_complete_matching_core_views=0)
        for r in refrows}
    if len(stats) != 22688:
        raise ValueError('Retain all22688 original frozen references')
    for r in joins:
        s = stats[key(r)]
        if not s['single_chain_supported'] or r['disagreeing_core_fine_rows']:
            raise ValueError('Unexpected unsupported join or original known-role disagreement')
        s['contained_original_views'] += 1
        s['exact_MAP_contained_views'] += r['exact_chain_in_unrestricted_view_MAP']
        s['complete_matching_core_views'] += r['whole_span_matching_fine_in_this_core']
        s['exact_MAP_complete_matching_core_views'] += (r['whole_span_matching_fine_in_this_core']
                                                      and r['exact_chain_in_unrestricted_view_MAP'])
    # Ownership13450203 has zero source mismatches on these original cores.
    # If an entire chain is in one core, deleting a negative-gain selected
    # chain for background is feasible under the unchanged gene-set grammar.
    # This is a diagnostic for replay review, not a new model or PASS gate.
    core_counterexamples = [dict(species=s['species'], seqid=s['seqid'], id=s['id'], strand=s['strand'],
        complete_matching_core_views=s['complete_matching_core_views'],
        exact_MAP_complete_matching_core_views=s['exact_MAP_complete_matching_core_views'])
        for s in stats.values()
        if s['whole_contig_category'] == 'supported_reference_negative_gain_vs_background'
        and s['exact_MAP_complete_matching_core_views'] > 0]
    with (OUT/'reference_view_coverage.jsonl').open('x') as handle:
        for s in stats.values():
            s['complete_core_view_MAP_pattern'] = complete_core_map_pattern(s)
            handle.write(json.dumps(s, allow_nan=False)+'\n')
    def groups(species=None):
        refs = [s for s in stats.values() if species is None or s['species'] == species]
        return {name:dict(reference=len(group), whole_contig_categories=dict(Counter(s['whole_contig_category'] for s in group)),
            references_with_any_exact_contained_view=sum(s['exact_MAP_contained_views'] > 0 for s in group),
            references_with_any_exact_complete_matching_core_view=sum(s['exact_MAP_complete_matching_core_views'] > 0 for s in group),
            complete_matching_core_view_exposures=sum(s['complete_matching_core_views'] for s in group),
            exact_MAP_complete_matching_core_view_exposures=sum(s['exact_MAP_complete_matching_core_views'] for s in group),
            complete_core_view_MAP_patterns={pattern:dict(reference=len(part),
                whole_contig_categories=dict(Counter(s['whole_contig_category'] for s in part)))
                for pattern in sorted({s['complete_core_view_MAP_pattern'] for s in group})
                for part in [[s for s in group if s['complete_core_view_MAP_pattern'] == pattern]]})
            for name, condition in (
                ('all_frozen_references', lambda s: True),
                ('supported_contained_in_any_original_view', lambda s: s['contained_original_views'] > 0),
                ('whole_span_matching_fine_in_at_least_one_original_core', lambda s: s['complete_matching_core_views'] > 0),
                ('union_matching_span_without_single_complete_matching_core', lambda s:
                 s['original_union_coverage'] == 'whole_reference_span_matching_fine_no_fine_disagreement'
                 and s['complete_matching_core_views'] == 0))
            for group in [[s for s in refs if condition(s)]]}
    report = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], source_jobs=sources,
        checkpoint=str(I.CHECKPOINT), scope='all_original_1491_unique_TRAIN_views_and1536_draw_exposures',
        pooled=dict(unique=totals(rows, False), per_pass_draw_exposure=totals(rows, True), references=groups()),
        per_species={species:dict(unique=totals([r for r in rows if r['species'] == species], False),
            per_pass_draw_exposure=totals([r for r in rows if r['species'] == species], True),
            references=groups(species)) for species in F.Q.ALLOW},
        original_targets_changed=False, predictions_modified=False, checkpoint_read=False, native_forward=False,
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False, seconds=time.perf_counter()-began,
        fine_coarse_MAP_agreement_is_not_observed_event_probability_or_primary_chain_metric=True,
        any_view_exact_is_exploratory_not_a_new_reference_selected_prediction_set=True,
        same_core_emission_background_dominance_counterexamples=core_counterexamples,
        mechanism_interpretation_requires_core_emission_counterexample_review=bool(core_counterexamples),
        interpretation='fixed_final_original_TRAIN_finite_view_fit_diagnosis_not_new_training_or_primary_evaluation')
    write(OUT/'summary.json', report)
    print(json.dumps(report, allow_nan=False), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('checks', 'task', 'aggregate'))
    parser.add_argument('--task', type=int)
    args = parser.parse_args()
    torch.set_num_threads(2)
    if args.mode == 'checks':
        checks()
    elif args.mode == 'task':
        if args.task is None or not 0 <= args.task < len(I.SCOPE):
            parser.error('Original six TRAIN task required')
        task(args.task)
    else:
        aggregate()
