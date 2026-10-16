#!/usr/bin/env python3
"""Reconstruct original TRAIN core observations and join fixed reference roles."""
from collections import Counter
import json
import os
import time

import numpy as np
import torch
import fit_full_grammar_coarse_channel as F
import diagnose_full_grammar_train_energy as E
from src.m28.labels import REGION_COLUMNS, build_targets

OUT = F.ROOT/'outputs/M28-FULL-GRAMMAR-TRAIN-OBSERVATION-DIAGNOSTIC'
GROUPS = dict(coding_body=(1, 2, 3), intron_body=(4, 5, 6),
              splice_coding_boundary=(7, 8, 9, 10, 11, 12), start=(13,), stop=(14,))
FINE = np.uint32((1 << 15)-1)
COARSE = np.uint32(7 << 15)
REGIONS = np.empty(15, dtype=np.uint32)
for label, columns in enumerate(REGION_COLUMNS):
    REGIONS[list(columns)] = label


def encode(fine, region):
    bits = np.zeros(len(fine), dtype=np.uint32)
    exact = fine >= 0
    coarse = (~exact) & (region >= 0)
    bits[exact] = np.left_shift(np.uint32(1), fine[exact].astype(np.uint32))
    bits[coarse] = np.left_shift(np.uint32(1), np.uint32(15)+region[coarse].astype(np.uint32))
    return bits


def observation_counts(bits, roles):
    expected = np.left_shift(np.uint32(1), roles.astype(np.uint32))
    expected_coarse = np.left_shift(np.uint32(1), np.uint32(15)+REGIONS[roles])
    return dict(base_rows=len(bits),
        any_observation_rows=int(np.count_nonzero(bits)),
        matching_fine_rows=int(np.count_nonzero(bits & expected)),
        any_fine_rows=int(np.count_nonzero(bits & FINE)),
        disagreeing_fine_rows=int(np.count_nonzero(bits & FINE & ~expected)),
        compatible_coarse_rows=int(np.count_nonzero(bits & expected_coarse)),
        any_coarse_rows=int(np.count_nonzero(bits & COARSE)),
        disagreeing_coarse_rows=int(np.count_nonzero(bits & COARSE & ~expected_coarse)))


def coverage_class(counts):
    if counts['disagreeing_fine_rows']:
        return 'some_fine_observation_disagrees_with_reference'
    if counts['matching_fine_rows'] == counts['base_rows']:
        return 'whole_reference_span_matching_fine_no_fine_disagreement'
    if counts['matching_fine_rows']:
        return 'partial_matching_fine_no_fine_disagreement'
    if counts['any_coarse_rows']:
        return 'coarse_only_span_observation'
    return 'no_direct_core_span_observation'


def summarize(rows):
    joined = [r for r in rows if r['single_chain_supported']]
    return dict(reference=len(rows), supported_role_join=len(joined),
        unsupported_reference_without_role_join=len(rows)-len(joined),
        categories=dict(Counter(r['category'] for r in rows)),
        by_coverage={name:dict(reference=len(group),
            categories=dict(Counter(r['category'] for r in group)),
            reference_chains_with_coarse_disagreement=sum(r['span']['disagreeing_coarse_rows'] > 0 for r in group),
            role_groups={role:{key:sum(r['role_groups'][role][key] for r in group)
                              for key in group[0]['role_groups'][role]} for role in GROUPS})
            for name in sorted({r['coverage_class'] for r in joined})
            for group in [[r for r in joined if r['coverage_class'] == name]]},
        by_score_category={category:dict(reference=len(group),
            coverage=dict(Counter(r['coverage_class'] for r in group)),
            role_groups={role:{key:sum(r['role_groups'][role][key] for r in group)
                              for key in group[0]['role_groups'][role]} for role in GROUPS})
            for category in sorted({r['category'] for r in joined})
            for group in [[r for r in joined if r['category'] == category]]})


def fixture():
    fine = np.array([1, 1, -100, 0, -100, -100])
    region = np.array([1, 1, 2, 0, -100, 2])
    bits = encode(fine, region)
    bits[1] |= np.uint32(1 << 2)
    roles = np.array([1, 1, 4, 13, 14, 1])
    actual = observation_counts(bits, roles)
    assert actual == dict(base_rows=6, any_observation_rows=5, matching_fine_rows=2,
        any_fine_rows=3, disagreeing_fine_rows=2, compatible_coarse_rows=1,
        any_coarse_rows=2, disagreeing_coarse_rows=1)
    assert coverage_class(actual) == 'some_fine_observation_disagrees_with_reference'
    # A reverse-strand core's local row0 is at the high genomic coordinate.
    genomic = np.zeros(10, dtype=np.uint32)
    genomic[2:8] |= bits[::-1]
    assert np.array_equal(genomic[2:8][::-1], bits)
    return 'fine_first_coarse_union_disagreement_and_minus_coordinate_fixture'


def main():
    began = time.perf_counter()
    checked = fixture()
    prep = json.loads((F.READY/'summary.json').read_text())
    fit = json.loads((F.FIT/'summary.json').read_text())
    energy = json.loads((E.OUT/'summary.json').read_text())
    if (any((root/'STATUS').read_text().strip() != 'COMPLETED' for root in (F.READY, F.FIT, E.OUT))
            or not prep['completed'] or not fit['completed'] or not energy['completed']
            or prep['configuration'] != F.CONFIG or fit['configuration'] != F.CONFIG
            or fit['optimizer_steps'] != 4608 or energy['DEV_used'] or energy['sealed_data_used']):
        raise ValueError('Completed unchanged original fit and score sources required')
    plan = [json.loads(line) for line in (F.READY/'training_plan.jsonl').read_text().splitlines()]
    draws, _, windows, catalogs = F.Q.load_train_inputs()
    if len(plan) != 1536 or len(windows) != 1491 or [d['window_id'] for d in draws] != [r['window_id'] for r in plan]:
        raise ValueError('Original draw order or core population differs')
    expected = {r['window_id']:r for r in plan}
    multiplicity = Counter(r['window_id'] for r in plan)
    sequences, descriptors = F.P.read_train_DNA()
    profiles = [json.loads(line) for line in (E.PROFILE/'reference_chain_profile.jsonl').read_text().splitlines()]
    key = lambda r: (r['seqid'], r['id'], r['strand'])
    references = {key(r):r for r in profiles}
    if len(references) != len(profiles) or len(profiles) != energy['pooled']['reference']:
        raise ValueError('Preserve all unique frozen primary references')
    exposure = {s:dict(fine_per_pass=np.zeros(15, dtype=np.int64), coarse_per_pass=np.zeros(3, dtype=np.int64),
                       unknown_per_pass=0, unique_windows=0, genomic_union_any_observation_rows=0,
                       genomic_union_any_fine_rows=0, genomic_union_any_coarse_rows=0) for s in F.Q.ALLOW}
    all_rows = []
    with (OUT/'reference_observations.jsonl').open('x') as handle:
        for species, seqid, length in E.I.SCOPE:
            chromosome = sequences[species][seqid]
            scores = [json.loads(line) for line in (E.OUT/seqid/'reference_energy.jsonl').read_text().splitlines()]
            if len(chromosome) != length or len(scores) != energy['chromosomes'][seqid]['reference']:
                raise ValueError('Original chromosome and frozen score population differ')
            for strand in ('+', '-'):
                observations = np.zeros(length, dtype=np.uint32)
                for w in windows.values():
                    if (w['species'], w['seqid'], w['strand']) != (species, seqid, strand):
                        continue
                    a, b, n = w['start'], w['end'], w['valid_bases']
                    raw = chromosome[a:b]
                    core = raw if strand == '+' else F.Q.M.rc(raw)
                    target = build_targets(w, catalogs[species], core)
                    fine, region = target['seg'][:n], target['region'][:n]
                    exact = fine >= 0
                    coarse = (~exact) & (region >= 0)
                    counts = (int(exact.sum()), int(coarse.sum()), int((~exact & ~coarse).sum()))
                    if counts != tuple(expected[w['id']][name] for name in ('fine_positions', 'coarse_positions', 'fully_unknown_positions')):
                        raise ValueError('Reconstructed original training labels differ from the executed fit')
                    bits = encode(fine, region)
                    observations[a:b] |= bits if strand == '+' else bits[::-1]
                    stats = exposure[species]
                    copies = multiplicity[w['id']]
                    stats['fine_per_pass'] += np.bincount(fine[exact], minlength=15)*copies
                    stats['coarse_per_pass'] += np.bincount(region[coarse], minlength=3)*copies
                    stats['unknown_per_pass'] += counts[2]*copies
                    stats['unique_windows'] += 1
                stats = exposure[species]
                stats['genomic_union_any_observation_rows'] += int(np.count_nonzero(observations))
                stats['genomic_union_any_fine_rows'] += int(np.count_nonzero(observations & FINE))
                stats['genomic_union_any_coarse_rows'] += int(np.count_nonzero(observations & COARSE))
                for score in scores:
                    if score['strand'] != strand:
                        continue
                    ref = references[key(score)]
                    row = {name:score[name] for name in ('species', 'seqid', 'id', 'strand', 'category',
                           'single_chain_supported', 'exact_primary_hit', 'original_draw_R4_positive_chain')}
                    row.update(coverage_class='unsupported_reference_without_role_join', span=None, role_groups=None)
                    if score['single_chain_supported']:
                        parts = E.oriented_parts(ref, length)
                        a, b = parts[0][0], parts[-1][1]
                        bits = observations[a:b] if strand == '+' else observations[length-b:length-a][::-1]
                        roles = E.role_path(parts)
                        counts = observation_counts(bits, roles)
                        row.update(coverage_class=coverage_class(counts), span=counts,
                            role_groups={name:observation_counts(bits[keep], roles[keep])
                                         for name, indices in GROUPS.items()
                                         for keep in [np.isin(roles, indices)]})
                    all_rows.append(row)
                    handle.write(json.dumps(row, allow_nan=False)+'\n')
                del observations
    for species, stats in exposure.items():
        source = prep['by_species'][species]
        if (int(stats['fine_per_pass'].sum()) != source['fine_exposures']
                or int(stats['coarse_per_pass'].sum()) != source['coarse_exposures']
                or stats['unknown_per_pass'] != source['unknown_exposures']
                or stats['unique_windows'] != source['unique_windows']):
            raise ValueError('Per-pass original observation census differs')
        stats['fine_per_pass'] = stats['fine_per_pass'].tolist()
        stats['coarse_per_pass'] = stats['coarse_per_pass'].tolist()
    pooled = summarize(all_rows)
    if (len(all_rows) != len(references) or len({key(r) for r in all_rows}) != len(references)
            or pooled['categories'] != energy['pooled']['categories']):
        raise ValueError('Frozen primary reference or energy category denominator differs')
    report = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], fixture=checked,
        source_preparation_job=prep['job_id'], source_fit_job=fit['job_id'], source_score_job=energy['job_id'],
        scope='original_1491_unique_TRAIN_core_observations_joined_to_all_frozen_primary_references',
        original_passes=3, exposure=exposure, pooled=pooled,
        per_species={s:summarize([r for r in all_rows if r['species'] == s]) for s in F.Q.ALLOW},
        original_TRAIN_DNA_sources=descriptors, original_targets_reconstructed=True,
        original_targets_changed=False, primary_reference_changed=False,
        feature_tensors_read=False, checkpoint_read=False, native_forward=False,
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False, predictions_modified=False,
        matching_and_disagreeing_observation_counts_can_overlap=True,
        union_core_observation_coverage_is_not_partition_gradient_or_encoder_exposure=True,
        whole_span_matching_fine_is_not_joint_gene_set_or_independent_test=True,
        interpretation='descriptive_original_supervision_coverage_and_reference_agreement_not_new_fit_or_causal_attribution',
        seconds=time.perf_counter()-began)
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('pooled', 'per_species')}, allow_nan=False), flush=True)


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
