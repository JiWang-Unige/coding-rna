#!/usr/bin/env python3
"""Eight original TRAIN witnesses: zero-potential observation support, no model."""
from collections import Counter
import json
import math
import os
import time

import torch

import run_query as Q
from src.m28.labels import build_targets, oriented_parts
from src.m28.full_grammar_observed_partition import (
    InfeasibleObservation, allowed_roles_from_local_targets,
    observed_full_grammar_log_partition)

ROOT = Q.ROOT
OUT = ROOT / 'outputs/M28-FULL-GRAMMAR-R4-SUPPORT-PROBE'
SOURCE = ROOT / 'outputs/M28-A44-LOCAL-STRUCTURED-PREFLIGHT-R1'
CHECKS = ROOT / 'outputs/M28-FULL-GRAMMAR-OBSERVED-CHECKS'
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_partition_preparation.md'
SOURCE_JOB, CHECKS_JOB = '13425772', '13427476'
CONFIG = dict(scope='original_eight_fixed_TRAIN_witnesses_not_a_new_training_subset',
    original_train_unique_windows=1491, original_draws=1536, max_exons=128,
    potentials='all_zero_CPU_FP32_15_roles_no_model', score_accumulation='FP64',
    targets='unchanged_original_R4_build_targets_on_oriented_cached_DNA',
    observation='fine_singleton_else_original_coarse_group_else_all_roles',
    boundary='original_complete_gene_grammar_no_repair_or_observation_deletion',
    result='existence_and_log_number_of_compatible_complete_gene_sets_not_NLL_fit',
    optimizer_steps=0)


def ready():
    for root, ident in ((SOURCE, SOURCE_JOB), (CHECKS, CHECKS_JOB)):
        if ((root / 'STATUS').read_text().strip() != 'COMPLETED'
                or (root / 'JOBID').read_text().strip() != ident):
            raise ValueError('Matching closed source witness and observation mathematics required')
    source = json.loads((SOURCE / 'summary.json').read_text())
    checks = json.loads((CHECKS / 'summary.json').read_text())
    if (not source['completed'] or source['windows'] != 1491 or source['original_draws'] != 1536
            or len(source['records']) != 8 or source['optimizer_steps_this_run'] != 0
            or not source['actual_original_R4_and_complete_support_geometry_cache_gold_reconciled']
            or not source['actual_local_unknown_callable_padding_masks_and_position_denominators_verified']
            or source['new_R4_policy_or_chain_labels'] or source['DEV_used']):
        raise ValueError('Original eight fixed, unchanged R4 TRAIN witness source required')
    if (not checks['completed'] or checks['job_id'] != CHECKS_JOB or checks['tests_run'] != 10
            or checks['failures'] or checks['errors'] or checks['optimizer_steps_this_run'] != 0
            or not checks['artificial_observed_mass_and_NLL_full_gradients_verified']
            or not checks['fine_coarse_unknown_padding_mapping_verified_on_artificial_targets']):
        raise ValueError('Actual closed observed partition mathematics required')
    return source['records']


def main():
    chosen = ready()
    # Reuse the original loader and its pre-parse TRAIN catalog filter. No DEV
    # arrays, labels, sequence or fitted weights are opened by this diagnostic.
    draws, index, windows, catalogs = Q.load_train_inputs()
    if len(draws) != 1536 or len(windows) != 1491:
        raise ValueError('Complete unchanged original TRAIN input population required')
    ids = [row['window_id'] for row in chosen]
    if (len(set(ids)) != 8 or not set(ids) <= set(windows)
            or Counter(row['species'] for row in chosen) != {s: 4 for s in Q.ALLOW}):
        raise ValueError('Original two-species four-category witness identity differs')
    began = time.perf_counter()
    records = []
    with (OUT / 'started.json').open('x') as handle:
        json.dump(dict(job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
            source_preflight_job=SOURCE_JOB, source_mathematical_checks_job=CHECKS_JOB,
            witness_ids=ids, original_TRAIN_allowlists={s: sorted(v) for s, v in Q.ALLOW.items()},
            optimizer_steps_this_run=0), handle, indent=2, allow_nan=False)
    with torch.inference_mode(), (OUT / 'windows.jsonl').open('x') as handle:
        for ordinal, witness in enumerate(chosen):
            started = time.perf_counter()
            ident = witness['window_id']; w = windows[ident]; n = w['valid_bases']
            if w['species'] != witness['species'] or w['split'] != 'train':
                raise ValueError('Original witness/window TRAIN species differs')
            cached = torch.load(ROOT / index[ident]['artifact'], map_location='cpu', weights_only=True)
            Q.validate_cached(cached, w)
            targets = build_targets(w, catalogs[w['species']], cached['sequence'])
            fine = torch.as_tensor(targets['seg']); region = torch.as_tensor(targets['region'])
            fine_count = int((fine[:n] >= 0).sum())
            coarse_count = int(((fine[:n] < 0) & (region[:n] >= 0) & (region[:n] < 3)).sum())
            if (fine_count != witness['fine_positions']
                    or coarse_count != witness['grouped_region_positions']):
                raise ValueError('Original witness target counts changed since source preflight')
            allowed = allowed_roles_from_local_targets(targets, n)
            free_rows = int(allowed.all(dim=1).sum())
            if n - free_rows != fine_count + coarse_count:
                raise ValueError('Actual fine/coarse/unknown observation row mapping differs')
            primary = [catalogs[w['species']][tid] for tid in w['overlapping_transcript_ids']
                       if catalogs[w['species']][tid]['training_primary']]
            primary_geometry = [(tx, oriented_parts(tx, w)) for tx in primary]
            boundary_or_incomplete = [dict(transcript_id=tx['id'], gene_id=tx['gene_id'],
                oriented_start=parts[0][0], oriented_end=parts[-1][1], exons=len(parts),
                crosses_left=parts[0][0] < 0, crosses_right=parts[-1][1] > n,
                CDS_complete_supported=tx['CDS_complete_supported'])
                for tx, parts in primary_geometry
                if parts[0][0] < 0 or parts[-1][1] > n or not tx['CDS_complete_supported']]
            # Necessary complete-set edge roles are background or start at the
            # first base, and background or stop at the last. Either obstruction
            # is sufficient for infeasibility, not an exhaustive cause taxonomy.
            left_obstruction = not bool(allowed[0, 0] or allowed[0, 13])
            right_obstruction = not bool(allowed[-1, 0] or allowed[-1, 14])
            potentials = torch.zeros((n, 15), dtype=torch.float32)
            forward_started = time.perf_counter()
            try:
                result = observed_full_grammar_log_partition(cached['sequence'], potentials, 128, allowed)
            except InfeasibleObservation:
                feasible, mass, counts = False, None, None
            else:
                feasible, mass, counts = True, float(result['log_partition']), result['counts']
                if not math.isfinite(mass) or mass < 0:
                    raise ValueError('Nonfinite or negative zero-potential observation log count')
                if left_obstruction or right_obstruction:
                    raise ValueError('Complete-set support contradicts a necessary boundary role')
                del result
            row = dict(ordinal=ordinal, window_id=ident, species=w['species'], category=witness['category'],
                original_witness_priority=witness['priority'], seqid=w['seqid'], strand=w['strand'],
                start=w['start'], end=w['end'], valid_bases=n, cache_artifact=index[ident]['artifact'],
                fine_positions=fine_count, coarse_positions=coarse_count, fully_unknown_positions=free_rows,
                callable_positions=int(targets['callable_bases'][:n].sum()),
                original_positive_chains=len(w['positive_ids']), original_overlapping_primaries=len(primary),
                original_chain_unknown_intervals=w['chain_unknown_intervals'],
                primary_boundary_or_incomplete_transcripts=boundary_or_incomplete,
                allowed_roles_left=torch.where(allowed[0])[0].tolist(),
                allowed_roles_right=torch.where(allowed[-1])[0].tolist(),
                left_boundary_role_obstruction=left_obstruction, right_boundary_role_obstruction=right_obstruction,
                observation_feasible=feasible, log_number_compatible_sets=mass, observed_counts=counts,
                observed_forward_seconds=time.perf_counter() - forward_started,
                load_targets_forward_seconds=time.perf_counter() - started)
            handle.write(json.dumps(row, allow_nan=False) + '\n'); handle.flush()
            records.append(row)
            print(json.dumps({k: row[k] for k in ('ordinal', 'window_id', 'observation_feasible',
                'left_boundary_role_obstruction', 'right_boundary_role_obstruction',
                'load_targets_forward_seconds')}, allow_nan=False), flush=True)
            del cached, targets, fine, region, allowed, potentials
    if len(records) != 8 or [row['window_id'] for row in records] != ids:
        raise ValueError('Incomplete or reordered original witness observation support diagnostic')
    species_counts = {}
    for species in Q.ALLOW:
        rows = [r for r in records if r['species'] == species]
        species_counts[species] = dict(windows=len(rows), feasible=sum(r['observation_feasible'] for r in rows),
            infeasible=sum(not r['observation_feasible'] for r in rows),
            left_boundary_role_obstructions=sum(r['left_boundary_role_obstruction'] for r in rows),
            right_boundary_role_obstructions=sum(r['right_boundary_role_obstruction'] for r in rows))
    summary = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL,
        configuration=CONFIG, source_preflight_job=SOURCE_JOB, source_mathematical_checks_job=CHECKS_JOB,
        original_train_unique_windows=1491, original_draws=1536, actual_witnesses=8,
        all_original_witnesses_checked=True, all_original_observations_feasible=all(r['observation_feasible'] for r in records),
        species_counts=species_counts, records=records, seconds=time.perf_counter() - began,
        original_R4_targets_regenerated=True, original_TRAIN_catalog_metadata_read=True,
        original_cached_GEN_and_DNA_artifacts_opened=8, cached_features_used=False,
        research_grammar_observation_support_evaluated=True, candidate_pool_or_prediction_written=False,
        optimizer_steps_this_run=0, optimizer_constructed=False, model_constructed=False,
        model_or_backbone_forward=False, research_weights_read=False, checkpoints_written=False,
        original_R4_policy_changed=False, observation_rows_deleted_or_boundary_repaired=False,
        old_cached_labels_used=False, loss_normalizer_or_fit_coefficients_selected=False,
        DEV_used=False, independent_species_used=False, A33_outputs_opened=False,
        training_or_inference_caller_changed=False,
        interpretation='eight_fixed_TRAIN_zero_potential_support_diagnostic_not_cohort_rate_fit_NLL_or_model_performance')
    with (OUT / 'summary.json').open('x') as handle:
        json.dump(summary, handle, indent=2, allow_nan=False)


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
