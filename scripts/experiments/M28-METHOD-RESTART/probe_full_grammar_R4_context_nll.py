#!/usr/bin/env python3
"""Same eight TRAIN witnesses: available real contexts, zero-potential NLL, no fit."""
from collections import Counter
import json
import math
import os
import resource
import time

import torch

import run_query as Q
import probe_full_grammar_R4_completion_witness as W
from src.m28.labels import build_targets
from src.m28.full_grammar_observed_partition import allowed_roles_from_local_targets, InfeasibleObservation
from src.m28.full_grammar_boundary_compiled import (
    _forward_messages, _reverse_messages, compiled_boundary_core_observation_nll)

ROOT = Q.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-R4-CONTEXT-NLL'
SOURCES = (('M28-FULL-GRAMMAR-R4-COMPLETION-WITNESS','13427817'),
           ('M28-FULL-GRAMMAR-R4-SUPPORT-PROBE','13427583'),
           ('M28-FULL-GRAMMAR-BOUNDARY-COMPILED-CHECKS','13428495'),
           ('M28-FULL-GRAMMAR-BOUNDARY-COMPILED-WINDOW-PROBE','13428516'))
CONFIG = dict(scope='all_original_eight_TRAIN_witnesses_not_a_new_training_subset',
    context='same_nine_previously_opened_original_TRAIN_DNA_cache_artifacts',
    bounds='same_annotation_derived_completion_witness_bounds_diagnostic_only',
    missing_context='NOT_EXECUTED_external_DNA_missing_never_fabricated',
    targets='unchanged_R4_build_targets_fine_else_coarse_else_unknown_all_core_rows',
    potentials='all_zero_CPU_FP32_all15_roles_no_model',role_weights=1,null=0,
    max_exons=128,checkpoint_bp=256,score_accumulation='FP64',
    loss='unweighted_unnormalized_global_logZ_all_minus_core_observed_logZ_diagnostic',
    closed_observed_logZ_comparison_rtol=0,closed_observed_logZ_comparison_atol=1e-8,
    NLL_and_zero_potential_logZ_nonnegative_atol=1e-8,
    NLL_common_row_shift_gradient_atol=1e-5,NLL_occupancy_difference_bound_atol=1e-6,
    optimizer_steps=0)


def ready():
    values = []
    for name,job in SOURCES:
        root = ROOT/'outputs'/name
        summary = json.loads((root/'summary.json').read_text())
        if ((root/'JOBID').read_text().strip() != job or (root/'STATUS').read_text().strip() != 'COMPLETED'
                or summary['job_id'] != job or not summary['completed'] or summary['optimizer_steps_this_run'] != 0):
            raise ValueError('Matching completed zero-update source required: '+name)
        values.append(summary)
    witness,support,checks,full = values
    if (witness['actual_witnesses'] != 8 or witness['DNA_verified_completion_witnesses'] != 5
            or witness['original_draws'] != 1536 or witness['original_train_unique_windows'] != 1491
            or witness['original_R4_policy_changed'] or witness['observation_rows_deleted'] or witness['DEV_used']
            or support['actual_witnesses'] != 8 or support['original_R4_policy_changed'] or support['DEV_used']
            or checks['tests_run'] != 10 or checks['failures'] or checks['errors']
            or not checks['core_NLL_and_outside_gradients_verified']
            or not checks['every_supported_output_and_incoming_message_gradient_verified']
            or not checks['checkpoint_recomputation_invariance_verified']
            or full['source_mathematical_checks_job'] != '13428495'
            or not full['full_window_forward_and_backward_completed'] or not full['all_input_gradients_finite']
            or full['length'] != 24576 or full['max_exons'] != 128):
        raise ValueError('Actual fixed completion, original observation, mathematical and full-resource evidence differs')
    ids = [r['window_id'] for r in witness['records']]
    if (len(ids) != 8 or len(set(ids)) != 8 or [r['ordinal'] for r in witness['records']] != list(range(8))
            or [r['window_id'] for r in support['records']] != ids):
        raise ValueError('All eight original ordered witness identities required')
    return witness,support


def main():
    began = time.perf_counter()
    witness,support = ready()
    rows = witness['records']; ids = [r['window_id'] for r in rows]
    with (OUT/'started.json').open('x') as handle:
        json.dump(dict(job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
            source_jobs=[job for _,job in SOURCES],witness_ids=ids,
            permitted_cached_artifacts=witness['cached_artifacts_opened'],
            TRAIN_allowlists={s:sorted(v) for s,v in Q.ALLOW.items()},optimizer_steps_this_run=0,
            forward_signatures_before=len(_forward_messages.nopython_signatures),
            reverse_signatures_before=len(_reverse_messages.nopython_signatures)),handle,indent=2,allow_nan=False)
    metadata_started = time.perf_counter()
    draws,index,windows,catalogs = Q.load_train_inputs()
    metadata_seconds = time.perf_counter()-metadata_started
    if len(draws) != 1536 or len(windows) != 1491:
        raise ValueError('Original TRAIN population changed')
    permitted = set(witness['cached_artifacts_opened'])
    permitted_windows = {i:w for i,w in windows.items() if index[i]['artifact'] in permitted}
    if len(permitted) != 9 or len(permitted_windows) != 9 or not set(ids) <= set(permitted_windows):
        raise ValueError('Exactly the same nine previously opened TRAIN caches required')
    opened,records = set(),[]
    with (OUT/'windows.jsonl').open('x') as handle:
        for original,closed in zip(rows,support['records']):
            row_started = time.perf_counter()
            ident = original['window_id']; w = windows[ident]; n = w['valid_bases']
            cached = W.load_cache(w,index,opened); core = cached['sequence']
            targets = build_targets(w,catalogs[w['species']],core)
            allowed = allowed_roles_from_local_targets(targets,n)
            fine = int((torch.as_tensor(targets['seg'][:n]) >= 0).sum())
            free = int(allowed.all(1).sum()); coarse = n-fine-free
            if (fine,coarse,free) != (original['fine_positions'],original['coarse_positions'],original['fully_unknown_positions']):
                raise ValueError('Original core fine/coarse/unknown counts changed: '+ident)
            primary = [catalogs[w['species']][i] for i in w['overlapping_transcript_ids']
                       if catalogs[w['species']][i]['training_primary']]
            if [tx['id'] for tx in primary] != original['original_primary_transcript_ids']:
                raise ValueError('Original complete primary set changed: '+ident)
            a = min([w['start']]+[tx['CDS'][0][0] for tx in primary])
            b = max([w['end']]+[tx['CDS'][-1][1] for tx in primary])
            offset = w['start']-a if w['strand'] == '+' else b-w['end']
            if [a,b] != original['diagnostic_genomic_context'] or offset != original['core_oriented_offset']:
                raise ValueError('Original diagnostic-only context geometry changed: '+ident)
            context,fragments,missing = W.cached_context(w,core,a,b,permitted_windows,index,opened)
            if fragments != original['cached_fragments'] or missing != original['missing_genomic_DNA_intervals']:
                raise ValueError('Actual original cache context or missing-DNA evidence changed: '+ident)
            row = dict(ordinal=original['ordinal'],window_id=ident,species=w['species'],category=original['category'],
                seqid=w['seqid'],strand=w['strand'],valid_core_bases=n,diagnostic_genomic_context=[a,b],
                context_bounds_annotation_derived_not_a_training_or_raw_inference_rule=True,
                core_oriented_interval=[offset,offset+n],original_primary_transcript_ids=original['original_primary_transcript_ids'],
                fine_positions=fine,coarse_positions=coarse,fully_unknown_positions=free,
                cached_fragments=fragments,missing_genomic_DNA_intervals=missing,
                source_completion_witness_status=original['completion_witness_status'],
                original_closed_window_observation_feasible=closed['observation_feasible'],
                original_closed_observed_logZ=closed['log_number_compatible_sets'],
                NLL_computed=False,observation_support_evaluated=False,gradients_computed=False)
            if context is None:
                if original['DNA_verified_completion_witness'] or not missing:
                    raise ValueError('Missing context contradicts completed source: '+ident)
                row.update(status='NOT_EXECUTED_EXTERNAL_DNA_MISSING',context_observation_feasible=None,
                    row_seconds=time.perf_counter()-row_started)
            else:
                if not original['DNA_verified_completion_witness'] or len(context) != b-a:
                    raise ValueError('Available context does not match actual DNA completion evidence: '+ident)
                cuts = tuple(x for x in (offset,offset+n) if 0 < x < len(context))
                logits = torch.zeros((len(context),15),dtype=torch.float32,requires_grad=True)
                print(json.dumps(dict(stage='actual_context_NLL_started',window_id=ident,context_length=len(context),cuts=cuts)),flush=True)
                forward_started = time.perf_counter()
                try:
                    result = compiled_boundary_core_observation_nll(context,logits,128,offset,offset+n,allowed,
                        cuts=cuts,checkpoint_bp=256)
                except InfeasibleObservation as error:
                    raise RuntimeError('Observed support contradicts the original actual DNA completion witness: '+ident) from error
                forward_seconds = time.perf_counter()-forward_started
                all_mass,obs_mass,loss = (float(result[k].detach()) for k in ('log_partition_all','log_partition_observed','loss'))
                if not all(math.isfinite(x) and x >= -1e-8 for x in (all_mass,obs_mass,loss)):
                    raise ValueError('Finite nonnegative zero-potential log counts and NLL required: '+ident)
                closed_match = None
                if len(context) == n and offset == 0:
                    closed_match = (closed['observation_feasible'] and math.isclose(obs_mass,
                        closed['log_number_compatible_sets'],rel_tol=0,abs_tol=1e-8))
                    if not closed_match: raise ValueError('Same closed-context observed logZ differs from the actual original result: '+ident)
                backward_started = time.perf_counter(); result['loss'].backward()
                backward_seconds = time.perf_counter()-backward_started
                gradient = logits.grad
                if (gradient is None or gradient.dtype != torch.float32 or tuple(gradient.shape) != (len(context),15)
                        or not bool(torch.isfinite(gradient).all())):
                    raise ValueError('Every actual-context input role must receive finite NLL derivatives: '+ident)
                row_shift = float(gradient.sum(1).abs().max())
                role_difference = float(gradient[:,1:].abs().max())
                covered_difference = float(gradient[:,1:].sum(1).abs().max())
                if row_shift > 1e-5 or role_difference > 1+1e-6 or covered_difference > 1+1e-6:
                    raise ValueError('NLL differences of role posteriors or row-shift derivatives exceed original bounds: '+ident)
                prefix = float(gradient[:offset].double().abs().sum())
                suffix = float(gradient[offset+n:].double().abs().sum())
                row.update(status='COMPLETED_AVAILABLE_CONTEXT_NLL',context_observation_feasible=True,
                    context_length=len(context),interior_message_cuts=list(cuts),log_partition_all=all_mass,
                    log_partition_core_observed=obs_mass,unnormalized_unweighted_NLL=loss,
                    same_closed_context_observed_logZ_matches_original=closed_match,
                    NLL_forward_seconds=forward_seconds,NLL_backward_seconds=backward_seconds,
                    all_context_NLL_gradients_finite=True,NLL_common_row_shift_gradient_max_abs=row_shift,
                    maximum_abs_nonBG_role_NLL_occupancy_difference=role_difference,
                    maximum_abs_nonBG_covered_NLL_probability_difference=covered_difference,
                    core_NLL_gradient_l1=float(gradient[offset:offset+n].double().abs().sum()),
                    outside_prefix_NLL_gradient_l1=prefix,outside_suffix_NLL_gradient_l1=suffix,
                    outside_NLL_gradient_nonzero=prefix+suffix > 0,
                    statistics_all_after_backward=[dict(v) for v in result['statistics_all']],
                    statistics_observed_after_backward=[dict(v) for v in result['statistics_observed']],
                    process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                    NLL_computed=True,observation_support_evaluated=True,gradients_computed=True,
                    row_seconds=time.perf_counter()-row_started)
                del result,logits,gradient
            records.append(row); handle.write(json.dumps(row,allow_nan=False)+'\n'); handle.flush()
            print(json.dumps({k:row[k] for k in ('ordinal','window_id','status','context_observation_feasible','NLL_computed','gradients_computed','row_seconds')}),flush=True)
            del cached,targets,allowed
    counts = Counter(r['status'] for r in records)
    if ([r['window_id'] for r in records] != ids or counts != {'COMPLETED_AVAILABLE_CONTEXT_NLL':5,'NOT_EXECUTED_EXTERNAL_DNA_MISSING':3}
            or opened != permitted):
        raise ValueError('All eight original witnesses and exactly the original nine cache artifacts must be accounted for')
    summary = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
        protocol='reports/M28-METHOD-RESTART/full_grammar_partition_preparation.md',source_jobs=[job for _,job in SOURCES],
        original_draws=1536,original_train_unique_windows=1491,actual_witnesses=8,all_original_witnesses_accounted=True,
        available_context_NLL_and_gradient_windows=5,missing_context_NLL_NOT_EXECUTED_windows=3,
        all_original_context_NLLs_executed=False,records=records,cached_artifacts_opened=sorted(opened),
        original_TRAIN_metadata_load_seconds=metadata_seconds,seconds_after_module_import=time.perf_counter()-began,
        process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        forward_signatures_after=len(_forward_messages.nopython_signatures),reverse_signatures_after=len(_reverse_messages.nopython_signatures),
        synthetic_only=False,original_TRAIN_DNA_and_catalog_read=True,original_R4_targets_regenerated=True,
        cached_features_used=False,old_cached_labels_used=False,research_potentials_or_weights_read=False,
        boundary_prior_chosen=False,training_context_extent_or_padding_rule_selected=False,
        loss_normalizer_importance_rule_or_fit_coefficients_selected=False,optimizer_steps_this_run=0,
        model_constructed=False,optimizer_constructed=False,model_or_backbone_forward=False,checkpoints_written=False,
        original_R4_policy_changed=False,observation_rows_deleted=False,primary_transcripts_dropped=False,
        external_DNA_fabricated=False,unindexed_FASTA_scanned=False,raw_FASTA_opened=False,
        DEV_used=False,independent_species_used=False,A33_outputs_opened=False,training_or_inference_caller_changed=False,
        interpretation='fixed_available_real_TRAIN_context_zero_potential_NLL_and_gradients_not_full_cohort_support_fit_padding_contract_or_performance')
    with (OUT/'summary.json').open('x') as handle: json.dump(summary,handle,indent=2,allow_nan=False)
    print(json.dumps({k:v for k,v in summary.items() if k != 'records'},allow_nan=False),flush=True)


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
