#!/usr/bin/env python3
"""All eight original TRAIN diagnostic contexts from original FASTA; no fit."""
import json
import math
import os
from pathlib import Path
import resource
import time

import torch

import probe_full_grammar_R4_context_nll as P
from src.m28.labels import oriented_parts

Q, W, ROOT = P.Q, P.W, P.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-R4-TRAIN-FASTA-CONTEXT'
PARENT = ROOT/'outputs/M28-FULL-GRAMMAR-R4-CONTEXT-NLL'
PARENT_JOB = '13428648'
CONFIG = dict(P.CONFIG,
    context='original_dataset_FASTA_only_exact_six_original_TRAIN_sequence_records',
    missing_context='STOP_if_original_FASTA_or_matching_core_DNA_unavailable',
    reader='unchanged_manifest.read_allowed_no_census.metadata_no_FASTA_index_creation',
    source_scope='original_TWO_species_TRAIN_only_nonTRAIN_FASTA_lines_discarded_not_parsed_or_retained',
    parent_zero_potential_logZ_and_NLL_comparison_rtol=0,
    parent_zero_potential_logZ_and_NLL_comparison_atol=1e-8)


def read_train_DNA():
    config = Q.M.C.yaml.safe_load(Q.M.C.CONFIG.read_text())
    # Select species before any data-path access. Do not call census.metadata:
    # that older helper reads split_val.jsonl.
    paths = {Path(rel).name:ROOT/rel for rel in config['data']['development_species']
             if Path(rel).name in Q.ALLOW}
    if set(paths) != set(Q.ALLOW):
        raise ValueError('Original two-species TRAIN dataset paths required')
    sequences, descriptors = {}, {}
    for species, path in paths.items():
        lengths = {}
        with (path/'reference.gff3').open() as handle:
            for line in handle:
                if not line.startswith('##sequence-region'):
                    continue
                fields = line.split()
                if fields[1] not in Q.ALLOW[species]:
                    continue
                if int(fields[2]) != 1:
                    raise ValueError('Original TRAIN sequence-region must start at one')
                lengths[fields[1]] = int(fields[3])
                if set(lengths) == set(Q.ALLOW[species]):
                    break
        if set(lengths) != set(Q.ALLOW[species]):
            raise ValueError('All three original TRAIN sequence lengths required: '+species)
        sequences[species] = Q.M.read_allowed(path/'genome.fa',lengths)
        descriptors[species] = dict(configured_directory=str(path.relative_to(ROOT)),
            resolved_FASTA=str((path/'genome.fa').resolve()),FASTA_bytes=(path/'genome.fa').stat().st_size,
            TRAIN_sequence_lengths=lengths,parsed_sequence_records=sorted(sequences[species]),
            GFF_sequence_region_header_metadata_only=True,annotation_features_parsed=False,
            reader='manifest.read_allowed',nonTRAIN_FASTA_lines_streamed_but_discarded=True,
            nonTRAIN_DNA_parsed_retained_or_used=False)
    return sequences,descriptors


def main():
    began = time.perf_counter()
    witness,support = P.ready()
    parent = json.loads((PARENT/'summary.json').read_text())
    rows = witness['records']; ids = [r['window_id'] for r in rows]
    if ((PARENT/'JOBID').read_text().strip() != PARENT_JOB
            or (PARENT/'STATUS').read_text().strip() != 'COMPLETED'
            or parent['job_id'] != PARENT_JOB or not parent['completed']
            or parent['configuration'] != P.CONFIG or parent['optimizer_steps_this_run'] != 0
            or parent['available_context_NLL_and_gradient_windows'] != 5
            or parent['missing_context_NLL_NOT_EXECUTED_windows'] != 3
            or [r['window_id'] for r in parent['records']] != ids):
        raise ValueError('Original completed five-available/three-unexecuted parent required')
    with (OUT/'started.json').open('x') as handle:
        json.dump(dict(job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
            source_jobs=[job for _,job in P.SOURCES]+[PARENT_JOB],witness_ids=ids,
            permitted_cached_artifacts=witness['cached_artifacts_opened'],
            TRAIN_allowlists={s:sorted(v) for s,v in Q.ALLOW.items()},optimizer_steps_this_run=0,
            forward_signatures_before=len(P._forward_messages.nopython_signatures),
            reverse_signatures_before=len(P._reverse_messages.nopython_signatures)),handle,indent=2,allow_nan=False)
    t = time.perf_counter()
    draws,index,windows,catalogs = Q.load_train_inputs()
    metadata_seconds = time.perf_counter()-t
    if len(draws) != 1536 or len(windows) != 1491:
        raise ValueError('Original TRAIN population changed')
    permitted = set(witness['cached_artifacts_opened'])
    permitted_windows = {i:w for i,w in windows.items() if index[i]['artifact'] in permitted}
    if len(permitted) != 9 or len(permitted_windows) != 9 or not set(ids) <= set(permitted_windows):
        raise ValueError('Exactly the original nine TRAIN caches required')
    t = time.perf_counter()
    sequences,descriptors = read_train_DNA()
    raw_read_seconds = time.perf_counter()-t
    opened,records,new_chains = set(),[],[]
    with (OUT/'windows.jsonl').open('x') as handle, (OUT/'contexts.jsonl').open('x') as context_handle:
        for original,closed,previous in zip(rows,support['records'],parent['records']):
            row_started = time.perf_counter()
            ident = original['window_id']; w = windows[ident]; n = w['valid_bases']
            cached = W.load_cache(w,index,opened); core = cached['sequence']
            targets = P.build_targets(w,catalogs[w['species']],core)
            allowed = P.allowed_roles_from_local_targets(targets,n)
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
            chromosome = sequences[w['species']][w['seqid']]
            if not 0 <= a < b <= len(chromosome):
                raise ValueError('Diagnostic context is outside the original TRAIN chromosome')
            genomic = chromosome[a:b]
            context = genomic if w['strand'] == '+' else Q.M.rc(genomic)
            if context[offset:offset+n] != core:
                raise ValueError('Raw original TRAIN FASTA disagrees with an original core DNA base: '+ident)
            old_context,fragments,missing = W.cached_context(w,core,a,b,permitted_windows,index,opened)
            if fragments != original['cached_fragments'] or missing != original['missing_genomic_DNA_intervals']:
                raise ValueError('Original cache fragments or missing-DNA evidence changed: '+ident)
            if old_context is not None and old_context != context:
                raise ValueError('Raw original TRAIN context disagrees with the complete cached context: '+ident)
            if (old_context is not None) != original['DNA_verified_completion_witness']:
                raise ValueError('Original available-context identity changed: '+ident)
            verified = []
            for tx,old_chain in zip(primary,original['primary_chains']):
                if tx['id'] != old_chain['transcript_id']:
                    raise ValueError('Original primary chain order changed')
                if old_chain['DNA_source'] != 'external_DNA_unavailable':
                    if old_chain['actual_full_grammar_legal'] is not True:
                        raise ValueError('Original cached primary-chain evidence is not legal')
                    continue
                parts = oriented_parts(tx,dict(start=a,end=b,strand=w['strand']))
                if len(parts) != old_chain['exons'] or [parts[0][0],parts[-1][1]] != old_chain['oriented_span']:
                    raise ValueError('Original unverified primary geometry changed')
                geometry = W.geometry_failures(parts,len(context))
                failures,actual = W.DNA_checks(context,parts)
                check = dict(window_id=ident,transcript_id=tx['id'],gene_id=tx['gene_id'],
                    oriented_parts=[list(p) for p in parts],geometry_failures=geometry,
                    direct_DNA_failures=failures,actual_carry_coding_verdict=actual,
                    DNA_source='original_TRAIN_FASTA_diagnostic_context',actual_full_grammar_legal=not geometry and not failures)
                print(json.dumps(dict(stage='previously_missing_primary_DNA_checked',**check)),flush=True)
                if geometry or failures or not actual:
                    raise ValueError('Actual original TRAIN boundary chain is rejected: '+tx['id'])
                verified.append(check); new_chains.append(check)
            context_handle.write(json.dumps(dict(window_id=ident,species=w['species'],seqid=w['seqid'],
                strand=w['strand'],diagnostic_genomic_context=[a,b],core_oriented_interval=[offset,offset+n],
                annotation_derived_bounds_diagnostic_only=True,sequence=context),allow_nan=False)+'\n')
            context_handle.flush()
            row = dict(ordinal=original['ordinal'],window_id=ident,species=w['species'],category=original['category'],
                seqid=w['seqid'],strand=w['strand'],valid_core_bases=n,diagnostic_genomic_context=[a,b],
                context_bounds_annotation_derived_not_a_training_or_raw_inference_rule=True,
                core_oriented_interval=[offset,offset+n],original_primary_transcript_ids=original['original_primary_transcript_ids'],
                fine_positions=fine,coarse_positions=coarse,fully_unknown_positions=free,
                original_cached_fragments=fragments,original_cache_missing_genomic_DNA_intervals=missing,
                source_completion_witness_status=original['completion_witness_status'],
                source_parent_context_NLL_status=previous['status'],raw_FASTA_core_matches_original_cache=True,
                raw_FASTA_full_context_matches_original_cached_context=True if old_context is not None else None,
                newly_verified_primary_chains=verified,original_closed_window_observation_feasible=closed['observation_feasible'],
                original_closed_observed_logZ=closed['log_number_compatible_sets'])
            cuts = tuple(x for x in (offset,offset+n) if 0 < x < len(context))
            logits = torch.zeros((len(context),15),dtype=torch.float32,requires_grad=True)
            print(json.dumps(dict(stage='original_FASTA_context_NLL_started',window_id=ident,context_length=len(context),cuts=cuts)),flush=True)
            t = time.perf_counter()
            try:
                result = P.compiled_boundary_core_observation_nll(context,logits,128,offset,offset+n,allowed,
                    cuts=cuts,checkpoint_bp=256)
            except P.InfeasibleObservation as error:
                raise RuntimeError('Core support contradicts the actual original TRAIN DNA completion: '+ident) from error
            forward_seconds = time.perf_counter()-t
            all_mass,obs_mass,loss = (float(result[k].detach()) for k in ('log_partition_all','log_partition_observed','loss'))
            if not all(math.isfinite(x) and x >= -1e-8 for x in (all_mass,obs_mass,loss)):
                raise ValueError('Finite nonnegative zero-potential log counts and NLL required: '+ident)
            parent_differences = None
            if previous['NLL_computed']:
                parent_differences = dict(zip(('log_partition_all','log_partition_core_observed','unnormalized_unweighted_NLL'),
                    (value-previous[key] for value,key in zip((all_mass,obs_mass,loss),
                    ('log_partition_all','log_partition_core_observed','unnormalized_unweighted_NLL')))))
                if any(abs(value) > 1e-8 for value in parent_differences.values()):
                    raise ValueError('Original available-context zero-potential control differs from parent: '+ident)
            closed_match = None
            if len(context) == n and offset == 0:
                closed_match = closed['observation_feasible'] and math.isclose(obs_mass,
                    closed['log_number_compatible_sets'],rel_tol=0,abs_tol=1e-8)
                if not closed_match:
                    raise ValueError('Same closed-context observed logZ differs from original result: '+ident)
            t = time.perf_counter(); result['loss'].backward()
            backward_seconds = time.perf_counter()-t
            gradient = logits.grad
            if (gradient is None or gradient.dtype != torch.float32 or tuple(gradient.shape) != (len(context),15)
                    or not bool(torch.isfinite(gradient).all())):
                raise ValueError('Every actual-context input role must receive finite NLL derivatives: '+ident)
            row_shift = float(gradient.sum(1).abs().max())
            role_difference = float(gradient[:,1:].abs().max())
            covered_difference = float(gradient[:,1:].sum(1).abs().max())
            if row_shift > 1e-5 or role_difference > 1+1e-6 or covered_difference > 1+1e-6:
                raise ValueError('NLL posterior differences or row-shift derivatives exceed original bounds: '+ident)
            prefix = float(gradient[:offset].double().abs().sum())
            suffix = float(gradient[offset+n:].double().abs().sum())
            row.update(status='COMPLETED_ORIGINAL_TRAIN_FASTA_CONTEXT_NLL',context_observation_feasible=True,
                context_length=len(context),interior_message_cuts=list(cuts),log_partition_all=all_mass,
                log_partition_core_observed=obs_mass,unnormalized_unweighted_NLL=loss,
                same_closed_context_observed_logZ_matches_original=closed_match,
                parent_available_context_logZ_and_NLL_differences=parent_differences,
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
            records.append(row); handle.write(json.dumps(row,allow_nan=False)+'\n'); handle.flush()
            print(json.dumps({k:row[k] for k in ('ordinal','window_id','status','context_observation_feasible','NLL_computed','gradients_computed','row_seconds')}),flush=True)
            del result,logits,gradient,cached,targets,allowed
    if len(records) != 8 or [r['window_id'] for r in records] != ids or opened != permitted or len(new_chains) != 3:
        raise ValueError('All original eight contexts, nine caches and three new DNA chains required')
    summary = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,
        protocol='reports/M28-METHOD-RESTART/full_grammar_partition_preparation.md',
        source_jobs=[job for _,job in P.SOURCES]+[PARENT_JOB],
        original_draws=1536,original_train_unique_windows=1491,actual_witnesses=8,all_original_witnesses_accounted=True,
        original_TRAIN_DNA_sources=descriptors,records=records,cached_artifacts_opened=sorted(opened),
        contexts_written=8,all_original_context_NLLs_executed=True,context_NLL_and_gradient_windows=8,
        newly_DNA_verified_primary_chains=new_chains,original_35_primary_DNA_verdicts_transferred_by_exact_DNA_comparison=True,
        original_TRAIN_metadata_load_seconds=metadata_seconds,original_TRAIN_FASTA_read_seconds=raw_read_seconds,
        seconds_after_module_import=time.perf_counter()-began,process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        forward_signatures_after=len(P._forward_messages.nopython_signatures),reverse_signatures_after=len(P._reverse_messages.nopython_signatures),
        synthetic_only=False,original_TRAIN_DNA_and_catalog_read=True,original_R4_targets_regenerated=True,
        cached_features_used=False,old_cached_labels_used=False,research_potentials_or_weights_read=False,
        boundary_prior_chosen=False,training_context_extent_or_padding_rule_selected=False,
        loss_normalizer_importance_rule_or_fit_coefficients_selected=False,optimizer_steps_this_run=0,
        model_constructed=False,optimizer_constructed=False,model_or_backbone_forward=False,checkpoints_written=False,
        original_R4_policy_changed=False,observation_rows_deleted=False,primary_transcripts_dropped=False,
        external_DNA_fabricated=False,unindexed_FASTA_streamed_with_TRAIN_record_gate=True,raw_FASTA_opened=True,
        nonTRAIN_FASTA_lines_streamed_but_discarded=True,nonTRAIN_DNA_parsed_retained_or_used=False,
        split_val_or_test_files_opened=False,census_metadata_called=False,FASTA_index_created=False,
        DEV_used=False,independent_species_used=False,A33_outputs_opened=False,training_or_inference_caller_changed=False,
        interpretation='all_eight_original_real_TRAIN_diagnostic_context_zero_potential_NLL_and_gradients_not_full_cohort_fit_padding_contract_or_performance')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(summary,handle,indent=2,allow_nan=False)
    print(json.dumps({k:v for k,v in summary.items() if k != 'records'},allow_nan=False),flush=True)


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
