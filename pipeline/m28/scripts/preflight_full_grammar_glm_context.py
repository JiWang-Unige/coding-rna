#!/usr/bin/env python3
"""Original TRAIN raw-view GLM -> fresh emissions -> full grammar; zero updates."""
from collections import Counter
import json
import math
import os
import re
import resource
import time

import torch

import probe_full_grammar_R4_train_fasta_context as P
from cache_full_train import load_native
from annotate_query import encode_window
from src.m28.core import M28Core
from src.m28.training import one_hot
from src.m28.contig_emission_lattice import context_segments, stitch_contig_emissions
from src.m28.full_grammar_context_view import context_geometry, owned_context_logits

Q,W,ROOT = P.Q,P.W,P.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-GLM-CONTEXT-PREFLIGHT'
SOURCE_JOB = '13428966'
CONFIG = dict(scope='same_original_eight_TRAIN_witnesses_all_core_observations',
    view='original_24576bp_core_plus_one_12288bp_stride_neighbor_each_side_raw_length_clipped_deduplicated',
    ownership='nearest_containing_genomic_window_center_per_base_lower_start_ties_same_as_context_segments',
    reference_boundaries_used_to_choose_view=False,view_candidate_for_preparation_not_approved_fit_contract=True,
    encoder_model=Q.F.MODEL,encoder_revision=Q.F.REVISION,encoder_frozen=True,
    encoder_blocks_bp=6144,encoder_token_bp=6,encoder_materialization='BF16_hidden_then_FP32_clone_outside_inference_mode',
    model='fresh_seed0_M28Core_C0_2048_128_366191_parameters_FP32',
    targets='unchanged_R4_build_targets_fine_else_coarse_else_unknown_all_core_rows',
    max_exons=128,checkpoint_bp=256,role_weights=1,null=0,score_accumulation='FP64',
    loss='unweighted_unnormalized_global_logZ_all_minus_core_observed_logZ_preflight',
    finite_NLL_nonnegative_atol=1e-8,NLL_common_row_shift_gradient_atol=1e-5,
    NLL_occupancy_difference_bound_atol=1e-6,owned_row_VJP_comparison_rtol=0,owned_row_VJP_comparison_atol=0,
    optimizer_steps=0)


def owned_stitch_checks():
    checks = []
    for windows in (((0,12),(6,18),(12,24)),((0,14),(6,19),(13,24))):
        geometry = dict(local_windows=windows,base_owned_segments=context_segments(24,windows))
        for strand in ('+','-'):
            sequence = 'ACGT'*6
            values = [(torch.arange((b-a)*15,dtype=torch.float32).reshape(b-a,15)+i*1000).requires_grad_()
                      for i,(a,b) in enumerate(windows)]
            old = stitch_contig_emissions(sequence,windows,strand,[dict(window_index=i,
                sequence=sequence[a:b] if strand == '+' else Q.M.rc(sequence[a:b]),local_potentials=v.detach())
                for i,((a,b),v) in enumerate(zip(windows,values))])['local_potentials']
            joined = owned_context_logits(values,geometry,strand)
            if not torch.equal(joined,old):
                raise ValueError('Differentiable ownership differs from original raw ownership')
            joined.sum().backward()
            for value,(a,b),(_,left,right) in zip(values,windows,geometry['base_owned_segments']):
                expected = torch.zeros_like(value)
                lo,hi = (left-a,right-a) if strand == '+' else (b-right,b-left)
                expected[lo:hi] = 1
                if not torch.equal(value.grad,expected):
                    raise ValueError('Selected owner rows, discarded rows or strand gradients differ')
            checks.append(dict(windows=windows,strand=strand,forward_exact=True,owned_VJP_exact=True))
    return checks


def main():
    began = time.perf_counter()
    source = json.loads((P.OUT/'summary.json').read_text())
    if ((P.OUT/'JOBID').read_text().strip() != SOURCE_JOB or (P.OUT/'STATUS').read_text().strip() != 'COMPLETED'
            or source['job_id'] != SOURCE_JOB or not source['completed'] or source['configuration'] != P.CONFIG
            or source['context_NLL_and_gradient_windows'] != 8 or source['optimizer_steps_this_run'] != 0
            or not source['all_original_context_NLLs_executed'] or source['DEV_used']):
        raise ValueError('Actual completed original eight real-DNA diagnostic source required')
    witness,_ = P.P.ready()
    ids = [r['window_id'] for r in source['records']]
    if ids != [r['window_id'] for r in witness['records']]:
        raise ValueError('All eight original ordered identities required')
    with (OUT/'started.json').open('x') as handle:
        json.dump(dict(job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,source_job=SOURCE_JOB,
            witness_ids=ids,TRAIN_allowlists={s:sorted(v) for s,v in Q.ALLOW.items()}),handle,indent=2,allow_nan=False)
    checks = owned_stitch_checks()
    t = time.perf_counter(); draws,index,windows,catalogs = Q.load_train_inputs()
    metadata_seconds = time.perf_counter()-t
    t = time.perf_counter(); sequences,descriptors = P.read_train_DNA()
    raw_seconds = time.perf_counter()-t
    census = {s:Counter() for s in Q.ALLOW}; outside_examples = []
    for w in windows.values():
        geometry = context_geometry(w['start'],w['end'],len(sequences[w['species']][w['seqid']]))
        a,b = geometry['genomic_context']; counts = census[w['species']]
        counts['windows'] += 1; counts['encoder_windows'] += len(geometry['genomic_windows'])
        outside = [tid for tid in w['overlapping_transcript_ids'] if catalogs[w['species']][tid]['training_primary']
                   and (catalogs[w['species']][tid]['CDS'][0][0] < a or catalogs[w['species']][tid]['CDS'][-1][1] > b)]
        counts['original_primary_window_pairs_outside_context'] += len(outside)
        counts['windows_with_original_primary_outside_context'] += bool(outside)
        if outside and len(outside_examples) < 12:
            outside_examples.append(dict(window_id=w['id'],genomic_context=[a,b],primary_ids=outside))
    t = time.perf_counter(); backbone,tokenizer = load_native(); torch.cuda.synchronize()
    encoder_load_seconds = time.perf_counter()-t
    torch.manual_seed(0); model = M28Core('C0',2048,128).to('cuda').train()
    if sum(p.numel() for p in model.parameters()) != 366191:
        raise ValueError('Original fresh C0 representation capacity changed')
    before = {k:v.detach().cpu().clone() for k,v in model.state_dict().items()}
    torch.cuda.reset_peak_memory_stats(); opened,records = set(),[]
    with (OUT/'windows.jsonl').open('x') as handle:
        for original in witness['records']:
            row_started = time.perf_counter(); ident = original['window_id']; w = windows[ident]; n = w['valid_bases']
            chromosome = sequences[w['species']][w['seqid']]
            geometry = context_geometry(w['start'],w['end'],len(chromosome)); a,b = geometry['genomic_context']
            raw_core = chromosome[w['start']:w['end']]
            core = raw_core if w['strand'] == '+' else Q.M.rc(raw_core)
            cached = W.load_cache(w,index,opened)
            if core != cached['sequence']:
                raise ValueError('Original raw core DNA differs from actual cached core: '+ident)
            del cached
            targets = P.P.build_targets(w,catalogs[w['species']],core)
            allowed = P.P.allowed_roles_from_local_targets(targets,n)
            fine = int((torch.as_tensor(targets['seg'][:n]) >= 0).sum())
            free = int(allowed.all(1).sum()); coarse = n-fine-free
            if (fine,coarse,free) != (original['fine_positions'],original['coarse_positions'],original['fully_unknown_positions']):
                raise ValueError('Original core observation counts changed: '+ident)
            genomic = re.sub('[^ACGT]','N',chromosome[a:b])
            context = genomic if w['strand'] == '+' else Q.M.rc(genomic)
            offset = w['start']-a if w['strand'] == '+' else b-w['end']
            if context[offset:offset+n] != re.sub('[^ACGT]','N',core):
                raise ValueError('Normalized oriented core or coordinates changed')
            model.zero_grad(set_to_none=True); tensors = []; encoding = []
            for left,right in geometry['genomic_windows']:
                raw = re.sub('[^ACGT]','N',chromosome[left:right])
                dna_string = raw if w['strand'] == '+' else Q.M.rc(raw)
                t = time.perf_counter()
                # Native encoder returns inference tensors. Clone AFTER that
                # scope so the trainable projection may save a normal input.
                features = encode_window(backbone,tokenizer,dna_string,'cuda').clone()
                dna = one_hot(dna_string,Q.M.WIDTH)[None].to('cuda')
                logits = model(features,dna,[right-left])['segmentation_logits'][0]
                logits.retain_grad(); tensors.append(logits)
                torch.cuda.synchronize()
                encoding.append(dict(genomic_window=[left,right],valid_bases=right-left,
                    feature_shape=list(features.shape),features_require_grad=features.requires_grad,
                    feature_is_inference_tensor=torch.is_inference(features),seconds=time.perf_counter()-t))
            values = [v[:right-left] for v,(left,right) in zip(tensors,geometry['genomic_windows'])]
            joined = owned_context_logits(values,geometry,w['strand']).to('cpu'); joined.retain_grad()
            if tuple(joined.shape) != (len(context),15) or joined.dtype != torch.float32 or not bool(torch.isfinite(joined).all()):
                raise ValueError('Finite actual-context FP32 neural potentials required')
            cuts = tuple(x for x in (offset,offset+n) if 0 < x < len(context))
            print(json.dumps(dict(stage='actual_GLM_neural_context_NLL_started',window_id=ident,
                context_length=len(context),core_oriented_interval=[offset,offset+n],encoder_windows=len(tensors))),flush=True)
            t = time.perf_counter()
            result = P.P.compiled_boundary_core_observation_nll(context,joined,128,offset,offset+n,allowed,cuts=cuts,checkpoint_bp=256)
            forward_seconds = time.perf_counter()-t
            masses = {k:float(result[k].detach()) for k in ('log_partition_all','log_partition_observed','loss')}
            if not all(math.isfinite(v) for v in masses.values()) or masses['loss'] < -1e-8:
                raise ValueError('Finite neural logZs and nonnegative NLL required; observed logZ may be negative')
            t = time.perf_counter(); result['loss'].backward(); torch.cuda.synchronize()
            backward_seconds = time.perf_counter()-t; gradient = joined.grad
            if gradient is None or not bool(torch.isfinite(gradient).all()):
                raise ValueError('Every assembled neural role derivative must be finite')
            row_shift = float(gradient.sum(1).abs().max())
            role_bound = float(gradient[:,1:].abs().max()); covered_bound = float(gradient[:,1:].sum(1).abs().max())
            if row_shift > 1e-5 or role_bound > 1+1e-6 or covered_bound > 1+1e-6:
                raise ValueError('Actual neural NLL posterior differences exceed original bounds')
            for value,(left,right),(_,lo,hi) in zip(tensors,geometry['local_windows'],geometry['base_owned_segments']):
                expected = torch.zeros(value.shape,dtype=torch.float32)
                x,y = (lo-left,hi-left) if w['strand'] == '+' else (right-hi,right-lo)
                g,h = (lo,hi) if w['strand'] == '+' else (len(context)-hi,len(context)-lo)
                expected[x:y] = gradient[g:h]
                if value.grad is None or not torch.equal(value.grad.cpu(),expected):
                    raise ValueError('Actual selected owner-row VJP or discarded/padding rows differ')
            groups = {}
            for name,module in (('project',model.shared.project),('context',model.shared.context),
                ('dna',model.shared.dna),('fuse',model.shared.fuse),('segmentation',model.segmentation)):
                params = list(module.parameters())
                if any(p.grad is None or not bool(torch.isfinite(p.grad).all()) for p in params):
                    raise ValueError('Missing/nonfinite actual neural parameter gradient: '+name)
                groups[name] = math.sqrt(sum(float(p.grad.double().square().sum()) for p in params))
                if groups[name] <= 0:
                    raise ValueError('Actual neural representation gradient is zero: '+name)
            row = dict(ordinal=original['ordinal'],window_id=ident,status='COMPLETED_NEURAL_CONTEXT_NLL_AND_SHARED_PARAMETER_GRADIENTS',
                geometry=geometry,strand=w['strand'],context_length=len(context),core_oriented_interval=[offset,offset+n],
                view_chosen_without_CDS_boundaries=True,raw_core_matches_original_cache=True,
                fine_positions=fine,coarse_positions=coarse,fully_unknown_positions=free,encoding=encoding,
                **masses,NLL_forward_seconds=forward_seconds,NLL_backward_seconds=backward_seconds,
                all_context_role_gradients_finite=True,all_shared_parameter_gradients_finite=True,
                shared_parameter_gradient_l2=groups,owned_window_row_VJP_exact=True,
                common_row_shift_gradient_max_abs=row_shift,maximum_abs_role_posterior_difference=role_bound,
                maximum_abs_covered_posterior_difference=covered_bound,
                prefix_role_gradient_l1=float(gradient[:offset].double().abs().sum()),
                core_role_gradient_l1=float(gradient[offset:offset+n].double().abs().sum()),
                suffix_role_gradient_l1=float(gradient[offset+n:].double().abs().sum()),
                row_seconds=time.perf_counter()-row_started)
            records.append(row); handle.write(json.dumps(row,allow_nan=False)+'\n'); handle.flush()
            print(json.dumps({k:row[k] for k in ('ordinal','window_id','status','loss','row_seconds')}),flush=True)
            del result,joined,gradient,tensors,values,features,dna,logits,targets,allowed
    if len(records) != 8 or [r['window_id'] for r in records] != ids:
        raise ValueError('All original eight real neural contexts required')
    if any(not torch.equal(v.cpu(),before[k]) for k,v in model.state_dict().items()):
        raise ValueError('Zero-update preparation changed the fresh neural state')
    if any(p.requires_grad or p.grad is not None for p in backbone.parameters()):
        raise ValueError('Released GLM must remain frozen with no parameter gradients')
    summary = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,source_job=SOURCE_JOB,
        protocol='reports/M28-METHOD-RESTART/full_grammar_partition_preparation.md',records=records,
        original_draws=len(draws),original_train_unique_windows=len(windows),actual_witnesses=8,
        context_geometry_population={s:dict(c) for s,c in census.items()},outside_primary_examples=outside_examples,
        complete_primary_extent_census_only_not_full_cohort_NLL_support=True,full_cohort_observation_support_proven=False,
        original_TRAIN_DNA_sources=descriptors,cached_artifacts_opened=sorted(opened),owned_stitch_checks=checks,
        pretrained_encoder_frozen_no_grad=True,actual_released_GLM_forward=True,actual_fresh_neural_forward=True,
        all_eight_NLL_role_and_shared_parameter_gradients_completed=True,new_model_weights_unchanged=True,
        optimizer_constructed=False,optimizer_steps_this_run=0,checkpoints_written=False,supervised_checkpoint_read=False,
        original_R4_policy_changed=False,observation_rows_deleted=False,primary_transcripts_dropped=False,
        boundary_prior_chosen=False,approved_fit_context_or_raw_chromosome_protocol_selected=False,
        loss_normalizer_importance_rule_or_fit_coefficients_selected=False,cached_features_used=False,
        old_cached_labels_used=False,prediction_candidate_pool_used=False,DEV_used=False,independent_species_used=False,
        A33_outputs_opened=False,A44_execution_chain_reopened=False,original_TRAIN_FASTA_read_seconds=raw_seconds,
        original_TRAIN_metadata_load_seconds=metadata_seconds,encoder_load_seconds=encoder_load_seconds,
        seconds_after_module_import=time.perf_counter()-began,process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved(),
        interpretation='raw_coordinate_view_released_GLM_and_full_grammar_NLL_training_wiring_not_fit_accuracy_or_generalization')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(summary,handle,indent=2,allow_nan=False)
    print(json.dumps({k:v for k,v in summary.items() if k!='records'},allow_nan=False),flush=True)


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
