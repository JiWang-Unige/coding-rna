#!/usr/bin/env python3
"""Six fixed raw-view extent warnings: actual core-event support, not a fit."""
import json
import math
import os
import re
import time
import torch

import preflight_full_grammar_glm_context as G
from src.m28.full_grammar_context_view import context_geometry
from src.m28.full_grammar_boundary_compiled import compiled_boundary_full_grammar_log_partition
from src.m28.full_grammar_observed_partition import InfeasibleObservation

Q,W,P,ROOT = G.Q,G.W,G.P,G.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-CONTEXT-EXTENT-WARNINGS'
SOURCE_JOB = '13429053'
CONFIG = dict(source_view_configuration=G.CONFIG,
    scope='all_six_exact_original_TRAIN_extent_warning_windows_no_resampling',
    question='core_role_event_support_not_complete_outside_primary_chain_representability',
    potentials='all_zero_CPU_FP32_support_independent_of_any_finite_neural_potentials',
    max_exons=128,checkpoint_bp=256,role_weights=1,null=0,score_accumulation='FP64',
    targets='unchanged_R4_build_targets_all_fine_coarse_unknown_core_rows',
    context='unchanged_original_stride_neighbor_view_no_CDS_based_expansion',optimizer_steps=0)


def main():
    began = time.perf_counter()
    source = json.loads((G.OUT/'summary.json').read_text())
    if ((G.OUT/'STATUS').read_text().strip() != 'COMPLETED' or (G.OUT/'JOBID').read_text().strip() != SOURCE_JOB
            or source['job_id'] != SOURCE_JOB or not source['completed'] or source['configuration'] != G.CONFIG
            or not source['all_eight_NLL_role_and_shared_parameter_gradients_completed']
            or source['optimizer_steps_this_run'] != 0 or source['DEV_used']):
        raise ValueError('Actual complete unchanged original neural-view preparation required')
    warnings = source['outside_primary_examples']; ids = [r['window_id'] for r in warnings]
    if len(ids) != 6 or len(set(ids)) != 6 or sum(c['windows_with_original_primary_outside_context']
            for c in source['context_geometry_population'].values()) != 6:
        raise ValueError('Exactly all six source extent warnings required')
    with (OUT/'started.json').open('x') as handle:
        json.dump(dict(job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,source_job=SOURCE_JOB,
            warning_window_ids=ids,TRAIN_allowlists={s:sorted(v) for s,v in Q.ALLOW.items()}),handle,indent=2)
    draws,index,windows,catalogs = Q.load_train_inputs()
    sequences,descriptors = P.read_train_DNA()
    actual_warnings = []
    for w in windows.values():
        geometry = context_geometry(w['start'],w['end'],len(sequences[w['species']][w['seqid']]))
        a,b = geometry['genomic_context']
        outside = [tid for tid in w['overlapping_transcript_ids'] if catalogs[w['species']][tid]['training_primary']
                   and (catalogs[w['species']][tid]['CDS'][0][0] < a or catalogs[w['species']][tid]['CDS'][-1][1] > b)]
        if outside:
            actual_warnings.append(dict(window_id=w['id'],genomic_context=[a,b],primary_ids=outside))
    if actual_warnings != warnings:
        raise ValueError('Original all-six warning census changed')
    opened,records = set(),[]
    with torch.inference_mode(), (OUT/'windows.jsonl').open('x') as handle:
        for ordinal,warning in enumerate(warnings):
            t = time.perf_counter(); ident = warning['window_id']; w = windows[ident]; n = w['valid_bases']
            chromosome = sequences[w['species']][w['seqid']]
            geometry = context_geometry(w['start'],w['end'],len(chromosome)); a,b = geometry['genomic_context']
            raw = chromosome[w['start']:w['end']]; core = raw if w['strand'] == '+' else Q.M.rc(raw)
            cached = W.load_cache(w,index,opened)
            if cached['sequence'] != core:
                raise ValueError('Original raw core differs from matching cache: '+ident)
            del cached
            targets = P.P.build_targets(w,catalogs[w['species']],core)
            allowed_core = P.P.allowed_roles_from_local_targets(targets,n)
            fine = int((torch.as_tensor(targets['seg'][:n]) >= 0).sum())
            free = int(allowed_core.all(1).sum()); coarse = n-fine-free
            genomic = re.sub('[^ACGT]','N',chromosome[a:b])
            context = genomic if w['strand'] == '+' else Q.M.rc(genomic)
            offset = w['start']-a if w['strand'] == '+' else b-w['end']
            if context[offset:offset+n] != re.sub('[^ACGT]','N',core):
                raise ValueError('Actual normalized oriented core does not match the fixed view')
            mask = torch.ones((len(context),15),dtype=torch.bool); mask[offset:offset+n] = allowed_core
            potentials = torch.zeros((len(context),15),dtype=torch.float32)
            cuts = tuple(x for x in (offset,offset+n) if 0 < x < len(context))
            print(json.dumps(dict(stage='fixed_warning_observed_support_started',window_id=ident)),flush=True)
            try:
                result = compiled_boundary_full_grammar_log_partition(context,potentials,128,
                    allowed=mask,cuts=cuts,checkpoint_bp=256)
                log_mass = float(result['log_partition'])
                if not math.isfinite(log_mass) or log_mass < -1e-8:
                    raise ValueError('Finite nonnegative zero-potential observed log count required')
                feasible,reason = True,None
            except InfeasibleObservation as error:
                feasible,log_mass,reason = False,None,str(error)
            row = dict(ordinal=ordinal,window_id=ident,species=w['species'],strand=w['strand'],
                diagnostic_genomic_context=[a,b],core_oriented_interval=[offset,offset+n],
                context_length=len(context),original_raw_core_matches_cache=True,
                fine_positions=fine,coarse_positions=coarse,fully_unknown_positions=free,
                original_outside_primary_ids=warning['primary_ids'],
                original_outside_primary_CDS_spans={tid:[catalogs[w['species']][tid]['CDS'][0][0],
                    catalogs[w['species']][tid]['CDS'][-1][1]] for tid in warning['primary_ids']},
                complete_outside_primary_chain_representable_in_view=False,
                core_observation_feasible=feasible,zero_potential_observed_log_count=log_mass,
                infeasible_reason=reason,all_core_observation_rows_preserved=True,
                status='FEASIBLE_CORE_EVENT' if feasible else 'INFEASIBLE_CORE_EVENT',
                seconds=time.perf_counter()-t)
            records.append(row); handle.write(json.dumps(row,allow_nan=False)+'\n'); handle.flush()
            print(json.dumps(row,allow_nan=False),flush=True)
            del targets,allowed_core,mask,potentials
    if [r['window_id'] for r in records] != ids or opened != {index[i]['artifact'] for i in ids}:
        raise ValueError('All six cases and only their original caches required')
    summary = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,source_job=SOURCE_JOB,
        records=records,actual_warning_windows=6,all_six_warning_windows_accounted=True,
        feasible_core_events=sum(r['core_observation_feasible'] for r in records),
        infeasible_core_events=sum(not r['core_observation_feasible'] for r in records),
        original_draws=len(draws),original_train_unique_windows=len(windows),
        original_TRAIN_DNA_sources=descriptors,cached_artifacts_opened=sorted(opened),
        full_cohort_observation_support_proven=False,complete_original_primary_coverage_proven=False,
        context_expanded_from_annotations=False,original_R4_policy_changed=False,observation_rows_deleted=False,
        primary_transcripts_dropped=False,model_or_backbone_forward=False,weights_read=False,
        optimizer_steps_this_run=0,checkpoints_written=False,cached_features_used=False,old_cached_labels_used=False,
        DEV_used=False,independent_species_used=False,A33_outputs_opened=False,A44_execution_chain_reopened=False,
        seconds_after_module_import=time.perf_counter()-began,
        interpretation='all_six_fixed_extent_warning_core_event_support_not_whole_primary_chain_or_full_cohort_or_model_performance')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(summary,handle,indent=2,allow_nan=False)
    print(json.dumps({k:v for k,v in summary.items() if k!='records'},allow_nan=False),flush=True)


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
