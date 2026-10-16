#!/usr/bin/env python3
"""A44 paired last emissions, label-free grammar, then original TRAIN labels."""
import argparse
from collections import Counter
import json
import resource
import time
import torch

import fit_local_structured_emission as F
import readout_natural_competitor_train as D

ROOT, R, PROTOCOL = F.ROOT, D.R, F.PROTOCOL
ARMS = F.ARMS
FIT_JOBS = {'A44_LOCAL': '13425946', 'A44_JOINT': '13425945'}
CONFIG = dict(source_checkpoints={a: str((F.output_root(a)/'step_004608.pt').relative_to(ROOT)) for a in ARMS},
    source_fit_jobs=FIT_JOBS, source_fit_configuration=F.CONFIG,
    frozen_parameters_per_arm=366191, source_optimizer_steps_per_arm=4608,
    population='original_1491_unique_TRAIN_windows_first_occurrence_order',
    model_precision='FP32', score_accumulation='FP64', max_exons=128, threshold=0., null=0.,
    all_role_weights=1., search=D.CONFIG['search'], prediction_candidate_pool_input=False,
    prediction_label_input=False, unique_TRAIN_windows=1491, original_draws=1536,
    arms=list(ARMS), shards_per_arm=32, tasks=64, task_assignment='arm_index_times32_plus_original_ordinal_modulo32',
    main_checkpoint='last_step_004608_only', optimizer_steps_this_run=0)
CHECKS = ROOT/'outputs/M28-A44-PAIRED-READOUT-CHECKS'
POTENTIALS = ROOT/'outputs/M28-A44-PAIRED-TRAIN-POTENTIALS'
DECODE = ROOT/'outputs/M28-A44-PAIRED-TRAIN-DECODE'
OUT = ROOT/'outputs/M28-A44-PAIRED-TRAIN-READOUT'
CHECK = ROOT/'outputs/M28-A44-PAIRED-RESULT-CHECK'
SCOPE = dict(optimizer_steps_this_run=0, new_backbone_forward=False, R4_policy_or_catalog_changed=False,
             DEV_used=False, test_setaria_used=False, independent_species_used=False, A33_outputs_opened=False)


def validate_fit(value, arm, state=False):
    if (arm not in ARMS or value['arm'] != arm or value['configuration'] != F.CONFIG
            or value['protocol'] != F.PROTOCOL or value['step' if state else 'optimizer_steps'] != 4608
            or value['revision' if state else 'feature_revision'] != R.C.I.F.REVISION):
        raise ValueError('Only the actual matching A44 paired last004608 allowed')
    if not state and (not value['completed'] or value['source_support_job'] != '13424089'
            or value['source_preflight_job'] != '13425772' or value['position_checks_job'] != '13425929'
            or value['final_checkpoint'] != CONFIG['source_checkpoints'][arm]
            or value['trainable_parameters'] != 366191
            or value['initialization'] != 'fresh_seed0_no_supervised_weight_inheritance'
            or not value['same_capacity_fresh_seed0_paired_initialization_verified']
            or value['loss_coefficients'] != F.CONFIG['loss_coefficients_by_arm'][arm]
            or value['structured_candidate_supervision_used'] != (arm == 'A44_JOINT')
            or not value['original_TRAIN_catalog_metadata_read']
            or not value['original_R4_local_targets_regenerated']
            or any(value[k] for k in ('supervised_checkpoint_read', 'old_cached_label_fields_used',
                'R4_policy_or_catalog_changed', 'generator_updated', 'new_backbone_forward',
                'DEV_used_for_selection', 'test_setaria_used', 'independent_species_used', 'A33_outputs_opened'))):
        raise ValueError('Actual complete original-R4 paired fit identity differs')


def require_fits():
    F.A.B.D.B.F.A.require_jobs(list(FIT_JOBS.values()))
    values = {a: R.completed(F.output_root(a), FIT_JOBS[a]) for a in ARMS}
    for arm, value in values.items(): validate_fit(value, arm)
    return values


def task_identity(task):
    if not 0 <= task < 64:
        raise ValueError('Only the fixed 64 paired tasks allowed')
    return ARMS[task//32], task % 32


def potential_path(arm, ordinal):
    return POTENTIALS/arm/('window_'+str(ordinal).zfill(6)+'.pt')


def shard_folder(task):
    return DECODE/('task_'+str(task).zfill(2))


def label_paired_record(rows, old, original, targets):
    result = None
    for arm in ARMS:
        if rows[arm]['arm'] != arm:
            raise ValueError('Paired label join arm differs')
        local = D.label_record(rows[arm], old, original, targets)
        if result is None:
            result = dict(window=local['window'], labels={'A31_A3': local['labels']['A31_A3']},
                counts={'A31_A3': local['counts']['A31_A3']},
                eligible_references=local['eligible_references'], eligible_five_or_more=local['eligible_five_or_more'])
        result['labels'][arm] = local['labels']['A43']
        result['counts'][arm] = local['counts']['A43']
    return result


def aggregate(rows):
    totals = {s: {a: Counter() for a in ('A31_A3', *ARMS)} for s in R.C.R.ALLOW}
    refs = {s: Counter(eligible=0, five_or_more=0) for s in R.C.R.ALLOW}
    for row in rows:
        s = row['window']['species']
        for arm in totals[s]: totals[s][arm].update(row['counts'][arm])
        refs[s].update(eligible=row['eligible_references'], five_or_more=row['eligible_five_or_more'])
    return ({s: {a: dict(v) for a, v in arms.items()} for s, arms in totals.items()},
            {s: dict(v) for s, v in refs.items()})


def comparison_gate(counts, candidate, baseline):
    checks = {}
    for s, values in counts.items():
        checks[s+'_correct_strictly_increased'] = values[candidate]['selected_correct'] > values[baseline]['selected_correct']
        checks[s+'_known_nonmatch_not_increased'] = values[candidate]['selected_known_nonmatch'] <= values[baseline]['selected_known_nonmatch']
    checks['pooled_five_or_more_correct_strictly_increased'] = (
        sum(v[candidate]['selected_five_or_more_correct'] for v in counts.values()) >
        sum(v[baseline]['selected_five_or_more_correct'] for v in counts.values()))
    return checks


def gates(counts):
    values = dict(joint_vs_A31_A3=comparison_gate(counts, 'A44_JOINT', 'A31_A3'),
        joint_vs_LOCAL=comparison_gate(counts, 'A44_JOINT', 'A44_LOCAL'),
        LOCAL_vs_A31_A3=comparison_gate(counts, 'A44_LOCAL', 'A31_A3'))
    return values, all(values['joint_vs_A31_A3'].values()) and all(values['joint_vs_LOCAL'].values())


def preparation():
    draws, windows, index = D.population()
    old_draws, free, prior, old_index, old = D.A42.B.population()
    if (draws != old_draws or [r['window'] for r in free] != [r['window'] for r in windows]
            or any(index[i]['artifact'] != old_index[i]['artifact'] for i in index)):
        raise ValueError('Original paired label-free prediction inputs differ')
    labels = R.C.R.rows(R.READOUT/'windows.jsonl')
    placeholders = [dict(window=r['window'], reference_used_in_scoring=False, **r['variants']['A31_A3']) for r in prior]
    rows, counts, refs = D.A42.B.paired_readout(placeholders, prior, labels)
    if (len(rows) != 1491 or sum(v['eligible'] for v in refs.values()) != 2916
            or sum(v['five_or_more'] for v in refs.values()) != 1086
            or any(counts[s]['A31_A3'] != old['selected_counts_per_species'][s]['A31_A3'] for s in counts)):
        raise ValueError('Complete original stronger baseline and denominators differ')
    R.C.R.write_summary(CHECKS, dict(completed=True, configuration=CONFIG, protocol=PROTOCOL,
        windows=1491, original_draws=1536, original_natural_candidates=16321,
        original_prediction_inputs_baseline_labels_and_denominators_reconciled=True,
        fitted_model_or_feature_tensors_opened=False, TRAIN_catalog_opened=False,
        model_or_backbone_forward=False, new_geometry_labels_computed=False, **SCOPE))


def require_ready(job):
    F.A.B.D.B.F.A.require_jobs([job]); value = R.completed(CHECKS, job)
    if (value['configuration'] != CONFIG or value['protocol'] != PROTOCOL
            or not value['original_prediction_inputs_baseline_labels_and_denominators_reconciled']):
        raise ValueError('Closed actual paired readout preparation required')


@torch.inference_mode()
def potentials(checks_job):
    require_ready(checks_job); fitted = require_fits(); models = {}; before = {}
    for arm in ARMS:
        state = torch.load(ROOT/fitted[arm]['final_checkpoint'], map_location='cpu', weights_only=True)
        validate_fit(state, arm, state=True)
        model = F.A.A36.new_model('cuda').eval().requires_grad_(False)
        model.load_state_dict(state['model'], strict=True); del state
        if sum(p.numel() for p in model.parameters()) != 366191:
            raise ValueError('Paired fixed source model capacity differs')
        models[arm] = model; before[arm] = {k: v.detach().clone() for k, v in model.state_dict().items()}
    draws, windows, index = D.population(); receipts = {a: [] for a in ARMS}
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    for ordinal, record in enumerate(windows):
        w = record['window']; cached = torch.load(ROOT/index[w['id']]['artifact'], map_location='cpu', weights_only=True)
        features, dna = R.Q.window_inputs(cached, w, 'cuda')
        for arm in ARMS:
            local = models[arm](features, dna, [w['valid_bases']])['segmentation_logits'][0][:w['valid_bases']]
            if local.dtype != torch.float32 or tuple(local.shape) != (w['valid_bases'], 15) or not bool(torch.isfinite(local).all()):
                raise ValueError('Finite paired last-model valid-base FP32 potentials required')
            path = potential_path(arm, ordinal)
            with path.open('xb') as handle:
                torch.save(dict(configuration=CONFIG, protocol=PROTOCOL, arm=arm, original_ordinal=ordinal,
                    fit_job=FIT_JOBS[arm], window=w, sequence=cached['sequence'], local_logits=local.cpu(),
                    revision=R.C.I.F.REVISION, reference_used_in_scoring=False), handle)
            receipts[arm].append(dict(original_ordinal=ordinal, window=w, artifact=str(path.relative_to(ROOT))))
        del cached, features, dna, local
        if (ordinal+1) % 64 == 0 or ordinal+1 == 1491: print(json.dumps(dict(paired_potential_windows=ordinal+1, total=1491)), flush=True)
    for arm, model in models.items():
        if (any(not torch.equal(v, before[arm][k]) for k, v in model.state_dict().items())
                or any(p.requires_grad or p.grad is not None for p in model.parameters())):
            raise ValueError('Frozen paired model changed or received gradients')
    R.C.R.write_summary(POTENTIALS, dict(completed=True, configuration=CONFIG, protocol=PROTOCOL,
        source_fit_jobs=FIT_JOBS, source_checks_job=checks_job, checkpoints=CONFIG['source_checkpoints'],
        windows_per_arm=1491, original_draws=1536, records=receipts, feature_revision=R.C.I.F.REVISION,
        all_potential_files_closed=True, frozen_source_weights_unchanged=True,
        all_source_parameters_frozen_and_gradients_absent=True, prediction_candidate_pool_opened=False,
        reference_labels_opened=False, TRAIN_catalog_opened=False, checkpoints_written=False,
        new_candidate_search=False, new_geometry_labels_computed=False, seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(), **SCOPE))


def require_potentials(job):
    F.A.B.D.B.F.A.require_jobs([job]); value = R.completed(POTENTIALS, job)
    if (value['configuration'] != CONFIG or value['protocol'] != PROTOCOL or value['source_fit_jobs'] != FIT_JOBS
            or value['checkpoints'] != CONFIG['source_checkpoints'] or value['windows_per_arm'] != 1491
            or any(len(value['records'][a]) != 1491 for a in ARMS)
            or not value['all_potential_files_closed'] or not value['frozen_source_weights_unchanged']
            or not value['all_source_parameters_frozen_and_gradients_absent']
            or value['prediction_candidate_pool_opened'] or value['reference_labels_opened']
            or value['TRAIN_catalog_opened'] or value['optimizer_steps_this_run'] != 0):
        raise ValueError('Complete matching paired label-free potentials required')
    return value


def load_potential(arm, ordinal, receipt, report):
    path = potential_path(arm, ordinal); value = torch.load(path, map_location='cpu', weights_only=True)
    w, local = value['window'], value['local_logits']
    if (receipt['artifact'] != str(path.relative_to(ROOT)) or receipt['original_ordinal'] != ordinal
            or value['original_ordinal'] != ordinal or value['arm'] != arm or value['fit_job'] != FIT_JOBS[arm]
            or value['configuration'] != CONFIG or value['protocol'] != PROTOCOL or w != receipt['window']
            or value['revision'] != report['feature_revision'] or value['reference_used_in_scoring']
            or w['split'] != 'train' or w['seqid'] not in R.C.R.ALLOW[w['species']]
            or len(value['sequence']) != w['valid_bases'] or local.dtype != torch.float32
            or tuple(local.shape) != (w['valid_bases'], 15) or not bool(torch.isfinite(local).all())):
        raise ValueError('Matching original paired prediction potential required')
    return value


@torch.inference_mode()
def decode(task, potential_job):
    arm, shard = task_identity(task); report = require_potentials(potential_job)
    ordinals = list(range(shard, 1491, 32)); out = shard_folder(task); began = time.perf_counter()
    with (out/'natural_scores.jsonl').open('x') as handle:
        for ordinal in ordinals:
            value = load_potential(arm, ordinal, report['records'][arm][ordinal], report); started = time.perf_counter()
            natural = D.A42.P.decode_emission_lattice(value['sequence'].upper(), value['local_logits'], CONFIG['max_exons'])
            D.validate_prediction(value, natural); seconds = time.perf_counter()-started
            row = dict(arm=arm, original_ordinal=ordinal, task=task, shard=shard, window=value['window'],
                chains=natural['chains'], scores=natural['chain_gains'], selected=list(range(len(natural['chains']))),
                gain=natural['gain'], potential_artifact=report['records'][arm][ordinal]['artifact'],
                feature_revision=value['revision'], potential_job=potential_job, source_fit_job=FIT_JOBS[arm],
                reference_used_in_scoring=False, decode_seconds=seconds)
            handle.write(json.dumps(row, allow_nan=False)+'\n'); handle.flush()
            print(json.dumps(dict(arm=arm, task=task, original_ordinal=ordinal, decode_seconds=seconds)), flush=True)
            del value, natural
    R.C.R.write_summary(out, dict(completed=True, configuration=CONFIG, protocol=PROTOCOL,
        arm=arm, task=task, shard=shard, potential_job=potential_job, source_fit_job=FIT_JOBS[arm],
        original_ordinals=ordinals, windows=len(ordinals), natural_score_file_closed=True,
        all_selected_ORF_motif_cap_span_and_direct_role_gains_verified=True,
        prediction_candidate_pool_opened=False, reference_labels_opened=False, TRAIN_catalog_opened=False,
        model_or_backbone_forward=False, new_candidate_search=True, new_geometry_labels_computed=False,
        seconds=time.perf_counter()-began, process_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, **SCOPE))


def closed_shards(decode_job, potential_job):
    F.A.B.D.B.F.A.require_jobs([decode_job+'_'+str(i) for i in range(64)]); rows = {a: [] for a in ARMS}
    for task in range(64):
        arm, shard = task_identity(task); out = shard_folder(task)
        value = R.completed(out, decode_job+'_'+str(task)); local = R.C.R.rows(out/'natural_scores.jsonl')
        expected = list(range(shard, 1491, 32))
        if (value['configuration'] != CONFIG or value['protocol'] != PROTOCOL or value['arm'] != arm
                or value['task'] != task or value['shard'] != shard or value['potential_job'] != potential_job
                or value['source_fit_job'] != FIT_JOBS[arm] or value['original_ordinals'] != expected
                or value['windows'] != len(expected) or [r['original_ordinal'] for r in local] != expected
                or any(r['arm'] != arm or r['task'] != task for r in local)
                or not value['natural_score_file_closed']
                or not value['all_selected_ORF_motif_cap_span_and_direct_role_gains_verified']
                or value['prediction_candidate_pool_opened'] or value['reference_labels_opened']
                or value['TRAIN_catalog_opened'] or value['new_geometry_labels_computed']):
            raise ValueError('All complete matching paired 64 tasks required')
        rows[arm].extend(local)
    return rows


@torch.inference_mode()
def readout(potential_job, decode_job):
    report = require_potentials(potential_job); draws, windows, _ = D.population()
    shards = closed_shards(decode_job, potential_job)
    saved = {a: D.A42.assemble_records(shards[a], windows) for a in ARMS}
    # Both arms' entire prediction files close before labels or old pools open.
    for arm in ARMS:
        with (OUT/(arm+'_natural_scores.jsonl')).open('x') as handle:
            for row in saved[arm]: handle.write(json.dumps(row, allow_nan=False)+'\n')
    old_draws, free, prior, _, old = D.A42.B.population()
    loaded_draws, _, target_windows, catalogs = R.Q.load_train_inputs()
    if (draws != old_draws or draws != loaded_draws or set(target_windows) != {r['window']['id'] for r in windows}
            or [r['window'] for r in free] != [r['window'] for r in windows]):
        raise ValueError('Original paired TRAIN-only target population differs')
    original = R.C.R.rows(R.READOUT/'windows.jsonl'); rows = []; began = time.perf_counter()
    if len(original) != 1491: raise ValueError('Complete original label population required')
    with (OUT/'windows.jsonl').open('x') as handle:
        for ordinal in range(1491):
            values = {a: load_potential(a, ordinal, report['records'][a][ordinal], report) for a in ARMS}
            if values[ARMS[0]]['sequence'] != values[ARMS[1]]['sequence']:
                raise ValueError('Paired original oriented DNA differs')
            w = target_windows[windows[ordinal]['window']['id']]; D.A42.validate_target_geometry(w, windows[ordinal]['window'])
            targets = R.Q.build_targets(w, catalogs[w['species']], values[ARMS[0]]['sequence'].upper())
            actual = label_paired_record({a: saved[a][ordinal] for a in ARMS}, prior[ordinal], original[ordinal], targets)
            rows.append(actual); handle.write(json.dumps(actual, allow_nan=False)+'\n')
            del values, targets
            if (ordinal+1) % 128 == 0 or ordinal+1 == 1491: print(json.dumps(dict(paired_labelled_windows=ordinal+1, total=1491)), flush=True)
    counts, refs = aggregate(rows); D.A42.check_baseline(counts, refs, old); checks, supported = gates(counts)
    R.C.R.write_summary(OUT, dict(completed=True, configuration=CONFIG, protocol=PROTOCOL,
        source_jobs=dict(fits=FIT_JOBS, potentials=potential_job, decode=decode_job),
        windows_per_arm=1491, original_draws=1536, original_candidates=16321, reference_injected_predictions=0,
        all_both_arm_natural_search_and_score_files_closed_before_labels=True, original_R4_labels_exact=True,
        new_geometry_labels_computed=True, selected_counts_per_species=counts, eligible_reference_pairs=refs,
        predeclared_TRAIN_gates=checks, TRAIN_preparation_supported=supported,
        model_or_backbone_forward=False, new_candidate_search=False, seconds=time.perf_counter()-began,
        limitation='Paired original TRAIN window-chain decisions, not raw assembly, independent generalization or novelty.', **SCOPE))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('preparation', 'potentials', 'decode', 'readout'), required=True)
    parser.add_argument('--checks-job'); parser.add_argument('--potential-job'); parser.add_argument('--decode-job')
    parser.add_argument('--task', type=int); args = parser.parse_args(); torch.set_num_threads(2)
    if args.mode == 'preparation': preparation()
    elif args.mode == 'potentials' and args.checks_job: potentials(args.checks_job)
    elif args.mode == 'decode' and args.task is not None and args.potential_job: decode(args.task, args.potential_job)
    elif args.mode == 'readout' and args.potential_job and args.decode_job: readout(args.potential_job, args.decode_job)
    else: parser.error('Actual source jobs and paired task required')
