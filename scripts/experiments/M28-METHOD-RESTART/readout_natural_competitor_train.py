#!/usr/bin/env python3
"""A43 last emissions -> label-free grammar searches -> frozen TRAIN labels."""
import argparse
from collections import Counter
import json
import resource
import time
import torch

import fit_natural_competitor as F
import run_emission_lattice_train as A42
from src.m28.labels import candidate_targets

R, ROOT = A42.R, F.ROOT
PROTOCOL = F.PROTOCOL
CONFIG = dict(source_checkpoint=str((F.OUT/'step_004608.pt').relative_to(ROOT)),
    source_arm='A43', frozen_source_parameters=366191, source_optimizer_steps=4608,
    source_fit_configuration=F.CONFIG, optimizer_steps_this_run=0, new_parameters=0,
    population='original_1491_unique_TRAIN_windows_first_occurrence_order',
    model_precision='FP32', score_accumulation='FP64', max_exons=128,
    threshold=0., null=0., all_role_weights=1.,
    search='all_motif_legal_oriented_window_emission_lattice_no_query_or_beam',
    prediction_candidate_pool_input=False, prediction_label_input=False,
    unique_TRAIN_windows=1491, original_draws=1536, shards=32,
    shard_assignment='original_ordinal_modulo_32')
CHECKS = ROOT/'outputs/M28-A43-NATURAL-COMPETITOR-READOUT-CHECKS'
POTENTIALS = ROOT/'outputs/M28-A43-NATURAL-COMPETITOR-TRAIN-POTENTIALS'
DECODE = ROOT/'outputs/M28-A43-NATURAL-COMPETITOR-TRAIN-DECODE'
OUT = ROOT/'outputs/M28-A43-NATURAL-COMPETITOR-TRAIN-READOUT'
CHECK = ROOT/'outputs/M28-A43-NATURAL-COMPETITOR-RESULT-CHECK'
FIT_JOB = '13424143'
ARMS = ('A31_A3', 'A43')


def ordered_prediction_inputs(draws, entries):
    """Original draw IDs and ten-field feature geometry, not supervised pools."""
    ids = list(dict.fromkeys(d['window_id'] for d in draws)); wanted = set(ids)
    local = [e for e in entries if e['window']['id'] in wanted]
    index = {e['window']['id']: e for e in local}
    if len(local) != len(ids) or set(index) != wanted:
        raise ValueError('Complete unique original prediction feature coverage required')
    for draw in draws:
        w = index[draw['window_id']]['window']
        if (set(w) != set(R.C.R.GEOMETRY) or w['split'] != 'train' or
                w['species'] != draw['species'] or w['seqid'] not in R.C.R.ALLOW[w['species']] or
                w['width'] != 24576 or not 0 < w['valid_bases'] <= w['width'] or
                w['end']-w['start'] != w['valid_bases'] or
                w['pad_right_after_orientation'] != w['width']-w['valid_bases']):
            raise ValueError('Original allowed label-free prediction window geometry differs')
    return [dict(window=index[i]['window']) for i in ids], index


def population():
    draws = R.C.R.rows(R.Q.DRAWS)
    windows, index = ordered_prediction_inputs(draws, R.C.R.rows(R.C.R.PLAN/'feature_plan.jsonl'))
    if len(draws) != 1536 or len(windows) != 1491:
        raise ValueError('Only unchanged original ordered TRAIN prediction population allowed')
    R.C.R.validate_feature_summary(json.loads((R.C.R.FEATURES/'summary.json').read_text()))
    return draws, windows, index


def validate_fit(value, state=False):
    if (value['arm'] != 'A43' or value['configuration'] != F.CONFIG or value['protocol'] != F.PROTOCOL or
            value['step' if state else 'optimizer_steps'] != 4608 or
            value['revision' if state else 'feature_revision'] != R.C.I.F.REVISION):
        raise ValueError('Only actual fixed A43 last004608 allowed')
    if not state and (not value['completed'] or value['source_support_job'] != '13424089' or
            value['preflight_job'] != '13424115' or value['final_checkpoint'] != CONFIG['source_checkpoint'] or
            value['trainable_parameters'] != 366191 or
            value['initialization'] != 'fresh_seed0_no_supervised_weight_inheritance' or
            any(value[k] for k in ('supervised_checkpoint_read', 'new_R4_labelling',
                'generator_updated', 'new_backbone_forward', 'DEV_used_for_selection',
                'test_setaria_used', 'independent_species_used', 'A33_outputs_opened'))):
        raise ValueError('Actual matching complete original TRAIN-only A43 fit required')


def require_fit():
    F.A.B.D.B.F.A.require_jobs([FIT_JOB])
    fitted = R.completed(F.OUT, FIT_JOB); validate_fit(fitted)
    return fitted


def preparation():
    draws, windows, index = population()
    old_draws, free, prior, old_index, old = A42.B.population()
    if (draws != old_draws or [r['window'] for r in free] != [r['window'] for r in windows] or
            any(index[i]['artifact'] != old_index[i]['artifact'] for i in index)):
        raise ValueError('Original natural baseline and A43 prediction inputs differ')
    saved_labels = R.C.R.rows(R.READOUT/'windows.jsonl')
    placeholders = [dict(window=r['window'], reference_used_in_scoring=False,
                         **r['variants']['A31_A3']) for r in prior]
    rows, counts, refs = A42.B.paired_readout(placeholders, prior, saved_labels)
    if (len(rows) != 1491 or sum(v['eligible'] for v in refs.values()) != 2916 or
            sum(v['five_or_more'] for v in refs.values()) != 1086 or
            any(counts[s]['A31_A3'] != old['selected_counts_per_species'][s]['A31_A3'] for s in counts)):
        raise ValueError('Original stronger baseline labels and complete reference denominators differ')
    R.C.R.write_summary(CHECKS, dict(completed=True, configuration=CONFIG, protocol=PROTOCOL,
        windows=1491, original_draws=1536, original_natural_candidates=16321,
        original_prediction_inputs_baseline_labels_and_denominators_reconciled=True,
        fitted_model_or_feature_tensors_opened=False, TRAIN_catalog_opened=False,
        optimizer_steps_this_run=0, model_or_backbone_forward=False, new_geometry_labels_computed=False,
        DEV_used=False, test_setaria_used=False, A33_outputs_opened=False,
        limitation='Readout adapter and original saved metadata readiness; no new predictions or performance.'))


def require_ready(checks_job):
    F.A.B.D.B.F.A.require_jobs([checks_job]); report = R.completed(CHECKS, checks_job)
    if (report['configuration'] != CONFIG or report['protocol'] != PROTOCOL or
            not report['original_prediction_inputs_baseline_labels_and_denominators_reconciled']):
        raise ValueError('Actual closed matching A43 readout readiness required')


def potential_path(ordinal):
    return POTENTIALS/('window_'+str(ordinal).zfill(6)+'.pt')


def shard_folder(shard):
    return DECODE/('shard_'+str(shard).zfill(2))


@torch.inference_mode()
def potentials(checks_job):
    require_ready(checks_job); fitted = require_fit()
    state = torch.load(ROOT/fitted['final_checkpoint'], map_location='cpu', weights_only=True)
    validate_fit(state, state=True)
    model = F.A.A36.new_model('cuda').eval().requires_grad_(False)
    model.load_state_dict(state['model'], strict=True); del state
    if sum(p.numel() for p in model.parameters()) != CONFIG['frozen_source_parameters']:
        raise ValueError('Fixed A43 parameter support differs')
    before = {k: v.detach().clone() for k, v in model.state_dict().items()}
    draws, windows, index = population()
    began = time.perf_counter(); receipts = []; torch.cuda.reset_peak_memory_stats()
    for ordinal, record in enumerate(windows):
        w = record['window']
        cached = torch.load(ROOT/index[w['id']]['artifact'], map_location='cpu', weights_only=True)
        features, dna = R.Q.window_inputs(cached, w, 'cuda')
        local = model(features, dna, [w['valid_bases']])['segmentation_logits'][0][:w['valid_bases']]
        if local.dtype != torch.float32 or tuple(local.shape) != (w['valid_bases'], 15) or not bool(torch.isfinite(local).all()):
            raise ValueError('Finite fixed last-model valid-base FP32 emission required')
        destination = potential_path(ordinal)
        with destination.open('xb') as handle:
            torch.save(dict(configuration=CONFIG, protocol=PROTOCOL, original_ordinal=ordinal,
                fit_job=FIT_JOB, window=w, sequence=cached['sequence'], local_logits=local.cpu(),
                revision=R.C.I.F.REVISION, reference_used_in_scoring=False), handle)
        receipts.append(dict(original_ordinal=ordinal, window=w, artifact=str(destination.relative_to(ROOT))))
        if (ordinal+1) % 64 == 0 or ordinal+1 == len(windows):
            print(json.dumps(dict(potential_windows=ordinal+1, total=1491)), flush=True)
        del cached, features, dna, local
    if (any(not torch.equal(v, before[k]) for k, v in model.state_dict().items()) or
            any(p.requires_grad or p.grad is not None for p in model.parameters())):
        raise ValueError('Frozen A43 last-model state changed or received gradient')
    R.C.R.write_summary(POTENTIALS, dict(completed=True, configuration=CONFIG, protocol=PROTOCOL,
        source_fit_job=FIT_JOB, source_checks_job=checks_job, checkpoint=fitted['final_checkpoint'],
        windows=1491, original_draws=1536, records=receipts, feature_revision=R.C.I.F.REVISION,
        all_potential_files_closed=True, frozen_source_weights_unchanged=True,
        all_source_parameters_frozen_and_gradients_absent=True,
        prediction_candidate_pool_opened=False, reference_labels_opened=False, TRAIN_catalog_opened=False,
        optimizer_steps_this_run=0, checkpoints_written=False, query_decoding=False, new_candidate_search=False,
        new_backbone_forward=False, new_geometry_labels_computed=False, R4_policy_or_catalog_changed=False,
        DEV_used=False, test_setaria_used=False, independent_species_used=False, A33_outputs_opened=False,
        seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='Complete frozen A43 emissions only; no search, labels, raw annotation or scientific result.'))


def require_potentials(job):
    F.A.B.D.B.F.A.require_jobs([job]); report = R.completed(POTENTIALS, job)
    if (report['configuration'] != CONFIG or report['protocol'] != PROTOCOL or
            report['source_fit_job'] != FIT_JOB or report['checkpoint'] != CONFIG['source_checkpoint'] or
            report['windows'] != 1491 or len(report['records']) != 1491 or
            not report['all_potential_files_closed'] or not report['frozen_source_weights_unchanged'] or
            not report['all_source_parameters_frozen_and_gradients_absent'] or
            report['prediction_candidate_pool_opened'] or report['reference_labels_opened'] or
            report['TRAIN_catalog_opened'] or report['optimizer_steps_this_run'] != 0):
        raise ValueError('Complete matching frozen label-free A43 potentials required')
    return report


def load_potential(ordinal, receipt, report):
    destination = potential_path(ordinal)
    value = torch.load(destination, map_location='cpu', weights_only=True)
    w = value['window']; local = value['local_logits']
    if (receipt['artifact'] != str(destination.relative_to(ROOT)) or receipt['original_ordinal'] != ordinal or
            value['original_ordinal'] != ordinal or value['configuration'] != CONFIG or value['protocol'] != PROTOCOL or
            value['fit_job'] != FIT_JOB or w != receipt['window'] or value['revision'] != report['feature_revision'] or
            value['reference_used_in_scoring'] or w['split'] != 'train' or
            w['seqid'] not in R.C.R.ALLOW[w['species']] or len(value['sequence']) != w['valid_bases'] or
            local.dtype != torch.float32 or tuple(local.shape) != (w['valid_bases'], 15) or
            not bool(torch.isfinite(local).all())):
        raise ValueError('Matching finite original A43 prediction potential required')
    return value


def validate_prediction(value, natural):
    """No training gold, labels or old candidate geometry is an input."""
    sequence = value['sequence'].upper(); previous_end = 0
    for chain in natural['chains']:
        A42.P.check_chain(sequence, chain)
        if chain[0][0] < previous_end:
            raise ValueError('Natural complete-gene spans overlap')
        previous_end = chain[-1][1]
    direct = A42.P.local_path_scores(value['local_logits'], natural['chains'])
    torch.testing.assert_close(direct, torch.tensor(natural['chain_gains'], dtype=torch.float64), rtol=1e-12, atol=1e-8)
    torch.testing.assert_close(direct.sum(), torch.tensor(natural['gain'], dtype=torch.float64), rtol=1e-12, atol=1e-8)
    if natural['reference_used'] or any(s <= 0 for s in natural['chain_gains']):
        raise ValueError('Reference-free natural positive-gain chains required')


@torch.inference_mode()
def decode(shard, potential_job):
    if not 0 <= shard < CONFIG['shards']:
        raise ValueError('Only fixed 32 original-ordinal shards allowed')
    report = require_potentials(potential_job); out = shard_folder(shard)
    ordinals = list(range(shard, 1491, CONFIG['shards'])); began = time.perf_counter()
    with (out/'started.json').open('x') as handle:
        json.dump(dict(configuration=CONFIG, protocol=PROTOCOL, shard=shard,
                       original_ordinals=ordinals, potential_job=potential_job), handle)
    with (out/'natural_scores.jsonl').open('x') as handle:
        for ordinal in ordinals:
            value = load_potential(ordinal, report['records'][ordinal], report); started = time.perf_counter()
            natural = A42.P.decode_emission_lattice(value['sequence'].upper(), value['local_logits'], CONFIG['max_exons'])
            validate_prediction(value, natural); seconds = time.perf_counter()-started
            row = dict(original_ordinal=ordinal, shard=shard, window=value['window'],
                chains=natural['chains'], scores=natural['chain_gains'], selected=list(range(len(natural['chains']))),
                gain=natural['gain'], potential_artifact=report['records'][ordinal]['artifact'],
                feature_revision=value['revision'], potential_job=potential_job, source_fit_job=FIT_JOB,
                reference_used_in_scoring=False, decode_seconds=seconds,
                search_counters={k: v for k, v in natural.items() if k not in ('chains', 'chain_gains', 'gain', 'reference_used')})
            handle.write(json.dumps(row, allow_nan=False)+'\n'); handle.flush()
            print(json.dumps(dict(shard=shard, original_ordinal=ordinal, decode_seconds=seconds)), flush=True)
    R.C.R.write_summary(out, dict(completed=True, configuration=CONFIG, protocol=PROTOCOL,
        shard=shard, potential_job=potential_job, source_fit_job=FIT_JOB, original_ordinals=ordinals, windows=len(ordinals),
        natural_score_file_closed=True, all_selected_ORF_motif_cap_span_and_direct_role_gains_verified=True,
        prediction_candidate_pool_opened=False, reference_labels_opened=False, TRAIN_catalog_opened=False,
        feature_revision=report['feature_revision'], optimizer_steps_this_run=0, model_or_backbone_forward=False,
        new_candidate_search=True, new_geometry_labels_computed=False, R4_policy_or_catalog_changed=False,
        DEV_used=False, test_setaria_used=False, independent_species_used=False, A33_outputs_opened=False,
        seconds=time.perf_counter()-began, process_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        limitation='Complete assigned label-free grammar predictions, not performance or global-optimality replay.'))


def closed_shards(decode_job, potential_job):
    F.A.B.D.B.F.A.require_jobs([decode_job+'_'+str(i) for i in range(CONFIG['shards'])]); rows = []
    for i in range(CONFIG['shards']):
        folder = shard_folder(i); report = R.completed(folder, decode_job+'_'+str(i))
        expected = list(range(i, 1491, CONFIG['shards'])); local = R.C.R.rows(folder/'natural_scores.jsonl')
        if (report['configuration'] != CONFIG or report['protocol'] != PROTOCOL or
                report['shard'] != i or report['potential_job'] != potential_job or report['source_fit_job'] != FIT_JOB or
                report['original_ordinals'] != expected or report['windows'] != len(expected) or
                [r['original_ordinal'] for r in local] != expected or not report['natural_score_file_closed'] or
                not report['all_selected_ORF_motif_cap_span_and_direct_role_gains_verified'] or
                report['prediction_candidate_pool_opened'] or report['reference_labels_opened'] or
                report['TRAIN_catalog_opened'] or report['new_geometry_labels_computed']):
            raise ValueError('All actual complete matching label-free A43 shards required')
        rows.extend(local)
    return rows


def label_record(row, old, original, targets):
    if row['window'] != old['window'] or row['window'] != original['window']:
        raise ValueError('Original paired R4 label window differs')
    baseline = old['variants']['A31_A3']
    old_labels = candidate_targets(baseline['chains'], targets); labels = candidate_targets(row['chains'], targets)
    old_counts = R.selected_counts(baseline, old_labels)
    if old_labels != original['labels']['A31_A3'] or old_counts != original['counts']['A31_A3']:
        raise ValueError('Original frozen R4 baseline labels or selected counts differ')
    eligible = len(targets['positive_chains']); long = sum(len(c) >= 5 for c in targets['positive_chains'])
    if eligible != original['eligible_references'] or long != original['eligible_five_or_more']:
        raise ValueError('Original complete reference denominators differ')
    return dict(window=row['window'], labels={'A31_A3': old_labels, 'A43': labels},
        counts={'A31_A3': A42.selected_only(old_counts), 'A43': A42.selected_only(R.selected_counts(row, labels))},
        eligible_references=eligible, eligible_five_or_more=long)


def aggregate(rows):
    totals = {s: {a: Counter() for a in ARMS} for s in R.C.R.ALLOW}
    refs = {s: Counter(eligible=0, five_or_more=0) for s in R.C.R.ALLOW}
    for row in rows:
        s = row['window']['species']
        for arm in ARMS: totals[s][arm].update(row['counts'][arm])
        refs[s].update(eligible=row['eligible_references'], five_or_more=row['eligible_five_or_more'])
    return ({s: {a: dict(v) for a, v in arms.items()} for s, arms in totals.items()}, {s: dict(v) for s, v in refs.items()})


def gate(counts):
    checks = {}
    for s, values in counts.items():
        checks[s+'_natural_selected_correct_strictly_increased'] = values['A43']['selected_correct'] > values['A31_A3']['selected_correct']
        checks[s+'_natural_selected_known_nonmatch_not_increased'] = values['A43']['selected_known_nonmatch'] <= values['A31_A3']['selected_known_nonmatch']
    checks['pooled_natural_selected_five_or_more_correct_strictly_increased'] = (
        sum(v['A43']['selected_five_or_more_correct'] for v in counts.values()) >
        sum(v['A31_A3']['selected_five_or_more_correct'] for v in counts.values()))
    return checks


@torch.inference_mode()
def readout(potential_job, decode_job):
    report = require_potentials(potential_job); draws, windows, _ = population()
    saved = A42.assemble_records(closed_shards(decode_job, potential_job), windows)
    with (OUT/'natural_scores.jsonl').open('x') as handle:
        for row in saved: handle.write(json.dumps(row, allow_nan=False)+'\n')
    # Only after every search and the merged prediction file have closed.
    old_draws, free, prior, _, old = A42.B.population()
    loaded_draws, _, target_windows, catalogs = R.Q.load_train_inputs()
    if draws != old_draws or draws != loaded_draws or set(target_windows) != {r['window']['id'] for r in windows}:
        raise ValueError('Original TRAIN-only target population differs')
    original_labels = R.C.R.rows(R.READOUT/'windows.jsonl'); rows = []; began = time.perf_counter()
    if len(original_labels) != 1491 or [r['window'] for r in free] != [r['window'] for r in windows]:
        raise ValueError('Complete original saved label/baseline windows required')
    with (OUT/'windows.jsonl').open('x') as handle:
        for ordinal, row in enumerate(saved):
            value = load_potential(ordinal, report['records'][ordinal], report)
            w = target_windows[row['window']['id']]; A42.validate_target_geometry(w, row['window'])
            targets = R.Q.build_targets(w, catalogs[w['species']], value['sequence'].upper())
            result = label_record(row, prior[ordinal], original_labels[ordinal], targets)
            rows.append(result); handle.write(json.dumps(result, allow_nan=False)+'\n')
            if (ordinal+1) % 128 == 0 or ordinal+1 == 1491:
                print(json.dumps(dict(labelled_windows=ordinal+1, total=1491)), flush=True)
            del value, targets
    counts, refs = aggregate(rows); A42.check_baseline(counts, refs, old); checks = gate(counts)
    R.C.R.write_summary(OUT, dict(completed=True, configuration=CONFIG, protocol=PROTOCOL,
        source_jobs=dict(fit=FIT_JOB, potentials=potential_job, decode=decode_job),
        windows=1491, original_draws=1536, original_candidates=16321, reference_injected_predictions=0,
        all_natural_search_and_score_files_closed_before_labels=True, original_R4_labels_exact=True,
        new_geometry_labels_computed=True, R4_policy_or_catalog_changed=False,
        selected_counts_per_species=counts, eligible_reference_pairs=refs,
        predeclared_TRAIN_gate=checks, TRAIN_preparation_supported=all(checks.values()),
        optimizer_steps_this_run=0, model_or_backbone_forward=False, new_candidate_search=False,
        DEV_used=False, test_setaria_used=False, independent_species_used=False, A33_outputs_opened=False,
        seconds=time.perf_counter()-began,
        limitation='Natural full-grammar decisions on original TRAIN windows; not raw assembly, independent generalization or causal ablation.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('preparation', 'potentials', 'decode', 'readout'), required=True)
    parser.add_argument('--checks-job'); parser.add_argument('--potential-job'); parser.add_argument('--decode-job')
    parser.add_argument('--shard', type=int); args = parser.parse_args(); torch.set_num_threads(2)
    if args.mode == 'preparation': preparation()
    elif args.mode == 'potentials' and args.checks_job: potentials(args.checks_job)
    elif args.mode == 'decode' and args.shard is not None and args.potential_job: decode(args.shard, args.potential_job)
    elif args.mode == 'readout' and args.potential_job and args.decode_job: readout(args.potential_job, args.decode_job)
    else: parser.error('Actual source jobs and shard required for chosen stage')
