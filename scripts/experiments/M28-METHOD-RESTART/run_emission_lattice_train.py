#!/usr/bin/env python3
"""A42 full original TRAIN lattice: frozen potentials, sharded search, then labels."""
import argparse
from collections import Counter
import json
import resource
import time

import torch

import preflight_emission_lattice_train as P
from src.m28.labels import candidate_targets

R, B, ROOT = P.R, P.B, P.ROOT
PROTOCOL = P.PROTOCOL
CONFIG = {k: v for k, v in P.CONFIG.items() if k not in ('cases', 'cases_count')}
CONFIG.update(unique_TRAIN_windows=1491, original_draws=1536,
              shards=32, shard_assignment='original_ordinal_modulo_32')
POTENTIALS = ROOT/'outputs/M28-A42-EMISSION-LATTICE-TRAIN-POTENTIALS'
DECODE = ROOT/'outputs/M28-A42-EMISSION-LATTICE-TRAIN-DECODE'
OUT = ROOT/'outputs/M28-A42-EMISSION-LATTICE-TRAIN-READOUT'
CHECKS = ROOT/'outputs/M28-A42-EMISSION-LATTICE-TRAIN-CHECKS'
CHECK = ROOT/'outputs/M28-A42-EMISSION-LATTICE-RESULT-CHECK'
ARMS = ('A31_A3', 'A42')
PREFLIGHT_JOBS = ['13421592', '13421848', '13421953'] + ['13421967_'+str(i) for i in range(8)]


def shard_folder(shard):
    return DECODE/('shard_'+str(shard).zfill(2))


def potential_path(ordinal):
    return POTENTIALS/('window_'+str(ordinal).zfill(6)+'.pt')


def selected_only(counts):
    return {k: v for k, v in counts.items() if k.startswith('selected')}


def validate_target_geometry(target_window, natural_window):
    # R4 windows additionally carry label metadata, never search inputs.
    if {k: target_window[k] for k in R.C.R.GEOMETRY} != natural_window:
        raise ValueError('Original TRAIN target and natural window geometry differ')


def assemble_records(records, free):
    """Restore original order; incomplete or duplicated shards must not form a result."""
    by_ordinal = {}
    for row in records:
        ordinal = row['original_ordinal']
        if (not isinstance(ordinal, int) or not 0 <= ordinal < len(free) or
                ordinal in by_ordinal or row['window'] != free[ordinal]['window'] or
                row['reference_used_in_scoring'] or row['shard'] != ordinal % CONFIG['shards']):
            raise ValueError('Unique corresponding reference-free original window/shard required')
        by_ordinal[ordinal] = row
    if set(by_ordinal) != set(range(len(free))):
        raise ValueError('Complete original window population required')
    return [by_ordinal[i] for i in range(len(free))]


def label_record(row, old, saved_labels, targets):
    if row['window'] != old['window'] or row['window'] != saved_labels['window']:
        raise ValueError('Original paired label window identity differs')
    baseline = old['variants']['A31_A3']
    old_labels = candidate_targets(baseline['chains'], targets)
    new_labels = candidate_targets(row['chains'], targets)
    full_old = R.selected_counts(baseline, old_labels)
    if old_labels != saved_labels['labels']['A31_A3'] or full_old != saved_labels['counts']['A31_A3']:
        raise ValueError('Original frozen R4 labels or baseline counts differ')
    refs = dict(eligible=len(targets['positive_chains']),
                five_or_more=sum(len(c) >= 5 for c in targets['positive_chains']))
    if (refs['eligible'] != saved_labels['eligible_references'] or
            refs['five_or_more'] != saved_labels['eligible_five_or_more']):
        raise ValueError('Original full reference denominators differ')
    return dict(window=row['window'], labels={'A31_A3': old_labels, 'A42': new_labels},
        counts={'A31_A3': selected_only(full_old),
                'A42': selected_only(R.selected_counts(row, new_labels))},
        eligible_references=refs['eligible'], eligible_five_or_more=refs['five_or_more'])


def gate(counts):
    checks = {}
    for species, values in counts.items():
        new, old = values['A42'], values['A31_A3']
        checks[species+'_natural_selected_correct_strictly_increased'] = new['selected_correct'] > old['selected_correct']
        checks[species+'_natural_selected_known_nonmatch_not_increased'] = new['selected_known_nonmatch'] <= old['selected_known_nonmatch']
    checks['pooled_natural_selected_five_or_more_correct_strictly_increased'] = (
        sum(v['A42']['selected_five_or_more_correct'] for v in counts.values()) >
        sum(v['A31_A3']['selected_five_or_more_correct'] for v in counts.values()))
    return checks


def aggregate(rows):
    totals = {s: {a: Counter() for a in ARMS} for s in R.C.R.ALLOW}
    references = {s: Counter(eligible=0, five_or_more=0) for s in R.C.R.ALLOW}
    for row in rows:
        species = row['window']['species']
        for arm in ARMS:
            totals[species][arm].update(row['counts'][arm])
        references[species].update(eligible=row['eligible_references'], five_or_more=row['eligible_five_or_more'])
    return ({s: {a: dict(c) for a, c in v.items()} for s, v in totals.items()},
            {s: dict(c) for s, c in references.items()})


def check_baseline(counts, references, old):
    if (sum(v['eligible'] for v in references.values()) != 2916 or
            sum(v['five_or_more'] for v in references.values()) != 1086 or
            any(counts[s]['A31_A3'] != selected_only(old['selected_counts_per_species'][s]['A31_A3'])
                for s in counts)):
        raise ValueError('Original full denominators or stronger baseline drifted')


def preparation():
    B.F.A.require_jobs(PREFLIGHT_JOBS)
    source = R.completed(P.POTENTIALS, '13421953')
    if source['configuration'] != P.CONFIG or not source['all_potential_files_closed']:
        raise ValueError('Actual unchanged eight-window preflight required')
    for i in range(8):
        report = R.completed(P.DECODE/('case_'+str(i)), '13421967_'+str(i))
        if (report['configuration'] != P.CONFIG or report['case'] != i or
                not report['selected_ORF_motif_cap_span_and_direct_role_gains_verified'] or
                not report['full_grammar_objective_not_lower_than_baseline']):
            raise ValueError('All actual eight decode preflight cases required')
    draws, free, prior, index, old = B.population()
    labelled = R.C.R.rows(R.READOUT/'windows.jsonl')
    placeholders = [dict(window=r['window'], reference_used_in_scoring=False,
                         **r['variants']['A31_A3']) for r in prior]
    # Existing original labels only; no catalog, tensors or new A42 geometry here.
    original_rows, counts, refs = B.paired_readout(placeholders, prior, labelled)
    if (len(original_rows) != 1491 or sum(v['five_or_more'] for v in refs.values()) != 1086 or
            sum(v['eligible'] for v in refs.values()) != 2916 or
            any(counts[s]['A31_A3'] != old['selected_counts_per_species'][s]['A31_A3'] for s in counts)):
        raise ValueError('Actual full original baseline metadata required')
    R.C.R.write_summary(CHECKS, dict(completed=True, configuration=CONFIG, protocol=PROTOCOL,
        source_preflight_jobs=PREFLIGHT_JOBS, windows=len(free), original_draws=len(draws),
        original_candidates=sum(len(p['chains']) for p in free), feature_index_records=len(index),
        original_labels_baseline_and_full_denominators_reconciled=True,
        eligible_reference_pairs=refs, optimizer_steps_this_run=0, feature_tensors_opened=False,
        model_or_backbone_forward=False, TRAIN_catalog_opened=False, new_geometry_labels_computed=False,
        DEV_used=False, test_setaria_used=False, A33_outputs_opened=False,
        limitation='Original metadata and new adapter checks only; no A42 search or scientific result.'))


def require_ready(job):
    B.F.A.require_jobs([job])
    value = R.completed(CHECKS, job)
    if (value['configuration'] != CONFIG or value['protocol'] != PROTOCOL or
            not value['original_labels_baseline_and_full_denominators_reconciled']):
        raise ValueError('Actual unchanged full TRAIN readiness required')


def require_potentials(job):
    B.F.A.require_jobs([job])
    report = R.completed(POTENTIALS, job)
    if (report['configuration'] != CONFIG or report['protocol'] != PROTOCOL or
            report['windows'] != 1491 or len(report['records']) != 1491 or
            not report['all_potential_files_closed'] or not report['frozen_source_weights_unchanged'] or
            not report['all_source_parameters_frozen_and_gradients_absent'] or
            report['reference_labels_opened'] or report['TRAIN_catalog_opened'] or
            report['optimizer_steps_this_run'] != 0):
        raise ValueError('Complete frozen label-free 1491-window potentials required')
    return report


def load_potential(ordinal, receipt, report):
    path = potential_path(ordinal)
    value = torch.load(path, map_location='cpu', weights_only=True)
    w = value['window']; sequence = value['sequence'].upper(); local = value['local_logits']
    if (receipt['artifact'] != str(path.relative_to(ROOT)) or
            receipt['original_ordinal'] != ordinal or value['original_ordinal'] != ordinal or
            value['configuration'] != CONFIG or value['protocol'] != PROTOCOL or
            w != receipt['window'] or value['revision'] != report['feature_revision'] or
            value['reference_used_in_scoring'] or w['split'] != 'train' or
            w['seqid'] not in R.C.R.ALLOW[w['species']] or len(sequence) != w['valid_bases'] or
            local.dtype != torch.float32 or tuple(local.shape) != (len(sequence), 15) or
            not bool(torch.isfinite(local).all())):
        raise ValueError('Matching finite original TRAIN potential identity required')
    return value


@torch.inference_mode()
def potentials(checks_job):
    require_ready(checks_job)
    draws, free, prior, index, _ = B.population()
    model = R.C.load_generator().eval().requires_grad_(False)
    if sum(p.numel() for p in model.parameters()) != CONFIG['frozen_source_parameters']:
        raise ValueError('Original frozen source parameter population differs')
    before = {k: v.detach().clone() for k, v in model.state_dict().items()}
    began = time.perf_counter(); receipts = []; torch.cuda.reset_peak_memory_stats()
    for ordinal, pool in enumerate(free):
        w = pool['window']; old = prior[ordinal]['variants']['A31_A3']
        cached = torch.load(ROOT/index[w['id']]['artifact'], map_location='cpu', weights_only=True)
        features, dna = R.Q.window_inputs(cached, w, 'cuda')
        hidden, _ = model.shared(features, dna, [w['valid_bases']])
        local = model.segmentation(hidden[0])[:w['valid_bases']]
        if local.dtype != torch.float32 or tuple(local.shape) != (w['valid_bases'], 15) or not bool(torch.isfinite(local).all()):
            raise ValueError('Finite frozen valid-base FP32 emissions required')
        scores = P.local_path_scores(local, pool['chains']).cpu()
        expected = torch.tensor(old['scores'], dtype=torch.float64)
        torch.testing.assert_close(scores, expected, rtol=1e-12, atol=1e-8)
        if P.select_nonoverlapping(pool['chains'], scores.tolist()) != old['selected']:
            raise ValueError('Original stronger baseline selection differs')
        for chain in pool['chains']:
            P.check_chain(cached['sequence'].upper(), chain)
        path = potential_path(ordinal)
        with path.open('xb') as handle:
            torch.save(dict(configuration=CONFIG, protocol=PROTOCOL, original_ordinal=ordinal,
                window=w, sequence=cached['sequence'], local_logits=local.cpu(),
                original_chains=pool['chains'], original_scores=old['scores'], original_selected=old['selected'],
                revision=R.C.I.F.REVISION, reference_used_in_scoring=False), handle)
        receipts.append(dict(original_ordinal=ordinal, window=w, artifact=str(path.relative_to(ROOT)),
            original_candidates=len(pool['chains']), original_A31_score_max_error=
            float((scores - expected).abs().max()) if len(scores) else 0.,
            original_selected_indices_exact=True, all_original_pool_chains_grammar_supported=True))
        if (ordinal+1) % 64 == 0 or ordinal+1 == len(free):
            print(json.dumps(dict(potential_windows=ordinal+1, total=1491)), flush=True)
        del cached, features, dna, hidden, local, scores, expected
    if (any(not torch.equal(v, before[k]) for k, v in model.state_dict().items()) or
            any(p.requires_grad or p.grad is not None for p in model.parameters())):
        raise ValueError('Frozen source changed or received gradients')
    R.C.R.write_summary(POTENTIALS, dict(completed=True, configuration=CONFIG, protocol=PROTOCOL,
        source_checks_job=checks_job, windows=len(free), original_draws=len(draws),
        original_natural_candidates=sum(len(p['chains']) for p in free), records=receipts,
        frozen_source_weights_unchanged=True, all_source_parameters_frozen_and_gradients_absent=True,
        all_potential_files_closed=True, feature_revision=R.C.I.F.REVISION,
        optimizer_steps_this_run=0, checkpoints_written=False, query_decoding=False, generator_updated=False,
        new_candidate_search=False, new_backbone_forward=False, reference_labels_opened=False,
        TRAIN_catalog_opened=False, new_geometry_labels_computed=False, R4_policy_or_catalog_changed=False,
        DEV_used=False, test_setaria_used=False, A33_outputs_opened=False,
        seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='All original TRAIN frozen emissions; no new search or scientific result yet.'))


def validate_prediction(value, natural):
    sequence = value['sequence'].upper(); local = value['local_logits']
    baseline = P.local_path_scores(local, value['original_chains'])
    torch.testing.assert_close(baseline, torch.tensor(value['original_scores'], dtype=torch.float64), rtol=1e-12, atol=1e-8)
    if P.select_nonoverlapping(value['original_chains'], baseline.tolist()) != value['original_selected']:
        raise ValueError('Actual original baseline selection differs')
    for chain in value['original_chains']:
        P.check_chain(sequence, chain)
    previous_end = 0
    for chain in natural['chains']:
        P.check_chain(sequence, chain)
        if chain[0][0] < previous_end:
            raise ValueError('Selected complete-gene spans overlap')
        previous_end = chain[-1][1]
    direct = P.local_path_scores(local, natural['chains'])
    torch.testing.assert_close(direct, torch.tensor(natural['chain_gains'], dtype=torch.float64), rtol=1e-12, atol=1e-8)
    torch.testing.assert_close(direct.sum(), torch.tensor(natural['gain'], dtype=torch.float64), rtol=1e-12, atol=1e-8)
    old_gain = float(baseline[value['original_selected']].sum())
    if (natural['reference_used'] or any(s <= 0 for s in natural['chain_gains']) or
            natural['gain'] + 1e-8 + 1e-12*abs(old_gain) < old_gain):
        raise ValueError('Reference-free full grammar must attain original supported pool objective')
    return old_gain


@torch.inference_mode()
def decode(shard, potential_job):
    if not 0 <= shard < CONFIG['shards']:
        raise ValueError('Only original fixed 32 shards permitted')
    report = require_potentials(potential_job)
    ordinals = list(range(shard, 1491, CONFIG['shards']))
    out = shard_folder(shard); began = time.perf_counter()
    with (out/'started.json').open('x') as handle:
        json.dump(dict(configuration=CONFIG, protocol=PROTOCOL, shard=shard,
                       original_ordinals=ordinals, potential_job=potential_job), handle)
    with (out/'natural_scores.jsonl').open('x') as handle:
        for ordinal in ordinals:
            value = load_potential(ordinal, report['records'][ordinal], report)
            started = time.perf_counter()
            natural = P.decode_emission_lattice(value['sequence'].upper(), value['local_logits'], CONFIG['max_exons'])
            seconds = time.perf_counter()-started
            old_gain = validate_prediction(value, natural)
            row = dict(original_ordinal=ordinal, shard=shard, window=value['window'],
                chains=natural['chains'], scores=natural['chain_gains'], selected=list(range(len(natural['chains']))),
                gain=natural['gain'], original_baseline_gain=old_gain,
                potential_artifact=report['records'][ordinal]['artifact'], feature_revision=value['revision'],
                potential_job=potential_job, reference_used_in_scoring=False,
                decode_seconds=seconds, search_counters={k: v for k, v in natural.items()
                    if k not in ('chains', 'chain_gains', 'gain', 'reference_used')})
            handle.write(json.dumps(row, allow_nan=False)+'\n'); handle.flush()
            print(json.dumps(dict(shard=shard, original_ordinal=ordinal, decode_seconds=seconds)), flush=True)
            del value, natural
    R.C.R.write_summary(out, dict(completed=True, configuration=CONFIG, protocol=PROTOCOL, shard=shard,
        potential_job=potential_job, original_ordinals=ordinals, windows=len(ordinals),
        natural_score_file_closed=True, all_original_baseline_scores_and_indices_reproduced=True,
        all_selected_ORF_motif_cap_span_and_direct_role_gains_verified=True,
        all_full_grammar_objectives_not_lower_than_baseline=True, feature_revision=report['feature_revision'],
        optimizer_steps_this_run=0, model_or_backbone_forward=False, new_candidate_search=True,
        reference_labels_opened=False, TRAIN_catalog_opened=False, new_geometry_labels_computed=False,
        R4_policy_or_catalog_changed=False, DEV_used=False, test_setaria_used=False, A33_outputs_opened=False,
        seconds=time.perf_counter()-began, process_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        limitation='Complete assigned natural TRAIN searches before labels; no scientific result yet.'))


def closed_shards(decode_job, potential_job):
    B.F.A.require_jobs([decode_job+'_'+str(i) for i in range(CONFIG['shards'])])
    rows = []
    for i in range(CONFIG['shards']):
        folder = shard_folder(i); report = R.completed(folder, decode_job+'_'+str(i))
        expected = list(range(i, 1491, CONFIG['shards']))
        shard_rows = R.C.R.rows(folder/'natural_scores.jsonl')
        if (report['configuration'] != CONFIG or report['protocol'] != PROTOCOL or
                report['shard'] != i or report['potential_job'] != potential_job or
                report['original_ordinals'] != expected or report['windows'] != len(expected) or
                [r['original_ordinal'] for r in shard_rows] != expected or
                not report['natural_score_file_closed'] or report['reference_labels_opened'] or
                report['TRAIN_catalog_opened'] or report['new_geometry_labels_computed'] or
                not report['all_selected_ORF_motif_cap_span_and_direct_role_gains_verified'] or
                not report['all_full_grammar_objectives_not_lower_than_baseline']):
            raise ValueError('All corresponding complete fixed label-free shards required')
        rows.extend(shard_rows)
    return rows


@torch.inference_mode()
def readout(potential_job, decode_job):
    report = require_potentials(potential_job)
    _, free, prior, _, old = B.population()
    saved = assemble_records(closed_shards(decode_job, potential_job), free)
    # All 1491 natural searches and merged predictions are closed before catalog/labels.
    with (OUT/'natural_scores.jsonl').open('x') as handle:
        for row in saved:
            handle.write(json.dumps(row, allow_nan=False)+'\n')
    draws, _, windows, catalogs = R.Q.load_train_inputs()
    if len(draws) != 1536 or set(windows) != {p['window']['id'] for p in free}:
        raise ValueError('Original TRAIN-only target population differs')
    labelled = R.C.R.rows(R.READOUT/'windows.jsonl')
    if len(labelled) != 1491:
        raise ValueError('Complete original saved labels required')
    rows = []; began = time.perf_counter()
    with (OUT/'windows.jsonl').open('x') as handle:
        for ordinal, row in enumerate(saved):
            value = load_potential(ordinal, report['records'][ordinal], report)
            w = windows[row['window']['id']]
            validate_target_geometry(w, row['window'])
            if value['window'] != row['window']:
                raise ValueError('Original TRAIN target and natural window geometry differ')
            targets = R.Q.build_targets(w, catalogs[w['species']], value['sequence'].upper())
            result = label_record(row, prior[ordinal], labelled[ordinal], targets)
            rows.append(result); handle.write(json.dumps(result, allow_nan=False)+'\n')
            del value, targets
            if (ordinal+1) % 128 == 0 or ordinal+1 == len(saved):
                print(json.dumps(dict(labelled_windows=ordinal+1, total=1491)), flush=True)
    counts, references = aggregate(rows); check_baseline(counts, references, old); checks = gate(counts)
    R.C.R.write_summary(OUT, dict(completed=True, configuration=CONFIG, protocol=PROTOCOL,
        source_jobs=dict(potentials=potential_job, decode=decode_job), windows=1491,
        original_draws=1536, original_candidates=16321, selected_counts_per_species=counts,
        eligible_reference_pairs=references, predeclared_TRAIN_gate=checks, TRAIN_preparation_supported=all(checks.values()),
        all_natural_search_and_score_files_closed_before_labels=True, original_R4_labels_exact=True,
        new_geometry_labels_computed=True, R4_policy_or_catalog_changed=False,
        feature_revision=report['feature_revision'], optimizer_steps_this_run=0, model_or_backbone_forward=False,
        new_candidate_search=False, reference_injected_predictions=0, DEV_used=False,
        test_setaria_used=False, A33_outputs_opened=False, seconds=time.perf_counter()-began,
        limitation='Selected gene sets on original TRAIN windows, not full candidate counts, raw assembly or generalization.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('preparation', 'potentials', 'decode', 'readout'), required=True)
    parser.add_argument('--checks-job'); parser.add_argument('--potential-job')
    parser.add_argument('--decode-job'); parser.add_argument('--shard', type=int)
    args = parser.parse_args(); torch.set_num_threads(2)
    if args.mode == 'preparation':
        preparation()
    elif args.mode == 'potentials' and args.checks_job:
        potentials(args.checks_job)
    elif args.mode == 'decode' and args.shard is not None and args.potential_job:
        decode(args.shard, args.potential_job)
    elif args.mode == 'readout' and args.potential_job and args.decode_job:
        readout(args.potential_job, args.decode_job)
    else:
        parser.error('Actual fixed source jobs and shard required for chosen stage')
