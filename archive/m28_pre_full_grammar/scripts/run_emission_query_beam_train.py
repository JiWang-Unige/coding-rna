#!/usr/bin/env python3
"""A33 frozen emission search: four natural shards, labels, independent CPU check."""
import argparse
from collections import Counter
import json
import os
import subprocess
import time
import torch
import run_local_path_caller as R
from preflight_emission_query_beam import CONFIG, PROTOCOL, OUT as PREFLIGHT
from query_search_score_replay import replay_node_score, STOP
from src.m28.candidates import Budget, ORFIndex
from src.m28.query_decode import _end_events, _legal_ends
from src.m28.emission_query_beam import decode_emission_query_beam, augment_pool

ROOT, C, Q = R.ROOT, R.C, R.Q
NATURAL = ROOT/'outputs/M28-A33-EMISSION-BEAM-TRAIN-NATURAL'
READOUT = ROOT/'outputs/M28-A33-EMISSION-BEAM-TRAIN-READOUT'
CHECK = ROOT/'outputs/M28-A33-EMISSION-BEAM-RESULT-CHECK'
ARMS = ('A31_A3', 'A33')
SHARDS = 4


def require_roots(jobs):
    result = subprocess.run(['sacct', '-X', '-j', ','.join(jobs), '-n', '-P',
        '--format=JobID%30,State%30,ExitCode'], check=True, capture_output=True, text=True)
    receipts = {f[0].strip(): [v.strip() for v in f[1:3]]
                for line in result.stdout.splitlines() if len(f := line.split('|')) >= 3}
    if any(receipts.get(j) != ['COMPLETED', '0:0'] for j in jobs):
        raise ValueError('Actual successful source jobs required: '+str(receipts))


def sources():
    require_roots(['13402927', '13398521', '13398524'])
    ready = R.completed(PREFLIGHT, 13402927)
    if (ready['configuration'] != CONFIG or ready['cases'] != 8 or
            not all(ready[k] for k in ('original_natural_pool_and_A31_scores_reproduced',
                'incremental_final_role_scores_verified', 'returned_raw_grammar_verified',
                'original_pool_geometry_preserved', 'model_weights_unchanged',
                'all_parameters_frozen_and_gradients_absent'))):
        raise ValueError('Closed matching eight-case A33 engineering preflight required')
    baseline = R.completed(R.READOUT, 13398521); checked = R.completed(R.CHECK, 13398524)
    if (baseline['configuration'] != R.CONFIG or checked['configuration'] != R.CONFIG or
            baseline['windows'] != 1491 or checked['windows'] != 1491 or
            not checked['actual_saved_natural_pool_and_greedy_subset_identity_verified'] or
            baseline['selected_counts_per_species'] != checked['selected_counts_per_species'] or
            baseline['eligible_reference_pairs'] != checked['eligible_reference_pairs']):
        raise ValueError('Closed complete original A31 baseline required')
    C.validate_head(R.completed(C.FIT))
    draws = C.R.rows(Q.DRAWS); ids = list(dict.fromkeys(d['window_id'] for d in draws))
    pools = C.R.rows(C.FIT/'free_pool.jsonl')
    saved = C.R.rows(R.READOUT/'natural_scores.jsonl')
    if (len(draws) != 1536 or len(ids) != 1491 or
            any([r['window']['id'] for r in rows] != ids for rows in (pools, saved)) or
            sum(len(r['chains']) for r in pools) != 16321):
        raise ValueError('Original ordered full TRAIN pool and score population required')
    index = {r['window_id']: r for r in Q.train_records(Q.FEATURES/'index.jsonl')
             if r['window_id'] in set(ids)}
    if set(index) != set(ids):
        raise ValueError('All original TRAIN-only cached features required')
    return draws, ids, index, pools, saved, baseline


def same_geometry(a, b):
    if any(a[k] != b[k] for k in C.R.GEOMETRY):
        raise ValueError('Paired TRAIN geometry differs')


def complete_order(records, ids):
    ordered = sorted(records, key=lambda r: r['ordinal'])
    if ([r['ordinal'] for r in ordered] != list(range(len(ids))) or
            [r['window']['id'] for r in ordered] != ids):
        raise ValueError('Complete unique original ordinal/window population required')
    return ordered


def gate(counts, references):
    checks = {}
    for species, values in counts.items():
        new, old = values['A33'], values['A31_A3']
        n, o = new['selected_correct'], old['selected_correct']
        nd = references[species]['eligible']+n+new['selected_known_nonmatch']
        od = references[species]['eligible']+o+old['selected_known_nonmatch']
        checks[species+'_proxy_F1_strictly_increased'] = n*od > o*nd
        checks[species+'_proxy_precision_not_decreased'] = R.precision_compare(new, old)
    checks['pooled_five_or_more_selected_correct_strictly_increased'] = (
        sum(v['A33']['selected_five_or_more_correct'] for v in counts.values()) >
        sum(v['A31_A3']['selected_five_or_more_correct'] for v in counts.values()))
    return checks


def raw_grammar(chains, sequence):
    sequence = sequence.upper(); orf = ORFIndex(sequence); events = _end_events(sequence)
    for chain in chains:
        if not 1 <= len(chain) <= CONFIG['max_exons']:
            raise ValueError('Completed chain outside frozen exon budget')
        carry, previous_end = '', -1
        for depth, (a, b) in enumerate(chain):
            lower = 0 if depth == 0 else previous_end+Budget().min_intron_bp
            if (not 0 <= a < b <= len(sequence) or a < lower or
                    (sequence[a:a+3] != 'ATG' if depth == 0 else sequence[a-2:a] != 'AG')):
                raise ValueError('Candidate entry outside original raw DNA grammar')
            kind = 'stop' if depth == len(chain)-1 else 'donor'
            choices = [(x, k, c) for x, k, c in _legal_ends(orf, events, a, carry, depth == 0)
                       if (x, k) == (b, kind)]
            if len(choices) != 1:
                raise ValueError('Candidate end outside original ORF/codon carry grammar')
            carry, previous_end = choices[0][2], b


def write(out, **values):
    C.R.write_summary(out, dict(completed=True, experiment='M28-A33-EMISSION-BEAM-TRAIN',
        configuration=CONFIG, protocol=PROTOCOL, optimizer_steps_this_run=0, new_parameters=0,
        new_candidate_search_this_phase=values.get('phase') == 'natural',
        new_backbone_forward=False, reference_injected_candidates=0,
        DEV_used=False, test_setaria_used=False, **values))


@torch.inference_mode()
def natural(shard):
    _, ids, index, pools, prior, _ = sources()
    out = NATURAL/('shard_'+str(shard)); (out/'logits').mkdir()
    model = C.load_generator()
    before = {k: v.detach().clone() for k, v in model.state_dict().items()}
    if sum(p.numel() for p in model.parameters()) != 649457:
        raise ValueError('Same frozen A3 capacity required')
    began = time.perf_counter(); counters = Counter(); depths = []; search_seconds = 0.
    error = 0.; replay_error = 0.; torch.cuda.reset_peak_memory_stats()
    ordinals = list(range(shard, len(ids), SHARDS))
    with (out/'natural_scores.jsonl').open('x') as handle:
        for position, ordinal in enumerate(ordinals):
            ident = ids[ordinal]; free = pools[ordinal]; w = free['window']
            same_geometry(w, prior[ordinal]['window'])
            old = prior[ordinal]['variants']['A31_A3']
            if (w['split'] != 'train' or w['seqid'] not in C.R.ALLOW[w['species']] or
                    free['reference_used_in_pool_generation'] or old['chains'] != free['chains']):
                raise ValueError('Closed reference-free original TRAIN-only pool required')
            cached = torch.load(ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
            features, dna = Q.window_inputs(cached, w, 'cuda')
            hidden = model(features, dna, [w['valid_bases']])['features'][0]
            logits = model.segmentation(hidden)
            if logits.dtype != torch.float32:
                raise ValueError('Frozen FP32 model required')
            scores = R.local_path_scores(logits, free['chains'], w['valid_bases'])
            torch.testing.assert_close(scores.cpu(), torch.tensor(old['scores'], dtype=torch.float64),
                                       rtol=1e-12, atol=1e-8)
            started = time.perf_counter()
            queries = decode_emission_query_beam(model.query, hidden, logits, cached['sequence'],
                                                128, 4, 4, w['valid_bases'])
            torch.cuda.synchronize(); search_seconds += time.perf_counter()-started
            emitted = [r for q in queries for r in q['completed']]
            direct = R.local_path_scores(logits, [r['chain'] for r in emitted], w['valid_bases']).cpu()
            incremental = torch.tensor([r['emission_gain'] for r in emitted], dtype=torch.float64)
            torch.testing.assert_close(incremental, direct, rtol=1e-12, atol=1e-8)
            error = max(error, max((abs(float(x-y)) for x, y in zip(incremental, direct)), default=0.))
            window = model.query.for_window(hidden, w['valid_bases'])
            for row in emitted:
                pointer = replay_node_score(window, cached['sequence'], row['query_index'],
                                            dict(chain=row['chain'], entry=None, event=STOP))
                torch.testing.assert_close(pointer, pointer.new_tensor(row['path_log_score']),
                                           rtol=1e-6, atol=1e-5)
                replay_error = max(replay_error, abs(float(pointer)-row['path_log_score']))
            chains, decoder_sources = augment_pool(free['chains'], free['decoder_sources'], queries)
            values = R.local_path_scores(logits, chains, w['valid_bases']).cpu().tolist()
            logits_path = out/'logits'/('%04d.pt' % ordinal)
            torch.save(dict(window=w, logits=logits[:w['valid_bases']].detach().cpu().clone()), logits_path)
            record = dict(ordinal=ordinal, window=w, variants=dict(A31_A3=old,
                A33=dict(chains=chains, scores=values, selected=R.select_nonoverlapping(chains, values))),
                decoder_sources=decoder_sources, emission_queries=queries,
                logits_artifact=str(logits_path.relative_to(ROOT)),
                reference_used_in_pool_generation=False, reference_used_in_scoring=False)
            handle.write(json.dumps(record, allow_nan=False)+'\n')
            counters.update(windows=1, old_candidates=len(free['chains']), candidates=len(chains),
                            new_unique=len(chains)-len(free['chains']), completed_query_paths=len(emitted))
            for q in queries:
                depths.append(q['counts']['depths_visited'])
                counters.update({'query_'+k: v for k, v in q['counts'].items() if k != 'depths_visited'})
            if (position+1) % 32 == 0 or position+1 == len(ordinals):
                handle.flush(); print(json.dumps(dict(shard=shard, windows=position+1,
                    total=len(ordinals), seconds=time.perf_counter()-began)), flush=True)
            del cached, features, dna, hidden, logits, window, queries, emitted
    if (any(not torch.equal(v, before[k]) for k, v in model.state_dict().items()) or
            any(p.requires_grad or p.grad is not None for p in model.parameters())):
        raise ValueError('Frozen A3 changed or gradients appeared')
    write(out, phase='natural', shard=shard, shards=SHARDS, windows=len(ordinals), ordinals=ordinals,
        counts=dict(counters), depth_range=[min(depths), max(depths)], model_weights_unchanged=True,
        all_parameters_frozen_and_gradients_absent=True, original_A31_A3_scores_reproduced=True,
        incremental_final_role_scores_verified=True, returned_raw_grammar_and_pointer_scores_verified=True,
        max_incremental_vs_direct_error=error, max_pointer_replay_error=replay_error,
        search_seconds=search_seconds, seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(), new_reference_labelling=False,
        limitation='Natural shard only; no effect labels, gate, gene F1 or generalization.')


def closed_natural(array_job, ids):
    jobs = [array_job+'_'+str(i) for i in range(SHARDS)]; require_roots(jobs)
    records, reports = [], []
    for shard, job in enumerate(jobs):
        out = NATURAL/('shard_'+str(shard)); report = R.completed(out, job)
        expected = list(range(shard, len(ids), SHARDS))
        if (report['configuration'] != CONFIG or report['phase'] != 'natural' or
                report['ordinals'] != expected or report['windows'] != len(expected) or
                not all(report[k] for k in ('model_weights_unchanged',
                    'all_parameters_frozen_and_gradients_absent', 'original_A31_A3_scores_reproduced',
                    'incremental_final_role_scores_verified', 'returned_raw_grammar_and_pointer_scores_verified'))):
            raise ValueError('Matching complete frozen natural shard required')
        rows = C.R.rows(out/'natural_scores.jsonl')
        if [r['ordinal'] for r in rows] != expected:
            raise ValueError('Complete actual natural shard file required')
        records.extend(rows); reports.append(report)
    return complete_order(records, ids), reports


def tabulate(records, draws, ids, index, pools, prior, baseline, independent=False):
    loaded, _, windows, catalogs = Q.load_train_inputs()
    if loaded != draws or set(windows) != set(ids):
        raise ValueError('Unchanged original TRAIN label population required')
    old_labels = C.R.rows(R.READOUT/'windows.jsonl')
    if [r['window']['id'] for r in old_labels] != ids:
        raise ValueError('Original complete validated per-window baseline labels required')
    counts = {s: {a: Counter() for a in ARMS} for s in C.R.ALLOW}
    references = {s: Counter(eligible=0, five_or_more=0) for s in C.R.ALLOW}; labelled = []
    max_score_error = 0.
    for ordinal, (record, pool, old, labels) in enumerate(zip(records, pools, prior, old_labels)):
        w = windows[ids[ordinal]]; species = w['species']; variants = record['variants']
        for other in (record['window'], pool['window'], old['window'], labels['window']):
            same_geometry(w, other)
        if (record['reference_used_in_pool_generation'] or record['reference_used_in_scoring'] or
                variants['A31_A3'] != old['variants']['A31_A3']):
            raise ValueError('Reference-free scoring and unchanged stronger baseline required')
        if len(record['emission_queries']) != 16 or any(len(q['completed']) > 4 for q in record['emission_queries']):
            raise ValueError('Original Q16/completed4 population required')
        chains, origins = augment_pool(pool['chains'], pool['decoder_sources'], record['emission_queries'])
        chains = [[list(e) for e in chain] for chain in chains]
        if chains != variants['A33']['chains'] or origins != record['decoder_sources']:
            raise ValueError('Actual old-pool-preserving augmented geometry/source identity differs')
        cached = torch.load(ROOT/index[w['id']]['artifact'], map_location='cpu', weights_only=True)
        Q.validate_cached(cached, w); sequence = cached['sequence']; del cached
        if independent:
            raw_grammar(chains, sequence)
            state = torch.load(ROOT/record['logits_artifact'], map_location='cpu', weights_only=True)
            same_geometry(w, state['window']); logits = state['logits']
            if logits.dtype != torch.float32 or tuple(logits.shape) != (w['valid_bases'], 15):
                raise ValueError('Saved valid FP32 emissions required for independent CPU sum')
            logp = logits.double().log_softmax(-1); positions = torch.arange(len(logp))
            for value in variants.values():
                scores = []
                for chain in value['chains']:
                    roles = R.candidate_roles(chain, len(logp))
                    scores.append(float((logp[positions, roles]-logp[:, 0]).sum()))
                torch.testing.assert_close(torch.tensor(scores, dtype=torch.float64),
                    torch.tensor(value['scores'], dtype=torch.float64), rtol=1e-12, atol=1e-8)
                max_score_error = max(max_score_error,
                    max((abs(x-y) for x, y in zip(scores, value['scores'])), default=0.))
            del state, logits, logp
        targets = Q.build_targets(w, catalogs[species], sequence)
        ref = dict(eligible=len(targets['positive_chains']),
                   five_or_more=sum(len(c) >= 5 for c in targets['positive_chains']))
        if (ref['eligible'] != labels['eligible_references'] or
                ref['five_or_more'] != labels['eligible_five_or_more']):
            raise ValueError('Original eligible R4 reference denominator differs')
        references[species].update(ref); actual_labels, totals = {}, {}
        for arm in ARMS:
            actual_labels[arm] = R.candidate_targets(variants[arm]['chains'], targets)
            totals[arm] = R.selected_counts(variants[arm], actual_labels[arm]); counts[species][arm].update(totals[arm])
        if (actual_labels['A31_A3'] != labels['labels']['A31_A3'] or
                totals['A31_A3'] != labels['counts']['A31_A3']):
            raise ValueError('Original per-window baseline labels and counts must reproduce')
        labelled.append(dict(window=record['window'], labels=actual_labels, counts=totals,
            eligible_references=ref['eligible'], eligible_five_or_more=ref['five_or_more']))
        if (ordinal+1) % 128 == 0:
            print(json.dumps(dict(phase='CPU_check' if independent else 'R4_labels', windows=ordinal+1)), flush=True)
    counts = {s: {a: dict(c) for a, c in v.items()} for s, v in counts.items()}
    references = {s: dict(c) for s, c in references.items()}
    if (sum(v['eligible'] for v in references.values()) != 2916 or
            references != baseline['eligible_reference_pairs'] or
            any(v['A31_A3'] != baseline['selected_counts_per_species'][s]['A31_A3'] for s, v in counts.items())):
        raise ValueError('All original eligible pairs and stronger-baseline totals must reproduce')
    return labelled, counts, references, max_score_error


def readout(array_job):
    began = time.perf_counter(); inputs = sources(); _, ids, _, _, _, _ = inputs
    records, reports = closed_natural(array_job, ids)
    with (READOUT/'natural_scores.jsonl').open('x') as handle:
        for record in records:
            handle.write(json.dumps(record, allow_nan=False)+'\n')
    # No TRAIN catalog or new effect labels opened until every natural file is closed.
    labelled, counts, references, _ = tabulate(records, *inputs)
    with (READOUT/'windows.jsonl').open('x') as handle:
        for record in labelled:
            handle.write(json.dumps(record, allow_nan=False)+'\n')
    checks = gate(counts, references)
    write(READOUT, phase='readout', windows=1491, natural_array_job=array_job,
        all_natural_files_closed_before_catalog_open=True, model_weights_unchanged=True,
        all_parameters_frozen_and_gradients_absent=True, original_pool_and_baseline_preserved=True,
        selected_counts_per_species=counts, eligible_reference_pairs=references,
        known_label_window_pair_proxies=R.proxies(counts, references), predeclared_TRAIN_gate=checks,
        TRAIN_preparation_supported=all(checks.values()), natural_shard_reports=reports,
        seconds=time.perf_counter()-began,
        limitation='TRAIN window/chain proxy only; requires independent CPU verification. '
            'Extra search compute is not matched. No native FASTA, gene F1, independent generalization or publication PASS.')


def check(array_job, source_job):
    began = time.perf_counter(); require_roots([source_job]); report = R.completed(READOUT, source_job)
    if (report['configuration'] != CONFIG or report['phase'] != 'readout' or report['windows'] != 1491 or
            report['natural_array_job'] != array_job or not report['all_natural_files_closed_before_catalog_open']):
        raise ValueError('Matching complete original A33 readout required')
    inputs = sources(); _, ids, _, _, _, _ = inputs
    records, _ = closed_natural(array_job, ids)
    if records != C.R.rows(READOUT/'natural_scores.jsonl'):
        raise ValueError('Actual whole natural source files must match merged file')
    labelled, counts, references, error = tabulate(records, *inputs, independent=True)
    checks = gate(counts, references)
    if (labelled != C.R.rows(READOUT/'windows.jsonl') or counts != report['selected_counts_per_species'] or
            references != report['eligible_reference_pairs'] or checks != report['predeclared_TRAIN_gate'] or
            R.proxies(counts, references) != report['known_label_window_pair_proxies'] or
            all(checks.values()) != report['TRAIN_preparation_supported']):
        raise ValueError('Actual saved labels, selection, denominators, metrics or fixed gates differ')
    write(CHECK, phase='check', windows=1491, natural_array_job=array_job, readout_job=source_job,
        actual_old_pool_augmented_geometry_and_sources_verified=True,
        actual_all_candidate_raw_grammar_verified=True,
        actual_saved_emission_direct_log_probability_scores_recomputed=True,
        actual_R4_labels_selection_reference_denominators_and_proxy_metrics_recomputed=True,
        actual_five_fixed_gate_integer_arithmetic_verified=True, max_direct_score_error=error,
        selected_counts_per_species=counts, eligible_reference_pairs=references,
        predeclared_TRAIN_gate=checks, TRAIN_preparation_supported=all(checks.values()),
        model_or_backbone_forward=False, seconds=time.perf_counter()-began,
        limitation='Independent CPU recalculates saved logits/grammar/R4 labels, not a new model forward. '
            'TRAIN proxy only, not full gene F1, native competitiveness, causal proof or generalization.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('natural', 'readout', 'check'), required=True)
    parser.add_argument('--shard', type=int, choices=range(SHARDS))
    parser.add_argument('--array-job')
    parser.add_argument('--source-job')
    args = parser.parse_args(); torch.set_num_threads(2)
    if args.mode == 'natural' and args.shard is not None:
        natural(args.shard)
    elif args.mode == 'readout' and args.array_job:
        readout(args.array_job)
    elif args.mode == 'check' and args.array_job and args.source_job:
        check(args.array_job, args.source_job)
    else:
        parser.error('Matching shard/actual source job arguments required')
