#!/usr/bin/env python3
"""A31 frozen local emission caller: zero-update preflight and TRAIN readout."""
import argparse
from collections import Counter
import json
import math
import time
import torch
import run_chain_quality as C
import run_query_relative_boundary as R23
import diagnose_quality_train as D
from src.m28.core import select_nonoverlapping
from src.m28.labels import candidate_targets
from src.m28.local_path_score import local_path_scores
from src.m28.phase_evidence_chain import candidate_roles
from src.m28.query_decode import decode_query_chains

ROOT, Q = C.ROOT, C.R.Q
PREFLIGHT = ROOT/'outputs/M28-A31-LOCAL-PATH-CALLER-PREFLIGHT'
READOUT = ROOT/'outputs/M28-A31-LOCAL-PATH-CALLER-TRAIN-READOUT'
CHECK = ROOT/'outputs/M28-A31-LOCAL-PATH-CALLER-RESULT-CHECK'
PROTOCOL = 'reports/M28-METHOD-RESTART/a31_local_path_caller_protocol.md'
CONFIG = dict(generators={'A3': C.GENERATOR,
    'A23': str((R23.FIT/'step_004608.pt').relative_to(ROOT))},
    score='sum_valid_candidate_span_role_logit_minus_intergenic_logit',
    all_role_weights=1., null=0., model_precision='FP32', score_accumulation='FP64',
    main_baseline='A3_completed_greedy_max_query_logit_per_geometry',
    pool='reuse_original_completed_greedy_plus_beam4_union', beam_width=4,
    completed_limit=4, max_exons=128, threshold=0., threshold_scan=False,
    new_parameters=0, optimizer_steps=0, unique_TRAIN_windows=1491)
ARMS = ('A3_greedy', 'A31_A3', 'A31_A23', 'A6_auxiliary')


def completed(folder, job=None):
    if (folder/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Completed source required: '+str(folder))
    if job is not None and (folder/'JOBID').read_text().strip() != str(job):
        raise ValueError('Source job identity differs')
    value = json.loads((folder/'summary.json').read_text())
    if not value['completed']:
        raise ValueError('Incomplete source summary')
    return value


def sources():
    C.validate_head(completed(C.FIT))
    R23.validate_fit(completed(R23.FIT, 13391608))
    report = completed(R23.READOUT, 13391609)
    check = completed(ROOT/'outputs/M28-A23-RELATIVE-BOUNDARY-RESULT-CHECK', 13391610)
    if not (check['actual_saved_pool_identity_and_three_fixed_gate_checks_verified'] and
            check['actual_new_and_original_natural_chain_labels_recomputed'] and
            report['configuration'] == R23.CONFIG and report['protocol'] == R23.PROTOCOL):
        raise ValueError('Actual original A23 source validation required')
    draws = C.R.rows(Q.DRAWS)
    ids = list(dict.fromkeys(d['window_id'] for d in draws))
    pools = {a: C.R.rows(p/'free_pool.jsonl') for a, p in
             (('A3', C.FIT), ('A23', R23.READOUT))}
    if len(draws) != 1536 or len(ids) != 1491:
        raise ValueError('Original TRAIN population required')
    for a, rows in pools.items():
        if [r['window']['id'] for r in rows] != ids:
            raise ValueError('Original ordered source windows required')
        if sum(len(r['chains']) for r in rows) != {'A3': 16321, 'A23': 18950}[a]:
            raise ValueError('Complete original natural pool required')
        for r in rows:
            w = r['window']
            if (r['reference_used_in_pool_generation'] or w['split'] != 'train' or
                    w['seqid'] not in C.R.ALLOW[w['species']] or
                    len(r['chains']) != len(r['decoder_sources']) or
                    len({tuple(map(tuple, c)) for c in r['chains']}) != len(r['chains'])):
                raise ValueError('Unique, reference-free allowed TRAIN pool required')
    if any(any(x['window'][k] != y['window'][k] for k in C.R.GEOMETRY)
           for x, y in zip(pools['A3'], pools['A23'])):
        raise ValueError('Paired source geometry differs')
    index = {r['window_id']: r for r in Q.train_records(Q.FEATURES/'index.jsonl')
             if r['window_id'] in set(ids)}
    if set(index) != set(ids):
        raise ValueError('Original TRAIN-only features required')
    return draws, ids, index, {a: {r['window']['id']: r for r in v} for a, v in pools.items()}


def load_models():
    a3 = C.load_generator()
    fitted = completed(R23.FIT, 13391608); R23.validate_fit(fitted)
    state = torch.load(ROOT/fitted['final_checkpoint'], map_location='cpu', weights_only=True)
    R23.validate_fit(state, state=True)
    a23 = R23.new_model().eval().requires_grad_(False)
    a23.load_state_dict(state['model'], strict=True)
    return {'A3': a3, 'A23': a23}


def greedy_record(model, hidden, sequence, free):
    queries = decode_query_chains(model.query, hidden, sequence, 128, len(sequence))
    geometry = {}
    for r in queries:
        if r['terminated']:
            if r['stop_reason'] != 'canonical_stop':
                raise ValueError('Only original completed canonical greedy is allowed')
            chain = tuple(map(tuple, r['chain']))
            geometry[chain] = max(geometry.get(chain, -math.inf), r['query_logit'])
    expected = {tuple(map(tuple, c)) for c, s in zip(free['chains'], free['decoder_sources'])
                if 'greedy' in s}
    if set(geometry) != expected:
        raise ValueError('Original A3 greedy source replay differs')
    chains = sorted(geometry)
    scores = [geometry[c] for c in chains]
    if not all(math.isfinite(v) for v in scores):
        raise ValueError('Nonfinite original greedy existence score')
    return dict(chains=chains, scores=scores, selected=select_nonoverlapping(chains, scores),
                source_geometry_replay_verified=True)


@torch.inference_mode()
def score_window(models, index, pools, ident, direct=False):
    w = pools['A3'][ident]['window']
    cached = torch.load(ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
    Q.validate_cached(cached, w)
    features, dna = Q.window_inputs(cached, w, 'cuda')
    variants, errors = {}, {}
    for a, model in models.items():
        hidden = model(features, dna, [w['valid_bases']])['features'][0]
        logits = model.segmentation(hidden)
        if logits.dtype != torch.float32:
            raise ValueError('Unchanged frozen FP32 forward required')
        chains = pools[a][ident]['chains']
        scores = local_path_scores(logits, chains, w['valid_bases']).cpu().tolist()
        variants['A31_'+a] = dict(chains=chains, scores=scores,
            selected=select_nonoverlapping(chains, scores))
        if direct:
            logp = logits[:w['valid_bases']].double().log_softmax(-1)
            values = []
            for chain in chains:
                roles = candidate_roles(chain, len(logp), logp.device)
                # Full-window reference-free direct log-probability difference.
                values.append(float((logp[torch.arange(len(logp), device=logp.device), roles]-logp[:, 0]).sum()))
            torch.testing.assert_close(torch.tensor(scores, dtype=torch.float64),
                torch.tensor(values, dtype=torch.float64), rtol=1e-12, atol=1e-8)
            errors[a] = max((abs(x-y) for x, y in zip(scores, values)), default=0.)
        if a == 'A3':
            variants['A3_greedy'] = greedy_record(model, hidden, cached['sequence'], pools[a][ident])
    return dict(window=w, variants=variants, reference_used_in_scoring=False), cached['sequence'], errors


def write(out, result):
    C.R.write_summary(out, dict(arm='A31', protocol=PROTOCOL, configuration=CONFIG,
        optimizer_steps_this_run=0, reference_injected_candidates=0,
        new_candidate_search=False, new_backbone_forward=False,
        DEV_forward_used=False, test_setaria_used=False, **result))


@torch.inference_mode()
def preflight():
    draws, ids, index, pools = sources()
    loaded, _, windows, catalogs = Q.load_train_inputs()
    if loaded != draws or set(windows) != set(ids):
        raise ValueError('Unchanged preflight metadata cases required')
    cases = Q.select_timing_draws(draws, windows, catalogs)
    if len(cases) != 8:
        raise ValueError('Original eight metadata cases required')
    models = load_models()
    before = {a: {k: v.detach().clone() for k, v in m.state_dict().items()} for a, m in models.items()}
    began = time.perf_counter(); records = []; torch.cuda.reset_peak_memory_stats()
    for case in cases:
        ident = case['draw']['window_id']; started = time.perf_counter()
        record, _, errors = score_window(models, index, pools, ident, direct=True)
        torch.cuda.synchronize()
        row = dict(window_id=ident, species=record['window']['species'], category=case['category'],
            draw_ordinal=case['ordinal'], max_direct_log_probability_error=errors,
            candidate_counts={a: len(v['chains']) for a, v in record['variants'].items()},
            original_greedy_source_replay_verified=True, two_model_seconds=time.perf_counter()-started)
        records.append(row); print(json.dumps(row, allow_nan=False), flush=True)
    for a, model in models.items():
        if (any(not torch.equal(v, before[a][k]) for k, v in model.state_dict().items()) or
                any(p.requires_grad or p.grad is not None for p in model.parameters())):
            raise ValueError('Frozen model changed or received gradient')
    write(PREFLIGHT, dict(completed=True, cases=8, records=records, model_weights_unchanged=True,
        original_greedy_source_replay_verified=True, direct_log_probability_identity_verified=True,
        all_parameters_frozen_and_gradients_absent=True,
        model_parameters={a: sum(p.numel() for p in m.parameters()) for a, m in models.items()},
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(), seconds=time.perf_counter()-began,
        limitation='Eight metadata cases establish engineering readiness only; no effect gate.'))


def selected_counts(value, labels):
    chains, scores, selected = value['chains'], value['scores'], value['selected']
    if (len(chains) != len(scores) or len(labels) != len(chains) or
            not all(math.isfinite(s) for s in scores) or
            selected != select_nonoverlapping(chains, scores)):
        raise ValueError('Saved scores and original selection do not reconcile')
    return dict(candidates=len(chains), pool_correct=labels.count(1),
        pool_known_nonmatch=labels.count(0), pool_unknown=labels.count(-1), selected=len(selected),
        selected_correct=sum(labels[i] == 1 for i in selected),
        selected_known_nonmatch=sum(labels[i] == 0 for i in selected),
        selected_unknown=sum(labels[i] == -1 for i in selected),
        selected_five_or_more_correct=sum(labels[i] == 1 and len(chains[i]) >= 5 for i in selected))


def precision_compare(new, old):
    n, o = new['selected_correct'], old['selected_correct']
    nd, od = n+new['selected_known_nonmatch'], o+old['selected_known_nonmatch']
    # Empty known selection has proxy precision zero, not a spurious perfect score.
    return (n if nd else 0)*(od or 1) >= (o if od else 0)*(nd or 1)


def gate(counts, references):
    checks = {}
    for arm in ('A31_A3', 'A31_A23'):
        local = {}
        for s, values in counts.items():
            new, old = values[arm], values['A3_greedy']
            n, o = new['selected_correct'], old['selected_correct']
            nd = references[s]['eligible']+n+new['selected_known_nonmatch']
            od = references[s]['eligible']+o+old['selected_known_nonmatch']
            local[s+'_proxy_F1_strictly_increased'] = n*od > o*nd
            local[s+'_proxy_precision_not_decreased'] = precision_compare(new, old)
        local['pooled_five_or_more_selected_correct_strictly_increased'] = (
            sum(v[arm]['selected_five_or_more_correct'] for v in counts.values()) >
            sum(v['A3_greedy']['selected_five_or_more_correct'] for v in counts.values()))
        checks[arm] = local
    return checks


def proxies(counts, references):
    result = {}
    for s, values in counts.items():
        result[s] = {}
        for a, c in values.items():
            n, fp, r = c['selected_correct'], c['selected_known_nonmatch'], references[s]['eligible']
            result[s][a] = dict(precision=n/(n+fp) if n+fp else 0.,
                recall=n/r if r else 0., F1=2*n/(r+n+fp) if r+n+fp else 0.)
    return result


def label_readout(saved, pools, windows, catalogs, sequences, auxiliary):
    counts = {s: {a: Counter() for a in ARMS} for s in C.R.ALLOW}
    references = {s: Counter(eligible=0, five_or_more=0) for s in C.R.ALLOW}
    labelled = []
    for record in saved:
        ident = record['window']['id']; w = windows[ident]; species = w['species']
        if record['reference_used_in_scoring'] or any(w[k] != record['window'][k] for k in C.R.GEOMETRY):
            raise ValueError('Natural scoring/label geometry mismatch')
        targets = Q.build_targets(w, catalogs[species], sequences[ident])
        references[species].update(eligible=len(targets['positive_chains']),
            five_or_more=sum(len(c) >= 5 for c in targets['positive_chains']))
        variants = record['variants']
        for a in ('A3', 'A23'):
            if variants['A31_'+a]['chains'] != pools[a][ident]['chains']:
                raise ValueError('Do not alter the original natural union pool')
        expected = {tuple(map(tuple, c)) for c, source in zip(pools['A3'][ident]['chains'],
                    pools['A3'][ident]['decoder_sources']) if 'greedy' in source}
        if {tuple(map(tuple, c)) for c in variants['A3_greedy']['chains']} != expected:
            raise ValueError('Saved original greedy subset differs')
        old = [r for r in auxiliary[ident]['candidates'] if r['origin'] == 'free']
        if [r['chain'] for r in old] != pools['A3'][ident]['chains']:
            raise ValueError('Auxiliary A6 must use the same original natural pool')
        chains, scores = [r['chain'] for r in old], [r['gain'] for r in old]
        variants = dict(variants, A6_auxiliary=dict(chains=chains, scores=scores,
            selected=select_nonoverlapping(chains, scores)))
        labels, totals = {}, {}
        for a, value in variants.items():
            labels[a] = candidate_targets(value['chains'], targets)
            if a == 'A6_auxiliary' and labels[a] != [r['label'] for r in old]:
                raise ValueError('Auxiliary original R4 labels differ')
            totals[a] = selected_counts(value, labels[a]); counts[species][a].update(totals[a])
        labelled.append(dict(window=record['window'], labels=labels, counts=totals,
            eligible_references=len(targets['positive_chains']),
            eligible_five_or_more=sum(len(c) >= 5 for c in targets['positive_chains'])))
    counts = {s: {a: dict(c) for a, c in v.items()} for s, v in counts.items()}
    references = {s: dict(c) for s, c in references.items()}
    return labelled, counts, references


@torch.inference_mode()
def readout():
    ready = completed(PREFLIGHT)
    if not (ready['configuration'] == CONFIG and ready['protocol'] == PROTOCOL and
            ready['model_weights_unchanged'] and ready['direct_log_probability_identity_verified'] and
            ready['original_greedy_source_replay_verified'] and ready['optimizer_steps_this_run'] == 0):
        raise ValueError('Matching completed A31 engineering preflight required')
    draws, ids, index, pools = sources(); models = load_models()
    before = {a: {k: v.detach().clone() for k, v in m.state_dict().items()} for a, m in models.items()}
    began = time.perf_counter(); saved = []; sequences = {}; torch.cuda.reset_peak_memory_stats()
    # No TRAIN annotation/catalog is opened until the full natural score file closes.
    with (READOUT/'natural_scores.jsonl').open('x') as handle:
        for ordinal, ident in enumerate(ids):
            record, sequence, _ = score_window(models, index, pools, ident)
            handle.write(json.dumps(record, allow_nan=False)+'\n'); saved.append(record); sequences[ident] = sequence
            if (ordinal+1) % 64 == 0:
                handle.flush(); print(json.dumps(dict(phase='natural_scores', windows=ordinal+1, total=1491)), flush=True)
    scoring_seconds = time.perf_counter()-began
    for a, model in models.items():
        if (any(not torch.equal(v, before[a][k]) for k, v in model.state_dict().items()) or
                any(p.requires_grad or p.grad is not None for p in model.parameters())):
            raise ValueError('Frozen model changed or received gradient')
    loaded, _, windows, catalogs = Q.load_train_inputs()
    if loaded != draws or set(windows) != set(ids):
        raise ValueError('Original TRAIN label population required')
    completed(D.OUT)
    auxiliary = {r['window']['id']: r for r in C.R.rows(D.OUT/'windows.jsonl')}
    labelled, counts, references = label_readout(saved, pools, windows, catalogs, sequences, auxiliary)
    with (READOUT/'windows.jsonl').open('x') as handle:
        for record in labelled:
            handle.write(json.dumps(record, allow_nan=False)+'\n')
    checks = gate(counts, references)
    write(READOUT, dict(completed=True, windows=1491, natural_score_file_closed_before_catalog_open=True,
        original_greedy_source_replay_verified=True, model_weights_unchanged=True,
        all_parameters_frozen_and_gradients_absent=True, original_natural_pools_unchanged=True,
        selected_counts_per_species=counts, eligible_reference_pairs=references,
        known_label_window_pair_proxies=proxies(counts, references), predeclared_TRAIN_gate=checks,
        TRAIN_preparation_supported={a: all(v.values()) for a, v in checks.items()},
        scoring_seconds=scoring_seconds, seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='In-sample window/chain readiness, not full gene F1, generalization or a publication result.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('preflight', 'readout'), required=True)
    args = parser.parse_args(); torch.set_num_threads(2)
    (preflight if args.mode == 'preflight' else readout)()
