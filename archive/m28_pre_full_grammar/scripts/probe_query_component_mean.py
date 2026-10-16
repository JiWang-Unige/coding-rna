#!/usr/bin/env python3
"""A11 zero-update score/component TRAIN probe; no final annotation or DEV."""
import json
import math
import time
from collections import Counter
import torch
import run_full_train as R
import run_chain_quality as C
import analyze_query_beam_train as B
from smoke_query_beam import check_chain
from src.m28.query_component_beam import ARMS, decode_query_component_pair
from src.m28.labels import build_targets

OUT = R.ROOT/'outputs/M28-A11-COMPONENT-MEAN-TRAIN-PROBE'
PROTOCOL = 'reports/M28-METHOD-RESTART/a11_component_mean_protocol.md'
LIMITS = dict(max_exons=128, beam_width=4, completed_limit=4)
ALL_ARMS = ('ORIGINAL',)+ARMS


def keys(rows):
    return {B.chain_key(r['chain']) for r in rows}


def metrics(rows):
    return {arm: sum(r['hits'][arm] for r in rows) for arm in ALL_ARMS}


@torch.inference_mode()
def main():
    source = json.loads((B.OUT/'summary.json').read_text())
    if ((B.OUT/'STATUS').read_text().strip() != 'COMPLETED' or not source['completed'] or
            source['windows'] != 8 or source['checkpoint'] != C.GENERATOR or
            source['pooled']['reference_window_pairs'] != 27):
        raise ValueError('Completed original eight-window/27-pair A5 source required')
    R.validate_feature_summary(json.loads((R.FEATURES/'summary.json').read_text()))
    predictions = R.rows(B.OUT/'predictions.jsonl')
    if len(predictions) != 8 or len({p['window']['id'] for p in predictions}) != 8:
        raise ValueError('Original eight TRAIN windows cannot be reselected or dropped')
    if Counter(p['window']['species'] for p in predictions) != Counter({s: 4 for s in R.ALLOW}):
        raise ValueError('Original four windows per species required')
    plan = {e['window']['id']: e for e in R.rows(R.PLAN/'feature_plan.jsonl')}
    generator = C.load_generator()
    if any(p.requires_grad for p in generator.parameters()):
        raise ValueError('A3 generator must remain frozen')
    began = time.perf_counter()
    torch.cuda.reset_peak_memory_stats()
    free, replay_equal, max_delta = {}, True, 0.
    with (OUT/'predictions.jsonl').open('x') as log:
        for pred in predictions:
            w = pred['window']
            if (w['split'] != 'train' or w['seqid'] not in R.ALLOW[w['species']] or
                    w['valid_bases'] != 24576 or pred['reference_used_in_decoding'] or
                    any(w[k] != plan[w['id']]['window'][k] for k in R.GEOMETRY)):
                raise ValueError('Unchanged allowed full-width TRAIN geometry required')
            cached = torch.load(R.ROOT/plan[w['id']]['artifact'], map_location='cpu', weights_only=True)
            if cached['sequence'] != pred['sequence']:
                raise ValueError('Original source DNA changed')
            features, dna = R.Q.window_inputs(cached, w, 'cuda')
            hidden = generator(features, dna, [w['valid_bases']])['features'][0]
            torch.cuda.synchronize()
            started = time.perf_counter()
            arms = decode_query_component_pair(generator.query, hidden, pred['sequence'], **LIMITS)
            torch.cuda.synchronize()
            for arm in ARMS:
                if [q['query_index'] for q in arms[arm]] != list(range(16)):
                    raise ValueError('Original Q16 identity changed')
                for q in arms[arm]:
                    if len(q['completed']) > 4 or q['counts']['score_bound_stopped_live'] != 0:
                        raise ValueError('Completion cap or no-bound protocol changed')
                    for row in q['completed']:
                        check_chain(row['chain'], pred['sequence'].upper())
                        if not math.isfinite(row['completed_score']):
                            raise ValueError('Nonfinite score receipt')
            for a, b in zip(arms['COMPONENT_SUM'], arms['COMPONENT_MEAN']):
                if a['counts'] != b['counts']:
                    raise ValueError('Weighted completed pools must share one actual search')
            matches = True
            for actual, old in zip(arms['RAW_SUM'], pred['beam_queries']):
                a = {B.chain_key(r['chain']): r['path_log_score'] for r in actual['completed']}
                b = {B.chain_key(r['chain']): r['path_log_score'] for r in old['completed']}
                matches &= a.keys() == b.keys()
                for chain in a.keys() & b.keys():
                    max_delta = max(max_delta, abs(a[chain]-b[chain]))
            replay_equal &= matches
            record = dict(window=w, arms=arms, weighted_search_shared=True,
                          raw_no_bound_matches_original=matches,
                          reference_used_in_decoding=False, seconds=time.perf_counter()-started)
            log.write(json.dumps(record, allow_nan=False)+'\n')
            log.flush()
            free[w['id']] = record
            print(json.dumps(dict(window_id=w['id'], seconds=record['seconds'],
                                  raw_no_bound_matches_original=matches)), flush=True)
            del cached, features, dna, hidden
    prediction_seconds = time.perf_counter()-began
    # Every new free prediction is closed before opening the TRAIN catalog.
    windows, catalogs = B.load_selected_labels(predictions)
    original = R.rows(B.OUT/'references.jsonl')
    if len(original) != 27 or len({(r['window_id'], r['reference_id']) for r in original}) != 27:
        raise ValueError('All original27 primary reference/window identities required')
    expected = {(w['id'], ident) for w in windows.values() for ident in w['positive_ids']}
    if expected != {(r['window_id'], r['reference_id']) for r in original}:
        raise ValueError('Original R4 primary identity set changed')
    for ref in original:
        w = windows[ref['window_id']]
        current = B.oriented_parts(catalogs[w['species']][ref['reference_id']], w)
        if B.chain_key(current) != B.chain_key(ref['chain']):
            raise ValueError('Original primary chain geometry changed')
    pools, candidate_records = {}, []
    for pred in predictions:
        w = windows[pred['window']['id']]
        greedy = [r for r in pred['greedy_queries'] if r['terminated']]
        old_beam = [r for q in pred['beam_queries'] for r in q['completed']]
        original_pool = keys(greedy+old_beam)
        per_arm = {'ORIGINAL': original_pool}
        beam_only = {}
        for arm in ARMS:
            rows = [r for q in free[w['id']]['arms'][arm] for r in q['completed']]
            beam_only[arm] = keys(rows)
            per_arm[arm] = beam_only[arm] | keys(greedy)
        targets = build_targets(w, catalogs[w['species']], pred['sequence'])
        pools[w['id']] = per_arm
        candidate_records.append(dict(window_id=w['id'], species=w['species'],
            arms={arm: dict(unique_geometry=len(chains),
                            reference_relative_labels=B.label_counts(chains, targets))
                  for arm, chains in per_arm.items()},
            beam_only_unique_geometry={arm: len(chains) for arm, chains in beam_only.items()}))
    references = []
    for ref in original:
        chain = B.chain_key(ref['chain'])
        references.append(dict(ref, hits={arm: chain in pools[ref['window_id']][arm] for arm in ALL_ARMS}))
    with (OUT/'references.jsonl').open('x') as log:
        for row in references:
            log.write(json.dumps(row, allow_nan=False)+'\n')
    counts = metrics(references)
    per_species = {s: metrics([r for r in references if r['species'] == s]) for s in R.ALLOW}
    new = [r for r in references if r['hits']['COMPONENT_MEAN'] and not r['hits']['COMPONENT_SUM']]
    lost = [r for r in references if r['hits']['COMPONENT_SUM'] and not r['hits']['COMPONENT_MEAN']]
    new_multi = sum(r['exons'] >= 5 for r in new)
    conditions = dict(mean_strictly_above_component_sum=counts['COMPONENT_MEAN'] > counts['COMPONENT_SUM'],
                      new_at_least_one_five_exon_pair=new_multi >= 1,
                      each_species_not_below_original=all(v['COMPONENT_MEAN'] >= v['ORIGINAL'] for v in per_species.values()),
                      raw_replay_identity_and_scores=replay_equal and max_delta <= 1e-3,
                      original27_retained=len(references) == 27)
    R.write_summary(OUT, dict(experiment='M28-A11-COMPONENT-MEAN-TRAIN-PROBE', completed=True,
        protocol=PROTOCOL, generator_checkpoint=C.GENERATOR, source_optimizer_steps=4608,
        windows=8, primary_reference_window_pairs=27, limits=LIMITS, arms=list(ALL_ARMS),
        counts=counts, per_species=per_species,
        by_exon_class={k: metrics([r for r in references if r['exon_class'] == k])
                       for k in ('one', 'two_to_four', 'five_or_more')},
        mean_vs_component_sum_new=len(new), mean_vs_component_sum_lost=len(lost),
        new_five_or_more_exon_pairs=new_multi, candidate_records=candidate_records,
        raw_no_bound_matches_original=replay_equal, common_raw_score_max_abs_difference=max_delta,
        weighted_search_shared=True, all_returned_chains_DNA_legal=True,
        gate_conditions=conditions, mechanism_gate_for_full_DEV_preparation=all(conditions.values()),
        optimizer_steps_this_run=0, new_head_parameters=0, new_backbone_forward=False,
        checkpoint_written=False, reference_used_in_decoding=False,
        TRAIN_reference_used_after_prediction=True, DEV_used=False, test_setaria_used=False,
        feature_spec=R.FEATURE_SPEC, prediction_seconds=prediction_seconds,
        seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitations='Repeated-use selected TRAIN pressure diagnostic; not final annotation, F1, background, independent generalization or novelty. Mean convention not exact teacher objective: retained A5 min-intron normalizer differs.'))


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
