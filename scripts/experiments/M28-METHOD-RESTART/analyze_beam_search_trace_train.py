#!/usr/bin/env python3
"""A28 actual frozen beam loss for every A27 locally-passing missing/control pair."""
import argparse
from collections import Counter, defaultdict
import json
import math
import time
import torch
import analyze_missing_true_prefix_train as A27
from trace_frozen_query_beam import trace_query, original_query, reconcile

R, C, ROOT = A27.R, A27.C, A27.ROOT
PREFLIGHT = ROOT/'outputs/M28-A28-BEAM-SEARCH-TRACE-PREFLIGHT-REPAIR-R1'
OUT = ROOT/'outputs/M28-A28-BEAM-SEARCH-TRACE-DIAGNOSTIC'
PROTOCOL = 'reports/M28-METHOD-RESTART/a28_beam_search_trace_protocol.md'


def cohort():
    if ((A27.OUT/'STATUS').read_text().strip() != 'COMPLETED' or
            (A27.OUT/'JOBID').read_text().strip() != '13392411'):
        raise ValueError('Closed original A27 required')
    summary = json.loads((A27.OUT/'summary.json').read_text())
    if (not summary['completed'] or summary['generator_checkpoint'] != C.GENERATOR or
            summary['query_profiles'] != 13760 or summary['windows'] != 396 or
            not summary['actual_saved_profiles_geometry_components_first_exclusions_and_summaries_reconciled'] or
            not summary['all_saved_natural_beam_controls_pass_necessary_local_conditions'] or
            summary['optimizer_steps_this_run'] != 0 or summary['DEV_forward_used'] or summary['test_setaria_used']):
        raise ValueError('Exact full A27 original-source and saved-profile verification required')
    index, windows, original = A27.cohort()
    rows = R.rows(A27.OUT/'references.jsonl')
    expected = [(ident, r['original_candidate_index']) for ident, refs in original.items() for r in refs]
    if [(r['window_id'], r['original_candidate_index']) for r in rows] != expected:
        raise ValueError('All original A27 reference/window identities and order required')
    sources = {(ident, r['original_candidate_index']):r for ident, refs in original.items() for r in refs}
    groups = defaultdict(list)
    for r in rows:
        source = sources[(r['window_id'], r['original_candidate_index'])]
        if any(r[k] != v for k,v in source.items()) or not r['grammar_supported']:
            raise ValueError('Original supported A27 source geometry changed')
        qs = r['local_survivor_queries']
        if qs != sorted(set(qs)) or any(q not in range(16) for q in qs):
            raise ValueError('Original local survivor queries required')
        if qs:
            groups[r['window_id']].append(r)
    selected = [r for refs in groups.values() for r in refs]
    kinds = Counter(r['cohort'] for r in selected)
    if kinds != Counter(missing=221, natural_hit_control=324) or sum(len(r['local_survivor_queries']) for r in selected) != 4247:
        raise ValueError('All 221 locally-passing missing pairs and 324 controls required')
    profiles = {}
    for p in R.rows(A27.OUT/'query_profiles.jsonl'):
        if p.get('local_width4_necessary_conditions_passed'):
            key = (p['window_id'], p['original_candidate_index'], p['query_index'])
            if key in profiles:
                raise ValueError('Duplicate original forced-path profile')
            profiles[key] = p
    keys = {(r['window_id'], r['original_candidate_index'], q)
            for r in selected for q in r['local_survivor_queries']}
    if set(profiles) != keys:
        raise ValueError('Saved A27 local survivor profiles and selected identities differ')
    free = {r['window']['id']:r for r in R.rows(C.FIT/'free_pool.jsonl')}
    return index, windows, dict(groups), profiles, free


def verify_loss(ref, forced, trace):
    chain = tuple(tuple(x) for x in ref['chain'])
    if trace['returned']:
        if trace['first_loss'] is not None:
            raise ValueError('Returned path has an alleged loss')
        return
    loss = trace['first_loss']; stage = loss['stage']; exon = loss['exon']
    if stage.startswith('local_'):
        raise ValueError('A27 locally-passing query failed the same original local rank during replay')
    if stage == 'global_entry_pool_pruned':
        expected_chain = chain[:exon-1]
        expected_score = sum(e['entry_logp']+e['end_logp']+e['eos_logp'] for e in forced['events'][:exon-1])
        expected_score += forced['events'][exon-1]['entry_logp']
        if loss['target']['entry'] != chain[exon-1][0]:
            raise ValueError('Wrong reference entry observed')
    elif stage == 'completion_pool_pruned':
        expected_chain = chain; expected_score = forced['path_log_score']
    elif stage in ('continuation_pool_pruned', 'retained_beam_score_bound', 'exon_cap'):
        expected_chain = chain[:exon]
        expected_score = sum(e['entry_logp']+e['end_logp']+e['eos_logp'] for e in forced['events'][:exon])
    else:
        raise ValueError('Unexpected first-loss category for a supported local survivor')
    if (tuple(tuple(x) for x in loss['target']['chain']) != expected_chain or
            not math.isclose(loss['target']['path_log_score'], expected_score, rel_tol=0., abs_tol=1e-3)):
        raise ValueError('Actual retained-prefix score or geometry differs from the frozen forced profile')
    if stage.endswith('_pruned') and loss['target_rank'] <= loss['retained_size']:
        raise ValueError('Allegedly pruned path was within the retained pool')


def aggregate(rows):
    stages = sum((Counter(r['query_first_loss_stages']) for r in rows), Counter())
    return dict(reference_window_pairs=len(rows), traced_reference_query_pairs=sum(r['traced_queries'] for r in rows),
        returned_reference_window_pairs=sum(bool(r['returned_queries']) for r in rows),
        query_first_loss_stages=dict(stages),
        reference_pairs_with_any_loss_stage={s:sum(s in r['query_first_loss_stages'] for r in rows) for s in stages})


@torch.inference_mode()
def main(mode, preflight_job):
    began = time.perf_counter(); index, windows, groups, forced, free = cohort()
    full_stats = dict(windows=len(groups), reference_window_pairs=545, missing_reference_window_pairs=221,
        control_reference_window_pairs=324, reference_query_pairs=4247,
        window_query_replays=sum(len({q for r in refs for q in r['local_survivor_queries']}) for refs in groups.values()))
    if mode == 'preflight':
        ids = []
        for species in A27.A24.SPECIES:
            local = [ident for ident in groups if windows[ident]['species'] == species]
            ids.extend(local[:4])
        if len(ids) != 8:
            raise ValueError('First four original eligible windows per species required')
        groups = {ident:groups[ident] for ident in ids}; out = PREFLIGHT
    else:
        if not preflight_job or (PREFLIGHT/'STATUS').read_text().strip() != 'COMPLETED' or (PREFLIGHT/'JOBID').read_text().strip() != preflight_job:
            raise ValueError('Actual completed preflight job identity required')
        report = json.loads((PREFLIGHT/'summary.json').read_text())
        if (not report['completed'] or report['windows'] != 8 or report['full_domain'] != full_stats or
                not report['actual_each_query_original_search_reconciled'] or
                not report['actual_saved_profiles_and_reference_summaries_reconciled']):
            raise ValueError('Actual eight-window replay and unchanged full-domain preflight required')
        out = OUT
    generator = C.load_generator()
    if any(p.requires_grad for p in generator.parameters()):
        raise ValueError('Frozen original A3 required')
    torch.cuda.reset_peak_memory_stats()
    profiles = []; natural_records = []; results = []; query_count = 0
    with (out/'query_profiles.jsonl').open('x') as plog, (out/'natural_queries.jsonl').open('x') as nlog:
        for ordinal, (ident, refs) in enumerate(groups.items()):
            w = windows[ident]
            cached = torch.load(ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
            features, dna = R.Q.window_inputs(cached, w, 'cuda')
            sequence = cached['sequence'].upper()
            hidden = generator(features, dna, [w['valid_bases']])['features'][0]
            window = generator.query.for_window(hidden, w['valid_bases'])
            original_beam = {tuple(tuple(x) for x in chain) for chain, sources in
                zip(free[ident]['chains'], free[ident]['decoder_sources']) if 'beam4' in sources}
            by_ref = defaultdict(list)
            queries = sorted({q for r in refs for q in r['local_survivor_queries']})
            for q in queries:
                qr = {r['original_candidate_index']:r['chain'] for r in refs if q in r['local_survivor_queries']}
                native, observed = trace_query(window, sequence, q, qr)
                original = original_query(window, sequence, q)
                reconcile(original, native)
                if any(r['chain'] not in original_beam for r in native['completed']):
                    raise ValueError('Actual replay completion not in the original immutable beam pool')
                nr = dict(window_id=ident, original_search_reconciled=True, **native)
                nlog.write(json.dumps(nr, allow_nan=False)+'\n'); natural_records.append(nr); query_count += 1
                for r in refs:
                    k = r['original_candidate_index']
                    if k not in observed:
                        continue
                    verify_loss(r, forced[(ident,k,q)], observed[k])
                    row = dict(window_id=ident, original_candidate_index=k, query_index=q,
                        cohort=r['cohort'], reference_used_for_observation_only=True,
                        original_search_reconciled=True, forced_prefix_score_and_geometry_verified=True, **observed[k])
                    plog.write(json.dumps(row, allow_nan=False)+'\n'); profiles.append(row); by_ref[k].append(row)
            for r in refs:
                rp = by_ref[r['original_candidate_index']]
                if [p['query_index'] for p in rp] != r['local_survivor_queries']:
                    raise ValueError('A local survivor was omitted or added')
                returned = [p['query_index'] for p in rp if p['returned']]
                if bool(returned) != (r['cohort'] == 'natural_hit_control'):
                    raise ValueError('Replay contradicts original natural missing/control identity')
                results.append(dict(window_id=ident, original_candidate_index=r['original_candidate_index'],
                    species=r['species'], chain=r['chain'], exon_class=r['exon_class'], cohort=r['cohort'],
                    A24_category=r['category'], local_survivor_queries=r['local_survivor_queries'],
                    traced_queries=len(rp), returned_queries=returned,
                    query_first_loss_stages=dict(Counter(p['first_loss']['stage'] for p in rp if not p['returned']))))
            plog.flush(); nlog.flush()
            print(json.dumps(dict(windows=ordinal+1, total_windows=len(groups),
                reference_query_pairs=len(profiles), window_query_replays=query_count,
                seconds=time.perf_counter()-began)), flush=True)
            del cached, features, dna, hidden, window
    if any(p.grad is not None for p in generator.parameters()):
        raise ValueError('Frozen generator received a gradient')
    with (out/'references.jsonl').open('x') as log:
        for r in results:
            log.write(json.dumps(r, allow_nan=False)+'\n')
    saved = R.rows(out/'query_profiles.jsonl'); saved_refs = R.rows(out/'references.jsonl')
    # JSON round-trip normalization preserves tuples as ordinary coordinate lists.
    if saved != json.loads(json.dumps(profiles)) or saved_refs != json.loads(json.dumps(results)) or R.rows(out/'natural_queries.jsonl') != json.loads(json.dumps(natural_records)):
        raise ValueError('Actual saved replay records differ')
    by_key = defaultdict(list)
    source_refs = {(ident,r['original_candidate_index']):r for ident, refs in groups.items() for r in refs}
    for p in saved:
        key = (p['window_id'],p['original_candidate_index'])
        verify_loss(source_refs[key], forced[(*key,p['query_index'])], p)
        by_key[key].append(p)
    if set(by_key) != set(source_refs):
        raise ValueError('Actual saved replay omitted a selected original reference')
    for r in saved_refs:
        ps = by_key[(r['window_id'],r['original_candidate_index'])]
        if (r['local_survivor_queries'] != [p['query_index'] for p in ps] or
                r['returned_queries'] != [p['query_index'] for p in ps if p['returned']] or
                r['query_first_loss_stages'] != dict(Counter(p['first_loss']['stage'] for p in ps if not p['returned']))):
            raise ValueError('Actual saved reference summaries and query observations differ')
    if mode == 'full' and (len(results) != 545 or len(profiles) != 4247 or query_count != full_stats['window_query_replays']):
        raise ValueError('Full original eligible reference/query coverage incomplete')
    missing = [r for r in saved_refs if r['cohort'] == 'missing']
    controls = [r for r in saved_refs if r['cohort'] == 'natural_hit_control']
    report = dict(completed=True, experiment='M28-A28-BEAM-SEARCH-TRACE-'+mode.upper(), protocol=PROTOCOL,
        source_A27_job='13392411', preflight_job=preflight_job, generator_checkpoint=C.GENERATOR,
        limits=A27.T.LIMITS, full_domain=full_stats, windows=len(groups),
        reference_window_pairs=len(results), query_profiles=len(profiles), window_query_replays=query_count,
        missing=aggregate(missing), natural_hit_controls=aggregate(controls),
        per_species={s:aggregate([r for r in missing if r['species'] == s]) for s in A27.A24.SPECIES},
        by_exon_class={c:aggregate([r for r in missing if r['exon_class'] == c]) for c in ('one','two_to_four','five_or_more')},
        actual_each_query_original_search_reconciled=True, actual_original_beam_pool_subset_verified=True,
        actual_forced_prefix_scores_and_geometry_verified=True, actual_natural_control_hits_verified=True,
        actual_saved_profiles_and_reference_summaries_reconciled=True,
        reference_used_in_search=False, reference_used_for_observation_only=True,
        new_decoder=False, search_budget_changed=False, optimizer_steps_this_run=0,
        checkpoint_written=False, new_backbone_forward=False, DEV_forward_used=False, test_setaria_used=False,
        seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='Observed first loss in the unchanged retained beam; not unpruned optimal search, causal training attribution, model improvement or independent gene recall.')
    R.write_summary(out, report)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('mode', choices=('preflight','full'))
    parser.add_argument('--preflight-job'); args = parser.parse_args()
    torch.set_num_threads(2); main(args.mode, args.preflight_job)
