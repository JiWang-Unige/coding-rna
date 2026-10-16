#!/usr/bin/env python3
"""A27 all original missing TRAIN chains, forced prefixes and natural-hit controls."""
from collections import Counter, defaultdict
import json
import math
import time
import torch
import run_chain_competition as P
import analyze_true_prefix_train as T
import audit_completed_prefix_coverage as A24

R, C, ROOT = P.C.R, P.C, P.ROOT
OUT = ROOT/'outputs/M28-A27-MISSING-TRUE-PREFIX-DIAGNOSTIC'
SOURCE = A24.OUT
PROTOCOL = 'reports/M28-METHOD-RESTART/a27_missing_true_prefix_protocol.md'
STATES = ('grammar_unsupported', 'all_queries_locally_excluded',
          'some_query_passes_local_necessary_conditions')


def summarize(ref, profiles, cohort):
    if len(profiles) != 16 or [p['query_index'] for p in profiles] != list(range(16)):
        raise ValueError('All ordered original Q16 required')
    grammar = {p['grammar_supported'] for p in profiles}
    if len(grammar) != 1:
        raise ValueError('Geometry/DNA grammar cannot differ by query')
    supported = grammar.pop()
    survivors = [p['query_index'] for p in profiles if p.get('local_width4_necessary_conditions_passed')]
    if not supported:
        if survivors:
            raise ValueError('Unsupported geometry is not locally passed')
        mechanism = 'grammar_unsupported'
    elif not survivors:
        mechanism = 'all_queries_locally_excluded'
    else:
        mechanism = 'some_query_passes_local_necessary_conditions'
    if cohort == 'natural_hit_control' and (not supported or not survivors):
        raise ValueError('Saved natural beam hit failed necessary local conditions')
    max_ranks = [max(max(e['entry_legal_rank'], e['end_same_kind_rank']) for e in p['events'])
                 for p in profiles] if supported else []
    return dict(**ref, cohort=cohort, grammar_supported=supported, mechanism=mechanism,
        local_survivor_queries=survivors,
        best_query_max_necessary_local_rank=min(max_ranks) if max_ranks else None,
        first_local_exclusion_by_query={str(p['query_index']): p.get('first_local_exclusion') for p in profiles},
        query_first_exclusion_stages=dict(Counter(
            p['first_local_exclusion']['stage'] for p in profiles
            if supported and p['first_local_exclusion'] is not None)),
        unsupported_reasons=sorted({p['unsupported_reason'] for p in profiles if not p['grammar_supported']}))


def aggregate(rows):
    counts = Counter(r['mechanism'] for r in rows)
    return dict(reference_window_pairs=len(rows), mechanisms={s: counts[s] for s in STATES},
        total_local_survivor_queries=sum(len(r['local_survivor_queries']) for r in rows),
        query_first_exclusion_stages=dict(sum((Counter(r['query_first_exclusion_stages']) for r in rows), Counter())))


def cohort():
    draws, index, windows, originals, _ = P.load_cohort()
    if (SOURCE/'STATUS').read_text().strip() != 'COMPLETED' or (SOURCE/'JOBID').read_text().strip() != '13391716':
        raise ValueError('Closed original A24 source required')
    report = json.loads((SOURCE/'summary.json').read_text())
    if (not report['completed'] or report['source_generator'] != C.GENERATOR or
        report['windows'] != 1491 or report['pooled']['reference_window_pairs'] != 2916 or
        report['pooled']['natural_hits'] != 2380 or report['pooled']['missing_pairs'] != 536):
        raise ValueError('Exact original A24 positive population required')
    refs = R.rows(SOURCE/'references.jsonl')
    free = {r['window']['id']: r for r in R.rows(C.FIT/'free_pool.jsonl')}
    expected = [(ident, i) for ident, rec in originals.items()
                for i, row in enumerate(rec['candidates']) if row['label'] == 1]
    if [(r['window_id'], r['original_candidate_index']) for r in refs] != expected:
        raise ValueError('Original ordered positive identities changed')
    by_id = defaultdict(list)
    for ref in refs:
        ident = ref['window_id']; source = originals[ident]['candidates'][ref['original_candidate_index']]
        if (ref['species'], ref['chain'], ref['origin'], ref['exon_class']) != (
                source['species'], source['chain'], source['origin'], source['exon_class']):
            raise ValueError('Original positive identity/geometry changed')
        classification = A24.classify(ref['chain'], free[ident]['chains'])
        if classification != {k: ref[k] for k in classification}:
            raise ValueError('Saved completed-pool classification differs')
        by_id[ident].append(ref)
    missing = [r for r in refs if r['category'] != 'natural_hit']
    missing_ids = {r['window_id'] for r in missing}
    controls = {}
    for ident in originals:
        if ident not in missing_ids:
            continue
        eligible = [r for r in by_id[ident] if r['category'] == 'natural_hit' and
                    'beam4' in free[ident]['decoder_sources'][r['exact_complete_candidate_indices'][0]]]
        if eligible:
            controls[ident] = eligible[0]
    if len(missing) != 536 or len(missing_ids) != 396 or len(controls) != 324:
        raise ValueError('Fixed all-missing cohort or metadata-first beam-hit controls changed')
    selected = {ident: ([r for r in by_id[ident] if r['category'] != 'natural_hit']+
                        ([controls[ident]] if ident in controls else []))
                for ident in originals if ident in missing_ids}
    for ident in selected:
        w = windows[ident]
        if (w['split'] != 'train' or w['seqid'] not in A24.ALLOW[w['species']] or
            w['width'] != 24576 or w['valid_bases'] != 24576):
            raise ValueError('Only original six TRAIN chromosomes/full-width oriented windows')
    return index, windows, selected


def verify_profiles(ref, profiles):
    for p in profiles:
        if p['grammar_supported']:
            if len(p['events']) != len(ref['chain']):
                raise ValueError('Complete forced chain was truncated')
            score = sum(e['entry_logp']+e['end_logp']+e['eos_logp'] for e in p['events'])
            if not math.isclose(score, p['path_log_score'], rel_tol=1e-7, abs_tol=1e-3):
                raise ValueError('Forced-path component sum differs')
            first = None
            for e in p['events']:
                if list(ref['chain'][e['exon']-1]) != [e['start'], e['end']]:
                    raise ValueError('Saved reference geometry changed during probe')
                if first is None and e['entry_legal_rank'] > 4:
                    first = dict(exon=e['exon'], stage='start' if e['exon'] == 1 else 'next_start',
                                 rank=e['entry_legal_rank'])
                elif first is None and e['end_same_kind_rank'] > 4:
                    first = dict(exon=e['exon'], stage='end_'+e['end_kind'], rank=e['end_same_kind_rank'])
            if first != p['first_local_exclusion'] or (first is None) != p['local_width4_necessary_conditions_passed']:
                raise ValueError('Saved first exclusion or local necessary-condition flag differs')


@torch.inference_mode()
def main():
    began = time.perf_counter(); index, windows, selected = cohort()
    generator = C.load_generator()
    if any(p.requires_grad for p in generator.parameters()):
        raise ValueError('Frozen original A3 required')
    torch.cuda.reset_peak_memory_stats()
    results = []; profile_count = 0; exon_events = 0
    with (OUT/'query_profiles.jsonl').open('x') as log:
        for ordinal, (ident, refs) in enumerate(selected.items()):
            w = windows[ident]
            cached = torch.load(ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
            features, dna = R.Q.window_inputs(cached, w, 'cuda')
            sequence = cached['sequence']
            if len(sequence) != w['valid_bases']:
                raise ValueError('Original oriented full-width DNA required')
            hidden = generator(features, dna, [w['valid_bases']])['features'][0]
            window = generator.query.for_window(hidden, w['valid_bases'])
            for ref in refs:
                kind = 'missing' if ref['category'] != 'natural_hit' else 'natural_hit_control'
                profiles = [T.probe_query(window, sequence, ref['chain'], q) for q in range(16)]
                verify_profiles(ref, profiles)
                results.append(summarize(ref, profiles, kind))
                for profile in profiles:
                    log.write(json.dumps(dict(window_id=ident,
                        original_candidate_index=ref['original_candidate_index'], cohort=kind,
                        reference_conditioned=True, **profile), allow_nan=False)+'\n')
                    profile_count += 1; exon_events += len(profile['events'])
            log.flush()
            if (ordinal+1) % 8 == 0:
                print(json.dumps(dict(windows=ordinal+1, total_windows=396, reference_pairs=len(results),
                                      query_profiles=profile_count, seconds=time.perf_counter()-began)), flush=True)
            del cached, features, dna, hidden, window
    if len(results) != 860 or profile_count != 13760 or any(p.grad is not None for p in generator.parameters()):
        raise ValueError('Incomplete original missing/control coverage or generator received gradient')
    with (OUT/'references.jsonl').open('x') as log:
        for rec in results:
            log.write(json.dumps(rec, allow_nan=False)+'\n')
    # Reconcile actual saved profiles independently of the in-memory summaries.
    grouped = defaultdict(list)
    for p in R.rows(OUT/'query_profiles.jsonl'):
        grouped[(p['window_id'], p['original_candidate_index'])].append(p)
    if set(grouped) != {(r['window_id'], r['original_candidate_index']) for r in results}:
        raise ValueError('Saved probe dropped or added an original identity')
    results_by_id = {(r['window_id'], r['original_candidate_index']): r for r in results}
    for ident, refs in selected.items():
        for ref in refs:
            kind = 'missing' if ref['category'] != 'natural_hit' else 'natural_hit_control'
            profiles = grouped[(ident, ref['original_candidate_index'])]
            verify_profiles(ref, profiles)
            recomputed = summarize(ref, profiles, kind)
            if recomputed != results_by_id[(ident, ref['original_candidate_index'])]:
                raise ValueError('Actual saved profiles and reference classification differ')
    if R.rows(OUT/'references.jsonl') != results:
        raise ValueError('Saved reference records differ')
    missing = [r for r in results if r['cohort'] == 'missing']
    controls = [r for r in results if r['cohort'] == 'natural_hit_control']
    report = dict(completed=True, experiment='M28-A27-MISSING-TRUE-PREFIX-DIAGNOSTIC', protocol=PROTOCOL,
        source_A24_job='13391716', generator_checkpoint=C.GENERATOR, limits=T.LIMITS,
        source_TRAIN_windows=1491, original_reference_window_pairs=2916,
        windows=396, missing_reference_window_pairs=536, control_reference_window_pairs=324,
        query_profiles=profile_count, exon_events=exon_events,
        missing=aggregate(missing), natural_hit_controls=aggregate(controls),
        per_species={s: aggregate([r for r in missing if r['species'] == s]) for s in A24.SPECIES},
        by_exon_class={c: aggregate([r for r in missing if r['exon_class'] == c])
                       for c in ('one', 'two_to_four', 'five_or_more')},
        by_completed_pool_category={c: aggregate([r for r in missing if r['category'] == c])
                                    for c in A24.CATEGORIES if c != 'natural_hit'},
        actual_saved_profiles_geometry_components_first_exclusions_and_summaries_reconciled=True,
        actual_original_missing_and_control_identity_verified=True,
        all_saved_natural_beam_controls_pass_necessary_local_conditions=True,
        reference_conditioned=True, natural_predictions_reused_without_change=True,
        new_candidate_search=False, optimizer_steps_this_run=0, checkpoint_written=False,
        new_backbone_forward=False, DEV_forward_used=False, test_setaria_used=False,
        seconds=time.perf_counter()-began, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='All original missing TRAIN candidate/window pairs, not independent gene recall, active-search trace, causal proof, new decoder, performance gate or novelty.')
    R.write_summary(OUT, report)


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
