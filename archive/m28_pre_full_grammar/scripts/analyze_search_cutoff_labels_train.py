#!/usr/bin/env python3
"""A29 annotate actual A28 cutoffs before defining a new search-training loss.

Set-level primary compatibility is deliberate: unlike A9's matched-slot label,
another eligible primary is positive here. This is diagnostic annotation only.
"""
from collections import Counter, defaultdict
import json
import time
import torch
import run_query as Q
from src.m28.labels import build_targets
from src.m28.query_prefix_value import START, NEXT_START, DONOR, STOP
from src.m28.query_prefix_value_training import compatible

ROOT = Q.ROOT
SOURCE = ROOT/'outputs/M28-A28-BEAM-SEARCH-TRACE-DIAGNOSTIC'
OUT = ROOT/'outputs/M28-A29-SEARCH-CUTOFF-LABELS'
PROTOCOL = 'reports/M28-METHOD-RESTART/a29_search_cutoff_protocol.md'
GROUPS = ('primary_compatible', 'unresolved_or_protected', 'known_nonmatch')
STAGES = ('global_entry_pool_pruned', 'continuation_pool_pruned',
          'completion_pool_pruned', 'retained_beam_score_bound')


def rows(path):
    with path.open() as handle:
        return [json.loads(line) for line in handle if line.strip()]


def node_at_loss(loss, target=False):
    stage = loss['stage']
    if stage not in STAGES:
        raise ValueError('Only original supported A28 first-loss stages required')
    row = loss['target'] if target else loss['completion_cutoff' if stage == STAGES[-1] else 'cutoff']
    chain = tuple(tuple(x) for x in row['chain'])
    if stage == 'global_entry_pool_pruned':
        entry = row['entry']; event = NEXT_START if chain else START
    else:
        entry = None
        event = DONOR if stage == 'continuation_pool_pruned' or (target and stage == STAGES[-1]) else STOP
    return dict(chain=chain, entry=entry, event=event)


def classify(node, targets, reference):
    primary = [i for i,c in enumerate(targets['positive_chains']) if compatible(node,c)]
    protected = any(compatible(node,c) for c in targets['protected_chains']
                    if tuple(tuple(x) for x in c) not in {tuple(tuple(x) for x in p) for p in targets['positive_chains']})
    chain = node['chain']; entry = node['entry']
    a = chain[0][0] if chain else entry
    b = entry+1 if entry is not None else chain[-1][1]
    if not 0 <= a < b <= len(targets['callable_bases']):
        raise ValueError('Actual cutoff span outside original valid oriented window')
    unknown = any(x < b and y > a for x,y in targets['unknown_chain_spans'])
    noncallable = not bool(targets['callable_bases'][a:b].all())
    label = 1 if primary else (-1 if protected or unknown or noncallable else 0)
    group = GROUPS[0] if label == 1 else (GROUPS[1] if label == -1 else GROUPS[2])
    return dict(label=label, group=group, matching_primary_indices=primary,
        target_primary_compatible=compatible(node,reference), protected_compatible=protected,
        overlaps_unknown_chain_span=unknown, contains_noncallable_span=noncallable,
        inspected_span=[a,b])


def aggregate(records):
    counts = Counter(r['cutoff_annotation']['group'] for r in records)
    return dict(cutoff_query_records=len(records), groups={g:counts[g] for g in GROUPS},
        other_primary_compatible=sum(r['cutoff_annotation']['label'] == 1 and
            not r['cutoff_annotation']['target_primary_compatible'] for r in records))


def source():
    if (SOURCE/'STATUS').read_text().strip() != 'COMPLETED' or (SOURCE/'JOBID').read_text().strip() != '13392879':
        raise ValueError('Only the completed full original A28 source is allowed')
    summary = json.loads((SOURCE/'summary.json').read_text())
    flags = ('actual_each_query_original_search_reconciled', 'actual_original_beam_pool_subset_verified',
        'actual_forced_prefix_scores_and_geometry_verified', 'actual_natural_control_hits_verified',
        'actual_saved_profiles_and_reference_summaries_reconciled')
    if (not summary['completed'] or summary['experiment'] != 'M28-A28-BEAM-SEARCH-TRACE-FULL' or
        summary['source_A27_job'] != '13392411' or summary['preflight_job'] != '13392831' or
        summary['generator_checkpoint'] != 'outputs/M28-A3-FIRST-FIT/step_004608.pt' or
        summary['limits'] != dict(beam_width=4, completed_limit_per_query=4, max_exons=128) or
        summary['full_domain'] != dict(windows=352, reference_window_pairs=545,
            missing_reference_window_pairs=221, control_reference_window_pairs=324,
            reference_query_pairs=4247, window_query_replays=3363) or
        any(not summary[k] for k in flags) or summary['optimizer_steps_this_run'] != 0 or
        summary['reference_used_in_search'] or summary['search_budget_changed'] or
        summary['DEV_forward_used'] or summary['test_setaria_used']):
        raise ValueError('Exact closed A28 original-search and saved-profile identity required')
    refs = rows(SOURCE/'references.jsonl'); profiles = rows(SOURCE/'query_profiles.jsonl')
    by_key = {(r['window_id'], r['original_candidate_index']):r for r in refs}
    expected = {(r['window_id'], r['original_candidate_index'], q) for r in refs for q in r['local_survivor_queries']}
    if (len(refs) != 545 or len(by_key) != 545 or len(profiles) != 4247 or
        len({r['window_id'] for r in refs}) != 352 or
        {(p['window_id'],p['original_candidate_index'],p['query_index']) for p in profiles} != expected):
        raise ValueError('Actual A28 reference/query identities or coverage differ')
    grouped = defaultdict(list)
    for ordinal,p in enumerate(profiles):
        key = (p['window_id'],p['original_candidate_index']); r = by_key[key]
        if p['cohort'] != r['cohort'] or (p['first_loss'] is None) != p['returned']:
            raise ValueError('Saved A28 profile and reference cohort/first-loss fields differ')
        if not p['returned']:
            grouped[p['window_id']].append((ordinal,p))
    if sum(p['cohort'] == 'missing' for values in grouped.values() for _,p in values) != 1708:
        raise ValueError('Every original locally-passing missing reference/query must remain a first-loss record')
    return refs, by_key, dict(grouped)


def main():
    began = time.perf_counter(); refs, by_key, grouped = source()
    _, index, windows, catalogs = Q.load_train_inputs()
    records = []; window_targets = {}; used = {r['window_id'] for r in refs}
    for ordinal,ident in enumerate(dict.fromkeys(r['window_id'] for r in refs)):
        w = windows[ident]
        if w['split'] != 'train' or w['seqid'] not in Q.ALLOW[w['species']] or w['valid_bases'] != 24576:
            raise ValueError('Only original six TRAIN chromosomes/full-width oriented windows')
        cached = torch.load(ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
        Q.validate_cached(cached,w)
        targets = build_targets(w,catalogs[w['species']],cached['sequence'])
        # Retain only small annotation objects for saved-record reconciliation.
        window_targets[ident] = {k:targets[k] for k in
            ('positive_chains','protected_chains','unknown_chain_spans','callable_bases')}
        for r in refs:
            if r['window_id'] == ident and tuple(tuple(x) for x in r['chain']) not in {
                    tuple(tuple(x) for x in c) for c in targets['positive_chains']}:
                raise ValueError('A28 source reference is not an original eligible R4 primary')
        for profile_index,p in grouped.get(ident,[]):
            r = by_key[(ident,p['original_candidate_index'])]
            loss = p['first_loss']; node = node_at_loss(loss)
            label = classify(node,targets,r['chain'])
            target = classify(node_at_loss(loss,True),targets,r['chain'])
            if target['label'] != 1 or not target['target_primary_compatible']:
                raise ValueError('Observed A28 target does not remain its original positive prefix')
            records.append(dict(source_profile_index=profile_index, window_id=ident,
                original_candidate_index=p['original_candidate_index'], query_index=p['query_index'],
                cohort=p['cohort'], species=r['species'], exon_class=r['exon_class'],
                first_loss_stage=loss['stage'], first_loss_exon=loss['exon'],
                cutoff_node=node, cutoff_annotation=label,
                target_primary_compatibility_verified=True))
        del cached,targets
        if (ordinal+1)%32 == 0:
            print(json.dumps(dict(windows=ordinal+1,total_windows=352,cutoff_records=len(records),
                                  seconds=time.perf_counter()-began)),flush=True)
    source_profiles = rows(SOURCE/'query_profiles.jsonl')
    if len(records) != sum(not p['returned'] for p in source_profiles) or len(used) != 352:
        raise ValueError('First-loss record coverage incomplete')
    with (OUT/'cutoff_annotations.jsonl').open('x') as log:
        for r in records:
            log.write(json.dumps(r,allow_nan=False)+'\n')
    saved = rows(OUT/'cutoff_annotations.jsonl')
    if saved != json.loads(json.dumps(records)):
        raise ValueError('Actual saved cutoff annotation records differ')
    grouped_refs = defaultdict(list)
    for r in saved:
        p = source_profiles[r['source_profile_index']]
        ref = by_key[(r['window_id'],r['original_candidate_index'])]
        node = node_at_loss(p['first_loss'])
        if ((r['window_id'],r['original_candidate_index'],r['query_index'],r['cohort']) != (
                p['window_id'],p['original_candidate_index'],p['query_index'],p['cohort']) or
            r['cutoff_node'] != json.loads(json.dumps(node)) or
            r['cutoff_annotation'] != classify(node,window_targets[r['window_id']],ref['chain'])):
            raise ValueError('Saved classification differs from actual original profile or R4 annotation')
        grouped_refs[(r['window_id'],r['original_candidate_index'])].append(r)
    summaries = []
    for ref in refs:
        rs = grouped_refs[(ref['window_id'],ref['original_candidate_index'])]
        negative = [r['query_index'] for r in rs if r['cutoff_annotation']['label'] == 0]
        summaries.append(dict(window_id=ref['window_id'], original_candidate_index=ref['original_candidate_index'],
            cohort=ref['cohort'], species=ref['species'], exon_class=ref['exon_class'],
            local_survivor_queries=ref['local_survivor_queries'], returned_queries=ref['returned_queries'],
            first_loss_queries=len(rs), known_nonmatch_cutoff_queries=negative,
            cutoff_groups=aggregate(rs)['groups']))
    with (OUT/'references.jsonl').open('x') as log:
        for r in summaries:
            log.write(json.dumps(r,allow_nan=False)+'\n')
    if rows(OUT/'references.jsonl') != summaries:
        raise ValueError('Saved reference coverage or cutoff summaries differ')
    missing = [r for r in saved if r['cohort'] == 'missing']
    controls = [r for r in saved if r['cohort'] == 'natural_hit_control']
    mr = [r for r in summaries if r['cohort'] == 'missing']
    report = dict(completed=True,experiment='M28-A29-SEARCH-CUTOFF-LABELS',protocol=PROTOCOL,
        source_A28_job='13392879',source_A27_job='13392411',source_query_profiles=4247,
        windows=352,reference_window_pairs=545,missing_reference_window_pairs=221,
        control_reference_window_pairs=324,cutoff_records=len(saved),missing_cutoff_records=len(missing),
        missing=aggregate(missing),natural_hit_controls=aggregate(controls),
        missing_reference_pairs_with_any_known_nonmatch_cutoff=sum(bool(r['known_nonmatch_cutoff_queries']) for r in mr),
        missing_reference_pairs_with_all_known_nonmatch_cutoffs=sum(
            len(r['known_nonmatch_cutoff_queries']) == len(r['local_survivor_queries']) for r in mr),
        per_species={s:aggregate([r for r in missing if r['species'] == s]) for s in Q.ALLOW},
        per_stage={s:aggregate([r for r in missing if r['first_loss_stage'] == s]) for s in STAGES},
        by_exon_class={c:aggregate([r for r in missing if r['exon_class'] == c])
                      for c in ('one','two_to_four','five_or_more')},
        actual_original_reference_and_profile_identities_verified=True,
        actual_original_R4_primary_geometry_and_protection_verified=True,
        actual_saved_cutoff_annotations_and_reference_summaries_reconciled=True,
        query_assignment_used=False,reference_used_for_annotation_only=True,
        model_or_backbone_forward=False,new_candidate_search=False,optimizer_steps_this_run=0,
        checkpoint_written=False,DEV_used=False,test_setaria_used=False,
        seconds=time.perf_counter()-began,
        limitation='One actual retained cutoff per first-loss query, not all outranking paths, query-slot labels, biological false-gene proof, a training objective, causal attribution or independent gene recall.')
    with (OUT/'summary.json').open('x') as log:
        json.dump(report,log,indent=2,allow_nan=False);log.write('\n')
    print(json.dumps(report,allow_nan=False),flush=True)


if __name__ == '__main__':
    torch.set_num_threads(2); main()
