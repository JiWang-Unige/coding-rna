#!/usr/bin/env python3
"""A24 saved complete-chain prefix coverage; no search or model inference."""
from collections import Counter
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT/'outputs/M28-A24-COMPLETED-PREFIX-COVERAGE'
PROTOCOL = 'reports/M28-METHOD-RESTART/a24_completed_prefix_protocol.md'
SOURCES = {
    'A6': ('M28-A6-CHAIN-QUALITY-FIT', '13358977'),
    'A12': ('M28-A12-QUALITY-TRAIN-DIAGNOSTIC', '13360650'),
    'A22': ('M28-A22-CHAIN-DISPLACEMENT-AUDIT', '13391139'),
}
SPECIES = ('arabidopsis_thaliana', 'oryza_sativa')
ALLOW = {'arabidopsis_thaliana': ('NC_003070.9', 'NC_003071.7', 'NC_003076.8'),
         'oryza_sativa': ('NC_089043.1', 'NC_089044.1', 'NC_089046.1')}
GEOMETRY = ('id', 'species', 'split', 'seqid', 'strand', 'start', 'end', 'width', 'valid_bases')
CATEGORIES = ('natural_hit', 'completed_start_absent', 'first_end_absent',
              'target_proper_prefix_of_completed', 'completed_prefix_only',
              'continuation_start_absent', 'continuation_end_absent')


def rows(path):
    with path.open() as handle:
        return [json.loads(line) for line in handle]


def key(chain):
    return tuple(tuple(exon) for exon in chain)


def common_prefix(a, b):
    count = 0
    for x, y in zip(a, b):
        if x != y:
            break
        count += 1
    return count


def classify(target, chains):
    """Classify only saved completed chains, never a decoder's active histories."""
    target, chains = key(target), [key(c) for c in chains]
    same_start = [i for i, c in enumerate(chains) if c[0][0] == target[0][0]]
    exact = [i for i, c in enumerate(chains) if c == target]
    prefix = [common_prefix(target, c) for c in chains]
    depth = max(prefix, default=0)
    frontier = [i for i in same_start if prefix[i] == depth]
    terminated = [i for i in frontier if len(chains[i]) == depth]
    continuing = [i for i in frontier if len(chains[i]) > depth]
    next_start = ([i for i in continuing if chains[i][depth][0] == target[depth][0]]
                  if depth < len(target) else [])
    if exact:
        category = 'natural_hit'
    elif not same_start:
        category = 'completed_start_absent'
    elif depth == 0:
        category = 'first_end_absent'
    elif depth == len(target):
        category = 'target_proper_prefix_of_completed'
    elif not continuing:
        category = 'completed_prefix_only'
    elif next_start:
        category = 'continuation_end_absent'
    else:
        category = 'continuation_start_absent'
    return dict(category=category, exact_complete_candidate_indices=exact,
        same_start_completed_candidate_indices=same_start, max_exact_prefix_exons=depth,
        deepest_prefix_candidate_indices=frontier,
        deepest_prefix_terminated_candidate_indices=terminated,
        deepest_prefix_continuing_candidate_indices=continuing,
        deepest_prefix_matching_next_start_candidate_indices=next_start,
        empty_completed_pool=not chains)


def aggregate(records):
    counts = Counter(r['category'] for r in records)
    return dict(reference_window_pairs=len(records), natural_hits=counts['natural_hit'],
        missing_pairs=len(records)-counts['natural_hit'],
        categories={c: counts[c] for c in CATEGORIES},
        missing_pairs_with_empty_completed_pool=sum(r['empty_completed_pool'] for r in records if r['category'] != 'natural_hit'),
        missing_pairs_with_both_terminated_and_continuing_deepest_prefix=sum(
            bool(r['deepest_prefix_terminated_candidate_indices']) and bool(r['deepest_prefix_continuing_candidate_indices'])
            for r in records if r['category'] != 'natural_hit'),
        missing_max_exact_prefix_exon_histogram=dict(sorted(Counter(
            r['max_exact_prefix_exons'] for r in records if r['category'] != 'natural_hit').items())))


def main():
    began = time.perf_counter(); folders = {}; summaries = {}
    for name, (directory, job) in SOURCES.items():
        folder = ROOT/'outputs'/directory; folders[name] = folder
        if (folder/'STATUS').read_text().strip() != 'COMPLETED' or (folder/'JOBID').read_text().strip() != job:
            raise ValueError('Incomplete or different original source '+name)
        summaries[name] = json.loads((folder/'summary.json').read_text())
        if not summaries[name]['completed']:
            raise ValueError('Source not completed '+name)
    if (summaries['A12']['source_generator'] != 'outputs/M28-A3-FIRST-FIT/step_004608.pt' or
        summaries['A12']['source_head'] != summaries['A6']['final_checkpoint'] or
        summaries['A6']['configuration']['generator_checkpoint'] != summaries['A12']['source_generator'] or
        not summaries['A22']['source_identity_and_actual_selection_reconciled'] or
        summaries['A22']['original_injected_positive_metadata'] != 536):
        raise ValueError('Only original closed A3 natural pool and A12 positive metadata required')
    raw = rows(folders['A6']/'free_pool.jsonl'); labelled = rows(folders['A12']/'windows.jsonl')
    ids = [r['window']['id'] for r in raw]
    if len(ids) != 1491 or len(set(ids)) != 1491 or [r['window']['id'] for r in labelled] != ids:
        raise ValueError('Original ordered 1491 unique TRAIN windows required')
    references, windows = [], []
    for free, original in zip(raw, labelled):
        w = original['window']; chains = free['chains']
        if (w['split'] != 'train' or w['seqid'] not in ALLOW[w['species']] or
            free['reference_used_in_pool_generation'] or
            any(free['window'][k] != w[k] for k in GEOMETRY)):
            raise ValueError('Allowed TRAIN geometry and reference-free original pool required')
        candidates = original['candidates']; natural = [r for r in candidates if r['origin'] == 'free']
        if ([r['chain'] for r in natural] != chains or len({key(c) for c in chains}) != len(chains) or
            len({key(r['chain']) for r in candidates}) != len(candidates) or
            any(r['origin'] not in ('free', 'reference_injected') or
                (r['origin'] == 'reference_injected' and r['label'] != 1) for r in candidates)):
            raise ValueError('Original natural chain identity or injected positive semantics changed')
        local = []
        for index, rec in enumerate(candidates):
            if rec['label'] != 1:
                continue
            classification = classify(rec['chain'], chains)
            if (classification['category'] == 'natural_hit') != (rec['origin'] == 'free'):
                raise ValueError('Completion identity differs from saved positive origin')
            expected_class = 'one' if len(rec['chain']) == 1 else 'two_to_four' if len(rec['chain']) <= 4 else 'five_or_more'
            if rec['exon_class'] != expected_class:
                raise ValueError('Original positive exon class changed')
            result = dict(window_id=w['id'], species=w['species'], original_candidate_index=index,
                chain=rec['chain'], origin=rec['origin'], exon_class=expected_class, **classification)
            references.append(result); local.append(result)
        windows.append(dict(window={k: w[k] for k in GEOMETRY}, natural_candidates=len(chains),
                            positive_metadata_pairs=len(local), counts=aggregate(local)))
    totals = aggregate(references)
    if (sum(len(r['chains']) for r in raw) != 16321 or
        sum(len(r['candidates']) for r in labelled) != summaries['A12']['unique_candidate_window_pairs'] or
        totals['reference_window_pairs'] != 2916 or totals['natural_hits'] != 2380 or totals['missing_pairs'] != 536):
        raise ValueError('Original frozen natural and positive population counts differ')
    for name, records in (('windows.jsonl', windows), ('references.jsonl', references)):
        with (OUT/name).open('x') as handle:
            for row in records:
                handle.write(json.dumps(row, allow_nan=False)+'\n')
    per_species = {s: aggregate([r for r in references if r['species'] == s]) for s in SPECIES}
    report = dict(completed=True, experiment='M28-A24-COMPLETED-PREFIX-COVERAGE', protocol=PROTOCOL,
        source_jobs={k: v[1] for k, v in SOURCES.items()}, source_generator=summaries['A12']['source_generator'],
        windows=1491, natural_candidate_window_pairs=16321,
        actual_original_pool_positive_origin_and_exon_class_identity_verified=True,
        pooled=totals, per_species=per_species,
        by_exon_class={c: aggregate([r for r in references if r['exon_class'] == c])
                       for c in ('one', 'two_to_four', 'five_or_more')},
        per_species_by_exon_class={s: {c: aggregate([r for r in references if r['species'] == s and r['exon_class'] == c])
            for c in ('one', 'two_to_four', 'five_or_more')} for s in SPECIES},
        all_chain_geometries_preserved=True, optimizer_steps_this_run=0, new_candidate_search=False,
        model_or_backbone_forward=False, checkpoint_read_or_written=False, catalog_or_DNA_read=False,
        original_predictions_labels_and_thresholds_unchanged=True, DEV_forward_used=False, test_setaria_used=False,
        seconds=time.perf_counter()-began,
        limitation='Saved completed-pool morphology, not active-prefix reachability, local ranks, causal attribution or independent gene F1.')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(dict(completed=True, windows=1491, pooled=totals, per_species=per_species,
                         seconds=report['seconds']), allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
