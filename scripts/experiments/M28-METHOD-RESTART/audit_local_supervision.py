#!/usr/bin/env python3
"""TRAIN metadata/record census; no candidate relabeling or model forward."""
from collections import Counter
import json
import time
import run_query_known_mean as A
from src.m28.labels import oriented_parts

Q = A.Q
OUT = Q.ROOT/'outputs/M28-A19-LOCAL-SUPERVISION-DIAGNOSTIC'
PROTOCOL = 'reports/M28-METHOD-RESTART/a19_local_supervision_diagnostic.md'


def union(intervals, length):
    merged = []
    for a, b in sorted(intervals):
        assert 0 <= a < b <= length
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return merged


def outside_bins(guard, length):
    return [i for i in range(16) if not any(a < max(i*length//16+1, (i+1)*length//16)
            and b > i*length//16 for a, b in guard)]


def aggregate(rows):
    labels = Counter(); reasons = Counter()
    guarded = [r for r in rows if r['protection_reasons']]
    regional = [r for r in guarded if not r['noncallable_positions_unavailable']]
    for r in rows:
        labels.update(r['query_label_counts']); reasons.update(r['protection_reasons'])
    bases = sum(r['valid_bases'] for r in regional)
    outside = sum(r['outside_guard_span_bases'] for r in regional)
    return dict(windows_or_draws=len(rows), query_label_counts=dict(labels), protection_reason_counts=dict(reasons),
        guarded_windows_or_draws=len(guarded),
        guarded_with_noncallable_positions_unavailable=sum(r['noncallable_positions_unavailable'] for r in guarded),
        regional_guarded_windows_or_draws=len(regional), regional_valid_bases=bases,
        regional_guard_span_bases=bases-outside, regional_outside_guard_span_bases=outside,
        regional_outside_guard_span_fraction=outside/bases if bases else None,
        regional_initial_bins_outside_guard_span=sum(len(r['initial_bins_outside_guard_span']) for r in regional),
        regional_initial_bins_total=16*len(regional))


def main():
    assert union([[5, 10], [0, 6], [10, 12], [20, 21]], 32) == [[0, 12], [20, 21]]
    assert outside_bins([[2, 4]], 32) == [i for i in range(16) if i != 1]
    print('interval union and half-open initial-bin checks PASS', flush=True)
    began = time.perf_counter(); source = A.STAGE_DIRS['fit']
    fitted = json.loads((source/'summary.json').read_text()); A.check_last_identity(fitted)
    assert (source/'STATUS').read_text().strip() == 'COMPLETED'
    source_job = (source/'JOBID').read_text().strip(); assert source_job.isdecimal()
    draws, _, windows, catalogs = Q.load_train_inputs()
    log = Q.read_rows(source/'training.jsonl'); assert len(log) == 4608
    records = {}
    for step, row in enumerate(log, 1):
        ordinal = (step-1) % 1536; draw = draws[ordinal]; w = windows[draw['window_id']]
        assert row['step'] == step and row['epoch'] == (step-1)//1536+1 and row['draw_ordinal'] == ordinal
        assert row['window_id'] == w['id'] and row['species'] == w['species'] == draw['species']
        assert row['importance_weight'] == draw['importance_weight_uniform_within_stratum']
        reasons = row['protection_reasons']; assert len(reasons) == len(set(reasons))
        assert set(reasons) <= {'unknown_chain_span', 'complete_protected_nonpositive_chain', 'noncallable_valid_base'}
        n = len(w['positive_ids']); assert 0 <= n <= 16
        expected = {'1': n, '0': 0 if reasons else 16-n, '-1': 16-n if reasons else 0}
        assert row['query_label_counts'] == expected
        identity = dict(protection_reasons=reasons, query_label_counts=expected)
        if w['id'] in records:
            assert records[w['id']] == identity
        else:
            records[w['id']] = identity
    assert set(records) == set(windows) and len(records) == 1491
    result = []
    for ident, w in windows.items():
        record = records[ident]; catalog = catalogs[w['species']]
        positives = {tuple(oriented_parts(catalog[t], w)) for t in w['positive_ids']}
        protected = {tuple(oriented_parts(catalog[t], w)) for t in w['overlapping_transcript_ids']}
        protected = {chain for chain in protected if chain and chain[0][0] >= 0 and chain[-1][1] <= w['valid_bases']}
        nonpositive = protected-positives; metadata_reasons = []
        if w['chain_unknown_intervals']:
            metadata_reasons.append('unknown_chain_span')
        if nonpositive:
            metadata_reasons.append('complete_protected_nonpositive_chain')
        assert set(record['protection_reasons'])-{'noncallable_valid_base'} == set(metadata_reasons)
        intervals = list(w['chain_unknown_intervals'])+[[chain[0][0], chain[-1][1]] for chain in nonpositive]
        guard = union(intervals, w['valid_bases']); unavailable = 'noncallable_valid_base' in record['protection_reasons']
        result.append(dict(window_id=ident, species=w['species'], valid_bases=w['valid_bases'],
            positive_chains=len(w['positive_ids']), **record, guard_span_union=guard,
            noncallable_positions_unavailable=unavailable,
            outside_guard_span_bases=None if unavailable else w['valid_bases']-sum(b-a for a, b in guard),
            initial_bins_outside_guard_span=None if unavailable else outside_bins(guard, w['valid_bases'])))
    by_id = {r['window_id']: r for r in result}; sampled = [by_id[d['window_id']] for d in draws]
    per_species = {s: aggregate([r for r in sampled if r['species'] == s]) for s in Q.ALLOW}
    for epoch in fitted['epochs']:
        for s in Q.ALLOW:
            expected = epoch['supervision_by_species'][s]; actual = per_species[s]
            assert all(actual['query_label_counts'][k] == expected['query_'+k] for k in ('-1', '0', '1'))
            assert actual['protection_reason_counts'] == {k: v for k, v in expected.items() if not k.startswith('query_')}
    with (OUT/'windows.jsonl').open('x') as handle:
        for row in result:
            handle.write(json.dumps(row, allow_nan=False)+'\n')
    report = dict(completed=True, protocol=PROTOCOL, source_fit_job=source_job,
        source_final_checkpoint=fitted['final_checkpoint'], original_training_records_verified=4608,
        unique_windows=1491, draws_per_original_pass=1536, original_passes=3,
        per_unique_window=aggregate(result), per_original_pass=aggregate(sampled),
        per_species_original_pass=per_species,
        original_three_pass_counts={s: {k: 3*v for k, v in r['query_label_counts'].items()} for s, r in per_species.items()},
        optimizer_steps_this_run=0, model_or_backbone_forward=False, checkpoint_loaded=False,
        feature_tensor_or_DNA_loaded=False, changes_to_label_policy=False, new_candidate_search=False,
        DEV_forward_used=False, test_setaria_used=False, seconds=time.perf_counter()-began,
        limitation='TRAIN metadata and actual query supervision counts. Outside guard spans or initial bins are not valid query-negative labels: pointers and continuation remain global. No causal attribution, model accuracy or independent generalization.')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(report, allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
