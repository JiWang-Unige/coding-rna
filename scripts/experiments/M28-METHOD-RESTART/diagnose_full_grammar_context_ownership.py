#!/usr/bin/env python3
"""Original TRAIN metadata-only fit-view versus whole-contig owner comparison."""
from bisect import bisect_right
from collections import Counter
import json
import os
import time

import infer_full_grammar_original_train as I
from src.m28.full_grammar_context_view import context_geometry
from src.m28.contig_emission_lattice import context_segments

OUT = I.ROOT/'outputs/M28-FULL-GRAMMAR-CONTEXT-OWNERSHIP-DIAGNOSTIC'
PREP = I.ROOT/'outputs/M28-FULL-GRAMMAR-COARSE-CHANNEL-FIT-PREPARATION'
FIT = I.ROOT/'outputs/M28-FULL-GRAMMAR-COARSE-CHANNEL-FIRST-FIT'


def labeled_segments(windows, length, offset=0):
    return [(a+offset, b+offset, (windows[i][0]+offset, windows[i][1]+offset))
            for i, a, b in context_segments(length, windows)]


def differences(local, global_segments, lo, hi):
    """Compare source window identity by interval intersections, not by DNA."""
    i = bisect_right([b for _, b, _ in local], lo)
    j = bisect_right([b for _, b, _ in global_segments], lo)
    cursor, result = lo, []
    while cursor < hi:
        a, b, owner = local[i]
        x, y, other = global_segments[j]
        if not a <= cursor < b or not x <= cursor < y:
            raise ValueError('Owner segment coverage gap')
        end = min(b, y, hi)
        if owner != other:
            result.append(dict(genomic_interval=[cursor, end],
                               fit_view_owner=list(owner), whole_contig_owner=list(other)))
        cursor = end
        if cursor == b:
            i += 1
        if cursor == y:
            j += 1
    return result


def known_mismatch_bounds(n, known, mismatch):
    # No per-base labels are read: counts alone give bounds, not exact overlap.
    return max(0, mismatch-(n-known)), min(mismatch, known)


def summarize(rows):
    result = dict(windows=len(rows), core_owner_mismatch_windows=0,
        off_whole_contig_grid_neighbor_windows=0, context_flank_owner_mismatch_windows=0,
        core_bp=0, context_bp=0, context_flank_bp=0, core_owner_mismatch_bp=0,
        context_flank_owner_mismatch_bp=0, known_core_rows=0,
        known_core_owner_mismatch_lower_bound=0, known_core_owner_mismatch_upper_bound=0)
    for row in rows:
        result['core_owner_mismatch_windows'] += row['core_owner_mismatch_bp'] > 0
        result['off_whole_contig_grid_neighbor_windows'] += bool(row['off_whole_contig_grid_neighbors'])
        result['context_flank_owner_mismatch_windows'] += row['context_flank_owner_mismatch_bp'] > 0
        for key in ('core_bp', 'context_bp', 'context_flank_bp', 'core_owner_mismatch_bp',
                    'context_flank_owner_mismatch_bp', 'known_core_rows',
                    'known_core_owner_mismatch_lower_bound', 'known_core_owner_mismatch_upper_bound'):
            result[key] += row[key]
    return result


def fixture():
    length = 40000
    windows = [(a, min(a+I.R.M.WIDTH, length)) for a in I.R.M.starts(length)]
    whole = labeled_segments(windows, length)
    for start in I.R.M.starts(length):
        end = start+I.R.M.WIDTH
        geometry = context_geometry(start, end, length)
        a, b = geometry['genomic_context']
        local = labeled_segments(geometry['local_windows'], b-a, a)
        actual = differences(local, whole, start, end)
        # Independent per-base nearest-containing-center calculation.
        expected = 0
        for p in range(start, end):
            choose = lambda ws: min((w for w in ws if w[0] <= p < w[1]),
                                   key=lambda w: (abs(2*p+1-sum(w)), w[0]))
            expected += choose(geometry['genomic_windows']) != choose(windows)
        assert sum(r['genomic_interval'][1]-r['genomic_interval'][0] for r in actual) == expected
        if start == length-I.R.M.WIDTH:
            assert expected == 10720
        else:
            assert expected == 0
        for strand in ('+', '-'):
            offset = start-a if strand == '+' else b-end
            inferred = (a+offset, a+offset+end-start) if strand == '+' else (b-offset-(end-start), b-offset)
            assert inferred == (start, end)
    assert known_mismatch_bounds(10, 8, 3) == (1, 3)
    assert known_mismatch_bounds(10, 0, 3) == (0, 0)
    return 'interval_sweep_matches_direct_base_owner_shifted_tail_strands_and_count_bounds'


def main():
    began = time.perf_counter()
    checked = fixture()
    prep = json.loads((PREP/'summary.json').read_text())
    fit = json.loads((FIT/'summary.json').read_text())
    if ((PREP/'STATUS').read_text().strip() != 'COMPLETED'
            or (FIT/'STATUS').read_text().strip() != 'COMPLETED'
            or not prep['completed'] or not fit['completed'] or fit['optimizer_steps'] != 4608
            or prep['configuration'] != fit['configuration'] or fit['DEV_used']):
        raise ValueError('Completed original fixed three-pass fit and preparation required')
    plan = [json.loads(line) for line in (PREP/'training_plan.jsonl').read_text().splitlines()]
    if (len(plan) != 1536 or len({r['window_id'] for r in plan}) != 1491
            or Counter(r['species'] for r in plan) != Counter({s:768 for s, _, _ in I.SCOPE})):
        raise ValueError('Original paired draw population differs')
    scopes = {seqid:(species, length) for species, seqid, length in I.SCOPE}
    grids = {seqid:[(a, min(a+I.R.M.WIDTH, length)) for a in I.R.M.starts(length)]
             for seqid, (_, length) in scopes.items()}
    segments = {seqid:labeled_segments(grid, scopes[seqid][1]) for seqid, grid in grids.items()}
    rows, unique = [], {}
    with (OUT/'windows.jsonl').open('x') as handle:
        for ordinal, original in enumerate(plan):
            seqid, strand = original['stratum']
            species, length = scopes[seqid]
            if original['draw_ordinal'] != ordinal or original['species'] != species or strand not in ('+', '-'):
                raise ValueError('Original draw scope or order differs')
            a, b = original['genomic_context']
            lo, hi = original['core_oriented_interval']
            n = hi-lo
            start, end = (a+lo, a+hi) if strand == '+' else (b-hi, b-lo)
            geometry = context_geometry(start, end, length)
            known = original['fine_positions']+original['coarse_positions']
            if (geometry['genomic_context'] != [a, b] or n != known+original['fully_unknown_positions']
                    or (start, end) not in grids[seqid]):
                raise ValueError('Stored fit-view geometry or observations differ')
            local = labeled_segments(geometry['local_windows'], b-a, a)
            core_diff = differences(local, segments[seqid], start, end)
            full_diff = differences(local, segments[seqid], a, b)
            core_bp = sum(r['genomic_interval'][1]-r['genomic_interval'][0] for r in core_diff)
            full_bp = sum(r['genomic_interval'][1]-r['genomic_interval'][0] for r in full_diff)
            lower, upper = known_mismatch_bounds(n, known, core_bp)
            row = dict(window_id=original['window_id'], species=species, seqid=seqid, strand=strand,
                genomic_core=[start, end], genomic_context=[a, b],
                core_bp=n, context_bp=b-a, context_flank_bp=b-a-n,
                core_owner_mismatch_bp=core_bp, context_flank_owner_mismatch_bp=full_bp-core_bp,
                known_core_rows=known, known_core_owner_mismatch_lower_bound=lower,
                known_core_owner_mismatch_upper_bound=upper,
                off_whole_contig_grid_neighbors=[list(w) for w in geometry['genomic_windows'] if w not in grids[seqid]],
                core_owner_differences=core_diff)
            if row['window_id'] in unique and unique[row['window_id']] != row:
                raise ValueError('Repeated original window has different geometry')
            unique[row['window_id']] = row
            rows.append(row)
        for row in unique.values():
            handle.write(json.dumps(row, allow_nan=False)+'\n')
    report = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], fixture=checked,
        source_preparation_job=prep['job_id'], source_fit_job=fit['job_id'],
        checkpoint=str(I.CHECKPOINT), original_passes=3,
        scope='original_1536_draws_1491_unique_TRAIN_metadata_owner_comparison_only',
        model_or_checkpoint_read=False, DNA_or_reference_read=False, labels_or_features_read=False,
        native_forward=False, optimizer_steps_this_run=0, predictions_modified=False,
        DEV_used=False, sealed_data_used=False,
        owner_comparison_is_geometric_not_logit_or_gradient_equivalence=True,
        known_mismatch_is_bound_not_per_base_observation_join=True,
        flank_rows_affect_finite_view_partition_even_without_direct_observations=True,
        counts_are_per_one_pass_not_unique_genomic_coverage=True,
        unique_windows=summarize(list(unique.values())), per_pass_draws=summarize(rows),
        per_species={s:dict(unique_windows=summarize([r for r in unique.values() if r['species'] == s]),
                           per_pass_draws=summarize([r for r in rows if r['species'] == s]))
                     for s in {r['species'] for r in rows}},
        affected_core_windows=[r for r in unique.values() if r['core_owner_mismatch_bp']],
        seconds=time.perf_counter()-began,
        interpretation='fit_view_vs_whole_contig_source_window_diagnosis_not_causal_error_attribution_or_protocol_change')
    if report['per_pass_draws']['known_core_rows'] != sum(
            r['fine_exposures']+r['coarse_exposures'] for r in prep['by_species'].values()):
        raise ValueError('Original known observation count reconciliation differs')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(report, allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
