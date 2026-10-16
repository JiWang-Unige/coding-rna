#!/usr/bin/env python3
"""Descriptive geometry errors on saved original-TRAIN chains; no model or DNA."""
from bisect import bisect_right
from collections import Counter, defaultdict
import json
import os
import time
import common_ruler as C
import infer_full_grammar_original_train as I

SOURCE = I.ROOT/'outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-CHAIN-PROFILE'
OUT = I.ROOT/'outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-ERROR-PROFILE-R1'


def introns(tx):
    return tuple((a[1], b[0]) for a, b in zip(tx['CDS'], tx['CDS'][1:]))


def ends(tx):
    lo, hi = tx['CDS'][0][0], tx['CDS'][-1][1]
    return (lo, hi) if tx['strand'] == '+' else (hi, lo)


def shared_bp(a, b):
    i = j = total = 0
    while i < len(a) and j < len(b):
        total += max(0, min(a[i][1], b[j][1])-max(a[i][0], b[j][0]))
        if a[i][1] <= b[j][1]:
            i += 1
        else:
            j += 1
    return total


def prediction_index(txs):
    groups = defaultdict(list)
    for tx in txs:
        groups[tx['seqid'], tx['strand']].append(tx)
    result = {}
    for key, rows in groups.items():
        rows.sort(key=lambda tx: (tx['CDS'][0][0], tx['CDS'][-1][1], tx['id']))
        if any(a['CDS'][-1][1] > b['CDS'][0][0] for a, b in zip(rows, rows[1:])):
            raise ValueError('Saved selected predictions no longer have disjoint same-strand spans')
        result[key] = (rows, [tx['CDS'][-1][1] for tx in rows])
    return result


def overlaps(ref, index, strand=None):
    rows, stop = index.get((ref['seqid'], strand or ref['strand']), ([], []))
    start, end = ref['CDS'][0][0], ref['CDS'][-1][1]
    found = []
    for at in range(bisect_right(stop, start), len(rows)):
        tx = rows[at]
        if tx['CDS'][0][0] >= end:
            break
        bp = shared_bp(ref['CDS'], tx['CDS'])
        if bp:
            found.append((tx, bp))
    return found


def category(ref, pred):
    if pred is None:
        return 'no_same_strand_CDS_overlap'
    if C.S.key(ref) == C.S.key(pred):
        return 'exact_primary_chain'
    start, stop = (a == b for a, b in zip(ends(ref), ends(pred)))
    if introns(ref) == introns(pred):
        if stop:
            return 'same_intron_chain_CDS_start_only_differs'
        if start:
            return 'same_intron_chain_CDS_stop_only_differs'
        return 'same_intron_chain_both_CDS_ends_differ'
    if start and stop:
        return 'same_CDS_ends_internal_structure_differs'
    return 'CDS_ends_and_internal_structure_differ'


def feature_scores(predicted, reference):
    def events(rows, name):
        result = set()
        for tx in rows:
            prefix = (tx['seqid'], tx['strand'])
            start, stop = ends(tx)
            values = {'CDS_start': (start,), 'CDS_stop': (stop,),
                      'CDS_segment': tuple(tuple(x[:2]) for x in tx['CDS']),
                      'intron_pair': introns(tx)}[name]
            result.update(prefix+(value,) for value in values)
        return result
    result = {}
    for name in ('CDS_start', 'CDS_stop', 'CDS_segment', 'intron_pair'):
        p, r = events(predicted, name), events(reference, name)
        result[name] = C.S.prf(len(p & r), len(p), len(r))
    return result


def describe(predicted, reference):
    index = prediction_index(predicted)
    pk, rk = ({C.S.key(tx) for tx in rows} for rows in (predicted, reference))
    bins = Counter()
    groups = {'single_CDS_segment': Counter(), 'multiple_CDS_segments': Counter()}
    prediction_degree = Counter()
    detailed = []
    union_bp = reference_bp = best_bp = opposite_only = multiple_preds = 0
    for ref in reference:
        neighbors = overlaps(ref, index)
        def rank(item):
            tx, bp = item
            return (-bp, -len(set(introns(ref)) & set(introns(tx))),
                    sum(abs(a-b) for a, b in zip(ends(ref), ends(tx))), C.S.key(tx), tx['id'])
        pred, overlap = min(neighbors, key=rank) if neighbors else (None, 0)
        label = category(ref, pred)
        if (label == 'exact_primary_chain') != (C.S.key(ref) in pk):
            raise ValueError('Descriptive matching disagrees with fixed exact-chain hit')
        bins[label] += 1
        group = 'single_CDS_segment' if len(ref['CDS']) == 1 else 'multiple_CDS_segments'
        groups[group][label] += 1
        # Selected prediction spans are disjoint, so their CDS overlaps cannot double-count bp.
        covered = sum(bp for _, bp in neighbors)
        size = sum(b-a for a, b, _ in ref['CDS'])
        if covered > size:
            raise ValueError('Same-strand union CDS overlap exceeds reference CDS length')
        opposite = bool(not neighbors and overlaps(ref, index, '-' if ref['strand'] == '+' else '+'))
        opposite_only += opposite
        multiple_preds += len(neighbors) > 1
        reference_bp += size; best_bp += overlap; union_bp += covered
        for tx, _ in neighbors:
            prediction_degree[C.S.key(tx)] += 1
        detailed.append(dict(species=ref.get('species'), seqid=ref['seqid'], id=ref['id'],
            strand=ref['strand'], CDS_segments=len(ref['CDS']), category=label,
            selected_prediction_id=pred['id'] if pred else None,
            selected_prediction_CDS_overlap_bp=overlap, reference_CDS_bp=size,
            same_strand_union_CDS_overlap_bp=covered,
            same_strand_CDS_overlap_prediction_count=len(neighbors),
            opposite_strand_CDS_overlap_only=opposite))
    if sum(bins.values()) != len(reference) or len(rk) != len(reference):
        raise ValueError('Reference taxonomy changed the frozen primary denominator')
    report = dict(exact_chain=C.S.prf(len(pk & rk), len(pk), len(rk)),
        reference_error_categories=dict(bins), by_CDS_segment_count={k:dict(v) for k,v in groups.items()},
        strand_sensitive_unique_geometry_events=feature_scores(predicted, reference),
        CDS_overlap_graph=dict(reference_CDS_bp_sum=reference_bp,
            best_prediction_overlap_bp_sum=best_bp, same_strand_union_overlap_bp_sum=union_bp,
            reference_weighted_best_prediction_CDS_coverage=best_bp/reference_bp if reference_bp else None,
            reference_weighted_same_strand_union_CDS_coverage=union_bp/reference_bp if reference_bp else None,
            references_with_multiple_CDS_overlapping_predictions=multiple_preds,
            predictions_with_multiple_CDS_overlapping_primary_references=sum(v > 1 for v in prediction_degree.values()),
            references_without_same_strand_but_with_opposite_strand_CDS_overlap=opposite_only,
            same_strand_CDS_overlap_edges=sum(prediction_degree.values())))
    return report, detailed


def fixture():
    def tx(ident, cds, strand='+', seqid='a'):
        return dict(id=ident, seqid=seqid, strand=strand, CDS=cds)
    ref = tx('r', [(0,9,0),(21,30,0)])
    assert category(ref, ref) == 'exact_primary_chain'
    assert category(ref, tx('p', [(3,9,0),(21,30,0)])) == 'same_intron_chain_CDS_start_only_differs'
    assert category(tx('r', ref['CDS'], '-'), tx('p', [(3,9,0),(21,30,0)], '-')) == 'same_intron_chain_CDS_stop_only_differs'
    assert category(ref, tx('p', [(0,12,0),(21,30,0)])) == 'same_CDS_ends_internal_structure_differs'
    assert category(ref, tx('p', [(3,12,0),(21,30,0)])) == 'CDS_ends_and_internal_structure_differ'
    assert not overlaps(ref, prediction_index([tx('intron', [(12,18,0)])]))
    result, rows = describe([tx('left', [(0,15,0)]), tx('right', [(15,30,0)])], [tx('r', [(0,30,0)])])
    assert result['CDS_overlap_graph']['same_strand_union_overlap_bp_sum'] == 30
    assert result['CDS_overlap_graph']['best_prediction_overlap_bp_sum'] == 15
    assert result['CDS_overlap_graph']['references_with_multiple_CDS_overlapping_predictions'] == 1
    assert rows[0]['selected_prediction_id'] == 'left'
    result, rows = describe([tx('other', [(0,9,0)], '-'), tx('b', [(0,9,0)], seqid='b')], [tx('r', [(0,9,0)])])
    assert result['reference_error_categories'] == {'no_same_strand_CDS_overlap': 1}
    assert rows[0]['opposite_strand_CDS_overlap_only']
    assert result['strand_sensitive_unique_geometry_events']['CDS_segment']['tp'] == 0
    # An exact shared junction does not make shifted CDS ends a complete gene hit.
    result, _ = describe([tx('p', [(3,9,0),(21,30,0)])], [ref])
    assert result['strand_sensitive_unique_geometry_events']['intron_pair']['tp'] == 1
    assert result['exact_chain']['tp'] == 0
    return 'passed_strand_contig_terminal_internal_intron_only_fragment_and_geometry_fixtures'


def main():
    began = time.perf_counter()
    checked = fixture()
    original = json.loads((SOURCE/'summary.json').read_text())
    if ((SOURCE/'STATUS').read_text().strip() != 'COMPLETED' or not original['completed'] or
            original['checkpoint'] != str(I.CHECKPOINT) or original['DEV_used'] or original['sealed_data_used']):
        raise ValueError('Completed frozen original-TRAIN profile required')
    scope = {seqid:(species,length) for species,seqid,length in I.SCOPE}
    refs = defaultdict(list)
    with (SOURCE/'reference_chain_profile.jsonl').open() as handle:
        for line in handle:
            row = json.loads(line)
            if row['seqid'] not in scope or row['species'] != scope[row['seqid']][0]:
                raise ValueError('Reference profile not confined to original TRAIN')
            refs[row['seqid']].append(row)
    report = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], fixture=checked,
        scope='fixed_original_six_TRAIN_saved_chain_geometry_diagnosis_only',
        checkpoint=str(I.CHECKPOINT), source_profile_job=original['job_id'],
        source_diagnostic_job=original['source_diagnostic_job'], reference_policy=original['reference_policy'],
        DEV_used=False, sealed_data_used=False, optimizer_steps_this_run=0,
        DNA_read=False, checkpoint_read=False, features_read=False, native_forward=False,
        predictions_modified=False, primary_reference_unchanged=True,
        association_rule='same_strand_positive_CDS_overlap_max_bp_then_exact_introns_then_CDS_end_distance_then_geometry_and_ID',
        association_is_descriptive_not_one_to_one_gene_matching_or_causal_module_diagnosis=True,
        overlaps_do_not_prove_fragmentation_or_fusion=True,
        geometry_events_are_not_gene_accuracy_or_phase_metrics=True,
        chromosomes={}, per_species={}, pooled={})
    species_rows = defaultdict(lambda: ([], []))
    with (OUT/'reference_error_profile.jsonl').open('x') as handle:
        for species,seqid,length in I.SCOPE:
            txs = C.S.E.primary_transcripts(C.S.E.parse_annotation(
                I.OUT/seqid/'prediction/full_grammar.gff3', {seqid:length}), complete_only=False)
            row, detailed = describe(txs, refs[seqid])
            if row['exact_chain'] != original['chromosomes'][seqid]['exact_chain']:
                raise ValueError('Error profile differs from fixed original-TRAIN score')
            report['chromosomes'][seqid] = row
            species_rows[species][0].extend(txs); species_rows[species][1].extend(refs[seqid])
            for item in detailed:
                handle.write(json.dumps(item)+'\n')
            print(json.dumps(dict(profiled_TRAIN=seqid, reference=len(refs[seqid]))), flush=True)
    for species,(predicted,reference) in species_rows.items():
        row, _ = describe(predicted, reference)
        if row['exact_chain'] != original['per_species'][species]['exact_chain']:
            raise ValueError('Species denominator or hit count changed')
        report['per_species'][species] = row
    pooled = tuple([item for rows in species_rows.values() for item in rows[i]] for i in (0,1))
    report['pooled'], _ = describe(*pooled)
    if report['pooled']['exact_chain'] != original['pooled']['exact_chain']:
        raise ValueError('Pooled denominator or hit count changed')
    report['seconds'] = time.perf_counter()-began
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(dict(completed=True, seconds=report['seconds'])), flush=True)


if __name__ == '__main__':
    main()
