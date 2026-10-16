#!/usr/bin/env python3
"""Decompose same-end internal chain transitions against frozen TRAIN primaries."""
from collections import Counter
import json
import os
import time
import diagnose_full_grammar_local_chain_transitions as T

OUT = T.ROOT/'outputs/M28-FULL-GRAMMAR-LOCAL-INTRON-EDITS'
INTERNAL = 'same_CDS_ends_internal_structure_differs'


def span_key(tx):
    return tx['seqid'], tx['strand'], tx['CDS'][0][0], tx['CDS'][-1][1]


def edits(reference, prediction):
    ref, pred = set(T.G.introns(reference)), set(T.G.introns(prediction))
    missing, extra = ref-pred, pred-ref
    if not missing and not extra:
        raise ValueError('Internal-structure transition must change an intron pair')
    label = 'missing_and_extra' if missing and extra else 'missing_only' if missing else 'extra_only'
    return dict(edit_class=label, reference_introns=len(ref), prediction_introns=len(pred),
        shared=len(ref & pred), missing=len(missing), extra=len(extra))


def summarize(rows):
    return dict(reference=len(rows),edit_classes=dict(Counter(r['edit_class'] for r in rows)),
        intron_pair_counts={k:sum(r[k] for r in rows) for k in
            ('reference_introns','prediction_introns','shared','missing','extra')})


def fixture():
    def tx(cds):
        return dict(seqid='a',strand='-',CDS=[(a,b,0) for a,b in cds])
    reference = tx([(0,9),(21,30)])
    assert edits(reference,tx([(0,3),(6,9),(21,30)])) == dict(
        edit_class='extra_only',reference_introns=1,prediction_introns=2,shared=1,missing=0,extra=1)
    assert edits(reference,tx([(0,30)]))['edit_class'] == 'missing_only'
    assert edits(reference,tx([(0,12),(21,30)]))['edit_class'] == 'missing_and_extra'
    assert span_key(reference) == ('a','-',0,30)
    try:
        edits(reference,reference)
    except ValueError:
        pass
    else:
        raise AssertionError('Unchanged introns cannot explain an internal transition')
    return 'extra_missing_mixed_and_unchanged_intron_fixtures_passed'


def main():
    began = time.perf_counter()
    checked = fixture()
    source = json.loads((T.OUT/'summary.json').read_text())
    if ((T.OUT/'STATUS').read_text().strip() != 'COMPLETED' or not source['completed']
            or source['job_id'] != '13453930' or source['DEV_used'] or source['sealed_data_used']
            or source['new_checkpoint'] != str(T.CHECKPOINT) or source['optimizer_steps_this_run'] != 0):
        raise ValueError('Completed unchanged TRAIN paired transitions required')
    references = T.read_rows(T.G.SOURCE/'reference_chain_profile.jsonl')
    selected = {'lost_same_ends':[], 'gained_same_ends':[]}
    with (T.OUT/'reference_transitions.jsonl').open() as handle:
        for line in handle:
            row = json.loads(line)
            if row['transition'] == 'lost' and row['new_category'] == INTERNAL:
                selected['lost_same_ends'].append(row)
            elif row['transition'] == 'gained' and row['old_category'] == INTERNAL:
                selected['gained_same_ends'].append(row)
    details = []
    for species,seqid,length in T.G.I.SCOPE:
        indexes = {}
        for direction, root in (('lost_same_ends',T.RAW),('gained_same_ends',T.G.I.OUT)):
            txs = T.G.C.S.E.primary_transcripts(T.G.C.S.E.parse_annotation(
                root/seqid/'prediction/full_grammar.gff3',{seqid:length}),complete_only=False)
            by_span = {span_key(tx):tx for tx in txs}
            if len(by_span) != len(txs):
                raise ValueError('Same-strand selected prediction spans must be unique')
            indexes[direction] = by_span
        for direction, rows in selected.items():
            for row in rows:
                if row['seqid'] != seqid:
                    continue
                reference = references[T.key(row)]
                prediction = indexes[direction].get(span_key(reference))
                if prediction is None:
                    raise ValueError('Saved same-end transition has no matching boundary pair')
                details.append(dict(species=species,seqid=seqid,id=row['id'],strand=row['strand'],
                    direction=direction,original_core_coverage=row['original_core_coverage'],
                    **edits(reference,prediction)))
    report = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],fixture=checked,
        source_transition_job=source['job_id'],old_checkpoint=source['old_checkpoint'],
        new_checkpoint=source['new_checkpoint'],reference_policy=source['reference_policy'],
        scope='saved_TRAIN_same_end_internal_transition_intron_pairs_against_fixed_primary_only',
        DEV_used=False,sealed_data_used=False,optimizer_steps_this_run=0,DNA_read=False,
        checkpoint_read=False,features_read=False,native_forward=False,predictions_modified=False,
        counts_are_reference_prediction_pairs_not_globally_unique_intron_events=True,
        extra_and_missing_are_relative_to_primary_not_biological_false_or_true_introns=True,
        alternate_assessable_isoform_membership_not_determined=True,
        old_anchor_has_fewer_updates_not_matched_update_objective_control=True,
        per_species={},pooled={})
    for label in (*source['per_species'], 'pooled'):
        target = report['pooled'] if label == 'pooled' else report['per_species'].setdefault(label,{})
        paired = source['pooled']['paired'] if label == 'pooled' else source['per_species'][label]['paired']
        for direction in selected:
            rows = [r for r in details if r['direction'] == direction and (label == 'pooled' or r['species'] == label)]
            target[direction] = summarize(rows)
            categories = paired['lost_destination_categories' if direction == 'lost_same_ends' else 'gained_source_categories']
            if len(rows) != categories.get(INTERNAL,0):
                raise ValueError('Same-end subgroup denominator differs from paired transition report')
    with (OUT/'reference_intron_edits.jsonl').open('x') as handle:
        for row in details:
            handle.write(json.dumps(row,allow_nan=False)+'\n')
    report['seconds'] = time.perf_counter()-began
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps(report['pooled']),flush=True)


if __name__ == '__main__':
    main()
