#!/usr/bin/env python3
"""Paired saved-chain errors after local-objective continuation, original TRAIN only."""
from collections import Counter, defaultdict
import json
import os
import time
import profile_full_grammar_original_train_errors as G

ROOT = G.I.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-LOCAL-CHAIN-TRANSITIONS'
OBS = ROOT/'outputs/M28-FULL-GRAMMAR-TRAIN-OBSERVATION-DIAGNOSTIC'
RAW = ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-TRAIN-RAW/LOCAL_CHANNEL'
EVAL = ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-TRAIN-DIAGNOSTIC/LOCAL_CHANNEL'
CHECKPOINT = ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST/LOCAL_CHANNEL/step_009216.pt'
EXACT = 'exact_primary_chain'


def key(row):
    return row['species'], row['seqid'], row['id'], row['strand']


def read_rows(filename):
    result = {}
    with filename.open() as handle:
        for line in handle:
            row = json.loads(line)
            ident = key(row)
            if ident in result:
                raise ValueError('Repeated primary reference identity')
            result[ident] = row
    return result


def completed(folder, checkpoint=None):
    report = json.loads((folder/'summary.json').read_text())
    if ((folder/'STATUS').read_text().strip() != 'COMPLETED' or not report['completed']
            or report['DEV_used'] or report['sealed_data_used'] or report['optimizer_steps_this_run'] != 0
            or (checkpoint is not None and report['checkpoint'] != str(checkpoint))):
        raise ValueError('Completed original TRAIN-only source required')
    return report


def transition(old, new):
    return { (True, True): 'retained', (True, False): 'lost',
             (False, True): 'gained', (False, False): 'still_missed'}[old == EXACT, new == EXACT]


def counts(rows):
    transitions = Counter(r['transition'] for r in rows)
    matrix = defaultdict(Counter)
    for row in rows:
        matrix[row['old_category']][row['new_category']] += 1
    old_tp = transitions['retained'] + transitions['lost']
    new_tp = transitions['retained'] + transitions['gained']
    return dict(reference=len(rows), transitions={k:transitions[k] for k in
        ('retained', 'lost', 'gained', 'still_missed')}, old_exact_TP=old_tp, new_exact_TP=new_tp,
        net_exact_TP=new_tp-old_tp, old_recall=old_tp/len(rows) if rows else None,
        new_recall=new_tp/len(rows) if rows else None,
        lost_destination_categories=dict(Counter(r['new_category'] for r in rows if r['transition'] == 'lost')),
        gained_source_categories=dict(Counter(r['old_category'] for r in rows if r['transition'] == 'gained')),
        category_transition_matrix={k:dict(v) for k,v in matrix.items()})


def summarize(rows):
    result = counts(rows)
    for field in ('CDS_segment_group', 'original_core_coverage', 'span_group', 'single_chain_supported'):
        result['by_'+field] = {str(value):counts([r for r in rows if r[field] == value])
            for value in sorted({r[field] for r in rows})}
    return result


def paired_rows(refs, old, new, observations):
    if not (refs.keys() == old.keys() == new.keys() == observations.keys()):
        raise ValueError('Paired reference populations differ; do not inner-join away records')
    result = []
    for ident, ref in refs.items():
        before, after, obs = old[ident], new[ident], observations[ident]
        if (ref['CDS_segments'] != before['CDS_segments'] or ref['CDS_segments'] != after['CDS_segments']
                or ref['exact_primary_hit'] != (before['category'] == EXACT)
                or obs['exact_primary_hit'] != ref['exact_primary_hit']):
            raise ValueError('Saved primary reference geometry or original hit identity differs')
        result.append(dict(species=ref['species'], seqid=ref['seqid'], id=ref['id'], strand=ref['strand'],
            CDS_segment_group='single' if ref['CDS_segments'] == 1 else 'multiple',
            span_group='over_24576' if ref['span_bp'] > 24576 else 'at_most_24576',
            original_core_coverage=obs['coverage_class'], single_chain_supported=obs['single_chain_supported'],
            old_category=before['category'], new_category=after['category'],
            transition=transition(before['category'], after['category']),
            old_CDS_overlapping_predictions=before['same_strand_CDS_overlap_prediction_count'],
            new_CDS_overlapping_predictions=after['same_strand_CDS_overlap_prediction_count']))
    return result


def fixture():
    rows = [dict(old_category=a, new_category=b, transition=transition(a,b)) for a,b in
        ((EXACT,EXACT), (EXACT,'no_same_strand_CDS_overlap'),
         ('same_CDS_ends_internal_structure_differs',EXACT),
         ('same_intron_chain_CDS_start_only_differs','CDS_ends_and_internal_structure_differ'))]
    result = counts(rows)
    assert result['transitions'] == dict(retained=1,lost=1,gained=1,still_missed=1)
    assert result['old_exact_TP'] == result['new_exact_TP'] == 2 and result['net_exact_TP'] == 0
    assert result['lost_destination_categories'] == {'no_same_strand_CDS_overlap':1}
    assert result['gained_source_categories'] == {'same_CDS_ends_internal_structure_differs':1}
    ref = dict(species='s',seqid='a',id='r',strand='-',CDS_segments=2,span_bp=30,exact_primary_hit=True)
    before = dict(ref,category=EXACT,same_strand_CDS_overlap_prediction_count=1)
    after = dict(ref,category='same_CDS_ends_internal_structure_differs',same_strand_CDS_overlap_prediction_count=1)
    obs = dict(ref,coverage_class='no_direct_core_span_observation',single_chain_supported=True)
    ident = key(ref)
    joined = paired_rows({ident:ref},{ident:before},{ident:after},{ident:obs})
    assert summarize(joined)['by_CDS_segment_group']['multiple']['transitions']['lost'] == 1
    for wrong in (dict(ref,strand='+'), dict(ref,seqid='b')):
        try:
            paired_rows({ident:ref},{ident:before},{key(wrong):after},{ident:obs})
        except ValueError:
            pass
        else:
            raise AssertionError('Strand or chromosome mismatch must reject the join')
    return 'four_paired_outcomes_category_direction_and_full_population_strand_contig_join_passed'


def main():
    began = time.perf_counter()
    tested = fixture()
    original = completed(G.SOURCE, G.I.CHECKPOINT)
    errors = completed(G.OUT, G.I.CHECKPOINT)
    observations = completed(OBS)
    candidate = completed(EVAL, CHECKPOINT)
    if (errors['source_profile_job'] != original['job_id'] or observations['source_fit_job'] != '13434603'
            or candidate['job_id'] != '13451623' or candidate['reference_policy'] != original['reference_policy']):
        raise ValueError('Completed declared fit, reference profile and LOCAL_CHANNEL evaluation required')
    refs = read_rows(G.SOURCE/'reference_chain_profile.jsonl')
    old = read_rows(G.OUT/'reference_error_profile.jsonl')
    observed = read_rows(OBS/'reference_observations.jsonl')
    scope = {seqid:species for species,seqid,_ in G.I.SCOPE}
    if (len(refs) != 22688 or any(scope.get(r['seqid']) != r['species'] for r in refs.values())):
        raise ValueError('Original six TRAIN chromosomes and entire frozen primary denominator required')
    report = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],fixture=tested,
        scope='saved_original_six_TRAIN_paired_chain_errors_after_LOCAL_CHANNEL_continuation',
        old_checkpoint=str(G.I.CHECKPOINT),new_checkpoint=str(CHECKPOINT),
        source_profile_job=original['job_id'],source_error_job=errors['job_id'],
        source_observation_job=observations['job_id'],source_candidate_evaluation_job=candidate['job_id'],
        reference_policy=original['reference_policy'],primary_reference_changed=False,
        DEV_used=False,sealed_data_used=False,optimizer_steps_this_run=0,DNA_read=False,
        checkpoint_read=False,features_read=False,native_forward=False,predictions_modified=False,
        old_anchor_has_fewer_updates_not_matched_update_objective_control=True,
        association_rule=errors['association_rule'],overlap_graph_is_not_biological_fusion_or_fragmentation=True,
        geometry_events_are_not_gene_accuracy_or_phase_metrics=True,
        original_core_coverage_is_not_all_gradient_or_encoder_exposure_or_independent_test=True,
        chromosomes={},per_species={})
    new, species_data = {}, defaultdict(lambda:([],[]))
    for species,seqid,length in G.I.SCOPE:
        reference = [r for r in refs.values() if r['seqid'] == seqid]
        predicted = G.C.S.E.primary_transcripts(G.C.S.E.parse_annotation(
            RAW/seqid/'prediction/full_grammar.gff3',{seqid:length}),complete_only=False)
        geometry, detailed = G.describe(predicted,reference)
        if geometry['exact_chain'] != candidate['chromosomes'][seqid]['exact_chain']:
            raise ValueError('Saved candidate geometry disagrees with completed unchanged TRAIN ruler')
        report['chromosomes'][seqid] = dict(old_geometry=errors['chromosomes'][seqid],new_geometry=geometry)
        new.update((key(r),r) for r in detailed)
        species_data[species][0].extend(predicted); species_data[species][1].extend(reference)
    rows = paired_rows(refs,old,new,observed)
    for seqid in scope:
        report['chromosomes'][seqid]['paired'] = summarize([r for r in rows if r['seqid'] == seqid])
    for species,(predicted,reference) in species_data.items():
        geometry, _ = G.describe(predicted,reference)
        if geometry['exact_chain'] != candidate['per_species'][species]['exact_chain']:
            raise ValueError('Species primary denominator or hit count differs')
        report['per_species'][species] = dict(old_geometry=errors['per_species'][species],new_geometry=geometry,
            paired=summarize([r for r in rows if r['species'] == species]))
    pooled = tuple([r for data in species_data.values() for r in data[i]] for i in (0,1))
    geometry, _ = G.describe(*pooled)
    if geometry['exact_chain'] != candidate['pooled']['exact_chain']:
        raise ValueError('Pooled primary denominator or hit count differs')
    report['pooled'] = dict(old_geometry=errors['pooled'],new_geometry=geometry,paired=summarize(rows))
    for label, data in [*report['chromosomes'].items(),*report['per_species'].items(),('pooled',report['pooled'])]:
        pairs = data['paired']
        if (pairs['old_exact_TP'] != data['old_geometry']['exact_chain']['tp']
                or pairs['new_exact_TP'] != data['new_geometry']['exact_chain']['tp']):
            raise ValueError('Paired retained/lost/gained counts disagree with fixed exact-chain measurements: '+label)
    with (OUT/'reference_transitions.jsonl').open('x') as handle:
        for row in rows:
            handle.write(json.dumps(row,allow_nan=False)+'\n')
    report['seconds'] = time.perf_counter()-began
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps(dict(completed=True,seconds=report['seconds'],paired=report['pooled']['paired']['transitions'])),flush=True)


if __name__ == '__main__':
    main()
