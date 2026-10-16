#!/usr/bin/env python3
"""Annotated assessable same-gene alternatives in paired original-TRAIN chain changes."""
from collections import defaultdict
import json
import os
import time
import diagnose_full_grammar_local_chain_transitions as T

C = T.G.C
OUT = T.ROOT/'outputs/M28-FULL-GRAMMAR-LOCAL-ISOFORM-TRANSITIONS'
OLD_EVAL = T.ROOT/'outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-DIAGNOSTIC'


def alternatives(ref, predicted, by_gene):
    primary = C.S.key(ref)
    if primary in predicted:
        raise ValueError('Lost/gained wrong-side primary must not be predicted')
    emitted = (by_gene.get(ref['gene_id'],set())-{primary}) & predicted
    bounds = primary[2][0][0], primary[2][-1][1]
    same_ends = {key for key in emitted if (key[2][0][0],key[2][-1][1]) == bounds}
    return dict(same_gene_other_assessable_CDS_chain_emitted=bool(emitted),
        emitted_other_assessable_chain_count=len(emitted),
        same_gene_other_assessable_CDS_chain_with_same_ends=bool(same_ends))


def summarize(rows):
    return dict(reference=len(rows),
        references_with_same_gene_other_assessable_chain=sum(r['same_gene_other_assessable_CDS_chain_emitted'] for r in rows),
        references_with_same_gene_other_assessable_chain_and_same_ends=sum(
            r['same_gene_other_assessable_CDS_chain_with_same_ends'] for r in rows),
        references_without_same_gene_other_assessable_chain=sum(not r['same_gene_other_assessable_CDS_chain_emitted'] for r in rows),
        emitted_other_assessable_chain_count_sum=sum(r['emitted_other_assessable_chain_count'] for r in rows))


def fixture():
    def tx(cds,seqid='a',strand='-'):
        return dict(seqid=seqid,strand=strand,gene_id='g',CDS=[(a,b,0) for a,b in cds])
    ref = tx([(0,30)])
    alt = tx([(0,9),(21,30)])
    shifted = tx([(3,9),(21,30)])
    catalog = {'g':{C.S.key(t) for t in (ref,alt,shifted)}}
    result = alternatives(ref,{C.S.key(alt)},catalog)
    assert result == dict(same_gene_other_assessable_CDS_chain_emitted=True,
        emitted_other_assessable_chain_count=1,same_gene_other_assessable_CDS_chain_with_same_ends=True)
    result = alternatives(ref,{C.S.key(shifted)},catalog)
    assert result['same_gene_other_assessable_CDS_chain_emitted'] and not result['same_gene_other_assessable_CDS_chain_with_same_ends']
    for wrong in (tx([(0,9),(21,30)],seqid='b'),
                  tx([(0,9),(21,30)],strand='+')):
        assert not alternatives(ref,{C.S.key(wrong)},catalog)['same_gene_other_assessable_CDS_chain_emitted']
    assert not alternatives(ref,{C.S.key(alt)},{'other':{C.S.key(alt)}})['same_gene_other_assessable_CDS_chain_emitted']
    try:
        alternatives(ref,{C.S.key(ref)},catalog)
    except ValueError:
        pass
    else:
        raise AssertionError('A correct primary cannot be counted as an alternative-only transition')
    return 'same_gene_same_ends_shifted_ends_strand_contig_and_primary_exclusion_fixtures_passed'


def main():
    began = time.perf_counter()
    checked = fixture()
    source = json.loads((T.OUT/'summary.json').read_text())
    old = T.completed(OLD_EVAL,T.G.I.CHECKPOINT)
    new = T.completed(T.EVAL,T.CHECKPOINT)
    if ((T.OUT/'STATUS').read_text().strip() != 'COMPLETED' or not source['completed']
            or source['job_id'] != '13453930' or source['DEV_used'] or source['sealed_data_used']
            or source['new_checkpoint'] != str(T.CHECKPOINT) or source['optimizer_steps_this_run'] != 0
            or old['job_id'] != '13443329' or new['job_id'] != '13451623'
            or old['reference_policy'] != source['reference_policy'] or new['reference_policy'] != source['reference_policy']):
        raise ValueError('Completed unchanged paired TRAIN sources and frozen ruler required')
    saved = T.read_rows(T.G.SOURCE/'reference_chain_profile.jsonl')
    transitions = T.read_rows(T.OUT/'reference_transitions.jsonl')
    if saved.keys() != transitions.keys() or len(saved) != 22688:
        raise ValueError('Entire frozen primary population must remain present')
    report = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],fixture=checked,
        source_transition_job=source['job_id'],source_old_evaluation_job=old['job_id'],
        source_new_evaluation_job=new['job_id'],reference_policy=source['reference_policy'],
        scope='original_six_TRAIN_saved_prediction_assessable_same_gene_isoform_membership_only',
        original_primary_reference_records=22688,primary_reference_changed=False,
        original_TRAIN_saved_input_FASTA_read_for_existing_assessability_definition=True,
        post_inference_reference_DNA_and_annotation_read=True,
        labels_or_references_provided_to_inference=False,predictions_modified=False,
        checkpoint_read=False,features_read=False,native_forward=False,optimizer_steps_this_run=0,
        DEV_used=False,sealed_data_used=False,
        membership_is_CDS_geometry_not_expression_evidence_UTR_or_all_isoform_recovery=True,
        nonmembership_is_not_proof_of_biologically_false_structure=True,
        pair_counts_are_not_globally_unique_predictions_or_additional_primary_TP=True,
        any_same_gene_alternative_need_not_be_the_maximum_CDS_overlap_prediction=True,
        lost_uses_LOCAL9216_predictions_gained_uses_original4608_predictions=True,
        old_anchor_has_fewer_updates_not_matched_update_objective_control=True,
        chromosomes={},per_species={},pooled={})
    details = []
    for species,seqid,length in T.G.I.SCOPE:
        reference = [r for r in saved.values() if r['seqid'] == seqid]
        sequences = C.M.read_allowed(T.G.I.OUT/seqid/'input.fasta',{seqid:length})
        annotation = T.ROOT/'data/m1_screen'/species/'reference.gff3'
        coding = C.S.E.parse_annotation(annotation,{seqid:length},protein_coding_only=True)
        with annotation.open() as handle:
            _,_,raw_cds,_ = C.P.raw_features(handle,seqid)
        primary, isoforms = C.select_reference(coding,raw_cds,sequences)
        if ({tx['id']:(tx['gene_id'],C.S.key(tx)) for tx in primary} !=
                {tx['id']:(tx['gene_id'],C.S.key(tx)) for tx in reference}
                or len(isoforms) != new['chromosomes'][seqid]['reference_assessable_isoform_records']):
            raise ValueError('Reconstructed primary or assessable isoform population differs from frozen evaluation')
        by_gene = defaultdict(set)
        for tx in isoforms:
            by_gene[tx['gene_id']].add(C.S.key(tx))
        reference_keys = {C.S.key(tx) for tx in reference}
        isoform_keys = {C.S.key(tx) for tx in isoforms}
        predictions = {}
        observed = {}
        for direction, root, evaluation in (('lost',T.RAW,new),('gained',T.G.I.OUT,old)):
            txs = C.S.E.primary_transcripts(C.S.E.parse_annotation(
                root/seqid/'prediction/full_grammar.gff3',{seqid:length}),complete_only=False)
            predicted = {C.S.key(tx) for tx in txs}
            score = C.S.prf(len(predicted & reference_keys),len(predicted),len(reference_keys))
            other = len((predicted & isoform_keys)-reference_keys)
            if (score != evaluation['chromosomes'][seqid]['exact_chain']
                    or other != evaluation['chromosomes'][seqid]['exact_other_assessable_isoform_not_primary_TP']):
                raise ValueError('Reconstructed exact primary or other-isoform score differs from completed ruler')
            predictions[direction] = predicted
            observed[direction] = dict(exact_chain=score,exact_other_assessable_isoform_not_primary_TP=other)
        for ref in reference:
            transition = transitions[T.key(ref)]
            direction = transition['transition']
            if direction not in predictions:
                continue
            details.append(dict(species=species,seqid=seqid,id=ref['id'],gene_id=ref['gene_id'],strand=ref['strand'],
                direction=direction,wrong_side_category=transition['new_category' if direction == 'lost' else 'old_category'],
                **alternatives(ref,predictions[direction],by_gene)))
        report['chromosomes'][seqid] = dict(species=species,reference_primary_records=len(primary),
            reference_assessable_isoform_records=len(isoforms),reconciled_existing_measurements=observed)
        del sequences, coding, raw_cds
    for label in (*source['per_species'],'pooled'):
        target = report['pooled'] if label == 'pooled' else report['per_species'].setdefault(label,{})
        paired = source['pooled']['paired'] if label == 'pooled' else source['per_species'][label]['paired']
        for direction in ('lost','gained'):
            rows = [r for r in details if r['direction'] == direction and (label == 'pooled' or r['species'] == label)]
            measured = summarize(rows)
            measured['by_wrong_side_error_category'] = {category:summarize([r for r in rows if r['wrong_side_category'] == category])
                for category in sorted({r['wrong_side_category'] for r in rows})}
            if len(rows) != paired['transitions'][direction]:
                raise ValueError('Lost or gained denominator differs from paired source')
            target[direction] = measured
    with (OUT/'reference_isoform_membership.jsonl').open('x') as handle:
        for row in details:
            handle.write(json.dumps(row,allow_nan=False)+'\n')
    report['seconds'] = time.perf_counter()-began
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps(dict(completed=True,seconds=report['seconds'],pooled=report['pooled'])),flush=True)


if __name__ == '__main__':
    main()
