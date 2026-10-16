#!/usr/bin/env python3
"""Apply the pre-existing TRAIN geometry definitions to the fixed four DEV outputs."""
import json
import os
import time
import common_ruler as C
import infer_full_grammar_fixed_dev as I
import evaluate_full_grammar_dev_baselines as E
from profile_full_grammar_original_train_errors import feature_scores

OUT = I.ROOT/'outputs/M28-FULL-GRAMMAR-FIXED-DEV-LOCAL-GEOMETRY'


def main():
    began = time.perf_counter()
    source = json.loads((E.OUT/'summary.json').read_text())
    if ((E.OUT/'STATUS').read_text().strip() != 'COMPLETED' or not source['completed'] or
            source['job_id'] != '13449926' or source['checkpoint'] != str(I.CHECKPOINT) or
            not source['input_records_directly_equal'] or not source['DEV_used'] or source['sealed_data_used']):
        raise ValueError('Completed fixed DEV common-ruler comparison required')
    report = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], source_evaluation_job=source['job_id'],
        status='fixed_DEV_descriptive_local_geometry_not_selection_or_independent_generalization',
        checkpoint=source['checkpoint'], reference_policy=source['reference_policy'],
        DEV_used=True, sealed_data_used=False, optimizer_steps_this_run=0,
        native_forward=False, checkpoint_read=False, features_read=False, predictions_modified=False,
        primary_reference_unchanged=True,
        definition='same_preexisting_TRAIN_strand_sensitive_unique_CDS_start_stop_segment_and_intron_pair_geometry',
        events_ignore_phase_and_are_not_gene_accuracy=True, no_bp_tolerance=True,
        species={}, pooled={})
    for species,seqid,length in I.SCOPE:
        data = I.ROOT/'data/m1_screen'/species
        sequences = C.M.read_allowed(I.OUT/seqid/'input.fasta', {seqid:length})
        coding = C.S.E.parse_annotation(data/'reference.gff3', {seqid:length}, protein_coding_only=True)
        with (data/'reference.gff3').open() as handle:
            _, _, raw_cds, _ = C.P.raw_features(handle, seqid)
        refs, _ = C.select_reference(coding, raw_cds, sequences)
        if len(refs) != E.EXPECTED[species] or len({C.S.key(t) for t in refs}) != E.EXPECTED[species]:
            raise ValueError('Fixed primary reference denominator differs')
        item = dict(seqid=seqid, length=length, reference_primary_chains=len(refs), methods={})
        for method in E.METHODS:
            prior = source['species'][species]['methods'][method]
            txs = C.S.E.primary_transcripts(C.S.E.parse_annotation(
                I.ROOT/prior['raw_prediction_path'], {seqid:length}), complete_only=False)
            p, r = ({C.S.key(t) for t in rows} for rows in (txs,refs))
            if C.S.prf(len(p & r), len(p), len(r)) != prior['exact_chain']:
                raise ValueError('Post-hoc geometry parse changed a frozen chain score')
            item['methods'][method] = feature_scores(txs, refs)
        report['species'][species] = item
        print(json.dumps(dict(profiled_DEV=seqid, reference=len(refs))), flush=True)
    # Fixed contigs are distinct, so event sets cannot collide across the two species.
    for method in E.METHODS:
        report['pooled'][method] = {}
        for event in ('CDS_start','CDS_stop','CDS_segment','intron_pair'):
            rows = [item['methods'][method][event] for item in report['species'].values()]
            report['pooled'][method][event] = C.S.prf(
                *(sum(row[key] for row in rows) for key in ('tp','predicted','reference')))
    report['seconds'] = time.perf_counter()-began
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(dict(completed=True, seconds=report['seconds'])), flush=True)


if __name__ == '__main__':
    main()
