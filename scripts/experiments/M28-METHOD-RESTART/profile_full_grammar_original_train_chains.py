#!/usr/bin/env python3
"""Describe fixed original-TRAIN chains, without rerunning or changing inference."""
from collections import Counter
import json
import os
import time
import common_ruler as C
import infer_full_grammar_original_train as I
import query_training_census as Q

OUT = I.ROOT/'outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-CHAIN-PROFILE'
DIAGNOSTIC = I.ROOT/'outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-DIAGNOSTIC'


def exposure_keys():
    """Global geometry of R4 positive_ids in the original draws; not all labels."""
    draws = [json.loads(line) for line in Q.DRAWS.read_text().splitlines() if line.strip()]
    wanted = {d['window_id'] for d in draws}
    if (len(draws) != 1536 or len(wanted) != 1491 or
            Counter(d['species'] for d in draws) != Counter({s:768 for s in Q.ALLOW})):
        raise ValueError('Original draw counts differ')
    windows = {}
    catalogs = {}
    counts = {}
    for species, allowed in Q.ALLOW.items():
        selected = {}
        for w in Q.train_records(Q.LABELS/species/'windows_train.jsonl'):
            if w['id'] not in wanted:
                continue
            if (w['species'] != species or w['seqid'] not in allowed or
                    w['label_policy'] != 'CDS_assessable_primary_R4' or w['id'] in windows):
                raise ValueError('Original drawn TRAIN window identity differs')
            selected[w['id']] = w
            windows[w['id']] = w
        needed = {ident for w in selected.values() for ident in w['positive_ids']}
        catalog = {}
        for tx in Q.train_records(Q.LABELS/species/'transcripts.jsonl'):
            if tx['id'] not in needed:
                continue
            if (tx['id'] in catalog or tx['seqid'] not in allowed or
                    not tx['eligible_positive'] or not tx['training_primary']):
                raise ValueError('Original R4 positive catalog identity differs')
            catalog[tx['id']] = tx
        if set(catalog) != needed:
            raise ValueError('Original R4 positive catalog join incomplete')
        for w in selected.values():
            for ident in w['positive_ids']:
                tx = catalog[ident]
                if (tx['seqid'], tx['strand']) != (w['seqid'], w['strand']):
                    raise ValueError('Original positive/window chromosome or strand differs')
        catalogs[species] = {C.S.key(tx) for tx in catalog.values()}
        counts[species] = dict(original_draws=768, unique_drawn_windows=len(selected),
                              unique_R4_positive_ids=len(needed),
                              unique_R4_positive_chains=len(catalogs[species]))
    if set(windows) != wanted:
        raise ValueError('Original draw/window join incomplete')
    for d in draws:
        w = windows[d['window_id']]
        if (w['species'], w['seqid'], w['strand']) != (d['species'], *d['stratum']):
            raise ValueError('Original draw/window identity differs')
    return catalogs, counts


def describe(predicted, reference, positive):
    result = {'exact_chain': C.S.prf(len(predicted & reference), len(predicted), len(reference))}
    def stratified(predicate):
        p = {key for key in predicted if predicate(key)}
        r = {key for key in reference if predicate(key)}
        return C.S.prf(len(p & r), len(p), len(r))
    result['CDS_segment_count'] = {
        'single': stratified(lambda key: len(key[2]) == 1),
        'multiple': stratified(lambda key: len(key[2]) > 1)}
    result['span_longer_than_native_window'] = stratified(
        lambda key: key[2][-1][1]-key[2][0][0] > C.M.WIDTH)
    result['R4_positive_ids_in_original_draws'] = {}
    for label, refs in (('present', reference & positive), ('absent', reference-positive)):
        tp = len(refs & predicted)
        result['R4_positive_ids_in_original_draws'][label] = dict(
            reference=len(refs), tp=tp, recall=tp/len(refs) if refs else None,
            by_CDS_segment_count={})
        for group, predicate in (('single', lambda key: len(key[2]) == 1),
                                 ('multiple', lambda key: len(key[2]) > 1)):
            subset = {key for key in refs if predicate(key)}
            hits = len(subset & predicted)
            result['R4_positive_ids_in_original_draws'][label]['by_CDS_segment_count'][group] = dict(
                reference=len(subset), tp=hits, recall=hits/len(subset) if subset else None)
    result['positive_metadata_chains_not_primary_reference'] = len(positive-reference)
    return result


def fixture():
    # Two contigs and a negative strand must not collide with equal CDS boundaries.
    a = ('a', '+', ((0,9),))
    b = ('a', '+', ((12,18),(30,36)))
    c = ('a', '-', ((0,9),))
    d = ('b', '+', ((0,9),))
    long = ('b', '-', ((100,106),(25000,25006)))
    result = describe({a,b,c,long}, {a,b,d,long}, {b,d})
    assert result['exact_chain'] == C.S.prf(3,4,4)
    assert result['CDS_segment_count']['single'] == C.S.prf(1,2,2)
    assert result['CDS_segment_count']['multiple'] == C.S.prf(2,2,2)
    assert result['span_longer_than_native_window'] == C.S.prf(1,1,1)
    assert result['R4_positive_ids_in_original_draws']['present']['recall'] == .5
    assert result['R4_positive_ids_in_original_draws']['absent']['recall'] == 1
    return 'passed_geometry_strand_contig_group_and_exposure_fixture'


def main():
    began = time.perf_counter()
    checked = fixture()
    original = json.loads((DIAGNOSTIC/'summary.json').read_text())
    if ((DIAGNOSTIC/'STATUS').read_text().strip() != 'COMPLETED' or
            not original['completed'] or original['checkpoint'] != str(I.CHECKPOINT) or
            original['DEV_used'] or original['sealed_data_used']):
        raise ValueError('Completed fixed original-TRAIN diagnostic required')
    positives, metadata = exposure_keys()
    report = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], fixture=checked,
        scope='fixed_original_six_TRAIN_chromosomes_output_diagnosis_only',
        checkpoint=str(I.CHECKPOINT), source_diagnostic_job=original['job_id'],
        reference_policy=original['reference_policy'],
        DEV_used=False, sealed_data_used=False, optimizer_steps_this_run=0,
        checkpoint_read=False, features_read=False, native_forward=False, predictions_modified=False,
        primary_reference_unchanged=True, native_window_bp=C.M.WIDTH,
        exposure_definition='chain_geometry_present_in_R4_positive_ids_of_original_1536_draws',
        exposure_is_not_complete_fine_label_or_all_coarse_context_encoder_exposure=True,
        absent_exposure_group_is_not_untrained_or_independent_test=True,
        long_span_is_not_automatically_unsupported_by_whole_contig_grammar=True,
        metadata=metadata, chromosomes={}, per_species={})
    species_sets = {}
    with (OUT/'reference_chain_profile.jsonl').open('x') as handle:
        for species, seqid, length in I.SCOPE:
            folder = I.OUT/seqid
            data = I.ROOT/'data/m1_screen'/species
            sequences = C.M.read_allowed(folder/'input.fasta', {seqid:length})
            coding = C.S.E.parse_annotation(data/'reference.gff3', {seqid:length}, protein_coding_only=True)
            with (data/'reference.gff3').open() as source:
                _, _, raw_cds, _ = C.P.raw_features(source, seqid)
            refs, _ = C.select_reference(coding, raw_cds, sequences)
            txs = C.S.E.primary_transcripts(C.S.E.parse_annotation(
                folder/'prediction/full_grammar.gff3', {seqid:length}), complete_only=False)
            pk, rk = {C.S.key(t) for t in txs}, {C.S.key(t) for t in refs}
            local_positive = {key for key in positives[species] if key[0] == seqid}
            row = describe(pk, rk, local_positive)
            if row['exact_chain'] != original['chromosomes'][seqid]['exact_chain']:
                raise ValueError('Profile disagrees with completed diagnostic; do not change reference')
            row.update(species=species, seqid=seqid)
            report['chromosomes'][seqid] = row
            p, r, e = species_sets.setdefault(species, (set(),set(),set()))
            p.update(pk); r.update(rk); e.update(local_positive)
            for tx in refs:
                key = C.S.key(tx)
                handle.write(json.dumps(dict(species=species, seqid=seqid, id=tx['id'],
                    gene_id=tx['gene_id'], strand=tx['strand'], CDS=tx['CDS'],
                    CDS_segments=len(tx['CDS']), span_bp=tx['CDS'][-1][1]-tx['CDS'][0][0],
                    exact_primary_hit=key in pk, original_draw_R4_positive_chain=key in local_positive))+'\n')
            print(json.dumps(dict(profiled_TRAIN=seqid, reference=len(rk))), flush=True)
    for species, sets in species_sets.items():
        report['per_species'][species] = describe(*sets)
        if report['per_species'][species]['exact_chain'] != original['per_species'][species]['exact_chain']:
            raise ValueError('Species profile count reconciliation failed')
    pooled = tuple(set().union(*(sets[i] for sets in species_sets.values())) for i in range(3))
    report['pooled'] = describe(*pooled)
    if report['pooled']['exact_chain'] != original['pooled']['exact_chain']:
        raise ValueError('Pooled profile count reconciliation failed')
    report['seconds'] = time.perf_counter()-began
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(dict(completed=True, seconds=report['seconds'])), flush=True)


if __name__ == '__main__':
    main()
