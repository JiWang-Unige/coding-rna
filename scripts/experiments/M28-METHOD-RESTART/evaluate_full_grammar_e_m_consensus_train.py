#!/usr/bin/env python3
"""Fixed post-outcome TRAIN CDS-geometry intersection; no ref-informed selection."""
import json
import os
import time

import common_ruler as C
import evaluate_full_grammar_original_train as F
import evaluate_native_train_baselines as N
import infer_full_grammar_joint_lora_train as E
import infer_full_grammar_dense_qv_train as M

ROOT = E.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-E-M-CONSENSUS-TRAIN'
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_e_m_consensus_train.md'
ARMS = {'E_RANK4': (E, '13499817', '13499823'), 'M_DENSE': (M, '13510324', '13510351')}
ROW_FIELDS = ('exact_chain', 'CDS_base_unstranded_union', 'background',
    'prediction_records', 'duplicate_chain_records', 'exact_other_assessable_isoform_not_primary_TP')


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')


def main():
    began = time.perf_counter()
    reports, raw_reports = {}, {}
    for name, (module, fit_job, evaluation_job) in ARMS.items():
        module.source(fit_job); raw_reports[name] = module.reconcile_inputs()
        report = read(module.EVAL/'summary.json')
        receipt = read(module.EVAL/'paired_input_match.json')
        if ((module.EVAL/'STATUS').read_text().strip() != 'COMPLETED' or not report['completed']
                or report['job_id'] != evaluation_job or report['checkpoint'] != str(module.CHECKPOINT)
                or report['reference_used_in_inference'] or not report['unsupported_reference_strata_retained']
                or report['DEV_used'] or report['sealed_data_used'] or report['optimizer_steps_this_run']
                or receipt['adapted_fit_job'] != fit_job or receipt['protocol'] != module.A.PROTOCOL
                or not receipt['original_six_TRAIN_inputs_directly_equal']
                or not receipt['comparison_before_reference_access']):
            raise ValueError('Completed fixed E/M original-TRAIN sources and prior DNA receipt required')
        reports[name] = report
    native = read(N.OUT/'summary.json')
    if ((N.OUT/'STATUS').read_text().strip() != 'COMPLETED' or not native['completed']
            or native['job_id'] != '13461307' or native['reference_primary_chains'] != 22688
            or native['reference_policy'] != reports['E_RANK4']['reference_policy']
            or native['DEV_used'] or native['sealed_data_used']):
        raise ValueError('Fixed native original-six TRAIN evaluation required')
    write(OUT/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL,
        source_fit_jobs={n:a[1] for n,a in ARMS.items()}, source_evaluation_jobs={n:a[2] for n,a in ARMS.items()},
        native_evaluation_job=native['job_id'], rule='exact_sequence_strand_all_CDS_intervals_intersection_only',
        native_forward=False, checkpoint_read=False, optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False))
    chromosomes, agreement = {}, {}
    for species, seqid, length in E.I.SCOPE:
        sequences = {name:C.M.read_allowed(module.OUT/seqid/'input.fasta', {seqid:length})[seqid]
                     for name,(module,_,_) in ARMS.items()}
        if sequences['E_RANK4'] != sequences['M_DENSE']:
            raise ValueError('Both original caller FASTAs must be directly identical before truth access')
        sequence = sequences['E_RANK4']; del sequences
        txs, keyed = {}, {}
        for name, (module, _, _) in ARMS.items():
            txs[name] = C.S.E.primary_transcripts(C.S.E.parse_annotation(
                module.OUT/seqid/'prediction/full_grammar.gff3', {seqid:length}), complete_only=False)
            keyed[name] = {C.S.key(t):t for t in txs[name]}
            if len(txs[name]) != len(keyed[name]) or len(txs[name]) != raw_reports[name][seqid]['primary_CDS_chain_predictions']:
                raise ValueError('Original unique prediction count must replay')
            for tx in txs[name]:
                failures, _ = C.M.chain_checks(tx, sequence)
                if failures or len(tx['CDS']) > 128:
                    raise ValueError('Original predictions must retain literal DNA, phase and cap validity')
        common_keys = set(keyed['E_RANK4']) & set(keyed['M_DENSE'])
        common = [tx for tx in txs['E_RANK4'] if C.S.key(tx) in common_keys]
        # The operation is complete before reference access; no truth/gain filter.
        for strand in ('+', '-'):
            spans = sorted((t['CDS'][0][0], t['CDS'][-1][1]) for t in common if t['strand'] == strand)
            if any(b > c for (_,b),(c,_) in zip(spans, spans[1:])):
                raise ValueError('Intersection must remain a compatible subset of each complete gene set')
        data = ROOT/'data/m1_screen'/species
        coding = C.S.E.parse_annotation(data/'reference.gff3', {seqid:length}, protein_coding_only=True)
        annotation = C.S.E.parse_annotation(data/'reference.gff3', {seqid:length})
        with (data/'reference.gff3').open() as handle:
            _, _, raw_cds, _ = C.P.raw_features(handle, seqid)
        refs, isoforms = C.select_reference(coding, raw_cds, {seqid:sequence})
        rk, ik = ({C.S.key(t) for t in block} for block in (refs, isoforms))
        if len(refs) != len(rk) or any(not C.R7.phase_chain_consistent(t) for t in refs):
            raise ValueError('Full assessable primary population and reference phases required')
        for name in ARMS:
            replay = C.score(txs[name], refs, isoforms, length, annotation)
            old = reports[name]['chromosomes'][seqid]
            if any(replay[field] != old[field] for field in ROW_FIELDS):
                raise ValueError('Both original full endpoint/base/background/isoform rows must replay')
            if N.reference_metadata(old) != N.reference_metadata(native['methods']['ANNEVO']['chromosomes'][seqid]):
                raise ValueError('Native and parent reference population/strata/background must agree')
        row = C.score(common, refs, isoforms, length, annotation)
        row['prediction_nonexact_primary'] = row.pop('reference_nonexact_primary')
        row.update(species=species, seqid=seqid, length=length,
            reference_primary_records=len(refs), reference_assessable_isoform_records=len(isoforms),
            retained_unsupported_reference_strata=reports['E_RANK4']['chromosomes'][seqid]['retained_unsupported_reference_strata'],
            phase_checked_prediction_records=len(common))
        chromosomes[seqid] = row
        agreement[seqid] = {group:dict(unique_predictions=len(keys), exact_primary=len(keys & rk),
            exact_other_assessable_isoform=len((keys & ik)-rk), neither_exact_primary_nor_any_assessable_isoform=len(keys-ik))
            for group,keys in (('common',common_keys), ('E_only',set(keyed['E_RANK4'])-common_keys),
                               ('M_only',set(keyed['M_DENSE'])-common_keys))}
        if row['exact_chain']['tp'] != agreement[seqid]['common']['exact_primary']:
            raise ValueError('Consensus matching and agreement partition must agree')
        print(json.dumps(dict(completed_consensus_TRAIN=seqid)), flush=True)
    per_species = {s:F.aggregate([r for r in chromosomes.values() if r['species'] == s])
                   for s in ('arabidopsis_thaliana','oryza_sativa')}
    pooled = F.aggregate(list(chromosomes.values()))
    f1s = [r['exact_chain']['f1'] for r in per_species.values()]
    consensus = dict(chromosomes=chromosomes, per_species=per_species, pooled=pooled,
        macro_species_exact_chain_F1=sum(f1s)/2, worst_species_exact_chain_F1=min(f1s))
    comparisons = {}
    targets = dict(reports, **{name:native['methods'][name] for name in ('ANNEVO','Tiberius','Helixer')})
    for name, target in targets.items():
        comparisons[name] = dict(
            per_species_F1_pp={s:100*(per_species[s]['exact_chain']['f1']-target['per_species'][s]['exact_chain']['f1']) for s in per_species},
            pooled_F1_pp=100*(pooled['exact_chain']['f1']-target['pooled']['exact_chain']['f1']),
            macro_species_F1_pp=100*(consensus['macro_species_exact_chain_F1']-target['macro_species_exact_chain_F1']),
            worst_species_F1_pp=100*(consensus['worst_species_exact_chain_F1']-target['worst_species_exact_chain_F1']))
    if pooled['exact_chain']['reference'] != 22688:
        raise ValueError('Retain all22688 primary references')
    report = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL,
        source_fit_jobs={n:a[1] for n,a in ARMS.items()}, source_evaluation_jobs={n:a[2] for n,a in ARMS.items()},
        native_evaluation_job=native['job_id'], consensus=consensus, agreement=agreement,
        consensus_minus_parent_or_native=comparisons, reference_policy=native['reference_policy'],
        direct_paired_DNA_equality_before_reference_access=True, full_parent_metrics_replayed=True,
        all_original_primary_references_retained=True, consensus_selected_without_reference_or_score=True,
        opposite_strand_overlaps_permitted=True, predictions_modified=False, new_GFF_written=False,
        native_forward=False, checkpoint_read=False, optimizer_steps_this_run=0,
        DEV_used=False, sealed_data_used=False, threshold_selected=False, frozen_DEV_family_changed=False,
        interpretation='post_outcome_TRAIN_two_model_intersection_diagnostic_not_independent_validation_or_new_raw_caller',
        seconds=time.perf_counter()-began)
    write(OUT/'summary.json', report)


if __name__ == '__main__': main()
