#!/usr/bin/env python3
"""Actual dual-native raw caller versus independent fixed mean-energy decode."""
import json
import os
import time

import torch
import annotate_full_grammar_e_m_energy_pool as A
import preflight_full_grammar_gene_posterior as N
import readout_full_grammar_gene_posterior_train as P
import run_full_grammar_e_m_energy_pool_train as Q
from src.m28.full_grammar_context_view import context_geometry

ROOT = P.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-E-M-RAW-NATIVE'


def main():
    began = time.perf_counter()
    P.write(OUT/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], protocol=A.PROTOCOL,
        identity=A.IDENTITY, source_fit_jobs={name:job for name,(_,job) in P.MODELS.items()},
        views=N.VIEWS, checkpoint_read=True, native_forward=True, optimizer_steps_this_run=0,
        reference_used=False, DEV_used=False, sealed_data_used=False))
    if A.IDENTITY != Q.IDENTITY:
        raise ValueError('The raw caller must preserve the completed fixed .5/.5 energy-pool rule')
    Q.qualification('13551105')
    for name in P.MODELS:
        P.source(name)
    sequences, source_views = {}, {}
    for species, seqid, length, start in N.VIEWS:
        geometry = context_geometry(start, start+24576, length)
        a, b = geometry['genomic_context']
        originals = {name:P.C.M.read_allowed(module.OUT/seqid/'input.fasta', {seqid:length})[seqid]
                     for name,(module,_) in P.MODELS.items()}
        if originals['E_RANK4'] != originals['M_DENSE']:
            raise ValueError('Literal original TRAIN DNA must agree before slicing the fixed raw view')
        view = originals['E_RANK4'][a:b]; del originals
        if len(view) != 49152 or list(A.R.M.starts(len(view))) != [0,12288,24576]:
            raise ValueError('Reuse the two declared original49152bp views, not new loci')
        sequences[seqid] = view
        source_views[seqid] = dict(species=species, seqid=seqid, chromosome_length=length,
            genomic_context=[a,b], core_start=start, native_windows=geometry['local_windows'])
    sequences['toy14'] = 'ATGCGTNAGTTTAA'
    fasta = OUT/'input.fasta'
    with fasta.open('x') as handle:
        for seqid, sequence in sequences.items():
            written = 'atgcgtragtttaa' if seqid == 'toy14' else sequence
            handle.write('>'+seqid+'\n'+written+'\n')
    if dict(A.R.fasta_records(fasta)) != sequences:
        raise ValueError('Existing FASTA normalization and multiple-record roundtrip must replay')
    checkpoints = {name:module.CHECKPOINT for name,(module,_) in P.MODELS.items()}
    prediction = OUT/'prediction'
    A.main(fasta, checkpoints, prediction)
    raw = P.read(prediction/'summary.json')
    txs = P.C.S.E.primary_transcripts(P.C.S.E.parse_annotation(prediction/'full_grammar.gff3',
        {seqid:len(sequence) for seqid,sequence in sequences.items()}), complete_only=False)
    scores = {}
    with (prediction/'full_grammar.gff3').open() as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if fields[2] == 'mRNA':
                scores[P.C.S.E.parse_attrs(fields[8])['ID']] = fields[5]
    if len(txs) != raw['primary_CDS_chain_predictions'] or len({P.C.S.key(t) for t in txs}) != len(txs):
        raise ValueError('Exported raw GFF must preserve every unique complete predicted chain')
    for tx in txs:
        failures, _ = A.R.M.chain_checks(tx, sequences[tx['seqid']])
        if failures:
            raise ValueError('Exported literal DNA and phase failure')
    models, _ = A.load_models(checkpoints)
    records = []
    for ordinal, (seqid, sequence) in enumerate(sequences.items()):
        length = len(sequence)
        windows = [(a,min(a+A.R.M.WIDTH,length)) for a in A.R.M.starts(length)]
        for strand in ('+', '-'):
            arrays = []
            for name, module in A.SOURCES.items():
                joined = module.J.contig_emissions(*models[name], sequence, windows, strand, 'cuda')
                arrays.append(joined['local_potentials']); del joined
            # Independent promotion and expression; no new caller/pool-helper invocation.
            averaged = (arrays[0].double()+arrays[1].double())*0.5
            del arrays
            dna = sequence if strand == '+' else A.R.M.rc(sequence)
            decoded = A.R.decode_full_grammar_skyline(dna, averaged, 128)
            del averaged
            expected = [A.R.as_transcript(A.R.genomic_chain(chain,0,length,strand),seqid,strand,
                'M28_POOL_'+str(ordinal)+'_'+('p' if strand=='+' else 'm')+'_'+str(i),gain)
                for i,(chain,gain) in enumerate(zip(decoded['chains'],decoded['chain_gains']))]
            actual = [tx for tx in txs if tx['seqid']==seqid and tx['strand']==strand]
            by_id = {tx['id']:tx for tx in actual}
            if (set(by_id) != {tx['id'] for tx in expected}
                    or any(by_id[tx['id']]['CDS'] != tx['CDS'] or scores[tx['id']] != format(tx['score'],'.8g')
                           for tx in expected)):
                raise ValueError('Raw pooled output coordinates, ordered phases and exported gains differ')
            scope = next(r for r in raw['scopes'] if r['seqid']==seqid and r['strand']==strand)
            segments = [list(x) for x in context_geometry(12288,36864,49152)['base_owned_segments']] if length==49152 else [[0,0,length]]
            if (scope['gain'] != decoded['gain'] or scope['selected_complete_CDS_chains'] != len(expected)
                    or scope['windows'] != [list(w) for w in windows]
                    or scope['genomic_context_segments'] != segments
                    or not scope['joint_whole_oriented_contig_decode'] or scope['reference_used']):
                raise ValueError('Full native mean-energy gain, ownership and oriented coverage must replay')
            spans = sorted((tx['CDS'][0][0],tx['CDS'][-1][1]) for tx in actual)
            if any(b>c for (_,b),(c,_) in zip(spans,spans[1:])):
                raise ValueError('Exported same-strand complete-gene spans must remain nonoverlapping')
            records.append(dict(seqid=seqid,strand=strand,length=length,
                predictions=len(expected),gain=decoded['gain'],source_view=source_views.get(seqid),
                exact_independent_mean_energy_MAP_and_GFF=True, exported_DNA_and_phase_valid=True,
                original_context_ownership_preserved=True, same_strand_spans_nonoverlapping=True))
            print(json.dumps(records[-1]),flush=True)
            del decoded
    if (len(raw['scopes']) != 6 or raw['identity'] != A.IDENTITY or not raw['inference_complete']
            or any(raw[k] for k in ('reference_used','labels_or_feature_cache_used','saved_role_arrays_read',
                'saved_prediction_GFF_read','observation_model_used_at_inference','optimizer_steps_this_run'))
            or any(raw['encoder_parameters_requires_grad'].values())
            or not raw['two_separate_native_backbones']):
        raise ValueError('Actual inference must remain separate-native, raw-only and gradient-free')
    P.write(OUT/'summary.json', dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],protocol=A.PROTOCOL,
        identity=A.IDENTITY,source_fit_jobs={name:job for name,(_,job) in P.MODELS.items()},
        original_energy_pool_checks_job='13551105',records=records,
        fixed_original_TRAIN_orientations=4,short_normalized_orientations=2,
        original_TRAIN_DNA_equality_before_slicing=True,multi_record_raw_FASTA_roundtrip_verified=True,
        all_six_native_raw_outputs_match_independent_mean_energy_decode=True,
        raw_caller_summary=str(prediction/'summary.json'),
        primary_CDS_chain_predictions=len(txs),
        raw_peak_GPU_allocated_bytes=raw['peak_GPU_allocated_bytes'],
        raw_peak_GPU_reserved_bytes=raw['peak_GPU_reserved_bytes'],
        raw_model_loading_seconds=raw['model_loading_seconds'],raw_wall_loop_seconds=raw['wall_loop_seconds'],
        checkpoint_read=True,native_forward=True,reference_used=False,
        labels_or_feature_cache_used=False,saved_role_arrays_read=False,saved_prediction_GFF_read=False,
        weights_modified=False,checkpoints_written=False,optimizer_steps_this_run=0,
        original_predictions_modified=False,DEV_used=False,sealed_data_used=False,
        threshold_selected=False,frozen_DEV_family_changed=False,
        full_original_six_TRAIN_end_to_end_reencoded=False,
        interpretation='native_raw_integration_and_mean_energy_parity_not_accuracy_or_generalization',
        seconds=time.perf_counter()-began))


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
