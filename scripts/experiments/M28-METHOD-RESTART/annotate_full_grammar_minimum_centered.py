#!/usr/bin/env python3
"""Strict raw-FASTA fixed local-disagreement caller; no arrays or reference input."""
import argparse
import json
from pathlib import Path
import time

import torch
import annotate_full_grammar_joint_lora as E
import annotate_full_grammar_dense_qv as M
from src.m28.contig_emission_lattice import context_segments
from src.m28.full_grammar_energy_pool import minimum_centered_role_energy_pool

R = E.R
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_minimum_centered_raw.md'
from run_full_grammar_minimum_centered_train import IDENTITY
SOURCES = {'E_RANK4':E, 'M_DENSE':M}


def load_models(checkpoints):
    models, configurations = {}, {}
    for name, module in SOURCES.items():
        state = torch.load(checkpoints[name], map_location='cpu', weights_only=True)
        model = module.model_from_state(state, 'cuda')
        backbone, tokenizer = module.load_native_adapted(state)
        configurations[name] = state['configuration']
        models[name] = model, backbone, tokenizer
        del state
    return models, configurations


@torch.inference_mode()
def annotate_contig(models, seqid, sequence, ordinal):
    length = len(sequence)
    windows = [(a, min(a+R.M.WIDTH, length)) for a in R.M.starts(length)]
    segments = context_segments(length, windows)
    transcripts, scopes = [], []
    for strand in ('+', '-'):
        dna = sequence if strand == '+' else R.M.rc(sequence)
        values, encoding_seconds = {}, {}
        for name, module in SOURCES.items():
            began = time.perf_counter()
            joined = module.J.contig_emissions(*models[name], sequence, windows, strand, 'cuda')
            if joined['sequence'] != dna or joined['genomic_context_segments'] != segments:
                raise ValueError('Both native branches must share literal oriented DNA and row ownership')
            values[name] = joined['local_potentials']
            encoding_seconds[name] = time.perf_counter()-began
            del joined
        began = time.perf_counter()
        pooled = minimum_centered_role_energy_pool(values['E_RANK4'], values['M_DENSE'])
        pooling_seconds = time.perf_counter()-began
        del values
        began = time.perf_counter()
        decoded = R.decode_full_grammar_skyline(dna, pooled, 128)
        decode_seconds = time.perf_counter()-began
        del pooled
        for i, (chain, gain) in enumerate(zip(decoded['chains'], decoded['chain_gains'])):
            genomic = R.genomic_chain(chain, 0, length, strand)
            identifier = 'M28_MINROLE_'+str(ordinal)+'_'+('p' if strand == '+' else 'm')+'_'+str(i)
            tx = R.as_transcript(genomic, seqid, strand, identifier, gain)
            failures, _ = R.M.chain_checks(tx, sequence)
            if failures or not 1 <= len(tx['CDS']) <= 128:
                raise ValueError('Raw pooled complete CDS DNA, phase or cap failure')
            transcripts.append(tx)
        scopes.append(dict(seqid=seqid, strand=strand, windows=windows,
            genomic_context_segments=segments, context_policy=R.CONTEXT,
            native_encoder_conditions={name:module.CONDITION for name,module in SOURCES.items()},
            selected_complete_CDS_chains=len(decoded['chains']), gain=decoded['gain'],
            counts=decoded['counts'], encoding_seconds=encoding_seconds,
            pooling_seconds=pooling_seconds, decode_seconds=decode_seconds,
            tie_policy=decoded['tie_policy'], score_accumulation=decoded['score_accumulation'],
            backend='exact_exon_budget_skyline', pruning=decoded['pruning'],
            joint_whole_oriented_contig_decode=True, reference_used=False))
        print(json.dumps(dict(stage='pooled_raw_strand_complete', seqid=seqid, strand=strand,
            selected_chains=len(decoded['chains']), encoding_seconds=encoding_seconds,
            decode_seconds=decode_seconds)), flush=True)
        del decoded
    return transcripts, scopes


def main(fasta, checkpoints, out):
    limit = 128*1024*1024
    out.mkdir(parents=True, exist_ok=False)
    began = time.perf_counter()
    models, configurations = load_models(checkpoints)
    if any(p.requires_grad for model,backbone,_ in models.values() for p in
           [*model.parameters(), *backbone.parameters()]):
        raise ValueError('Both restored inference branches must be gradient-free')
    model_loading_seconds = time.perf_counter()-began
    torch.cuda.reset_peak_memory_stats()
    began = time.perf_counter()
    transcripts, scopes, lengths = [], [], {}
    for ordinal, (seqid, sequence) in enumerate(R.fasta_records(fasta)):
        txs, records = annotate_contig(models, seqid, sequence, ordinal)
        transcripts.extend(txs); scopes.extend(records); lengths[seqid] = len(sequence)
    R.write_gff3(out/'full_grammar.gff3', transcripts, lengths, 'M28_MINROLE')
    report = dict(inference_complete=True, input_fasta=str(fasta),
        checkpoints={name:str(path) for name,path in checkpoints.items()},
        configurations=configurations, protocol=PROTOCOL, identity=IDENTITY,
        encoder_model=R.F.MODEL, encoder_revision=R.F.REVISION,
        lengths=lengths, scopes=scopes, primary_CDS_chain_predictions=len(transcripts),
        score_kind='fixed_minimum_centered_fifteen_role_whole_contig_gain',
        source_head_steps=9216, source_native_steps=4608, two_separate_native_backbones=True,
        encoder_parameters_requires_grad={name:sum(p.numel() for p in backbone.parameters() if p.requires_grad)
            for name,(_,backbone,_) in models.items()},
        reference_used=False, labels_or_feature_cache_used=False, saved_role_arrays_read=False,
        saved_prediction_GFF_read=False, observation_model_used_at_inference=False,
        backend='exact_exon_budget_skyline', optimizer_steps_this_run=0,
        output_gib_limit=.125, model_loading_seconds=model_loading_seconds,
        wall_loop_seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved(),
        limitations='primary_complete_CDS_only_no_UTR_all_isoforms_cross_strand_nonoverlap_or_accuracy_claim')
    report['artifact_bytes_before_summary'] = R.output_bytes(out, limit)
    with (out/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False); handle.write('\n')
    R.output_bytes(out, limit)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta', type=Path, required=True)
    parser.add_argument('--e-checkpoint', type=Path, required=True)
    parser.add_argument('--m-checkpoint', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(); torch.set_num_threads(2)
    main(args.fasta.resolve(), {'E_RANK4':args.e_checkpoint.resolve(), 'M_DENSE':args.m_checkpoint.resolve()},
         args.output_dir.resolve())
