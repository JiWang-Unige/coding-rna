#!/usr/bin/env python3
"""Strict raw-FASTA caller with trained, centered complete-CDS codon scores."""
import json
import time
import torch
import annotate_full_grammar_fresh_joint24 as F
import annotate_full_grammar_joint24 as J
from src.m28.full_grammar_codon import CODONS, centered_codon_scores
from src.m28.full_grammar_codon_skyline import decode_full_grammar_codon

R = J.R
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_codon_factor_preparation.md'
IDENTITY = dict(F.IDENTITY,
    model='source4608_C0_fresh_joint24576_with_64_raw_complete_CDS_codon_parameters',
    codon_order=list(CODONS), codon_parameter_count=64,
    codon_centering='separate_61_sense_and_3_stop_means_removed',
    codon_energy='one_score_per_complete_spliced_CDS_codon_including_start_and_stop',
    codon_initialization='zero', codon_optimizer_steps=4608,
    optimizer_state='restore_original_head_AdamW_step4608_new_codon_group_step0',
    clip_parameters='joint_C0_and_raw_codon_coefficients_norm1')


def validate_checkpoint(state):
    if (state['step'] != 9216 or state['epoch'] != 6 or state['protocol'] != PROTOCOL
            or state['encoder_revision'] != R.F.REVISION
            or any(state['configuration'].get(k) != v for k, v in IDENTITY.items())):
        raise ValueError('Matching final9216 codon checkpoint required')
    raw = state['raw_codon_scores']
    centered_codon_scores(raw)
    if raw.dtype != torch.float32:
        raise ValueError('Saved raw codon parameters must retain their FP32 training dtype')
    return state['configuration']


def model_from_state(state, device):
    validate_checkpoint(state)
    model = R.M28Core('C0', 2048, 128)
    model.load_state_dict(state['model'], strict=True)
    return model.eval().requires_grad_(False).to(device), state['raw_codon_scores'].detach().cpu().clone()


@torch.inference_mode()
def annotate_contig(model, raw_scores, backbone, tokenizer, seqid, sequence, ordinal, device):
    locations = [(a, min(a + R.M.WIDTH, len(sequence))) for a in R.M.starts(len(sequence))]
    transcripts, scopes = [], []
    for strand in ('+', '-'):
        began = time.perf_counter()
        joined = J.contig_emissions(model, backbone, tokenizer, sequence, locations, strand, device)
        encoding_seconds = time.perf_counter() - began
        began = time.perf_counter()
        result = decode_full_grammar_codon(joined['sequence'], joined['local_potentials'], 128, raw_scores)
        decode_seconds = time.perf_counter() - began
        for i, (chain, gain) in enumerate(zip(result['chains'], result['chain_gains'])):
            genomic = R.genomic_chain(chain, 0, len(sequence), strand)
            identifier = 'M28_CODON_' + str(ordinal) + '_' + ('p' if strand == '+' else 'm') + '_' + str(i)
            tx = R.as_transcript(genomic, seqid, strand, identifier, gain)
            if (not 1 <= len(genomic) <= 128 or any(not 0 <= a < b <= len(sequence) for a, b in genomic)
                    or any(c - b < 5 for (_, b), (c, _) in zip(genomic, genomic[1:]))):
                raise ValueError('Original complete-chain bounds, intron minimum and cap required')
            failures, _ = R.M.chain_checks(tx, sequence)
            if failures:
                raise ValueError('Raw codon DNA/phase failure: ' + ','.join(failures))
            transcripts.append(tx)
        scopes.append(dict(seqid=seqid, strand=strand, windows=locations,
            genomic_context_segments=joined['genomic_context_segments'], context_policy=R.CONTEXT,
            encoder_condition='JOINT24576', selected_complete_CDS_chains=len(result['chains']),
            gain=result['gain'], counts=result['counts'], encoding_seconds=encoding_seconds,
            decode_seconds=decode_seconds, tie_policy=result['tie_policy'],
            score_accumulation=result['score_accumulation'], backend='exact_codon_exon_budget_skyline',
            pruning=result['pruning'], joint_whole_oriented_contig_decode=True, reference_used=False))
        del joined, result
    return transcripts, scopes


def main(fasta, checkpoint, out):
    limit = 128*1024*1024
    out.mkdir(parents=True, exist_ok=False)
    state = torch.load(checkpoint, map_location='cpu', weights_only=True)
    model, raw_scores = model_from_state(state, 'cuda')
    configuration = state['configuration']
    del state
    backbone, tokenizer = R.load_native()
    torch.cuda.reset_peak_memory_stats()
    began = time.perf_counter()
    transcripts, scopes, lengths = [], [], {}
    for ordinal, (seqid, sequence) in enumerate(R.fasta_records(fasta)):
        txs, records = annotate_contig(model, raw_scores, backbone, tokenizer, seqid, sequence, ordinal, 'cuda')
        transcripts.extend(txs)
        scopes.extend(records)
        lengths[seqid] = len(sequence)
    R.write_gff3(out/'full_grammar.gff3', transcripts, lengths, 'M28_CODON')
    report = dict(inference_complete=True, input_fasta=str(fasta), checkpoint=str(checkpoint),
        configuration=configuration, protocol=PROTOCOL, encoder_model=R.F.MODEL,
        encoder_revision=R.F.REVISION, lengths=lengths, scopes=scopes,
        primary_CDS_chain_predictions=len(transcripts),
        raw_codon_scores=raw_scores.tolist(), centered_codon_scores=centered_codon_scores(raw_scores).tolist(),
        score_kind='unit_fifteen_role_logit_minus_background_plus_complete_spliced_CDS_codon_energy',
        reference_used=False, labels_or_feature_cache_used=False, observation_model_used_at_inference=False,
        backend='exact_codon_exon_budget_skyline', optimizer_steps_this_run=0,
        output_gib_limit=.125, wall_loop_seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved(),
        limitations='primary_complete_CDS_only_no_UTR_all_isoforms_cross_strand_nonoverlap_or_accuracy_claim')
    report['artifact_bytes_before_summary'] = R.output_bytes(out, limit)
    with (out/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write('\n')
    R.output_bytes(out, limit)
