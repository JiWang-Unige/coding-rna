#!/usr/bin/env python3
"""Strict final9216 raw caller; no training module or annotation access."""
import json
import math
import time
import torch
import annotate_full_grammar_coarse_channel as R

ARMS = ('STRUCTURED_CONTINUE', 'LOCAL_CHANNEL', 'HALF_STRUCTURED_HALF_LOCAL')
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_objective_contrast_contract.md'
IDENTITY = dict(R.IDENTITY,
    model='final4608_M28Core_C0_2048_128_366191_parameters_FP32_continuation',
    optimizer_steps=9216, additional_optimizer_steps=4608, starting_optimizer_steps=4608,
    starting_checkpoint='outputs/M28-FULL-GRAMMAR-COARSE-CHANNEL-FIRST-FIT/step_004608.pt',
    starting_fit_job='13434603', optimizer_state='restore_original_final4608_AdamW_moments_and_step',
    loss='arm_specific_same_fine_coarse_channel_same_original_importance_and_known_denominator',
    formal_checkpoint='last_cumulative_step_009216_only_no_DEV_selection',
    arms=list(ARMS), mixed_structured_coefficient=0.5, mixed_local_coefficient=0.5,
    context='original_24576bp_core_plus_one_12288bp_stride_neighbor_each_side_raw_length_clipped_deduplicated',
    ownership='nearest_containing_genomic_window_center_per_base_lower_start_ties_same_as_context_segments',
    reference_view_expansion=False, optimizer='AdamW', lr=3e-4, weight_decay=.01,
    betas=[.9, .999], eps=1e-8, clip_norm=1.)


def validate_checkpoint(state, arm):
    if (arm not in ARMS or state['step'] != 9216 or state['epoch'] != 6
            or state['protocol'] != PROTOCOL or state['encoder_revision'] != R.F.REVISION
            or state['configuration'].get('arm') != arm
            or any(state['configuration'].get(k) != v for k, v in IDENTITY.items())):
        raise ValueError('Matching declared arm and final cumulative9216 checkpoint required')
    return state['configuration']


def main(fasta, checkpoint, out, arm, output_gib=0.125):
    if not math.isfinite(output_gib) or output_gib <= 0:
        raise ValueError('Positive output limit required')
    limit = int(output_gib*1024**3)
    out.mkdir(parents=True, exist_ok=False)
    state = torch.load(checkpoint, map_location='cpu', weights_only=True)
    configuration = validate_checkpoint(state, arm)
    model = R.M28Core('C0', 2048, 128)
    model.load_state_dict(state['model'], strict=True)
    del state
    model.eval().requires_grad_(False).to('cuda')
    backbone, tokenizer = R.load_native()
    torch.cuda.reset_peak_memory_stats()
    began = time.perf_counter()
    transcripts, scopes, lengths = [], [], {}
    for ordinal, (seqid, sequence) in enumerate(R.fasta_records(fasta)):
        txs, records = R.annotate_contig(model, backbone, tokenizer, seqid, sequence, ordinal, 'cuda')
        transcripts.extend(txs)
        scopes.extend(records)
        lengths[seqid] = len(sequence)
    R.write_gff3(out/'full_grammar.gff3', transcripts, lengths, 'M28_FULL')
    report = dict(inference_complete=True, input_fasta=str(fasta), checkpoint=str(checkpoint),
        configuration=configuration, protocol=PROTOCOL, arm=arm,
        encoder_model=R.F.MODEL, encoder_revision=R.F.REVISION, lengths=lengths, scopes=scopes,
        primary_CDS_chain_predictions=len(transcripts),
        score_kind='unit_fifteen_role_logit_minus_background_joint_oriented_contig_gain',
        reference_used=False, labels_or_feature_cache_used=False, observation_model_used_at_inference=False,
        backend='exact_exon_budget_skyline', optimizer_steps_this_run=0, output_gib_limit=output_gib,
        wall_loop_seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved(),
        limitations='primary_complete_CDS_only_no_UTR_all_isoforms_cross_strand_nonoverlap_or_accuracy_claim')
    report['artifact_bytes_before_summary'] = R.output_bytes(out, limit)
    with (out/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write('\n')
    R.output_bytes(out, limit)
