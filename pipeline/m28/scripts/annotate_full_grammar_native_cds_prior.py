#!/usr/bin/env python3
"""Strict final9216 native-CDS-prior raw caller; no labels or fit imports."""
import json
import math
import time

import torch
from torch import nn
from transformers import AutoModelForTokenClassification, AutoTokenizer
import annotate_full_grammar_coarse_channel as R
from native_cds_role_prior import NativeCDSRolePrior, native_cds_log_odds

PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_native_cds_prior_fit.md'
IDENTITY = dict(R.IDENTITY,
    model='final4608_C0_plus_zero_initialized_28_parameter_two_native_CDS_role_prior',
    starting_checkpoint='outputs/M28-FULL-GRAMMAR-COARSE-CHANNEL-FIRST-FIT/step_004608.pt',
    starting_fit_job='13434603', starting_optimizer_steps=4608,
    optimizer_steps=9216, additional_optimizer_steps=4608,
    optimizer_state='restore_old_C0_AdamW_step4608_then_append_new_prior_group_step0',
    native_head_indices=[0, 1], native_heads_frozen=True, head_selected_from_labels=False,
    prior='zero_initialized_Linear_2_to14_no_bias_pad_unchanged_background_role0',
    evidence='native_BF16_two_head_6x2_CDS_minus_NON_CODING_FP32_per_base',
    additional_parameters=28, total_trainable_parameters=366219,
    encoder_blocks_bp=6144, formal_checkpoint='last_cumulative_step_009216_only_no_DEV_selection',
    context='original_24576bp_core_plus_one_12288bp_stride_neighbor_each_side_raw_length_clipped_deduplicated',
    ownership='nearest_containing_genomic_window_center_per_base_lower_start_ties_same_as_context_segments',
    reference_view_expansion=False, optimizer='AdamW', lr=3e-4, weight_decay=.01,
    betas=[.9, .999], eps=1e-8, clip_norm=1.)


def load_native():
    tokenizer = AutoTokenizer.from_pretrained(R.F.MODEL, revision=R.F.REVISION,
        trust_remote_code=True, local_files_only=True)
    full = AutoModelForTokenClassification.from_pretrained(R.F.MODEL, revision=R.F.REVISION,
        trust_remote_code=True, local_files_only=True, attn_implementation='sdpa', torch_dtype=torch.bfloat16)
    if (full.k != 6 or full.num_labels != 2 or full.num_prediction_heads != 2 or full.feature_layer != -1
            or full.config.id2label != {0:'NON_CODING', 1:'CDS'} or not isinstance(full.score, nn.ModuleList)
            or len(full.score) != 2 or any(not isinstance(h, nn.Linear) or h.in_features != 2048
                or h.out_features != 12 or h.bias is not None for h in full.score)):
        raise ValueError('Both unchanged released six-slot native CDS heads required')
    return full.eval().requires_grad_(False).to('cuda'), tokenizer


class NativeCDSHeadAdapter(nn.Module):
    """Keep the original three-input neural caller, without serializing heads."""
    def __init__(self, core, native_heads):
        super().__init__()
        self.core = core
        self.native_heads = native_heads

    def forward(self, features, dna, valid_bases):
        return self.core(features, dna, valid_bases, native_cds_log_odds(self.native_heads, features))


def validate_checkpoint(state):
    if (state['step'] != 9216 or state['epoch'] != 6 or state['protocol'] != PROTOCOL
            or state['encoder_revision'] != R.F.REVISION
            or any(state['configuration'].get(k) != v for k, v in IDENTITY.items())):
        raise ValueError('Declared native-prior final cumulative9216 checkpoint required')
    return state['configuration']


def main(fasta, checkpoint, out, output_gib=0.125):
    if not math.isfinite(output_gib) or output_gib <= 0:
        raise ValueError('Positive finite output limit required')
    limit = int(output_gib*1024**3)
    out.mkdir(parents=True, exist_ok=False)
    state = torch.load(checkpoint, map_location='cpu', weights_only=True)
    configuration = validate_checkpoint(state)
    core = NativeCDSRolePrior()
    core.load_state_dict(state['model'], strict=True)
    del state
    core.eval().requires_grad_(False).to('cuda')
    full, tokenizer = load_native()
    model = NativeCDSHeadAdapter(core, full.score).eval()
    torch.cuda.reset_peak_memory_stats()
    began = time.perf_counter()
    transcripts, scopes, lengths = [], [], {}
    for ordinal, (seqid, sequence) in enumerate(R.fasta_records(fasta)):
        txs, records = R.annotate_contig(model, full.model, tokenizer, seqid, sequence, ordinal, 'cuda')
        transcripts.extend(txs); scopes.extend(records); lengths[seqid] = len(sequence)
    R.write_gff3(out/'full_grammar.gff3', transcripts, lengths, 'M28_FULL')
    report = dict(inference_complete=True, input_fasta=str(fasta), checkpoint=str(checkpoint),
        configuration=configuration, protocol=PROTOCOL, encoder_model=R.F.MODEL, encoder_revision=R.F.REVISION,
        lengths=lengths, scopes=scopes, primary_CDS_chain_predictions=len(transcripts),
        score_kind='unit_fifteen_role_logit_minus_background_joint_oriented_contig_gain',
        reference_used=False, labels_or_feature_cache_used=False, observation_model_used_at_inference=False,
        both_frozen_native_heads_used=True, native_prior_weights='strict_final_checkpoint',
        backend='exact_exon_budget_skyline', optimizer_steps_this_run=0, output_gib_limit=output_gib,
        wall_loop_seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved(),
        limitations='primary_complete_CDS_only_no_UTR_all_isoforms_cross_strand_nonoverlap_or_accuracy_claim')
    report['artifact_bytes_before_summary'] = R.output_bytes(out, limit)
    with (out/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False); handle.write('\n')
    R.output_bytes(out, limit)
