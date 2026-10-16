#!/usr/bin/env python3
"""Strict final9216 raw FASTA: original C0 plus the fixed native encoder adapters."""
import argparse
import json
from pathlib import Path
import time
import torch
import annotate_full_grammar_fresh_joint24 as B
import annotate_full_grammar_joint24 as J
from src.m28 import joint_encoder_lora as L

R = J.R
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_joint_lora_fit.md'
CONDITION = 'JOINT24576_LAST4_QV_R4_LORA'
IDENTITY = dict(B.IDENTITY, encoder_frozen=False, encoder_base_weights_frozen=True,
    model='source4608_C0_plus102400_joint_encoder_adapter_parameters_final9216',
    encoder_adaptation='last_four_layers22_25_q_proj_v_proj_only', adapter_rank=4,
    adapter_alpha=4, adapter_scale=1., adapter_seed=0, adapter_dtype='FP32',
    adapter_A='seed0_uniform_plus_minus_inverse_sqrt_input_dim', adapter_B='zero_initialization',
    adapter_parameters=102400, adapter_optimizer_steps=4608, adapter_dropout=0.,
    gradient_checkpointing=False,
    adapter_output='FP32_residual_rounded_to_BF16_added_to_BF16_projection',
    encoder_materialization='differentiable_joint24576_BF16_then_FP32_no_detach',
    optimizer_state='original_head_AdamW4608_restored_new_adapter_AdamW0_same_hyperparameters',
    formal_checkpoint='last_cumulative_head_step009216_adapter_step004608_only_no_DEV_selection')


def adapter_shapes():
    return {f'layers.{index}.self_attn.{projection}.{key}':
        (4, 2048) if key == 'A' else (2048 if projection == 'q_proj' else 256, 4)
        for index in L.LAYERS for projection in L.PROJECTIONS for key in ('A', 'B')}


def validate_adapter_state(values):
    shapes = adapter_shapes()
    if set(values) != set(shapes):
        raise ValueError('All sixteen fixed native adapter tensors required')
    if any(tuple(values[k].shape) != shape or values[k].dtype != torch.float32
            or not bool(torch.isfinite(values[k]).all()) for k, shape in shapes.items()):
        raise ValueError('Matching finite FP32 native adapter geometry required')


def validate_checkpoint(state):
    if (state['step'] != 9216 or state['epoch'] != 6 or state['adapter_optimizer_step'] != 4608
            or state['protocol'] != PROTOCOL or state['encoder_revision'] != R.F.REVISION
            or any(state['configuration'].get(k) != value for k, value in IDENTITY.items())):
        raise ValueError('Matching final head9216 and adapter4608 checkpoint required')
    validate_adapter_state(state['adapters'])
    return state['configuration']


def model_from_state(state, device):
    validate_checkpoint(state)
    model = R.M28Core('C0', 2048, 128)
    model.load_state_dict(state['model'], strict=True)
    return model.eval().requires_grad_(False).to(device)


def load_native_adapted(state):
    validate_checkpoint(state)
    backbone, tokenizer = R.load_native()
    adapters = L.attach_adapters(backbone)
    L.restore_adapters(adapters, state['adapters'])
    backbone.eval().requires_grad_(False)
    return backbone, tokenizer


def main(fasta, checkpoint, out):
    limit = 128*1024*1024
    out.mkdir(parents=True, exist_ok=False)
    state = torch.load(checkpoint, map_location='cpu', weights_only=True)
    model = model_from_state(state, 'cuda')
    configuration = state['configuration']
    backbone, tokenizer = load_native_adapted(state)
    del state
    torch.cuda.reset_peak_memory_stats()
    began = time.perf_counter()
    transcripts, scopes, lengths = [], [], {}
    for ordinal, (seqid, sequence) in enumerate(R.fasta_records(fasta)):
        txs, records = J.annotate_contig(model, backbone, tokenizer, seqid, sequence, ordinal, 'cuda')
        for record in records:
            record['encoder_condition'] = CONDITION
        transcripts.extend(txs); scopes.extend(records); lengths[seqid] = len(sequence)
    R.write_gff3(out/'full_grammar.gff3', transcripts, lengths, 'M28_FULL')
    report = dict(inference_complete=True, input_fasta=str(fasta), checkpoint=str(checkpoint),
        configuration=configuration, protocol=PROTOCOL, encoder_model=R.F.MODEL,
        encoder_revision=R.F.REVISION, lengths=lengths, scopes=scopes,
        primary_CDS_chain_predictions=len(transcripts), encoder_condition=CONDITION,
        native_adapter_tensors_restored=16, native_adapter_parameters=L.PARAMETERS,
        encoder_parameters_requires_grad=sum(p.numel() for p in backbone.parameters() if p.requires_grad),
        score_kind='unit_fifteen_role_logit_minus_background_joint_oriented_contig_gain',
        reference_used=False, labels_or_feature_cache_used=False, observation_model_used_at_inference=False,
        backend='exact_exon_budget_skyline', optimizer_steps_this_run=0, output_gib_limit=.125,
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
    parser.add_argument('--checkpoint', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    torch.set_num_threads(2)
    main(args.fasta.resolve(), args.checkpoint.resolve(), args.output_dir.resolve())
