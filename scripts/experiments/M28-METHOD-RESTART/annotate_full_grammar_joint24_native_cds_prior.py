#!/usr/bin/env python3
"""Joint24576 encoding plus both frozen native CDS heads, strict final9216 raw caller."""
import json
import time
import torch
import annotate_full_grammar_joint24 as J
import annotate_full_grammar_native_cds_prior as N
from native_cds_role_prior import NativeCDSRolePrior

R = J.R
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_joint24_native_prior_plan.md'
IDENTITY = dict(N.IDENTITY,
    model='final4608_C0_plus_zero_initialized_28_parameter_two_native_CDS_role_prior_joint24576',
    encoder_blocks_bp=24576, encoder_materialization=J.IDENTITY['encoder_materialization'],
    native_attention=J.IDENTITY['native_attention'],short_window_padding=J.IDENTITY['short_window_padding'])


def validate_checkpoint(state):
    if (state['step'] != 9216 or state['epoch'] != 6 or state['protocol'] != PROTOCOL
            or state['encoder_revision'] != R.F.REVISION
            or any(state['configuration'].get(k) != v for k,v in IDENTITY.items())):
        raise ValueError('Declared joint24576 dual-native-prior final9216 checkpoint required')
    return state['configuration']


def core_from_state(state,device):
    validate_checkpoint(state)
    core = NativeCDSRolePrior()
    core.load_state_dict(state['model'],strict=True)
    return core.eval().requires_grad_(False).to(device)


def main(fasta,checkpoint,out):
    limit = 128*1024*1024
    out.mkdir(parents=True,exist_ok=False)
    state = torch.load(checkpoint,map_location='cpu',weights_only=True)
    core = core_from_state(state,'cuda')
    configuration = state['configuration']
    del state
    full,tokenizer = N.load_native()
    adapter = N.NativeCDSHeadAdapter(core,full.score).eval()
    torch.cuda.reset_peak_memory_stats()
    began = time.perf_counter()
    transcripts,scopes,lengths = [],[],{}
    for ordinal,(seqid,sequence) in enumerate(R.fasta_records(fasta)):
        txs,records = J.annotate_contig(adapter,full.model,tokenizer,seqid,sequence,ordinal,'cuda')
        transcripts.extend(txs); scopes.extend(records); lengths[seqid] = len(sequence)
    R.write_gff3(out/'full_grammar.gff3',transcripts,lengths,'M28_FULL')
    report = dict(inference_complete=True,input_fasta=str(fasta),checkpoint=str(checkpoint),
        configuration=configuration,protocol=PROTOCOL,encoder_model=R.F.MODEL,encoder_revision=R.F.REVISION,
        lengths=lengths,scopes=scopes,primary_CDS_chain_predictions=len(transcripts),
        score_kind='unit_fifteen_role_logit_minus_background_joint_oriented_contig_gain',
        reference_used=False,labels_or_feature_cache_used=False,observation_model_used_at_inference=False,
        both_frozen_native_heads_used=True,native_prior_weights='strict_final_checkpoint',
        backend='exact_exon_budget_skyline',optimizer_steps_this_run=0,output_gib_limit=.125,
        wall_loop_seconds=time.perf_counter()-began,peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved(),
        limitations='primary_complete_CDS_only_no_UTR_all_isoforms_cross_strand_nonoverlap_or_accuracy_claim')
    report['artifact_bytes_before_summary'] = R.output_bytes(out,limit)
    with (out/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    R.output_bytes(out,limit)
