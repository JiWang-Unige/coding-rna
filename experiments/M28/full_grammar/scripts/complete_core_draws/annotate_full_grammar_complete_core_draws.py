#!/usr/bin/env python3
"""Final importance-corrected complete-core9216 raw caller; no fit readiness or reference read."""
import argparse
import json
from pathlib import Path
import time
import torch
import annotate_full_grammar_coarse_channel as R

PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_complete_core_draws_plan.md'
IDENTITY = dict(R.IDENTITY,
    scope='original_two_species_six_TRAIN_chromosomes_complete_core4608_draws_4107_unique',
    model='final4608_M28Core_C0_2048_128_366191_parameters_FP32_complete_core_importance_continuation',
    unique_windows=4107,optimizer_steps=9216,additional_optimizer_steps=4608,
    starting_optimizer_steps=4608,starting_fit_job='13434603',
    starting_checkpoint='outputs/M28-FULL-GRAMMAR-COARSE-CHANNEL-FIRST-FIT/step_004608.pt',
    optimizer_state='restore_original_final4608_AdamW_moments_and_step',
    sampling='half_uniform_half_complete_fine_chain_incidence_with_replacement',
    uniform_mixture_probability=.5,sampling_seed=1,preparation_job='13462399',support_job='13462462',
    training_order='fixed_complete_core_preparation_plan4608_rows_in_original_stratum_order',
    importance='inverse_stratum_size_times_proposal_probability',
    formal_checkpoint='last_cumulative_step_009216_only_no_DEV_selection')


def validate_checkpoint(state):
    if (state['step'] != 9216 or state['epoch'] != 6 or state['protocol'] != PROTOCOL
            or state['encoder_revision'] != R.F.REVISION
            or any(state['configuration'].get(k) != v for k,v in IDENTITY.items())):
        raise ValueError('Matching complete-core-proposal final9216 checkpoint required')
    return state['configuration']


def model_from_state(state,device):
    validate_checkpoint(state)
    model = R.M28Core('C0',2048,128)
    model.load_state_dict(state['model'],strict=True)
    return model.eval().requires_grad_(False).to(device)


def main(fasta,checkpoint,out):
    limit = 128*1024*1024
    out.mkdir(parents=True,exist_ok=False)
    state = torch.load(checkpoint,map_location='cpu',weights_only=True)
    model = model_from_state(state,'cuda')
    configuration = state['configuration']
    del state
    backbone,tokenizer = R.load_native()
    torch.cuda.reset_peak_memory_stats()
    began = time.perf_counter()
    transcripts,scopes,lengths = [],[],{}
    for ordinal,(seqid,sequence) in enumerate(R.fasta_records(fasta)):
        txs,records = R.annotate_contig(model,backbone,tokenizer,seqid,sequence,ordinal,'cuda')
        transcripts.extend(txs)
        scopes.extend(records)
        lengths[seqid] = len(sequence)
    R.write_gff3(out/'full_grammar.gff3',transcripts,lengths,'M28_FULL')
    report = dict(inference_complete=True,input_fasta=str(fasta),checkpoint=str(checkpoint),
        configuration=configuration,protocol=PROTOCOL,encoder_model=R.F.MODEL,encoder_revision=R.F.REVISION,
        lengths=lengths,scopes=scopes,primary_CDS_chain_predictions=len(transcripts),
        score_kind='unit_fifteen_role_logit_minus_background_joint_oriented_contig_gain',
        reference_used=False,labels_or_feature_cache_used=False,observation_model_used_at_inference=False,
        backend='exact_exon_budget_skyline',optimizer_steps_this_run=0,output_gib_limit=.125,
        wall_loop_seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved(),
        limitations='primary_complete_CDS_only_no_UTR_all_isoforms_cross_strand_nonoverlap_or_accuracy_claim')
    report['artifact_bytes_before_summary'] = R.output_bytes(out,limit)
    with (out/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False)
        handle.write('\n')
    R.output_bytes(out,limit)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta',type=Path,required=True)
    parser.add_argument('--checkpoint',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    main(args.fasta,args.checkpoint,args.out)
