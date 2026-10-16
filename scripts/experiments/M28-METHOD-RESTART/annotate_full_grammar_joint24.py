#!/usr/bin/env python3
"""Joint24576 native encoding and the unchanged whole-contig gene-set caller."""
import json
import time
import torch
from torch.nn.attention import SDPBackend, sdpa_kernel
import annotate_full_grammar_coarse_channel as R

PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_joint24_fit.md'
IDENTITY = dict(R.IDENTITY,
    model='final4608_M28Core_C0_2048_128_366191_parameters_FP32_joint24576_continuation',
    optimizer_steps=9216, additional_optimizer_steps=4608, starting_optimizer_steps=4608,
    starting_fit_job='13434603',
    starting_checkpoint='outputs/M28-FULL-GRAMMAR-COARSE-CHANNEL-FIRST-FIT/step_004608.pt',
    optimizer_state='restore_original_final4608_AdamW_moments_and_step',
    encoder_blocks_bp=24576, encoder_token_bp=6,
    encoder_materialization='fresh_joint_native_BF16_then_FP32_clone_each_read_window',
    native_attention='original_bidirectional_SDPA_flash_only_no_RoPE_change',
    short_window_padding='same_original_right_N_to24576_after_orientation',
    formal_checkpoint='last_cumulative_step_009216_only_no_DEV_selection')


def validate_checkpoint(state):
    if (state['step'] != 9216 or state['epoch'] != 6 or state['protocol'] != PROTOCOL
            or state['encoder_revision'] != R.F.REVISION
            or any(state['configuration'].get(k) != v for k,v in IDENTITY.items())):
        raise ValueError('Matching declared joint24576 final9216 checkpoint required')
    return state['configuration']


@torch.inference_mode()
def encode_window(backbone, tokenizer, sequence, device):
    if not 0 < len(sequence) <= R.M.WIDTH:
        raise ValueError('Nonempty original-width oriented window required')
    padded = sequence+'N'*(R.M.WIDTH-len(sequence))
    ids = tokenizer(padded, add_special_tokens=False, return_tensors='pt', truncation=False)['input_ids']
    if tuple(ids.shape) != (1,4096):
        raise ValueError('All4096 six-mer tokens without specials or truncation required')
    ids = ids.to(device)
    with sdpa_kernel(SDPBackend.FLASH_ATTENTION):
        hidden = backbone(input_ids=ids, attention_mask=torch.ones_like(ids)).last_hidden_state
    if tuple(hidden.shape) != (1,4096,2048) or not bool(torch.isfinite(hidden).all()):
        raise ValueError('Finite full-window native features required')
    return hidden.to(torch.bfloat16).float()


def window_logits(model, backbone, tokenizer, sequence, device):
    # Clone outside native inference mode so C0 training can save its input.
    features = encode_window(backbone,tokenizer,sequence,device).clone()
    dna = R.one_hot(sequence,R.M.WIDTH)[None].to(device)
    return model(features,dna,[len(sequence)])['segmentation_logits'][0,:len(sequence)]


@torch.inference_mode()
def contig_emissions(model,backbone,tokenizer,sequence,locations,strand,device):
    def records():
        for i,(a,b) in enumerate(locations):
            oriented = sequence[a:b] if strand == '+' else R.M.rc(sequence[a:b])
            values = window_logits(model,backbone,tokenizer,oriented,device).cpu()
            yield dict(window_index=i,sequence=oriented,local_potentials=values)
            if (i+1)%128 == 0 or i+1 == len(locations):
                print(json.dumps(dict(stage='joint24576_encoding',strand=strand,
                    completed_windows=i+1,total_windows=len(locations))),flush=True)
    return R.stitch_contig_emissions(sequence,locations,strand,records())


@torch.inference_mode()
def annotate_contig(model,backbone,tokenizer,seqid,sequence,ordinal,device):
    locations = [(a,min(a+R.M.WIDTH,len(sequence))) for a in R.M.starts(len(sequence))]
    transcripts,scopes = [],[]
    for strand in ('+','-'):
        began = time.perf_counter()
        joined = contig_emissions(model,backbone,tokenizer,sequence,locations,strand,device)
        encoding_seconds = time.perf_counter()-began
        began = time.perf_counter()
        result = R.decode_full_grammar_skyline(joined['sequence'],joined['local_potentials'],128)
        decode_seconds = time.perf_counter()-began
        for i,(chain,gain) in enumerate(zip(result['chains'],result['chain_gains'])):
            genomic = R.genomic_chain(chain,0,len(sequence),strand)
            identifier = 'M28_FULL_'+str(ordinal)+'_'+('p' if strand == '+' else 'm')+'_'+str(i)
            tx = R.as_transcript(genomic,seqid,strand,identifier,gain)
            if (not 1 <= len(genomic) <= 128 or any(not 0 <= a < b <= len(sequence) for a,b in genomic)
                    or any(c-b < 5 for (_,b),(c,_) in zip(genomic,genomic[1:]))):
                raise ValueError('Original complete-chain bounds, intron minimum and cap required')
            failures,_ = R.M.chain_checks(tx,sequence)
            if failures:
                raise ValueError('Raw joint24576 DNA/phase failure: '+','.join(failures))
            transcripts.append(tx)
        scopes.append(dict(seqid=seqid,strand=strand,windows=locations,
            genomic_context_segments=joined['genomic_context_segments'],context_policy=R.CONTEXT,
            encoder_condition='JOINT24576',selected_complete_CDS_chains=len(result['chains']),
            gain=result['gain'],counts=result['counts'],encoding_seconds=encoding_seconds,
            decode_seconds=decode_seconds,tie_policy=result['tie_policy'],
            score_accumulation=result['score_accumulation'],backend='exact_exon_budget_skyline',
            pruning=result['pruning'],joint_whole_oriented_contig_decode=True,reference_used=False))
        print(json.dumps(dict(stage='strand_complete',seqid=seqid,strand=strand,
            decode_seconds=decode_seconds,selected_chains=len(result['chains']))),flush=True)
        del joined,result
    return transcripts,scopes


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
        txs,records = annotate_contig(model,backbone,tokenizer,seqid,sequence,ordinal,'cuda')
        transcripts.extend(txs)
        scopes.extend(records)
        lengths[seqid] = len(sequence)
    R.write_gff3(out/'full_grammar.gff3',transcripts,lengths,'M28_FULL')
    report = dict(inference_complete=True,input_fasta=str(fasta),checkpoint=str(checkpoint),
        configuration=configuration,protocol=PROTOCOL,encoder_model=R.F.MODEL,
        encoder_revision=R.F.REVISION,lengths=lengths,scopes=scopes,
        primary_CDS_chain_predictions=len(transcripts),
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
