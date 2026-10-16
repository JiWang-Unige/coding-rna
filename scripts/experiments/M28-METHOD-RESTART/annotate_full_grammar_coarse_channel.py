#!/usr/bin/env python3
"""Strict raw-FASTA caller for the fixed final coarse-channel checkpoint.

No training readiness, catalog, observation likelihood, or old execution chain
is used at inference. Actual input/resource execution still needs its scope.
"""
import argparse
import json
import math
from pathlib import Path
import time

import torch
import manifest as M
import feature_smoke as F
from annotate_query import encode_window, fasta_records, output_bytes
from cache_full_train import load_native
from src.m28.core import M28Core
from src.m28.training import one_hot
from src.m28.contig_emission_lattice import stitch_contig_emissions
from src.m28.full_grammar_skyline import decode_full_grammar_skyline
from src.m28.c0_decode import genomic_chain
from src.m28.export import as_transcript, write_gff3

PROTOCOL='reports/M28-METHOD-RESTART/full_grammar_first_fit_contract.md'
CONTEXT='nearest_containing_genomic_window_center_per_base_lower_start_ties'
IDENTITY=dict(encoder_model=F.MODEL,encoder_revision=F.REVISION,encoder_frozen=True,
    model='fresh_seed0_M28Core_C0_2048_128_366191_parameters_FP32',seed=0,
    draws_per_pass=1536,unique_windows=1491,passes=3,optimizer_steps=4608,
    fine='unchanged_exact_role',
    coarse='normalized_three_class_correct_1_minus_epsilon_other_each_epsilon_over_two',
    unknown='likelihood_one',epsilon=0.1,epsilon_trainable=False,
    loss='importance_times_complete_gene_set_NLL_divided_by_max_1_fine_plus_coarse_core_rows',
    auxiliary_losses=False,max_exons=128,role_weights=1,null=0,score_accumulation='FP64',
    formal_checkpoint='last_step_004608_only_no_DEV_selection')


def validate_checkpoint(state):
    if (state['step']!=4608 or state['epoch']!=3 or state['protocol']!=PROTOCOL
            or state['encoder_revision']!=F.REVISION
            or any(state['configuration'].get(k)!=v for k,v in IDENTITY.items())):
        raise ValueError('Matching final step4608 fixed coarse-channel checkpoint required')
    return state['configuration']


@torch.inference_mode()
def contig_emissions(model,backbone,tokenizer,sequence,locations,strand,device):
    def records():
        for i,(a,b) in enumerate(locations):
            dna_string=sequence[a:b] if strand=='+' else M.rc(sequence[a:b])
            features=encode_window(backbone,tokenizer,dna_string,device)
            dna=one_hot(dna_string,M.WIDTH)[None].to(device)
            logits=model(features,dna,[b-a])['segmentation_logits'][0,:b-a].cpu()
            yield dict(window_index=i,sequence=dna_string,local_potentials=logits)
            if (i+1)%128==0 or i+1==len(locations):
                print(json.dumps(dict(stage='encoding',strand=strand,
                    completed_windows=i+1,total_windows=len(locations))),flush=True)
            del features,dna,logits
    return stitch_contig_emissions(sequence,locations,strand,records())


@torch.inference_mode()
def annotate_contig(model,backbone,tokenizer,seqid,sequence,ordinal,device):
    locations=[(a,min(a+M.WIDTH,len(sequence))) for a in M.starts(len(sequence))]
    transcripts,scopes=[],[]
    for strand in ('+','-'):
        began=time.perf_counter()
        joined=contig_emissions(model,backbone,tokenizer,sequence,locations,strand,device)
        encoding_seconds=time.perf_counter()-began; began=time.perf_counter()
        print(json.dumps(dict(stage='decoding',seqid=seqid,strand=strand,
            encoding_seconds=encoding_seconds)),flush=True)
        result=decode_full_grammar_skyline(joined['sequence'],joined['local_potentials'],128)
        decode_seconds=time.perf_counter()-began
        for i,(chain,gain) in enumerate(zip(result['chains'],result['chain_gains'])):
            genomic=genomic_chain(chain,0,len(sequence),strand)
            identifier='M28_FULL_'+str(ordinal)+'_'+('p' if strand=='+' else 'm')+'_'+str(i)
            tx=as_transcript(genomic,seqid,strand,identifier,gain)
            if (not 1<=len(genomic)<=128 or any(not 0<=a<b<=len(sequence) for a,b in genomic)
                    or any(c-b<5 for (_,b),(c,_) in zip(genomic,genomic[1:]))):
                raise ValueError('Complete in-bounds output chain and original intron/cap required')
            failures,_=M.chain_checks(tx,sequence)
            if failures: raise ValueError('Raw output DNA/phase failure: '+','.join(failures))
            transcripts.append(tx)
        scopes.append(dict(seqid=seqid,strand=strand,windows=locations,
            genomic_context_segments=joined['genomic_context_segments'],context_policy=CONTEXT,
            selected_complete_CDS_chains=len(result['chains']),gain=result['gain'],counts=result['counts'],
            encoding_seconds=encoding_seconds,decode_seconds=decode_seconds,
            tie_policy=result['tie_policy'],score_accumulation=result['score_accumulation'],
            backend='exact_exon_budget_skyline',pruning=result['pruning'],
            joint_whole_oriented_contig_decode=True,reference_used=False))
        print(json.dumps(dict(stage='strand_complete',seqid=seqid,strand=strand,
            decode_seconds=decode_seconds,selected_chains=len(result['chains']))),flush=True)
        del joined,result
    return transcripts,scopes


def main(fasta,checkpoint,out,output_gib):
    if not math.isfinite(output_gib) or output_gib<=0:
        raise ValueError('Positive finite explicit output budget required')
    limit=int(output_gib*1024**3); out.mkdir(parents=True,exist_ok=False)
    state=torch.load(checkpoint,map_location='cpu',weights_only=True)
    configuration=validate_checkpoint(state)
    model=M28Core('C0',2048,128)
    model.load_state_dict(state['model'],strict=True); del state
    model.eval().requires_grad_(False).to('cuda')
    backbone,tokenizer=load_native()
    torch.cuda.reset_peak_memory_stats(); began=time.perf_counter()
    transcripts,scopes,lengths=[],[],{}
    for ordinal,(seqid,sequence) in enumerate(fasta_records(fasta)):
        txs,records=annotate_contig(model,backbone,tokenizer,seqid,sequence,ordinal,'cuda')
        transcripts.extend(txs); scopes.extend(records); lengths[seqid]=len(sequence)
        print(json.dumps(dict(completed_contig=seqid,length=len(sequence),predictions=len(txs))),flush=True)
    write_gff3(out/'full_grammar.gff3',transcripts,lengths,'M28_FULL')
    report=dict(inference_complete=True,input_fasta=str(fasta),checkpoint=str(checkpoint),
        configuration=configuration,protocol=PROTOCOL,encoder_model=F.MODEL,encoder_revision=F.REVISION,
        lengths=lengths,scopes=scopes,primary_CDS_chain_predictions=len(transcripts),
        score_kind='unit_fifteen_role_logit_minus_background_joint_oriented_contig_gain',
        reference_used=False,labels_or_feature_cache_used=False,observation_model_used_at_inference=False,
        backend='exact_exon_budget_skyline',
        optimizer_steps_this_run=0,output_gib_limit=output_gib,
        wall_loop_seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved(),
        limitations='primary_complete_CDS_only_no_UTR_all_isoforms_cross_strand_nonoverlap_or_accuracy_claim')
    report['artifact_bytes_before_summary']=output_bytes(out,limit)
    with (out/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    output_bytes(out,limit)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta',type=Path,required=True)
    parser.add_argument('--checkpoint',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--output-gib-limit',type=float,required=True)
    args=parser.parse_args(); torch.set_num_threads(2)
    main(args.fasta.resolve(),args.checkpoint.resolve(),args.output_dir.resolve(),args.output_gib_limit)
