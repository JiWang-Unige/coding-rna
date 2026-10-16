#!/usr/bin/env python3
"""Conditional A41 raw caller; actual TRAIN PASS and separate execution scope required."""
import argparse
from collections import Counter
import json
import math
from pathlib import Path
import re
import time
import torch
import readout_emission_residual_train as D
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript, write_gff3
from src.m28.training import one_hot

ROOT,C,M=D.ROOT,D.R.C,D.R.C.I.M
FIT_JOB,READOUT_JOB,CHECK_JOB='13420820','13420821','13420824'
DEV_PROTOCOL='reports/M28-METHOD-RESTART/a41_emission_residual_dev_protocol.md'
SCORE_KIND='A41_FP64_frozen_A31_gain_plus_separate_residual_role_gain'
CPU_PROOFS=(
    'actual_original_draw_weights_counts_and_4608_updates_verified',
    'actual_saved_fixed_composite_loss_and_importance_arithmetic_verified',
    'actual_three_affine_checkpoints_and_optimizer_exposure_verified',
    'actual_three_frozen_source_states_equal_original_A3_and_excluded_from_optimizer_verified',
    'actual_saved_FP32_local_potentials_and_FP64_gain_arithmetic_verified',
    'actual_saved_labels_DP_baseline_denominators_and_five_conditions_recomputed',
)


def validate_preparation(fitted,prepared,checked):
    D.validate_fit(fitted)
    for report in (prepared,checked):
        if (report['completed'] is not True or report['configuration']!=D.CONFIG or
                report['protocol']!=D.PROTOCOL or report['windows']!=1491 or report['natural_candidates']!=16321 or
                report['potential_receipts']!=1491 or report['optimizer_steps_this_run']!=0 or
                any(report[k] is not False for k in ('new_candidate_search','new_R4_labelling',
                    'DEV_used','test_setaria_used','A33_outputs_opened'))):
            raise ValueError('Complete original TRAIN-only A41 preparation required')
    if (prepared['source_fit_job']!=FIT_JOB or prepared['source_checks_job']!='13420720' or
            prepared['source_preflight_job']!=D.F.PREFLIGHT_JOB or prepared['checkpoint']!=fitted['final_checkpoint'] or
            prepared['original_draws']!=1536 or prepared['reference_injected_predictions']!=0 or
            checked['source_jobs']!=dict(fit=FIT_JOB,readout=READOUT_JOB) or
            any(checked[k] is not True for k in CPU_PROOFS) or
            checked['model_or_backbone_forward'] is not False or checked['independent_neural_forward_replay'] is not False or
            prepared['label_catalog_opened'] is not False or prepared['new_backbone_forward'] is not False or
            any(prepared[k] is not True for k in ('all_scores_and_potential_files_closed_before_saved_labels_open',
                'original_natural_geometry_order_sources_unchanged','original_saved_labels_reused',
                'model_weights_unchanged','frozen_source_weights_unchanged','all_parameters_frozen_and_gradients_absent'))):
        raise ValueError('Actual matching residual, saved potentials and six CPU reconciliation proofs required')
    counts,refs=prepared['selected_counts_per_species'],prepared['eligible_reference_pairs']
    expected_refs={'arabidopsis_thaliana':{'eligible':2172,'five_or_more':854},
                   'oryza_sativa':{'eligible':744,'five_or_more':232}}
    old_counts={'arabidopsis_thaliana':(1531,347,548),'oryza_sativa':(478,251,152)}
    if (set(counts)!=set(old_counts) or refs!=expected_refs or
            checked['selected_counts_per_species']!=counts or checked['eligible_reference_pairs']!=refs or
            any(tuple(counts[s]['A31_A3'][k] for k in ('selected_correct','selected_known_nonmatch',
                'selected_five_or_more_correct'))!=v for s,v in old_counts.items())):
        raise ValueError('Original stronger A31 baseline and full species denominators required')
    gates=D.gate(counts)
    if (prepared['predeclared_TRAIN_gate']!=gates or checked['predeclared_TRAIN_gate']!=gates or
            prepared['TRAIN_preparation_supported'] is not True or checked['TRAIN_preparation_supported'] is not True or
            not all(gates.values())):
        raise ValueError('Original five-condition A41 TRAIN gate must actually pass')


def readiness():
    D.F.A.A.require_jobs([FIT_JOB,READOUT_JOB,CHECK_JOB,'13420720'])
    fitted=D.R.completed(D.F.OUT,FIT_JOB)
    validate_preparation(fitted,D.R.completed(D.OUT,READOUT_JOB),D.R.completed(D.CHECK,CHECK_JOB))
    return fitted


def load_models(fitted):
    generator=C.load_generator()
    state=torch.load(ROOT/fitted['final_checkpoint'],map_location='cpu',weights_only=True)
    D.validate_fit(state,state=True)
    if (state['frozen_source_checkpoint']!=D.CONFIG['frozen_source'] or
            set(state['frozen_source_model'])!=set(generator.state_dict()) or
            any(not torch.equal(v.cpu(),state['frozen_source_model'][k]) for k,v in generator.state_dict().items())):
        raise ValueError('Raw generator must equal the actual original frozen A41 fit source')
    scorer=D.F.A.new_model('cuda').eval().requires_grad_(False)
    scorer.load_state_dict(state['model'],strict=True)
    return generator,scorer


def verify_score_components(row):
    scores,base,residual=row['score'],row['frozen_path_gain'],row['residual_path_gain']
    if (row['null']!=0. or len(base)!=len(scores) or len(residual)!=len(scores) or
            any(not math.isfinite(v) for values in (scores,base,residual) for v in values) or
            any(a+b!=s for a,b,s in zip(base,residual,scores))):
        raise ValueError('Exact separate FP64 frozen plus residual potential arithmetic required')


@torch.inference_mode()
def local_record(generator,hidden,scorer,sequence,index,ident):
    chains,sources=C.free_pool(generator,hidden,sequence)
    local=generator.segmentation(hidden)
    result=scorer(hidden,local,chains,len(sequence))
    if (any(result[k].dtype!=torch.float64 or not bool(torch.isfinite(result[k]).all())
            for k in ('logits','frozen_path_gain','residual_path_gain')) or
            result['residual_local_logits'].dtype!=torch.float32 or
            not bool(torch.isfinite(result['residual_local_logits']).all()) or float(result['null'])!=0.):
        raise ValueError('Finite FP32 role corrections and separate FP64 potentials required')
    scores=result['logits'].cpu().tolist()
    row=dict(window_index=index,window_id=ident,chains=chains,score=scores,decoder_sources=sources,
        score_kind=SCORE_KIND,variant='A41',score_threshold=0.,
        frozen_path_gain=result['frozen_path_gain'].cpu().tolist(),
        residual_path_gain=result['residual_path_gain'].cpu().tolist(),null=float(result['null']),
        frozen_source_checkpoint=D.CONFIG['frozen_source'],reference_used=False,
        counts=dict(free_candidates=len(chains),nonpositive_candidates=sum(v<=0 for v in scores)))
    verify_score_components(row)
    return row


@torch.inference_mode()
def annotate_contig(generator,scorer,backbone,tokenizer,species,seqid,sequence,ordinal,out,device):
    sequence=re.sub('[^ACGT]','N',sequence.upper())
    if not sequence: raise ValueError('Empty contig')
    locations=[[a,min(a+M.WIDTH,len(sequence))] for a in M.starts(len(sequence))]
    transcripts,scopes=[],[]
    for strand in ('+','-'):
        path=out/('contig_'+str(ordinal).zfill(4)+'_'+('p' if strand=='+' else 'm')+'.jsonl')
        started=time.perf_counter(); counts=Counter()
        with path.open('x') as handle:
            for index,(a,b) in enumerate(locations):
                raw=sequence[a:b]; oriented=raw if strand=='+' else M.rc(raw)
                features=C.I.encode_window(backbone,tokenizer,oriented,device)
                dna=one_hot(oriented,M.WIDTH)[None].to(device)
                shared,_=generator.shared(features,dna,[b-a]); hidden=shared[0]
                ident=species+'|'+seqid+'|'+str(a)+'|'+strand
                row=local_record(generator,hidden,scorer,oriented,index,ident)
                handle.write(json.dumps(row,allow_nan=False)+'\n'); counts.update(row['counts'])
                if (index+1)%32==0:
                    handle.flush(); print(json.dumps(dict(variant='A41',species=species,strand=strand,
                        windows=index+1,total=len(locations))),flush=True)
                del features,dna,shared,hidden,row
        assembled=assemble_predictions(C.I.saved_records(path),locations,strand)
        for index,row in enumerate(assembled['predictions']):
            ident='M28_A41_'+str(ordinal)+'_'+('p' if strand=='+' else 'm')+'_'+str(index)
            transcripts.append(as_transcript(row['chain'],seqid,strand,ident,row['score']))
        scopes.append(dict(species=species,seqid=seqid,strand=strand,ordinal=ordinal,windows=locations,
            file=path.name,candidate_counts=dict(counts),assembly_counts=assembled['counts'],
            seconds=time.perf_counter()-started))
    return transcripts,scopes


@torch.inference_mode()
def main(out,fasta=None):
    fitted=readiness()
    # This file is a scientific contract, not evidence of execution permission.
    if not (ROOT/DEV_PROTOCOL).is_file():
        raise ValueError('Prepare the separate fixed A41 raw-FASTA protocol before native execution')
    if out.exists(): raise FileExistsError('A41 raw inference requires a new output directory')
    generator,scorer=load_models(fitted); backbone,tokenizer=C.I.load_native()
    models=dict(generator=generator,scorer=scorer)
    before={a:{k:v.detach().clone() for k,v in m.state_dict().items()} for a,m in models.items()}
    out.mkdir(parents=True,exist_ok=False); began=time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    transcripts,scopes,lengths=[],[],{}
    for ordinal,(species,seqid,sequence) in enumerate(C.I.input_contigs(fasta)):
        if seqid in lengths: raise ValueError('Raw FASTA contig identifiers must be unique for GFF coordinates')
        txs,local=annotate_contig(generator,scorer,backbone,tokenizer,species,seqid,sequence,ordinal,out,'cuda')
        transcripts.extend(txs); scopes.extend(local); lengths[seqid]=len(sequence)
    windows=sum(len(s['windows']) for s in scopes)
    if not lengths: raise ValueError('Raw FASTA must contain at least one nonempty contig')
    if fasta is None and windows!=8690: raise ValueError('Incomplete natural DEV grid')
    for arm,model in models.items():
        if (any(not torch.equal(v,before[arm][k]) for k,v in model.state_dict().items()) or
                any(p.requires_grad or p.grad is not None for p in model.parameters())):
            raise ValueError('Frozen raw caller altered parameters or received gradients')
    write_gff3(out/'A41.gff3',transcripts,lengths,'M28_A41')
    C.R.write_summary(out,dict(arm='A41',completed=True,inference_complete=True,
        configuration=D.CONFIG,protocol=DEV_PROTOCOL,TRAIN_protocol=D.PROTOCOL,fit_protocol=D.F.PROTOCOL,
        TRAIN_preparation_gate_passed=True,source_TRAIN_jobs=dict(fit=FIT_JOB,readout=READOUT_JOB,CPU=CHECK_JOB),
        generator_checkpoint=C.GENERATOR,scorer_checkpoint=fitted['final_checkpoint'],
        frozen_source_checkpoint=D.CONFIG['frozen_source'],frozen_source_matches_actual_fit=True,source_updated=False,
        generator_optimizer_steps=4608,scorer_optimizer_steps=4608,
        model_parameters={a:sum(p.numel() for p in m.parameters()) for a,m in models.items()},
        feature_revision=C.I.F.REVISION,feature_spec=C.R.FEATURE_SPEC,
        input_mode='fixed_DEV_raw_FASTA' if fasta is None else 'explicit_raw_FASTA',
        input_fasta=str(fasta) if fasta is not None else
            {s:str(ROOT/'data/m1_screen'/s/'genome.fa') for s in C.I.DEV_LENGTHS},
        score_kind=SCORE_KIND,reference_used=False,labels_or_feature_cache_used=False,
        GEN_encoded_once_per_oriented_window=True,A3_shared_encoded_once_per_oriented_window=True,
        natural_proposal_generator='unchanged_A3_greedy_plus_beam4',
        model_weights_unchanged=True,all_parameters_frozen_and_gradients_absent=True,
        optimizer_steps_this_run=0,DEV_used_for_selection=False,test_setaria_used=False,A33_outputs_opened=False,
        windows=windows,scopes=scopes,lengths=lengths,primary_CDS_chain_predictions=len(transcripts),
        wall_seconds=time.perf_counter()-began,peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='Primary CDS caller with frozen supervised A3 and affine residual; not enlarged affine function '
            'class, complete alternative isoforms, matched exposure, accuracy, novelty or independent generalization.'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--dev',action='store_true'); group.add_argument('--fasta',type=Path)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args(); torch.set_num_threads(2)
    main(args.output_dir.resolve(),args.fasta.resolve() if args.fasta else None)
