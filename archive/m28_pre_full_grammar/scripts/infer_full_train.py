#!/usr/bin/env python3
"""Fixed last full-TRAIN C0/Q checkpoint caller: raw FASTA, no annotation/cache input."""
import argparse
import json
import math
import re
import time
from collections import Counter
from pathlib import Path
import torch
import run_full_train as R
from annotate_query import encode_window, fasta_records, saved_records
from cache_full_train import load_native
from infer_pilot import DEV_LENGTHS
from src.m28.training import one_hot
from src.m28.query_output import query_window_candidates
from src.m28.c0_decode import decode_c0
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript, write_gff3

ROOT = R.ROOT
M = R.Q.M
F = R.Q.F
FIT = {arm:ROOT/'outputs'/('M28-FULL-TRAIN-'+arm+'-FIT') for arm in R.CONFIGS}
INFER = {arm:ROOT/'outputs'/('M28-FULL-TRAIN-'+arm+'-INFER') for arm in R.CONFIGS}


def validate_fit(s, arm):
    expected = str((FIT[arm]/'step_076536.pt').relative_to(ROOT))
    if (s['completed'] is not True or s['arm'] != arm or
        s['configuration'] != R.CONFIGS[arm] or s['optimizer_steps'] != R.STEPS or
        s['feature_revision'] != F.REVISION or s['feature_spec'] != R.FEATURE_SPEC or
        s['dataset']['windows'] != R.N or s['final_checkpoint'] != expected or
        s['main_evaluation'] != 'last_step_076536_only' or s['protocol'] != R.PROTOCOL or
        s['DEV_used_for_selection'] is not False or s['test_setaria_used'] is not False or
        [e['steps'] for e in s['epochs']] != [R.N]*R.PASSES):
        raise ValueError('Complete matching full-TRAIN final fit required; no DEV checkpoint selection')


def validate_checkpoint(s, arm):
    if (s['arm'] != arm or s['configuration'] != R.CONFIGS[arm] or
        s['step'] != R.STEPS or s['revision'] != F.REVISION or s['protocol'] != R.PROTOCOL):
        raise ValueError('Only the matching full-TRAIN step 076536 checkpoint is allowed')


def prediction_id(arm, ordinal, strand, index):
    return 'M28_FULL_'+arm+'_'+str(ordinal)+'_'+('p' if strand=='+' else 'm')+'_'+str(index)


@torch.inference_mode()
def window_record(arm, model, outputs, sequence, index, ident):
    if arm == 'Q':
        c = R.CONFIGS['Q']
        return query_window_candidates(model.query,outputs['features'][0],sequence,index,ident,
                                       c['max_exons'],c['existence_logit_threshold'],len(sequence))
    if arm != 'C0': raise ValueError('Fixed arm required')
    p = outputs['segmentation_logits'][0,:len(sequence)].softmax(-1).cpu().numpy()
    decoded = decode_c0(sequence,p)
    if any(not math.isfinite(x) for x in decoded['gains']): raise ValueError('Nonfinite C0 chain gain')
    return {'window_index':index,'window_id':ident,'chains':decoded['chains'],'score':decoded['gains'],
            'counts':{'partial_paths':decoded['partial_paths']},'reference_used':False,
            'score_kind':'native_HMM_chain_removal_gain'}


@torch.inference_mode()
def annotate_contig(arm, model, backbone, tokenizer, species, seqid, sequence, ordinal, out, device):
    if not sequence: raise ValueError('Empty contig')
    length = len(sequence)
    locations = [(a,min(a+M.WIDTH,length)) for a in M.starts(length)]
    transcripts,scopes = [],[]
    for strand in ('+','-'):
        path = out/('contig_'+str(ordinal).zfill(4)+'_'+('p' if strand=='+' else 'm')+'.jsonl')
        began = time.perf_counter(); counts = Counter()
        with path.open('x') as handle:
            for i,(a,b) in enumerate(locations):
                raw = sequence[a:b]
                oriented = raw if strand=='+' else M.rc(raw)
                # Same normalized encoder tokens and original oriented DNA as TRAIN cache.
                features = encode_window(backbone,tokenizer,re.sub('[^ACGT]','N',oriented),device)
                dna = one_hot(oriented,M.WIDTH)[None].to(device)
                outputs = model(features,dna,[b-a])
                ident = species+'|'+seqid+'|'+str(a)+'|'+strand
                free = window_record(arm,model,outputs,oriented,i,ident)
                handle.write(json.dumps(free,allow_nan=False)+'\n')
                counts.update(free['counts'])
                if (i+1)%64 == 0:
                    handle.flush()
                    print(json.dumps({'arm':arm,'seqid':seqid,'strand':strand,
                                      'windows_completed':i+1,'windows':len(locations)}),flush=True)
                del features,dna,outputs,free
        assembled = assemble_predictions(saved_records(path),locations,strand)
        for i,row in enumerate(assembled['predictions']):
            transcripts.append(as_transcript(row['chain'],seqid,strand,
                                prediction_id(arm,ordinal,strand,i),row['score']))
        scopes.append({'species':species,'seqid':seqid,'strand':strand,'ordinal':ordinal,
                       'windows':locations,'file':path.name,
                       'query_counts' if arm=='Q' else 'decoder_counts':dict(counts),
                       'assembly_counts':assembled['counts'],'seconds':time.perf_counter()-began})
    return transcripts,scopes


def input_contigs(fasta):
    if fasta is not None:
        for seqid,sequence in fasta_records(fasta): yield 'user_input',seqid,sequence
    else:
        # Materialize only the two permitted DEV chromosomes, not other FASTA records.
        for species,lengths in DEV_LENGTHS.items():
            sequences = M.read_allowed(ROOT/'data/m1_screen'/species/'genome.fa',lengths)
            for seqid in lengths: yield species,seqid,sequences[seqid]


def main(arm, out, fasta=None):
    fitted = json.loads((FIT[arm]/'summary.json').read_text()); validate_fit(fitted,arm)
    if (FIT[arm]/'STATUS').read_text().strip() != 'COMPLETED': raise ValueError('Terminal fit required')
    checkpoint = ROOT/fitted['final_checkpoint']
    state = torch.load(checkpoint,map_location='cpu',weights_only=True); validate_checkpoint(state,arm)
    model = R.new_model(arm,'cpu'); model.load_state_dict(state['model'],strict=True); del state
    model.eval().requires_grad_(False).to('cuda')
    backbone,tokenizer = load_native()
    out.mkdir(parents=True,exist_ok=True)
    if (out/'summary.json').exists(): raise FileExistsError('Completed inference exists')
    torch.cuda.reset_peak_memory_stats(); began = time.perf_counter()
    transcripts,scopes,lengths = [],[],{}
    for ordinal,(species,seqid,sequence) in enumerate(input_contigs(fasta)):
        txs,ss = annotate_contig(arm,model,backbone,tokenizer,species,seqid,sequence,ordinal,out,'cuda')
        transcripts.extend(txs); scopes.extend(ss); lengths[seqid] = len(sequence)
    write_gff3(out/(arm+'.gff3'),transcripts,lengths,'M28_FULL_'+arm)
    report = {'inference_complete':True,'arm':arm,'configuration':R.CONFIGS[arm],
              'checkpoint':fitted['final_checkpoint'],'optimizer_steps':R.STEPS,
              'encoder_model':F.MODEL,'encoder_revision':F.REVISION,'feature_spec':R.FEATURE_SPEC,
              'input_mode':'fixed_DEV_raw_FASTA' if fasta is None else 'explicit_raw_FASTA',
              'input_fasta':str(fasta) if fasta is not None else
                            {s:str(ROOT/'data/m1_screen'/s/'genome.fa') for s in DEV_LENGTHS},
              'lengths':lengths,'scopes':scopes,'windows':sum(len(s['windows']) for s in scopes),
              'primary_CDS_chain_predictions':len(transcripts),'reference_used':False,
              'labels_or_feature_cache_used':False,'optimizer_steps_this_run':0,
              'test_setaria_used':False,'DEV_used_for_selection':False,
              'wall_loop_seconds':time.perf_counter()-began,
              'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
              'peak_GPU_reserved_bytes':torch.cuda.max_memory_reserved(),
              'assembly':'unchanged nearest owner, global deduplication, positive-gain nonoverlap',
              'dependency':'unchanged ANNEVO native HMM' if arm=='C0' else 'no ANNEVO HMM',
              'limitations':'primary canonical CDS only; no UTR/all-isoform/accuracy/generalization claim'}
    with (out/'summary.json').open('x') as f:
        json.dump(report,f,indent=2,allow_nan=False); f.write('\n')
    print(json.dumps({k:v for k,v in report.items() if k!='scopes'},allow_nan=False),flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--arm',choices=tuple(R.CONFIGS),required=True)
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument('--dev',action='store_true')
    group.add_argument('--fasta',type=Path)
    p.add_argument('--output-dir',type=Path,required=True)
    a = p.parse_args(); torch.set_num_threads(2)
    main(a.arm,a.output_dir.resolve(),a.fasta.resolve() if a.fasta else None)
