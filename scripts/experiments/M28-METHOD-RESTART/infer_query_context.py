#!/usr/bin/env python3
"""A10 last matched checkpoint: raw FASTA and one frozen 4096-token forward/window."""
import argparse
import json
import re
import time
from collections import Counter
from pathlib import Path
import torch
import run_query_context as R
import cache_draw_context as C
import infer_full_train as I
from src.m28.training import one_hot
from src.m28.query_output import query_window_candidates
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript, write_gff3

ROOT, M, F = R.ROOT, I.M, I.F
DEV_LENGTHS = I.DEV_LENGTHS
ARM = 'A10_CONTEXT'
OUT = ROOT/'outputs/M28-A10-CONTEXT-INFER'


def validate_fit(fitted):
    R.check_last_identity(fitted)
    if (fitted['final_checkpoint'] != str((R.FIT/'step_004608.pt').relative_to(ROOT)) or
        fitted['main_evaluation'] != 'last_step_004608_only' or
        fitted['trainable_parameters'] != 649457 or
        [e['steps'] for e in fitted['epochs']] != [1536]*3 or
        fitted['DEV_used_for_selection'] is not False or fitted['test_setaria_used'] is not False):
        raise ValueError('Complete fixed matched A10 fit required, without checkpoint selection')


def prediction_id(ordinal, strand, index):
    return 'M28_A10_'+str(ordinal)+'_'+('p' if strand == '+' else 'm')+'_'+str(index)


@torch.inference_mode()
def encode_window(backbone, tokenizer, sequence, device):
    # Reuse the exact TRAIN BF16 materialization before the original FP32 head.
    return C.encode_full_window(backbone,tokenizer,sequence).float()[None].to(device)


@torch.inference_mode()
def annotate_contig(model, backbone, tokenizer, species, seqid, sequence, ordinal, out, device):
    sequence = re.sub('[^ACGT]','N',sequence.upper())
    if not sequence:
        raise ValueError('Empty raw contig')
    length = len(sequence)
    locations = [(a,min(a+M.WIDTH,length)) for a in M.starts(length)]
    transcripts,scopes = [],[]
    for strand in ('+','-'):
        path = out/('contig_'+str(ordinal).zfill(4)+'_'+('p' if strand == '+' else 'm')+'.jsonl')
        began = time.perf_counter(); counts = Counter()
        with path.open('x') as handle:
            for index,(a,b) in enumerate(locations):
                raw = sequence[a:b]; oriented = raw if strand == '+' else M.rc(raw)
                features = encode_window(backbone,tokenizer,oriented,device)
                dna = one_hot(oriented,M.WIDTH)[None].to(device)
                outputs = model(features,dna,[b-a])
                ident = species+'|'+seqid+'|'+str(a)+'|'+strand
                free = query_window_candidates(model.query,outputs['features'][0],oriented,index,ident,
                    R.CONFIG['max_exons'],R.CONFIG['existence_logit_threshold'],b-a)
                handle.write(json.dumps(free,allow_nan=False)+'\n'); counts.update(free['counts'])
                if (index+1)%64 == 0:
                    handle.flush()
                    print(json.dumps(dict(arm=ARM,seqid=seqid,strand=strand,
                        windows_completed=index+1,windows=len(locations))),flush=True)
                del features,dna,outputs,free
        assembled = assemble_predictions(I.saved_records(path),locations,strand)
        for index,row in enumerate(assembled['predictions']):
            transcripts.append(as_transcript(row['chain'],seqid,strand,
                prediction_id(ordinal,strand,index),row['score']))
        scopes.append(dict(species=species,seqid=seqid,strand=strand,ordinal=ordinal,
            windows=locations,file=path.name,query_counts=dict(counts),assembly_counts=assembled['counts'],
            seconds=time.perf_counter()-began,reference_used=False,labels_or_feature_cache_used=False))
    return transcripts,scopes


def main(out, fasta=None):
    if (R.FIT/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Terminal matched fit required before raw inference')
    fitted = json.loads((R.FIT/'summary.json').read_text()); validate_fit(fitted)
    state = torch.load(ROOT/fitted['final_checkpoint'],map_location='cpu',weights_only=True)
    R.check_last_identity(state,state=True)
    model = R.A.new_model(); model.load_state_dict(state['model'],strict=True); del state
    model.eval().requires_grad_(False)
    backbone,tokenizer = I.load_native()
    if any(p.requires_grad for p in backbone.parameters()) or int(backbone.config.max_position_embeddings) < 4096:
        raise ValueError('Same frozen context-capable native encoder required')
    out.mkdir(parents=True,exist_ok=True)
    if (out/'summary.json').exists():
        raise FileExistsError('Completed inference exists')
    torch.cuda.reset_peak_memory_stats(); began = time.perf_counter()
    transcripts,scopes,lengths = [],[],{}
    for ordinal,(species,seqid,sequence) in enumerate(I.input_contigs(fasta)):
        if seqid in lengths:
            raise ValueError('Raw FASTA sequence IDs must be unique')
        txs,ss = annotate_contig(model,backbone,tokenizer,species,seqid,sequence,ordinal,out,'cuda')
        transcripts.extend(txs); scopes.extend(ss); lengths[seqid] = len(sequence)
    if not scopes:
        raise ValueError('Raw FASTA contains no contigs')
    write_gff3(out/(ARM+'.gff3'),transcripts,lengths,'M28_A10')
    report = dict(inference_complete=True,arm=ARM,configuration=R.CONFIG,
        checkpoint=fitted['final_checkpoint'],optimizer_steps=4608,
        encoder_model=F.MODEL,encoder_revision=F.REVISION,feature_spec=C.FULL_SPEC,
        input_mode='fixed_DEV_raw_FASTA' if fasta is None else 'explicit_raw_FASTA',
        input_fasta=str(fasta) if fasta is not None else
            {s:str(ROOT/'data/m1_screen'/s/'genome.fa') for s in DEV_LENGTHS},
        lengths=lengths,scopes=scopes,windows=sum(len(s['windows']) for s in scopes),
        primary_CDS_chain_predictions=len(transcripts),reference_used=False,
        labels_or_feature_cache_used=False,optimizer_steps_this_run=0,
        test_setaria_used=False,DEV_used_for_selection=False,
        wall_loop_seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved(),
        assembly='unchanged nearest owner, global deduplication, positive-gain nonoverlap',
        dependency='no ANNEVO HMM',
        limitations='primary canonical CDS only; no UTR/all-isoform/accuracy/generalization claim')
    with (out/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps({k:v for k,v in report.items() if k != 'scopes'},allow_nan=False),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--dev',action='store_true'); group.add_argument('--fasta',type=Path)
    parser.add_argument('--output-dir',type=Path,required=True)
    args = parser.parse_args(); torch.set_num_threads(2)
    main(args.output_dir.resolve(),args.fasta.resolve() if args.fasta else None)
