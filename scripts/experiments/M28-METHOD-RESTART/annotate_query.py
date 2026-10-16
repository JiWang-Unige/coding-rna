#!/usr/bin/env python3
"""A raw-FASTA caller for the pending fixed A checkpoint.

Production use requires its own approved input/resource scope; it is not part
of the cached-feature first-fit GPU budget. No annotation or label input.
Frozen official GLM encoding follows the existing four-block feature cache.
"""
import argparse
import json
import math
import re
import time
from pathlib import Path
import torch
import run_query as Q
from src.m28.query_model import GeneQueryModel
from src.m28.query_output import query_window_candidates
from src.m28.training import one_hot
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript, write_gff3


def fasta_records(path):
    seen = set(); ident = None; parts = []
    with path.open() as handle:
        for line in handle:
            if not line.strip(): continue
            if line.startswith('>'):
                if ident is not None:
                    sequence = re.sub('[^ACGT]','N',''.join(parts).upper())
                    if not sequence: raise ValueError('Empty FASTA sequence')
                    yield ident,sequence
                words = line[1:].split()
                if not words or words[0] in seen:
                    raise ValueError('Expected nonempty unique FASTA sequence IDs')
                ident = words[0]; seen.add(ident); parts = []
            else:
                if ident is None: raise ValueError('Sequence before FASTA header')
                parts.append(line.strip())
    if ident is None: raise ValueError('No FASTA records')
    sequence = re.sub('[^ACGT]','N',''.join(parts).upper())
    if not sequence: raise ValueError('Empty FASTA sequence')
    yield ident,sequence


@torch.inference_mode()
def encode_window(backbone, tokenizer, sequence, device):
    if not 0 < len(sequence) <= Q.M.WIDTH:
        raise ValueError('Expected a nonempty valid oriented window')
    padded = sequence+'N'*(Q.M.WIDTH-len(sequence))
    whole = tokenizer(padded,add_special_tokens=False,return_tensors='pt',truncation=False)['input_ids']
    blocks = [tokenizer(padded[a:a+Q.F.BLOCK],add_special_tokens=False,
                        return_tensors='pt',truncation=False)['input_ids']
              for a in range(0,len(padded),Q.F.BLOCK)]
    if any(tuple(ids.shape) != (1,1024) for ids in blocks) or not torch.equal(torch.cat(blocks,dim=1),whole):
        raise ValueError('Whole/four-block six-mer token alignment mismatch')
    features = []
    for ids in blocks:
        ids = ids.to(device)
        hidden = backbone(input_ids=ids,attention_mask=torch.ones_like(ids)).last_hidden_state
        if tuple(hidden.shape) != (1,1024,2048) or not torch.isfinite(hidden).all():
            raise ValueError('Invalid frozen GLM hidden features')
        # Reproduce the existing cache's BF16 materialization before FP32 heads.
        features.append(hidden[0].to(dtype=torch.bfloat16))
    return torch.cat(features).contiguous().float()[None]


def output_bytes(directory, limit):
    size = sum(p.stat().st_size for p in directory.rglob('*') if p.is_file())
    if size > limit: raise RuntimeError('Raw-FASTA output budget exceeded')
    return size


def saved_records(path):
    with path.open() as handle:
        for line in handle: yield json.loads(line)


@torch.inference_mode()
def annotate_contig(model, backbone, tokenizer, seqid, sequence, ordinal, out, device, byte_limit):
    length = len(sequence); locations = [(a,min(a+Q.M.WIDTH,length)) for a in Q.M.starts(length)]
    transcripts,scopes = [],[]
    for strand in ('+','-'):
        letter = 'p' if strand == '+' else 'm'
        path = out/('contig_'+str(ordinal).zfill(4)+'_'+letter+'.jsonl')
        counts = {}
        with path.open('x') as handle:
            for i,(a,b) in enumerate(locations):
                oriented = sequence[a:b] if strand == '+' else Q.M.rc(sequence[a:b])
                features = encode_window(backbone,tokenizer,oriented,device)
                dna = one_hot(oriented,Q.M.WIDTH)[None].to(device)
                outputs = model(features,dna,[b-a])
                free = query_window_candidates(model.query,outputs['features'][0],oriented,i,
                                               seqid+'|'+str(a)+'|'+strand,Q.CONFIG['max_exons'],
                                               Q.CONFIG['existence_logit_threshold'],b-a)
                handle.write(json.dumps(free,allow_nan=False)+'\n'); handle.flush()
                for key,value in free['counts'].items(): counts[key] = counts.get(key,0)+value
                output_bytes(out,byte_limit)
                del features,dna,outputs,free
        assembled = assemble_predictions(saved_records(path),locations,strand)
        for i,row in enumerate(assembled['predictions']):
            identifier = 'M28_A_'+str(ordinal)+'_'+letter+'_'+str(i)
            transcripts.append(as_transcript(row['chain'],seqid,strand,identifier,row['score']))
        scopes.append({'seqid':seqid,'strand':strand,'windows':locations,'file':path.name,
                       'query_counts':counts,'assembly_counts':assembled['counts']})
    return transcripts,scopes


def main(fasta, checkpoint, out, output_gib):
    if not math.isfinite(output_gib) or output_gib <= 0:
        raise ValueError('A positive finite explicit output budget is required')
    byte_limit = int(output_gib*1024**3)
    out.mkdir(parents=True,exist_ok=False)
    state = torch.load(checkpoint,map_location='cpu',weights_only=True)
    Q.check_last_identity(state,state=True)
    model = GeneQueryModel(Q.CONFIG['feature_dim'],Q.CONFIG['num_queries'],Q.CONFIG['hidden'])
    model.load_state_dict(state['model'],strict=True); del state
    model.eval().requires_grad_(False).to('cuda')
    # Native cached implementation only; no downloads or replacement revision.
    from transformers import AutoTokenizer, AutoModelForTokenClassification
    tokenizer = AutoTokenizer.from_pretrained(Q.F.MODEL,revision=Q.F.REVISION,
                                              trust_remote_code=True,local_files_only=True)
    full = AutoModelForTokenClassification.from_pretrained(
        Q.F.MODEL,revision=Q.F.REVISION,trust_remote_code=True,local_files_only=True,
        attn_implementation='sdpa',torch_dtype=torch.bfloat16)
    if int(getattr(full,'k',getattr(full.config,'k',0))) != 6:
        raise ValueError('Released encoder six-mer identity mismatch')
    backbone = full.model; del full
    backbone.eval().requires_grad_(False).to('cuda')
    torch.cuda.reset_peak_memory_stats(); began = time.perf_counter()
    transcripts,scopes,lengths = [],[],{}
    for ordinal,(seqid,sequence) in enumerate(fasta_records(fasta)):
        txs,record_scopes = annotate_contig(model,backbone,tokenizer,seqid,sequence,ordinal,out,'cuda',byte_limit)
        transcripts.extend(txs); scopes.extend(record_scopes); lengths[seqid] = len(sequence)
        print(json.dumps({'completed_contig':seqid,'length':len(sequence),'predictions':len(txs)}),flush=True)
    write_gff3(out/'A.gff3',transcripts,lengths,'M28_A')
    report = {'inference_complete':True,'input_fasta':str(fasta),'checkpoint':str(checkpoint),
              'configuration':Q.CONFIG,'encoder_model':Q.F.MODEL,'encoder_revision':Q.F.REVISION,
              'lengths':lengths,'scopes':scopes,'primary_CDS_chain_predictions':len(transcripts),
              'reference_used':False,'labels_or_feature_cache_used':False,'optimizer_steps_this_run':0,
              'score_kind':'query_existence_logit_minus_threshold','output_gib_limit':output_gib,
              'wall_loop_seconds':time.perf_counter()-began,
              'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
              'peak_GPU_reserved_bytes':torch.cuda.max_memory_reserved(),
              'limitations':'greedy canonical CDS only; no UTR/all-isoform/global-optimum or accuracy claim'}
    report['artifact_bytes_before_summary'] = output_bytes(out,byte_limit)
    with (out/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    output_bytes(out,byte_limit)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta',type=Path,required=True)
    parser.add_argument('--checkpoint',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--output-gib-limit',type=float,required=True)
    args = parser.parse_args(); torch.set_num_threads(2)
    main(args.fasta.resolve(),args.checkpoint.resolve(),args.output_dir.resolve(),args.output_gib_limit)
