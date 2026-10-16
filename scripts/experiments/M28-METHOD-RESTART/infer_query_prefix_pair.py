#!/usr/bin/env python3
"""Prepared A9 raw-contig paired caller; no experiment/DEV execution entry yet.

Imports reusable native encoding/FASTA geometry, never reads an annotation or
feature artifact. Caller supplies already validated frozen model/head. An actual
DEV launch still requires the fixed TRAIN probe gate and a new full protocol.
"""
import json,re,time
from collections import Counter
import torch
import infer_full_train as I
from src.m28.training import one_hot
from src.m28.query_prefix_pair import ARMS,paired_query_records
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript


def arm_records(path,arm):
    for pair in I.saved_records(path):yield pair['arms'][arm]


@torch.no_grad()
def annotate_contig(generator,head,backbone,tokenizer,species,seqid,sequence,ordinal,out,device,limits):
    """One native forward per oriented window, unchanged owner/selector per arm.

    This returns transcripts/scopes, not an accuracy/complete-grid receipt.
    It does not load checkpoints, choose a dataset, or submit any work.
    """
    sequence=re.sub('[^ACGT]','N',sequence.upper())
    if not sequence:raise ValueError('Empty raw contig')
    length=len(sequence);locations=[(a,min(a+I.M.WIDTH,length)) for a in I.M.starts(length)]
    transcripts={a:[] for a in ARMS};scopes=[]
    for strand in ('+','-'):
        letter='p' if strand=='+' else 'm';path=out/('contig_'+str(ordinal).zfill(4)+'_'+letter+'.jsonl')
        began=time.perf_counter();counts={a:Counter() for a in ARMS};common=0
        with path.open('x') as handle:
            for index,(a,b) in enumerate(locations):
                raw=sequence[a:b];oriented=raw if strand=='+' else I.M.rc(raw)
                features=I.encode_window(backbone,tokenizer,oriented,device)
                dna=one_hot(oriented,I.M.WIDTH)[None].to(device)
                hidden=generator(features,dna,[b-a])['features'][0]
                ident=species+'|'+seqid+'|'+str(a)+'|'+strand
                pair=paired_query_records(generator.query,hidden,oriented,head,index,ident,**limits,valid_bases=b-a)
                handle.write(json.dumps(pair,allow_nan=False)+'\n');common+=pair['common_query_chain_records']
                for arm in ARMS:counts[arm].update(pair['arms'][arm]['counts'])
                if (index+1)%32==0:
                    handle.flush();print(json.dumps(dict(phase='paired_raw_infer',species=species,strand=strand,
                        windows=index+1,total=len(locations))),flush=True)
                del features,dna,hidden,pair
        assembly_counts={}
        for arm in ARMS:
            assembled=assemble_predictions(arm_records(path,arm),locations,strand)
            assembly_counts[arm]=assembled['counts']
            for i,row in enumerate(assembled['predictions']):
                identifier='M28_A9_'+arm+'_'+str(ordinal)+'_'+letter+'_'+str(i)
                transcripts[arm].append(as_transcript(row['chain'],seqid,strand,identifier,row['score']))
        scopes.append(dict(species=species,seqid=seqid,strand=strand,ordinal=ordinal,windows=locations,
            file=path.name,decoder_counts={a:dict(c) for a,c in counts.items()},assembly_counts=assembly_counts,
            common_query_chain_records=common,common_query_chain_values_identical=True,
            seconds=time.perf_counter()-began,reference_used=False,labels_or_feature_cache_used=False))
    return transcripts,scopes


if __name__=='__main__':
    raise SystemExit('Preparation helper only: no full DEV protocol/resources or execution entry is fixed yet.')
