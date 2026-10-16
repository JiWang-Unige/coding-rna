#!/usr/bin/env python3
"""A13 fixed last-head raw-FASTA caller, conditional on the frozen TRAIN gate."""
import argparse
import json
import math
import re
import time
from collections import Counter
from pathlib import Path
import torch
import run_chain_competition as R
from src.m28.training import one_hot
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript, write_gff3

C = R.C
ROOT = R.ROOT
M = C.I.M
OUT = ROOT/'outputs/M28-A13-CHAIN-COMPETITION-INFER'
DEV_PROTOCOL = 'reports/M28-METHOD-RESTART/a13_chain_competition_dev_protocol.md'
SCORE_KIND = 'A13_complete_chain_logit_minus_window_null'


def validate_fit(fitted):
    R.validate_fit(fitted)
    if (fitted['feature_spec'] != C.R.FEATURE_SPEC or
        fitted['quality_trainable_parameters'] != 165635 or fitted['generator_trainable_parameters'] != 0 or
        [e['steps'] for e in fitted['epochs']] != [1536]*3 or
        fitted['source_free_pool'] != 'outputs/M28-A6-CHAIN-QUALITY-FIT/free_pool.jsonl' or
        fitted['source_TRAIN_cohort'] != 'outputs/M28-A12-QUALITY-TRAIN-DIAGNOSTIC'):
        raise ValueError('Matched A6 architecture, original exposure and four-block features required')


def validate_train_preparation(fitted, prepared):
    validate_fit(fitted)
    checks = {}
    for species in C.R.ALLOW:
        a, b = (prepared['per_species'][species][arm] for arm in ('A6','A13'))
        if any(a[k]['candidates'] != b[k]['candidates'] for k in a):
            raise ValueError('TRAIN preparation denominators changed')
        checks[species+'_free_correct_positive_gain_increased'] = b['free_positive']['positive_gain'] > a['free_positive']['positive_gain']
        checks[species+'_all_correct_positive_gain_increased'] = b['positive']['positive_gain'] > a['positive']['positive_gain']
        checks[species+'_known_nonmatch_positive_gain_not_increased'] = b['known_nonmatch']['positive_gain'] <= a['known_nonmatch']['positive_gain']
    a, b = (prepared['pooled'][arm] for arm in ('A6','A13'))
    if any(a[k]['candidates'] != b[k]['candidates'] for k in a):
        raise ValueError('Pooled TRAIN preparation denominators changed')
    checks['pooled_five_or_more_correct_positive_gain_increased'] = b['five_or_more_positive']['positive_gain'] > a['five_or_more_positive']['positive_gain']
    checks['pooled_correct_not_above_known_conflict_not_increased'] = b['positive']['positive_candidates_not_above_max_known_overlap'] <= a['positive']['positive_candidates_not_above_max_known_overlap']
    if (prepared['experiment'] != 'M28-A13-CHAIN-COMPETITION-TRAIN-READOUT' or
        prepared['completed'] is not True or prepared['configuration'] != R.CONFIG or
        prepared['protocol'] != R.PROTOCOL or prepared['checkpoint'] != fitted['final_checkpoint'] or
        prepared['source_generator'] != C.GENERATOR or
        prepared['source_A6'] != 'outputs/M28-A12-QUALITY-TRAIN-DIAGNOSTIC' or
        prepared['windows'] != 1491 or prepared['original_draws'] != 1536 or
        prepared['candidate_window_pairs'] != 16857 or prepared['optimizer_steps_this_run'] != 0 or
        prepared['candidate_geometry_labels_and_origins_unchanged'] is not True or
        any(prepared[k] is not False for k in ('new_candidate_search','new_backbone_forward',
            'DEV_forward_used','test_setaria_used','threshold_scan','head_weights_updated','checkpoints_written')) or
        prepared['threshold'] != 0. or prepared['predeclared_TRAIN_preparation_gate'] != checks or
        prepared['TRAIN_preparation_supported'] is not True or not all(checks.values())):
        raise ValueError('Complete fixed TRAIN preparation gate must pass before A13 DEV execution')


def prediction_id(ordinal, strand, index):
    return 'M28_A13_'+str(ordinal)+'_'+strand+'_'+str(index)


@torch.inference_mode()
def quality_record(generator, head, hidden, sequence, index, ident):
    chains, sources = C.free_pool(generator, hidden, sequence)
    scores = head(hidden, chains, len(sequence))['logits'].float().cpu().tolist()
    if any(not math.isfinite(v) for v in scores):
        raise ValueError('Nonfinite A13 chain gain')
    return dict(window_index=index, window_id=ident, chains=chains, score=scores,
        decoder_sources=sources, score_kind=SCORE_KIND, quality_logit_threshold=0., reference_used=False,
        counts=dict(free_candidates=len(chains), nonpositive_candidates=sum(v <= 0 for v in scores)))


@torch.inference_mode()
def annotate_contig(generator, head, backbone, tokenizer, species, seqid, sequence, ordinal, out, device):
    sequence = re.sub('[^ACGT]','N',sequence.upper())
    if not sequence:
        raise ValueError('Empty contig')
    length = len(sequence)
    locations = [[a,min(a+M.WIDTH,length)] for a in M.starts(length)]
    transcripts, scopes = [],[]
    for strand in ('+','-'):
        path = out/('contig_'+str(ordinal).zfill(4)+'_'+('p' if strand == '+' else 'm')+'.jsonl')
        started = time.perf_counter(); counts = Counter()
        with path.open('x') as handle:
            for index,(a,b) in enumerate(locations):
                raw = sequence[a:b]; oriented = raw if strand == '+' else M.rc(raw)
                features = C.I.encode_window(backbone,tokenizer,oriented,device)
                dna = one_hot(oriented,M.WIDTH)[None].to(device)
                hidden = generator(features,dna,[b-a])['features'][0]
                ident = species+'|'+seqid+'|'+str(a)+'|'+strand
                free = quality_record(generator,head,hidden,oriented,index,ident)
                handle.write(json.dumps(free,allow_nan=False)+'\n'); counts.update(free['counts'])
                if (index+1)%32 == 0:
                    handle.flush(); print(json.dumps(dict(species=species,strand=strand,
                        windows=index+1,total=len(locations))),flush=True)
                del features,dna,hidden,free
        assembled = assemble_predictions(C.I.saved_records(path),locations,strand)
        for index,row in enumerate(assembled['predictions']):
            transcripts.append(as_transcript(row['chain'],seqid,strand,
                prediction_id(ordinal,strand,index),row['score']))
        scopes.append(dict(species=species,seqid=seqid,strand=strand,ordinal=ordinal,
            windows=locations,file=path.name,quality_counts=dict(counts),
            assembly_counts=assembled['counts'],seconds=time.perf_counter()-started))
    return transcripts,scopes


def main(out, fasta=None):
    if any((p/'STATUS').read_text().strip() != 'COMPLETED' for p in (R.FIT,R.READOUT)):
        raise ValueError('Terminal A13 fit and TRAIN readout required')
    fitted = json.loads((R.FIT/'summary.json').read_text())
    prepared = json.loads((R.READOUT/'summary.json').read_text()); validate_train_preparation(fitted,prepared)
    state = torch.load(ROOT/fitted['final_checkpoint'],map_location='cpu',weights_only=True)
    R.validate_fit(state,state=True)
    head = C.new_quality('cuda').eval().requires_grad_(False)
    head.load_state_dict(state['head'],strict=True); del state
    generator = C.load_generator()
    backbone,tokenizer = C.I.load_native()
    out.mkdir(parents=True,exist_ok=True)
    if (out/'summary.json').exists():
        raise FileExistsError('Completed A13 inference exists')
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    transcripts,scopes,lengths = [],[],{}
    for ordinal,(species,seqid,sequence) in enumerate(C.I.input_contigs(fasta)):
        txs,ss = annotate_contig(generator,head,backbone,tokenizer,species,seqid,sequence,ordinal,out,'cuda')
        transcripts.extend(txs); scopes.extend(ss); lengths[seqid] = len(sequence)
    windows = sum(len(s['windows']) for s in scopes)
    if fasta is None and windows != 8690:
        raise ValueError('Incomplete natural DEV grid')
    write_gff3(out/'A13.gff3',transcripts,lengths,'M28_A13')
    C.R.write_summary(out,dict(arm='A13',completed=True,inference_complete=True,
        configuration=R.CONFIG,protocol=DEV_PROTOCOL,TRAIN_protocol=R.PROTOCOL,
        TRAIN_preparation_gate_passed=True,checkpoint=fitted['final_checkpoint'],optimizer_steps=4608,
        generator_checkpoint=C.GENERATOR,feature_revision=C.I.F.REVISION,feature_spec=C.R.FEATURE_SPEC,
        input_mode='fixed_DEV_raw_FASTA' if fasta is None else 'explicit_raw_FASTA',
        input_fasta=str(fasta) if fasta is not None else
            {s:str(ROOT/'data/m1_screen'/s/'genome.fa') for s in C.I.DEV_LENGTHS},
        reference_used=False,labels_or_feature_cache_used=False,optimizer_steps_this_run=0,
        DEV_used_for_selection=False,test_setaria_used=False,windows=windows,scopes=scopes,lengths=lengths,
        primary_CDS_chain_predictions=len(transcripts),wall_seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved(),
        limitations='Complete raw FASTA output, not accuracy or independent generalization before separate evaluation.'))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument('--dev',action='store_true'); group.add_argument('--fasta',type=Path)
    p.add_argument('--output-dir',type=Path,required=True)
    a = p.parse_args(); torch.set_num_threads(2)
    main(a.output_dir.resolve(),a.fasta.resolve() if a.fasta else None)
