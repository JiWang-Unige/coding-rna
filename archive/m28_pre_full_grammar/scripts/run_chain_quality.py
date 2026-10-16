#!/usr/bin/env python3
"""Fixed A6 TRAIN chain-quality fit and reference/cache-free raw-FASTA caller."""
import argparse
from collections import Counter
import json
from pathlib import Path
import re
import time
import torch
import run_full_train as R
import run_query_known_mean as A
import infer_full_train as I
from src.m28.core import ChainHead, chain_loss
from src.m28.labels import build_targets, candidate_targets
from src.m28.training import inject_training, one_hot
from src.m28.query_beam_pool import completed_union
from src.m28.query_beam import decode_query_beam
from src.m28.query_decode import decode_query_chains
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript, write_gff3

ROOT = R.ROOT
FIT = ROOT/'outputs/M28-A6-CHAIN-QUALITY-FIT'
INFER = ROOT/'outputs/M28-A6-CHAIN-QUALITY-INFER'
PROTOCOL = 'reports/M28-METHOD-RESTART/a6_chain_quality_fit_protocol.md'
GENERATOR = 'outputs/M28-A3-FIRST-FIT/step_004608.pt'
CONFIG = {'generator_checkpoint':GENERATOR,'generator_optimizer_steps':4608,
    'generator_frozen':True,'hidden':128,'quality_seed':0,'precision':'FP32',
    'beam_width':4,'completed_limit':4,'max_exons':128,
    'pool':'all_completed_greedy_plus_beam4_geometry_union',
    'quality_head':'existing_ChainHead_without_link_loss','quality_trainable_parameters':165635,
    'loss':'known_candidate_mean_BCE_with_TRAIN_positive_injection',
    'draws_per_pass':1536,'unique_TRAIN_windows':1491,'passes':3,'optimizer_steps':4608,
    'optimizer':'AdamW','lr':3e-4,'weight_decay':.01,'betas':[.9,.999],
    'eps':1e-8,'clip_norm':1.,'quality_logit_threshold':0.,
    'importance':'unchanged_original_uniform_within_stratum_draw_weight',
    'main_evaluation':'last_head_step_004608_only'}


def load_generator():
    fitted = json.loads((A.STAGE_DIRS['fit']/'summary.json').read_text())
    A.check_last_identity(fitted)
    if ((A.STAGE_DIRS['fit']/'STATUS').read_text().strip()!='COMPLETED' or
        fitted['final_checkpoint']!=GENERATOR):
        raise ValueError('Only completed fixed A3 last generator is allowed')
    state = torch.load(ROOT/GENERATOR,map_location='cpu',weights_only=True)
    A.check_last_identity(state,state=True)
    model = A.new_model().eval().requires_grad_(False)
    model.load_state_dict(state['model'],strict=True)
    return model


def new_quality(device):
    torch.manual_seed(CONFIG['quality_seed'])
    head = ChainHead(CONFIG['hidden']).to(device)
    head.link.requires_grad_(False)
    return head


def validate_head(s, state=False):
    if (s['arm']!='A6' or s['configuration']!=CONFIG or
        s['revision' if state else 'feature_revision']!=I.F.REVISION or
        s['step' if state else 'optimizer_steps']!=4608 or s['protocol']!=PROTOCOL):
        raise ValueError('Only fixed A6 last head checkpoint is allowed')
    if not state and (not s['completed'] or s['final_checkpoint']!=
                     str((FIT/'step_004608.pt').relative_to(ROOT)) or
                     s['DEV_used_for_selection'] or s['test_setaria_used']):
        raise ValueError('Complete matching A6 fixed TRAIN fit required')


@torch.no_grad()
def free_pool(generator, hidden, sequence):
    c = CONFIG
    greedy = decode_query_chains(generator.query,hidden,sequence,c['max_exons'],len(sequence))
    beam = decode_query_beam(generator.query,hidden,sequence,c['max_exons'],
                             c['beam_width'],c['completed_limit'],len(sequence))
    chains,sources = completed_union(greedy,beam)
    return chains,sources


def fit():
    preflight = ROOT/'outputs/M28-A6-CHAIN-QUALITY-PREFLIGHT'
    prepared = json.loads((preflight/'summary.json').read_text())
    if ((preflight/'STATUS').read_text().strip()!='COMPLETED' or not prepared['completed'] or
        prepared['generator_checkpoint']!=GENERATOR or prepared['optimizer_steps_this_run']!=0 or
        not prepared['head_weights_unchanged'] or prepared['quality_trainable_parameters']!=165635):
        raise ValueError('Completed A6 zero-update gradient preflight required')
    R.validate_feature_summary(json.loads((R.FEATURES/'summary.json').read_text()))
    draws = R.rows(R.Q.DRAWS)
    ids = list(dict.fromkeys(d['window_id'] for d in draws))
    if len(draws)!=1536 or len(ids)!=1491:
        raise ValueError('Original fixed TRAIN draw population required')
    plan = {e['window']['id']:e for e in R.rows(R.PLAN/'feature_plan.jsonl')}
    generator = load_generator()
    head = new_quality('cuda').train()
    parameters = [p for p in head.parameters() if p.requires_grad]
    if sum(p.numel() for p in parameters)!=165635:
        raise ValueError('Fixed quality parameter count differs')
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    pools = {}; pool_counts = Counter()
    # No TRAIN label/catalog is opened until every free candidate pool is closed.
    with (FIT/'free_pool.jsonl').open('x') as handle:
        for ordinal,ident in enumerate(ids):
            entry = plan[ident]; w = entry['window']
            if w['split']!='train' or w['seqid'] not in R.ALLOW[w['species']]:
                raise ValueError('TRAIN-only geometry required')
            cached = torch.load(ROOT/entry['artifact'],map_location='cpu',weights_only=True)
            features,dna = R.Q.window_inputs(cached,w,'cuda')
            with torch.no_grad(): hidden = generator(features,dna,[w['valid_bases']])['features'][0]
            chains,sources = free_pool(generator,hidden,cached['sequence'])
            record = {'window':w,'chains':chains,'decoder_sources':sources,
                      'reference_used_in_pool_generation':False}
            handle.write(json.dumps(record,allow_nan=False)+'\n')
            pools[ident] = record
            pool_counts.update({'windows':1,'free_candidates':len(chains)})
            if (ordinal+1)%32==0:
                handle.flush(); print(json.dumps({'phase':'free_pool','windows':ordinal+1,'total':1491}),flush=True)
            del cached,features,dna,hidden
    pool_seconds = time.perf_counter()-began
    loaded_draws,index,windows,catalogs = R.Q.load_train_inputs()
    if loaded_draws!=draws or set(windows)!=set(pools):
        raise ValueError('Do not change original TRAIN draw order or population')
    for ident,w in windows.items():
        if any(w[k]!=pools[ident]['window'][k] for k in R.GEOMETRY):
            raise ValueError('Free pool/label window geometry mismatch')
    optimizer = torch.optim.AdamW(parameters,lr=CONFIG['lr'],weight_decay=CONFIG['weight_decay'],
                                 betas=tuple(CONFIG['betas']),eps=CONFIG['eps'])
    epochs = []; step = 0
    with (FIT/'training.jsonl').open('x') as log:
        for epoch in range(3):
            counts = Counter(); total_loss = 0.; epoch_began = time.perf_counter()
            for ordinal,d in enumerate(draws):
                w = windows[d['window_id']]; started = time.perf_counter()
                cached = torch.load(ROOT/index[w['id']]['artifact'],map_location='cpu',weights_only=True)
                features,dna = R.Q.window_inputs(cached,w,'cuda')
                with torch.no_grad(): hidden = generator(features,dna,[w['valid_bases']])['features'][0]
                targets = build_targets(w,catalogs[w['species']],cached['sequence'])
                pool,origins = inject_training(pools[w['id']]['chains'],targets['positive_chains'])
                labels = candidate_targets(pool,targets)
                optimizer.zero_grad(set_to_none=True)
                scores = head(hidden,pool,w['valid_bases'])
                loss = chain_loss(scores,labels)*d['importance_weight_uniform_within_stratum']
                if not bool(torch.isfinite(loss)) or not bool(torch.isfinite(scores['logits']).all()):
                    raise ValueError('Nonfinite A6 score/loss')
                loss.backward()
                norm = torch.nn.utils.clip_grad_norm_(parameters,CONFIG['clip_norm'],error_if_nonfinite=True)
                if any(p.grad is not None for p in generator.parameters()):
                    raise ValueError('Frozen generator received gradient')
                optimizer.step(); step += 1; torch.cuda.synchronize()
                record = {'step':step,'epoch':epoch+1,'draw_ordinal':ordinal,'window_id':w['id'],
                    'species':w['species'],'importance_weight':d['importance_weight_uniform_within_stratum'],
                    'weighted_loss':float(loss.detach()),'gradient_norm_before_clip':float(norm),
                    'free_candidates':len(pools[w['id']]['chains']),
                    'reference_injected_candidates':origins.count('reference_injected'),
                    'label_counts':{str(y):labels.count(y) for y in (1,0,-1)},
                    'known_candidates':sum(y>=0 for y in labels),'seconds':time.perf_counter()-started}
                log.write(json.dumps(record,allow_nan=False)+'\n')
                counts.update({'windows':1,'free_candidates':record['free_candidates'],
                    'reference_injected_candidates':record['reference_injected_candidates'],
                    'zero_known_windows':int(record['known_candidates']==0),
                    **{'label_'+str(y):labels.count(y) for y in (1,0,-1)}})
                total_loss += record['weighted_loss']
                if step%64==0:
                    log.flush(); print(json.dumps({'phase':'fit','step':step,'total':4608}),flush=True)
                del cached,features,dna,hidden,targets,pool,origins,labels,scores,loss
            checkpoint = FIT/('step_'+str(step).zfill(6)+'.pt')
            if checkpoint.exists(): raise FileExistsError(checkpoint)
            torch.save({'arm':'A6','configuration':CONFIG,'step':step,'revision':I.F.REVISION,
                        'head':head.state_dict(),'optimizer':optimizer.state_dict(),'protocol':PROTOCOL},checkpoint)
            epochs.append({'epoch':epoch+1,'steps':1536,'counts':dict(counts),
                           'mean_weighted_loss':total_loss/1536,'seconds':time.perf_counter()-epoch_began})
    if step!=4608: raise ValueError('Incomplete fixed A6 fit')
    R.write_summary(FIT,{'arm':'A6','completed':True,'configuration':CONFIG,'protocol':PROTOCOL,
        'feature_revision':I.F.REVISION,'feature_spec':R.FEATURE_SPEC,'optimizer_steps':step,
        'epochs':epochs,'pool_counts':dict(pool_counts),'pool_generation_seconds':pool_seconds,
        'generator_trainable_parameters':0,'quality_trainable_parameters':165635,
        'final_checkpoint':str((FIT/'step_004608.pt').relative_to(ROOT)),
        'DEV_used_for_selection':False,'test_setaria_used':False,
        'wall_seconds':time.perf_counter()-began,'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
        'limitation':'additional head exposure on original TRAIN draws; not a matched-exposure causal comparison'})


@torch.inference_mode()
def quality_record(generator,head,hidden,sequence,index,ident):
    chains,sources = free_pool(generator,hidden,sequence)
    scores = head(hidden,chains,len(sequence))['logits'].float().cpu().tolist()
    if any(not torch.isfinite(torch.tensor(v)) for v in scores):
        raise ValueError('Nonfinite A6 inference score')
    return {'window_index':index,'window_id':ident,'chains':chains,'score':scores,
            'decoder_sources':sources,'score_kind':'A6_complete_chain_logit_minus_window_null',
            'quality_logit_threshold':0.,'reference_used':False,
            'counts':{'free_candidates':len(chains),'nonpositive_candidates':sum(v<=0 for v in scores)}}


@torch.inference_mode()
def infer():
    fitted = json.loads((FIT/'summary.json').read_text()); validate_head(fitted)
    if (FIT/'STATUS').read_text().strip()!='COMPLETED': raise ValueError('Terminal A6 fit required')
    state = torch.load(ROOT/fitted['final_checkpoint'],map_location='cpu',weights_only=True)
    validate_head(state,state=True)
    head = new_quality('cuda'); head.load_state_dict(state['head'],strict=True)
    head.eval().requires_grad_(False); del state
    generator = load_generator()
    backbone,tokenizer = I.load_native()
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    transcripts,scopes,lengths = [],[],{}
    for ordinal,(species,seqid,sequence) in enumerate(I.input_contigs(None)):
        length = len(sequence); lengths[seqid] = length
        locations = [[a,min(a+I.M.WIDTH,length)] for a in I.M.starts(length)]
        for strand in ('+','-'):
            path = INFER/('contig_'+str(ordinal).zfill(4)+'_'+('p' if strand=='+' else 'm')+'.jsonl')
            started = time.perf_counter(); counts = Counter()
            with path.open('x') as handle:
                for index,(a,b) in enumerate(locations):
                    raw = sequence[a:b]; oriented = raw if strand=='+' else I.M.rc(raw)
                    features = I.encode_window(backbone,tokenizer,re.sub('[^ACGT]','N',oriented),'cuda')
                    dna = one_hot(oriented,I.M.WIDTH)[None].to('cuda')
                    hidden = generator(features,dna,[b-a])['features'][0]
                    ident = species+'|'+seqid+'|'+str(a)+'|'+strand
                    free = quality_record(generator,head,hidden,oriented,index,ident)
                    handle.write(json.dumps(free,allow_nan=False)+'\n'); counts.update(free['counts'])
                    if (index+1)%32==0:
                        handle.flush(); print(json.dumps({'phase':'infer','species':species,'strand':strand,
                            'windows':index+1,'total':len(locations)}),flush=True)
                    del features,dna,hidden,free
            assembled = assemble_predictions(I.saved_records(path),locations,strand)
            for index,row in enumerate(assembled['predictions']):
                transcripts.append(as_transcript(row['chain'],seqid,strand,
                    'M28_A6_'+str(ordinal)+'_'+strand+'_'+str(index),row['score']))
            scopes.append({'species':species,'seqid':seqid,'strand':strand,'ordinal':ordinal,
                'windows':locations,'file':path.name,'quality_counts':dict(counts),
                'assembly_counts':assembled['counts'],'seconds':time.perf_counter()-started})
    if sum(len(s['windows']) for s in scopes)!=8690: raise ValueError('Incomplete natural DEV grid')
    write_gff3(INFER/'A6.gff3',transcripts,lengths,'M28_A6')
    R.write_summary(INFER,{'arm':'A6','completed':True,'inference_complete':True,
        'configuration':CONFIG,'protocol':PROTOCOL,'checkpoint':fitted['final_checkpoint'],
        'optimizer_steps':4608,'feature_revision':I.F.REVISION,'feature_spec':R.FEATURE_SPEC,
        'input_mode':'fixed_DEV_raw_FASTA','reference_used':False,'labels_or_feature_cache_used':False,
        'optimizer_steps_this_run':0,'DEV_used_for_selection':False,'test_setaria_used':False,
        'windows':8690,'scopes':scopes,'lengths':lengths,'primary_CDS_chain_predictions':len(transcripts),
        'wall_seconds':time.perf_counter()-began,'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
        'limitations':'complete raw FASTA DEV output, not independent test or accuracy until separate evaluation'})


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('fit','infer'),required=True)
    args = parser.parse_args(); torch.set_num_threads(2)
    {'fit':fit,'infer':infer}[args.mode]()
