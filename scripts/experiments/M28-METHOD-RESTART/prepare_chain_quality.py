#!/usr/bin/env python3
"""A6 fixed TRAIN gradient preflight, no optimizer or checkpoint writing."""
import json
import math
import time
from collections import Counter
import torch
import run_full_train as R
import run_query_known_mean as A
from analyze_query_beam_train import load_selected_labels
from smoke_query_beam import check_chain, LIMITS
from src.m28.core import ChainHead, chain_loss
from src.m28.labels import build_targets, candidate_targets
from src.m28.training import inject_training
from src.m28.query_beam_pool import completed_union

SOURCE = R.ROOT/'outputs/M28-A5-TRAINED-TRAIN-DIAGNOSTIC'
OUT = R.ROOT/'outputs/M28-A6-CHAIN-QUALITY-PREFLIGHT'
PROTOCOL = 'reports/M28-METHOD-RESTART/a6_chain_quality_preparation.md'
GROUPS = ('exon','sequence','additive','nonadditive','null')


def gradient_groups(head):
    norms = {}
    for name in GROUPS:
        parameters = list(getattr(head,name).parameters())
        if any(p.grad is None or not bool(torch.isfinite(p.grad).all()) for p in parameters):
            raise ValueError('Missing/nonfinite quality gradient: '+name)
        norm = math.sqrt(sum(float(p.grad.detach().double().square().sum()) for p in parameters))
        if not math.isfinite(norm) or norm <= 0:
            raise ValueError('Zero/nonfinite quality branch gradient: '+name)
        norms[name] = norm
    return norms


def main():
    R.require_readiness()
    R.validate_feature_summary(json.loads((R.FEATURES/'summary.json').read_text()))
    source = json.loads((SOURCE/'summary.json').read_text())
    if ((SOURCE/'STATUS').read_text().strip()!='COMPLETED' or not source['completed'] or
        source['experiment']!='M28-A5-TRAINED-TRAIN-DIAGNOSTIC' or
        source['checkpoint']!='outputs/M28-A3-FIRST-FIT/step_004608.pt' or
        source['source_optimizer_steps']!=4608 or source['limits']!=LIMITS or
        source['feature_spec']!=R.FEATURE_SPEC or source['windows']!=8 or
        source['reference_used_during_decoding'] is not False or
        source['optimizer_steps_this_run']!=0 or source['checkpoints_written'] is not False or
        source['DEV_forward_used'] is not False or source['test_setaria_used'] is not False):
        raise ValueError('Completed fixed TRAIN no-reference A5 snapshot required')
    predictions = R.rows(SOURCE/'predictions.jsonl')
    witnesses = {r['window_id']:r for r in R.rows(R.READY/'witnesses.jsonl')}
    source_rows = {r['window']['id']:r for r in source['records']}
    if (len(predictions)!=8 or len(witnesses)!=8 or
        {r['window']['id'] for r in predictions}!=set(witnesses) or
        Counter(r['window']['species'] for r in predictions)!=Counter({s:4 for s in R.ALLOW})):
        raise ValueError('Exactly the eight existing TRAIN pressure sentinels required')
    plan = {r['window']['id']:r for r in R.rows(R.PLAN/'feature_plan.jsonl')}
    free_records = []
    # Close every free pool record before loading any TRAIN reference targets.
    with (OUT/'free_pool.jsonl').open('x') as log:
        for pred in predictions:
            w = pred['window']; witness = witnesses[w['id']]
            if (w['split']!='train' or w['seqid'] not in R.ALLOW[w['species']] or
                w['valid_bases']!=24576 or w!=plan[w['id']]['window'] or
                w!=source_rows[w['id']]['window'] or pred['category']!=witness['category'] or
                pred['reference_used_in_decoding'] is not False or
                len(pred['greedy_queries'])!=16 or len(pred['beam_queries'])!=16):
                raise ValueError('TRAIN snapshot geometry or source identity changed')
            chains,sources = completed_union(pred['greedy_queries'],pred['beam_queries'])
            expected = source_rows[w['id']]
            if len(chains)!=expected['greedy_unique_candidates']+expected['candidate_unique_only_beam']:
                raise ValueError('Greedy/beam union count differs from completed source')
            for chain in chains:
                check_chain(chain,pred['sequence'].upper())
            record = {'window':w,'category':pred['category'],'chains':chains,
                      'decoder_sources':sources,'reference_used_in_pool_generation':False}
            log.write(json.dumps(record,allow_nan=False)+'\n')
            free_records.append(record)
    windows,catalogs = load_selected_labels(predictions)
    fitted = json.loads((A.STAGE_DIRS['fit']/'summary.json').read_text())
    if (A.STAGE_DIRS['fit']/'STATUS').read_text().strip()!='COMPLETED':
        raise ValueError('Completed A3 generator required')
    A.check_last_identity(fitted)
    checkpoint = R.ROOT/fitted['final_checkpoint']
    state = torch.load(checkpoint,map_location='cpu',weights_only=True)
    A.check_last_identity(state,state=True)
    generator = A.new_model().eval().requires_grad_(False)
    generator.load_state_dict(state['model'],strict=True); del state
    torch.manual_seed(0)
    quality = ChainHead(A.CONFIG['hidden']).to('cuda').train()
    quality.link.requires_grad_(False)  # No link objective in completed-chain-only scorer.
    initial = {name:p.detach().clone() for name,p in quality.named_parameters()}
    torch.cuda.reset_peak_memory_stats(); began = time.perf_counter()
    records = []
    with (OUT/'records.jsonl').open('x') as log:
        for pred,free in zip(predictions,free_records):
            started = time.perf_counter(); w = windows[pred['window']['id']]
            cached = torch.load(R.ROOT/plan[w['id']]['artifact'],map_location='cpu',weights_only=True)
            features,dna = R.Q.window_inputs(cached,w,'cuda')
            if cached['sequence']!=pred['sequence']:
                raise ValueError('Frozen A5 DNA and scoring cache sequence differ')
            # Ordinary detached tensors, not inference tensors: the head needs backward.
            with torch.no_grad():
                hidden = generator(features,dna,[w['valid_bases']])['features'][0]
            if hidden.requires_grad or hidden.is_inference():
                raise ValueError('Quality head must receive normal detached generator features')
            targets = build_targets(w,catalogs[w['species']],cached['sequence'])
            if (len(targets['positive_chains'])!=pred['expected_positive_chains'] or
                sum(map(len,targets['positive_chains']))!=pred['expected_positive_exons']):
                raise ValueError('Existing TRAIN positive denominator changed')
            pool,origins = inject_training(free['chains'],targets['positive_chains'])
            labels = candidate_targets(pool,targets)
            if (origins[:len(free['chains'])]!=['free']*len(free['chains']) or
                any(label!=1 for label,origin in zip(labels,origins) if origin=='reference_injected')):
                raise ValueError('TRAIN reference injection origin or eligibility changed')
            for chain in pool:
                check_chain(chain,cached['sequence'].upper())
            quality.zero_grad(set_to_none=True)
            scores = quality(hidden,pool,w['valid_bases'])
            if not bool(torch.isfinite(scores['logits']).all()):
                raise ValueError('Nonfinite quality score')
            loss = chain_loss(scores,labels)
            if not bool(torch.isfinite(loss)):
                raise ValueError('Nonfinite known-candidate-mean BCE')
            loss.backward()
            gradients = gradient_groups(quality)
            if (any(p.requires_grad or p.grad is not None for p in generator.parameters()) or
                any(p.requires_grad or p.grad is not None for p in quality.link.parameters())):
                raise ValueError('Frozen generator/unused link received gradients')
            torch.cuda.synchronize()
            free_sources = {tuple(map(tuple,c)):s for c,s in zip(free['chains'],free['decoder_sources'])}
            record = {'window':pred['window'],'category':pred['category'],
                      'free_candidates':len(free['chains']),'training_candidates':len(pool),
                      'reference_injected_candidates':origins.count('reference_injected'),
                      'label_counts':{str(k):labels.count(k) for k in (1,0,-1)},
                      'free_label_counts':{str(k):labels[:len(free['chains'])].count(k) for k in (1,0,-1)},
                      'known_candidate_mean_BCE':float(loss.detach()),
                      'quality_gradient_norms':gradients,'generator_gradients':False,
                      'hidden_requires_grad':hidden.requires_grad,'hidden_is_inference':hidden.is_inference(),
                      'seconds':time.perf_counter()-started,
                      'candidates':[{'chain':c,'origin':o,'reference_relative_label':y,
                                     'decoder_sources':free_sources.get(tuple(map(tuple,c)),[]),
                                     'initial_quality_logit':float(v)}
                                    for c,o,y,v in zip(pool,origins,labels,scores['logits'].detach().cpu())]}
            log.write(json.dumps(record,allow_nan=False)+'\n'); log.flush()
            records.append(record)
            print(json.dumps({k:v for k,v in record.items() if k!='candidates'},allow_nan=False),flush=True)
            del cached,features,dna,hidden,targets,pool,origins,labels,scores,loss
    if any(not torch.equal(p.detach(),initial[name]) for name,p in quality.named_parameters()):
        raise ValueError('Zero-update preflight changed head weights')
    report = {'experiment':'M28-A6-CHAIN-QUALITY-PREFLIGHT','completed':True,
              'scope':'eight selected TRAIN windows; gradient/identity preparation only',
              'protocol':PROTOCOL,'source_prediction_snapshot':str(SOURCE.relative_to(R.ROOT)),
              'generator_checkpoint':str(checkpoint.relative_to(R.ROOT)),'generator_optimizer_steps':4608,
              'generator_trainable_parameters':sum(p.numel() for p in generator.parameters() if p.requires_grad),
              'quality_head':'existing ChainHead: exon/additive/ordered-GRU/nonadditive/window-null',
              'quality_trainable_parameters':sum(p.numel() for p in quality.parameters() if p.requires_grad),
              'unused_link_frozen_parameters':sum(p.numel() for p in quality.link.parameters()),
              'quality_seed':0,'precision':'FP32','loss':'mean BCE over known candidate labels; -1 excluded',
              'decoder_limits':LIMITS,'windows':len(records),
              'free_candidates_per_window_sum':sum(r['free_candidates'] for r in records),
              'reference_injected_per_window_sum':sum(r['reference_injected_candidates'] for r in records),
              'training_candidates_per_window_sum':sum(r['training_candidates'] for r in records),
              'records':[{k:v for k,v in r.items() if k!='candidates'} for r in records],
              'optimizer_steps_this_run':0,'head_weights_unchanged':True,'checkpoints_written':False,
              'generator_gradients':False,'reference_used_in_pool_generation':False,
              'TRAIN_reference_used_for_loss':True,'new_backbone_forward':False,
              'DEV_forward_used':False,'test_setaria_used':False,'feature_spec':R.FEATURE_SPEC,
              'seconds':time.perf_counter()-began,'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
              'limitations':'not trained calibration, raw-FASTA performance, precision/F1/background or scientific gate'}
    R.write_summary(OUT,report)


if __name__=='__main__':
    torch.set_num_threads(2)
    main()
