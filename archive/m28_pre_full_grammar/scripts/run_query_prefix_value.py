#!/usr/bin/env python3
"""Fixed A9 critic-only fit and TRAIN paired generation probe; no DEV execution."""
import argparse,json,time
from collections import Counter
import torch
import run_full_train as R
import run_chain_quality as C
import run_query_known_mean as A
import analyze_query_beam_train as B
from smoke_query_beam import check_chain
from src.m28.labels import build_targets,oriented_parts,candidate_targets
from src.m28.query_model import unmatched_query_policy
from src.m28.query_known_mean import known_query_components
from src.m28.query_decode import decode_query_chains
from src.m28.query_beam import decode_query_beam
from src.m28.query_prefix_value import PrefixValueHead,EVENT_NAMES,feature_sums
from src.m28.query_prefix_value_training import collect_free,collect_teacher,prefix_labels,prefix_loss
from src.m28.query_prefix_value_decode import decode_prefix_values,complete_query_values

FIT=R.ROOT/'outputs/M28-A9-QUERY-PREFIX-FIT'
PROBE=R.ROOT/'outputs/M28-A9-QUERY-PREFIX-TRAIN-PROBE'
PROTOCOL='reports/M28-METHOD-RESTART/a9_query_prefix_fit_protocol.md'
CONFIG=dict(generator=C.GENERATOR,generator_steps=4608,generator_frozen=True,hidden=128,head_parameters=67745,
    seed=0,precision='FP32',draws=1536,unique_windows=1491,passes=3,steps=4608,
    loss='mean known BCE per present event then equal event means',teacher='frozen_A3_original_Hungarian',
    free_sampling='all_initial_ATG_plus_original_greedy_beam4_prefixes',
    teacher_sampling='all_true_prefixes_raw_top4_entry_counterfactuals_all_legal_end_counterfactuals',
    max_exons=128,beam_width=4,completed_limit=4,lr=3e-4,weight_decay=.01,betas=[.9,.999],eps=1e-8,clip_norm=1.,
    checkpoint='last_step_004608_only',importance='original_uniform_within_stratum')
LIMITS={k:CONFIG[k] for k in ('max_exons','beam_width','completed_limit')}


def new_head():
    torch.manual_seed(0);return PrefixValueHead(128).to('cuda')


def identity(record,state=False):
    if record['configuration']!=CONFIG or record['protocol']!=PROTOCOL or record['step' if state else 'optimizer_steps']!=4608:
        raise ValueError('Fixed A9 last checkpoint required')
    if not state and (not record['completed'] or record['DEV_used_for_selection'] or record['test_setaria_used']):
        raise ValueError('Complete TRAIN-only fit required')


def fit():
    pre=R.ROOT/'outputs/M28-A9-QUERY-PREFIX-PREFLIGHT';prepared=json.loads((pre/'summary.json').read_text())
    if ((pre/'STATUS').read_text().strip()!='COMPLETED' or not prepared['completed'] or
        prepared['generator_checkpoint']!=C.GENERATOR or prepared['prefix_trainable_parameters']!=67745 or
        prepared['optimizer_steps_this_run']!=0 or not prepared['head_weights_unchanged'] or prepared['windows']!=8):
        raise ValueError('Completed zero-update preflight required')
    R.validate_feature_summary(json.loads((R.FEATURES/'summary.json').read_text()))
    draws=R.rows(R.Q.DRAWS);ids=list(dict.fromkeys(d['window_id'] for d in draws))
    if len(draws)!=1536 or len(ids)!=1491:raise ValueError('Original TRAIN population required')
    plan={e['window']['id']:e for e in R.rows(R.PLAN/'feature_plan.jsonl')}
    generator=C.load_generator();head=new_head().train();parameters=list(head.parameters())
    if sum(p.numel() for p in parameters)!=67745:raise ValueError('Unexpected new-head capacity')
    began=time.perf_counter();torch.cuda.reset_peak_memory_stats();pools={}
    # Every original free trace is closed before current TRAIN catalogs are read.
    with (FIT/'free_traces.jsonl').open('x') as log:
        for ordinal,ident in enumerate(ids):
            e=plan[ident];w=e['window']
            if w['split']!='train' or w['seqid'] not in R.ALLOW[w['species']]:raise ValueError('TRAIN-only free traces required')
            cached=torch.load(R.ROOT/e['artifact'],map_location='cpu',weights_only=True)
            features,dna=R.Q.window_inputs(cached,w,'cuda')
            with torch.no_grad():hidden=generator(features,dna,[w['valid_bases']])['features'][0]
            greedy=decode_query_chains(generator.query,hidden,cached['sequence'],128,w['valid_bases'])
            beam=decode_query_beam(generator.query,hidden,cached['sequence'],128,4,4,w['valid_bases'])
            record=dict(window=w,greedy_queries=greedy,beam_queries=beam,reference_used_in_decoding=False)
            log.write(json.dumps(record,allow_nan=False)+'\n');pools[ident]=record
            if (ordinal+1)%32==0:log.flush();print(json.dumps(dict(phase='fixed_free_trace',windows=ordinal+1,total=1491)),flush=True)
            del cached,features,dna,hidden
    pool_seconds=time.perf_counter()-began
    loaded,index,windows,catalogs=R.Q.load_train_inputs()
    if loaded!=draws or set(windows)!=set(pools):raise ValueError('Original draw order/identity changed')
    for ident,w in windows.items():
        if any(w[k]!=pools[ident]['window'][k] for k in R.GEOMETRY):raise ValueError('Free trace/label mismatch')
    optimizer=torch.optim.AdamW(parameters,lr=CONFIG['lr'],weight_decay=CONFIG['weight_decay'],betas=tuple(CONFIG['betas']),eps=CONFIG['eps'])
    step=0;epochs=[]
    with (FIT/'training.jsonl').open('x') as log:
        for epoch in range(3):
            counts=Counter();sums=0.;started=time.perf_counter()
            for ordinal,d in enumerate(draws):
                w=windows[d['window_id']];t=time.perf_counter();pool=pools[w['id']]
                cached=torch.load(R.ROOT/index[w['id']]['artifact'],map_location='cpu',weights_only=True)
                features,dna=R.Q.window_inputs(cached,w,'cuda')
                with torch.no_grad():
                    hidden=generator(features,dna,[w['valid_bases']])['features'][0]
                    window=generator.query.for_window(hidden,w['valid_bases'])
                samples=collect_free(window,feature_sums(window),cached['sequence'],pool['greedy_queries'],pool['beam_queries'])
                free_count=len(samples.records());targets=build_targets(w,catalogs[w['species']],cached['sequence'])
                policy=unmatched_query_policy(targets,w['valid_bases'])
                with torch.no_grad():matched=known_query_components(window,targets['positive_chains'],A.CONFIG['loss_coefficients'][1:],policy['unmatched_label'])
                unavailable=collect_teacher(samples,feature_sums(window),cached['sequence'],matched['positive_chains'])
                nodes=samples.records();labels=prefix_labels(nodes,matched['positive_chains'],matched['query_labels'],targets)
                if labels.count(1)!=2*sum(map(len,targets['positive_chains'])):raise ValueError('Primary positive events dropped')
                optimizer.zero_grad(set_to_none=True);inputs=samples.inputs();values=head(*inputs)
                loss=prefix_loss(values,labels,inputs[-1])*d['importance_weight_uniform_within_stratum']
                if not torch.isfinite(values).all() or not torch.isfinite(loss):raise ValueError('Nonfinite critic fit')
                loss.backward();norm=torch.nn.utils.clip_grad_norm_(parameters,1.,error_if_nonfinite=True)
                if any(p.grad is not None for p in generator.parameters()):raise ValueError('Frozen generator gradient')
                optimizer.step();step+=1;torch.cuda.synchronize()
                record=dict(step=step,epoch=epoch+1,draw_ordinal=ordinal,window_id=w['id'],species=w['species'],
                    importance_weight=d['importance_weight_uniform_within_stratum'],weighted_loss=float(loss.detach()),
                    gradient_norm_before_clip=float(norm),fixed_free_events=free_count,total_events=len(nodes),
                    label_counts={str(y):labels.count(y) for y in (1,0,-1)},
                    grammar_unavailable_teacher_events=len(unavailable),generator_gradients=False,seconds=time.perf_counter()-t)
                log.write(json.dumps(record,allow_nan=False)+'\n');sums+=record['weighted_loss']
                counts.update({str(y):labels.count(y) for y in (1,0,-1)});counts.update(dict(windows=1,zero_known_windows=int(not any(y>=0 for y in labels)),grammar_unavailable_teacher_events=len(unavailable)))
                if step%64==0:log.flush();print(json.dumps(dict(phase='fit',step=step,total=4608)),flush=True)
                del cached,features,dna,hidden,window,samples,targets,matched,nodes,labels,inputs,values,loss
            torch.save(dict(configuration=CONFIG,protocol=PROTOCOL,step=step,head=head.state_dict(),optimizer=optimizer.state_dict()),FIT/('step_'+str(step).zfill(6)+'.pt'))
            epochs.append(dict(epoch=epoch+1,steps=1536,counts=dict(counts),mean_weighted_loss=sums/1536,seconds=time.perf_counter()-started))
    if step!=4608:raise ValueError('Incomplete fixed fit')
    R.write_summary(FIT,dict(completed=True,configuration=CONFIG,protocol=PROTOCOL,optimizer_steps=step,
        final_checkpoint=str((FIT/'step_004608.pt').relative_to(R.ROOT)),epochs=epochs,free_trace_windows=len(pools),
        free_trace_seconds=pool_seconds,seconds=time.perf_counter()-began,peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        generator_trainable_parameters=0,prefix_trainable_parameters=67745,DEV_used_for_selection=False,test_setaria_used=False,
        limitations='Additional critic exposure on existing TRAIN, off-policy; not matched total exposure or evidence of improved free annotation.'))


def candidate_readout_report(arms,targets):
    """Reference-relative pool burden, not annotation or biological precision.

    The same geometry may occur at several query slots. The zero-value view
    explicitly maxes over those slots and is only a diagnostic of this pool.
    """
    gains={};records={};query_maps={}
    for arm,rows in arms.items():
        scores={};query_maps[arm]={}
        for row in rows:
            chain=B.chain_key(row['chain']);value=row['shared_complete_value']
            query_maps[arm][(row['query_index'],chain)]=value
            scores[chain]=max(value,scores.get(chain,float('-inf')))
        pool=sorted(scores);positive=[c for c in pool if scores[c]>0]
        labels=candidate_targets(pool,targets);positive_labels=candidate_targets(positive,targets)
        records[arm]=dict(completed_records=len(rows),unique_query_geometry=len(query_maps[arm]),
            unique_geometry=len(pool),positive_value_geometry=len(positive),
            reference_relative_labels={str(y):labels.count(y) for y in (1,0,-1)},
            positive_value_reference_relative_labels={str(y):positive_labels.count(y) for y in (1,0,-1)})
        gains[arm]=scores
    common=query_maps['early'].keys()&query_maps['late'].keys()
    equal=all(query_maps['early'][k]==query_maps['late'][k] for k in common)
    return records,gains,equal


@torch.no_grad()
def probe():
    fitted=json.loads((FIT/'summary.json').read_text());identity(fitted)
    if (FIT/'STATUS').read_text().strip()!='COMPLETED':raise ValueError('Terminal fit required')
    state=torch.load(R.ROOT/fitted['final_checkpoint'],map_location='cpu',weights_only=True);identity(state,True)
    head=new_head().eval().requires_grad_(False);head.load_state_dict(state['head'],strict=True);del state
    generator=C.load_generator();predictions=R.rows(B.OUT/'predictions.jsonl')
    if len(predictions)!=8 or len({p['window']['id'] for p in predictions})!=8:raise ValueError('Fixed eight-window TRAIN probe required')
    plan={e['window']['id']:e for e in R.rows(R.PLAN/'feature_plan.jsonl')};free={}
    began=time.perf_counter();torch.cuda.reset_peak_memory_stats()
    with (PROBE/'predictions.jsonl').open('x') as log:
        for pred in predictions:
            w=pred['window']
            if w['split']!='train' or w['seqid'] not in R.ALLOW[w['species']] or pred['reference_used_in_decoding']:
                raise ValueError('Original allowed TRAIN baseline required')
            cached=torch.load(R.ROOT/plan[w['id']]['artifact'],map_location='cpu',weights_only=True)
            if cached['sequence']!=pred['sequence']:raise ValueError('Frozen source DNA changed')
            features,dna=R.Q.window_inputs(cached,w,'cuda');hidden=generator(features,dna,[w['valid_bases']])['features'][0]
            early=decode_prefix_values(generator.query,hidden,pred['sequence'],head,**LIMITS)
            late=[r for r in pred['greedy_queries'] if r['terminated']]+[r for q in pred['beam_queries'] for r in q['completed']]
            early_rows=[r for q in early for r in q['completed']]
            identities=[(r['query_index'],B.chain_key(r['chain'])) for r in late+early_rows]
            values=complete_query_values(generator.query,hidden,pred['sequence'],head,identities)
            arms={arm:[dict(query_index=r['query_index'],chain=r['chain'],shared_complete_value=values[(r['query_index'],B.chain_key(r['chain']))]) for r in rows]
                  for arm,rows in (('early',early_rows),('late',late))}
            for q,c in values:check_chain(c,pred['sequence'].upper())
            record=dict(window=w,arms=arms,early_counts=[q['counts'] for q in early],reference_used_in_decoding=False)
            log.write(json.dumps(record,allow_nan=False)+'\n');log.flush();free[w['id']]=record
            del cached,features,dna,hidden
    windows,catalogs=B.load_selected_labels(predictions);references=[];candidate_records=[];gain_maps={};same_values=True
    for pred in predictions:
        w=windows[pred['window']['id']]
        targets=build_targets(w,catalogs[w['species']],pred['sequence'])
        report,gains,equal=candidate_readout_report(free[w['id']]['arms'],targets)
        candidate_records.append(dict(window_id=w['id'],species=w['species'],arms=report))
        gain_maps[w['id']]=gains;same_values &= equal
    original=R.rows(B.OUT/'references.jsonl')
    expected={(w['id'],ident) for w in windows.values() for ident in w['positive_ids']}
    if len(expected)!=27 or len(original)!=27 or {(r['window_id'],r['reference_id']) for r in original}!=expected:
        raise ValueError('All original 27 primary reference/window pairs required')
    with (PROBE/'references.jsonl').open('x') as log:
        for ref in original:
            w=windows[ref['window_id']]
            if B.chain_key(oriented_parts(catalogs[w['species']][ref['reference_id']],w))!=B.chain_key(ref['chain']):
                raise ValueError('Original primary geometry changed')
            scores=gain_maps[ref['window_id']];chain=B.chain_key(ref['chain'])
            hits={a:chain in scores[a] for a in scores}
            state='both' if all(hits.values()) else ('early_only' if hits['early'] else ('late_only' if hits['late'] else 'neither'))
            record=dict(ref,hits=hits,free_recovery_state=state,
                        positive_value_hits={a:scores[a].get(chain,float('-inf'))>0 for a in scores})
            log.write(json.dumps(record,allow_nan=False)+'\n');references.append(record)
    counts={a:sum(r['hits'][a] for r in references) for a in ('early','late')}
    new=[r for r in references if r['hits']['early'] and not r['hits']['late']]
    lost=[r for r in references if r['hits']['late'] and not r['hits']['early']]
    gate=counts['early']>counts['late'] and any(r['exons']>=5 for r in new) and same_values
    R.write_summary(PROBE,dict(completed=True,configuration=CONFIG,protocol=PROTOCOL,windows=8,primary_reference_window_pairs=len(references),
        counts=counts,new_free_reference_pairs=len(new),lost_free_reference_pairs=len(lost),
        new_multiexon_reference_pairs=sum(r['exons']>=5 for r in new),common_query_chain_values_identical=same_values,
        per_species={s:{a:sum(r['hits'][a] for r in references if r['species']==s) for a in counts} for s in R.ALLOW},
        free_recovery_states={state:sum(r['free_recovery_state']==state for r in references)
                              for state in ('both','early_only','late_only','neither')},
        per_species_recovery_states={s:{state:sum(r['species']==s and r['free_recovery_state']==state for r in references)
            for state in ('both','early_only','late_only','neither')} for s in R.ALLOW},
        by_exon_class={c:{a:sum(r['hits'][a] for r in references if r['exon_class']==c) for a in counts}
                       for c in ('one','two_to_four','five_or_more')},
        positive_value_reference_pairs={a:sum(r['positive_value_hits'][a] for r in references) for a in counts},
        candidate_records=candidate_records,
        candidate_counts_window_sum={a:{k:sum(r['arms'][a][k] for r in candidate_records)
            for k in ('completed_records','unique_query_geometry','unique_geometry','positive_value_geometry')} for a in counts},
        candidate_reference_relative_labels_window_sum={a:{str(y):sum(r['arms'][a]['reference_relative_labels'][str(y)] for r in candidate_records)
            for y in (1,0,-1)} for a in counts},
        preparation_gate_for_full_DEV=gate,optimizer_steps_this_run=0,checkpoint=fitted['final_checkpoint'],
        reference_used_in_decoding=False,TRAIN_reference_used_after_prediction=True,DEV_used=False,test_setaria_used=False,
        seconds=time.perf_counter()-began,peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitations='Metadata-selected TRAIN diagnostic, not independent accuracy. Late is original G+beam4, early critic-only beam4; shared STOP value does not equalize candidate exposure or prove F1/background.'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=('fit','probe'));args=p.parse_args()
    torch.set_num_threads(2);globals()[args.stage]()
