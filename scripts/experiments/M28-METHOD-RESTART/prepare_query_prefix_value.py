#!/usr/bin/env python3
"""A9 TRAIN-only all-site query critic/gradient preparation, zero updates."""
import json
import math
import time
from collections import Counter
import torch
import run_full_train as R
import run_chain_quality as C
import run_query_known_mean as A
import analyze_query_beam_train as B
from smoke_query_beam import check_chain
from src.m28.labels import build_targets
from src.m28.query_model import unmatched_query_policy
from src.m28.query_known_mean import known_query_components
from src.m28.query_prefix_value import PrefixValueHead,EVENT_NAMES,feature_sums
from src.m28.query_prefix_value_decode import decode_prefix_values
from src.m28.query_prefix_value_training import collect_free,collect_teacher,prefix_labels,prefix_loss

OUT=R.ROOT/'outputs/M28-A9-QUERY-PREFIX-PREFLIGHT'
PROTOCOL='reports/M28-METHOD-RESTART/a9_query_prefix_preparation.md'
LIMITS=dict(max_exons=128,beam_width=4,completed_limit=4)


def main():
    source=json.loads((B.OUT/'summary.json').read_text())
    if ((B.OUT/'STATUS').read_text().strip()!='COMPLETED' or not source['completed'] or source['windows']!=8 or
        source['checkpoint']!=C.GENERATOR or source['pooled']['reference_window_pairs']!=27):
        raise ValueError('Fixed completed A5 TRAIN cohort required')
    predictions=R.rows(B.OUT/'predictions.jsonl')
    witnesses=R.rows(R.READY/'witnesses.jsonl');ids={p['window']['id'] for p in predictions}
    if len(ids)!=8 or ids!={w['window_id'] for w in witnesses} or Counter(p['window']['species'] for p in predictions)!=Counter({s:4 for s in R.ALLOW}):
        raise ValueError('Exactly existing eight TRAIN windows, no reselection')
    plan={e['window']['id']:e for e in R.rows(R.PLAN/'feature_plan.jsonl') if e['window']['id'] in ids}
    generator=C.load_generator();torch.manual_seed(0)
    head=PrefixValueHead(128).to('cuda').train()
    if sum(p.numel() for p in head.parameters())!=67745:raise ValueError('Fixed 128-dimensional critic required')
    initial={n:p.detach().clone() for n,p in head.named_parameters()}
    began=time.perf_counter();torch.cuda.reset_peak_memory_stats();free_records=[]
    # Finish all actual no-reference decoding before opening current TRAIN catalogs.
    with (OUT/'predictions.jsonl').open('x') as log:
        for pred in predictions:
            w=pred['window'];entry=plan[w['id']]
            if (w['split']!='train' or w['seqid'] not in R.ALLOW[w['species']] or w['valid_bases']!=24576 or
                any(w[k]!=entry['window'][k] for k in R.GEOMETRY) or pred['reference_used_in_decoding']):
                raise ValueError('Allowed original full-width TRAIN geometry required')
            cached=torch.load(R.ROOT/entry['artifact'],map_location='cpu',weights_only=True)
            if cached['sequence']!=pred['sequence']:raise ValueError('Frozen source DNA changed')
            features,dna=R.Q.window_inputs(cached,w,'cuda')
            with torch.no_grad():hidden=generator(features,dna,[w['valid_bases']])['features'][0]
            torch.cuda.synchronize();started=time.perf_counter()
            decoded=decode_prefix_values(generator.query,hidden,pred['sequence'],head,**LIMITS)
            if len(decoded)!=16 or [q['query_index'] for q in decoded]!=list(range(16)):
                raise ValueError('Frozen Q16 changed')
            for q in decoded:
                if len(q['completed'])>4:raise ValueError('Completion cap changed')
                for row in q['completed']:check_chain(row['chain'],pred['sequence'].upper())
            torch.cuda.synchronize()
            record=dict(window=w,queries=decoded,seconds=time.perf_counter()-started,reference_used_in_decoding=False)
            log.write(json.dumps(record,allow_nan=False)+'\n');log.flush();free_records.append(record)
            print(json.dumps(dict(phase='free_untrained_critic',window_id=w['id'],seconds=record['seconds'])),flush=True)
            del cached,features,dna,hidden
    windows,catalogs=B.load_selected_labels(predictions);records=[]
    with (OUT/'supervision.jsonl').open('x') as log:
        for pred,free in zip(predictions,free_records):
            w=windows[pred['window']['id']];started=time.perf_counter()
            if any(w[k]!=free['window'][k] for k in R.GEOMETRY):raise ValueError('No-reference/target geometry mismatch')
            cached=torch.load(R.ROOT/plan[w['id']]['artifact'],map_location='cpu',weights_only=True)
            features,dna=R.Q.window_inputs(cached,w,'cuda')
            with torch.no_grad():
                hidden=generator(features,dna,[w['valid_bases']])['features'][0]
                window=generator.query.for_window(hidden,w['valid_bases'])
            if hidden.requires_grad or hidden.is_inference():raise ValueError('Normal detached tensors required for head backward')
            sums=feature_sums(window)
            samples=collect_free(window,sums,pred['sequence'],pred['greedy_queries'],pred['beam_queries'])
            free_count=len(samples.records())
            targets=build_targets(w,catalogs[w['species']],pred['sequence'])
            if len(targets['positive_chains'])!=pred['expected_positive_chains'] or sum(map(len,targets['positive_chains']))!=pred['expected_positive_exons']:
                raise ValueError('All original R4 primary positives required')
            policy=unmatched_query_policy(targets,w['valid_bases'])
            with torch.no_grad():matched=known_query_components(window,targets['positive_chains'],A.CONFIG['loss_coefficients'][1:],policy['unmatched_label'])
            unavailable=collect_teacher(samples,sums,pred['sequence'],matched['positive_chains'])
            nodes=samples.records();labels=prefix_labels(nodes,matched['positive_chains'],matched['query_labels'],targets)
            if labels.count(1)!=2*pred['expected_positive_exons']:
                raise ValueError('Every assigned positive exon must supply its entry and end, no target dropping')
            inputs=samples.inputs();head.zero_grad(set_to_none=True);values=head(*inputs)
            loss=prefix_loss(values,labels,inputs[-1])
            if not torch.isfinite(values).all() or not torch.isfinite(loss):raise ValueError('Nonfinite critic score/loss')
            loss.backward()
            if any(p.requires_grad or p.grad is not None for p in generator.parameters()):raise ValueError('Frozen generator received gradient')
            if any(p.grad is None or not torch.isfinite(p.grad).all() for p in head.parameters()):raise ValueError('Missing/nonfinite head gradient')
            norm=math.sqrt(sum(float(p.grad.double().square().sum()) for p in head.parameters()))
            known=sum(y>=0 for y in labels)
            if (known and norm<=0) or (not known and (float(loss)!=0 or norm!=0)):raise ValueError('Known/unknown gradient policy mismatch')
            event_counts={event:{str(y):sum(n['event']==k and label==y for n,label in zip(nodes,labels)) for y in (1,0,-1)}
                          for k,event in enumerate(EVENT_NAMES)}
            teacher_only_positive=sum(y==1 and n['origins']=={'TRAIN_teacher_or_counterfactual'} for n,y in zip(nodes,labels))
            torch.cuda.synchronize()
            record=dict(window=w,assigned={str(q):chain for q,chain in matched['positive_chains'].items()},
                query_labels=matched['query_labels'],protection_reasons=policy['protection_reasons'],
                total_events=len(nodes),fixed_free_events=free_count,TRAIN_added_events=len(nodes)-free_count,
                TRAIN_only_positive_events=teacher_only_positive,grammar_unavailable_teacher_events=unavailable,
                label_counts={str(y):labels.count(y) for y in (1,0,-1)},event_label_counts=event_counts,
                known_event_mean_BCE=float(loss.detach()),head_gradient_norm=norm,generator_gradients=False,
                seconds=time.perf_counter()-started,
                samples=[dict(n,reference_relative_slot_label=y,initial_prefix_logit=float(v))
                         for n,y,v in zip(samples.public_records(),labels,values.detach().cpu())])
            log.write(json.dumps(record,allow_nan=False)+'\n');log.flush()
            records.append({k:v for k,v in record.items() if k!='samples'})
            print(json.dumps(dict(phase='supervision_backward',**records[-1]),allow_nan=False),flush=True)
            del cached,features,dna,hidden,window,sums,samples,targets,matched,nodes,labels,inputs,values,loss
    if any(not torch.equal(p.detach(),initial[n]) for n,p in head.named_parameters()):raise ValueError('Zero-update preparation changed head weights')
    R.write_summary(OUT,dict(experiment='M28-A9-QUERY-PREFIX-PREFLIGHT',completed=True,protocol=PROTOCOL,windows=len(records),
        generator_checkpoint=C.GENERATOR,generator_optimizer_steps=4608,generator_trainable_parameters=0,
        prefix_trainable_parameters=67745,critic_seed=0,precision='FP32',limits=LIMITS,
        teacher_assignment='frozen A3 known-query mean Hungarian with original coefficients',
        loss='mean known BCE per present event type, then equal event-type means; unknown -1 excluded',
        records=records,free_decode_seconds=sum(r['seconds'] for r in free_records),
        label_counts_window_sum={str(y):sum(r['label_counts'][str(y)] for r in records) for y in (1,0,-1)},
        primary_reference_window_pairs=sum(p['expected_positive_chains'] for p in predictions),
        optimizer_steps_this_run=0,head_weights_unchanged=True,checkpoints_written=False,generator_gradients=False,
        reference_used_in_decoding=False,TRAIN_reference_used_for_loss=True,new_backbone_forward=False,
        DEV_forward_used=False,test_setaria_used=False,feature_spec=R.FEATURE_SPEC,
        seconds=time.perf_counter()-began,peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitations='Untrained critic, engineering/supervision readiness only; no accuracy/novelty/independent evidence or scientific PASS. Off-policy teacher samples do not resolve B2 distribution-shift risk.'))


if __name__=='__main__':
    torch.set_num_threads(2);main()
