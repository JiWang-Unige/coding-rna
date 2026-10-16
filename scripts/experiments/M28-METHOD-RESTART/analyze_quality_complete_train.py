#!/usr/bin/env python3
"""Fixed eight-window TRAIN test of paired complete-pool retention, zero updates."""
import json
import time
from collections import Counter
import torch
import run_full_train as R
import run_chain_quality as C
import analyze_query_beam_train as B
from smoke_query_beam import check_chain
from src.m28.labels import build_targets, oriented_parts
from src.m28.quality_complete import decode_query_complete_pool

OUT = R.ROOT/'outputs/M28-A8-QUALITY-COMPLETE-TRAIN'
MODES = ('baseline','continued_path','quality','expanded_quality','all_examined')
LIMITS = dict(max_exons=128, beam_width=4, completed_limit=4)


def pools(pred, decoded):
    greedy = {B.chain_key(r['chain']) for r in pred['greedy_queries'] if r['terminated']}
    beam = {B.chain_key(r['chain']) for q in pred['beam_queries'] for r in q['completed']}
    path = {B.chain_key(r['chain']) for q in decoded['queries'] for r in q['path_completed']}
    quality = {B.chain_key(r['chain']) for q in decoded['queries'] for r in q['quality_completed']}
    all_examined = {B.chain_key(r['chain']) for r in decoded['all_completed']}
    return dict(baseline=greedy|beam, continued_path=greedy|path, quality=greedy|quality,
                expanded_quality=greedy|beam|quality, all_examined=greedy|all_examined)


def aggregate(references):
    return dict(reference_window_pairs=len(references),
                hits={m:sum(r['hits'][m] for r in references) for m in MODES},
                expanded_quality_new=sum(r['hits']['expanded_quality'] and not r['hits']['baseline'] for r in references),
                expanded_quality_lost=sum(r['hits']['baseline'] and not r['hits']['expanded_quality'] for r in references),
                quality_vs_continued_path_new=sum(r['hits']['quality'] and not r['hits']['continued_path'] for r in references),
                quality_vs_continued_path_lost=sum(r['hits']['continued_path'] and not r['hits']['quality'] for r in references))


@torch.inference_mode()
def main():
    old = json.loads((B.OUT/'summary.json').read_text())
    if ((B.OUT/'STATUS').read_text().strip()!='COMPLETED' or not old['completed'] or old['windows']!=8 or
        old['pooled']['reference_window_pairs']!=27 or old['limits']!=dict(beam_width=4,completed_limit_per_query=4,max_exons=128)):
        raise ValueError('Fixed completed A5 eight-window natural TRAIN cohort required')
    fitted = json.loads((C.FIT/'summary.json').read_text()); C.validate_head(fitted)
    if (C.FIT/'STATUS').read_text().strip()!='COMPLETED':
        raise ValueError('Completed fixed A6 last fit required')
    state = torch.load(R.ROOT/fitted['final_checkpoint'],map_location='cpu',weights_only=True)
    C.validate_head(state,state=True)
    head = C.new_quality('cuda').eval().requires_grad_(False)
    head.load_state_dict(state['head'],strict=True); del state
    generator = C.load_generator()
    predictions = R.rows(B.OUT/'predictions.jsonl')
    witnesses = R.rows(R.READY/'witnesses.jsonl')
    ids = {p['window']['id'] for p in predictions}
    if (len(predictions)!=8 or len(ids)!=8 or ids!={w['window_id'] for w in witnesses} or
        Counter(p['window']['species'] for p in predictions)!=Counter({s:4 for s in R.ALLOW})):
        raise ValueError('Do not change original TRAIN pressure sentinels')
    plan = {e['window']['id']:e for e in R.rows(R.PLAN/'feature_plan.jsonl') if e['window']['id'] in ids}
    began = time.perf_counter(); torch.cuda.reset_peak_memory_stats()
    decoded_windows, timing, control_differences = {}, [], []
    # No current reference geometry/catalog read until all free records are saved.
    with (OUT/'predictions.jsonl').open('x') as log:
        for pred in predictions:
            w = pred['window']; entry = plan[w['id']]
            if (w['split']!='train' or w['seqid'] not in R.ALLOW[w['species']] or w['valid_bases']!=24576 or
                any(w[k]!=entry['window'][k] for k in R.GEOMETRY) or pred['reference_used_in_decoding']):
                raise ValueError('Fixed allowed full-width TRAIN geometry required')
            cached = torch.load(R.ROOT/entry['artifact'],map_location='cpu',weights_only=True)
            if cached['sequence']!=pred['sequence']:
                raise ValueError('Unchanged natural source DNA required')
            features,dna = R.Q.window_inputs(cached,w,'cuda')
            hidden = generator(features,dna,[w['valid_bases']])['features'][0]
            def quality(chains):
                return head(hidden,chains,w['valid_bases'])['logits'].float().cpu().tolist()
            torch.cuda.synchronize(); started=time.perf_counter()
            decoded = decode_query_complete_pool(generator.query,hidden,pred['sequence'],quality,**LIMITS)
            torch.cuda.synchronize(); seconds=time.perf_counter()-started
            actual = pools(pred,decoded)
            for chain in actual['all_examined']|actual['baseline']:
                check_chain(chain,pred['sequence'].upper())
            # Removing a valid cumulative-score bound should not change path top4.
            for q, original in zip(decoded['queries'],pred['beam_queries']):
                previous={B.chain_key(c['chain']):c['path_log_score'] for c in original['completed']}
                current={B.chain_key(c['chain']):c['path_log_score'] for c in q['path_completed']}
                if set(previous)!=set(current) or any(abs(previous[c]-current[c])>.001 for c in previous):
                    control_differences.append(dict(window_id=w['id'],query_index=q['query_index']))
            if len(decoded['queries'])!=16 or [q['query_index'] for q in decoded['queries']]!=list(range(16)):
                raise ValueError('Frozen Q16 changed')
            record=dict(window=w,decoded=decoded,seconds=seconds,reference_used_in_decoding=False)
            log.write(json.dumps(record,allow_nan=False)+'\n');log.flush()
            decoded_windows[w['id']]=decoded; timing.append(dict(window_id=w['id'],seconds=seconds))
            print(json.dumps(dict(window_id=w['id'],seconds=seconds,unique_candidates={k:len(v) for k,v in actual.items()})),flush=True)
            del cached,features,dna,hidden
    windows,catalogs = B.load_selected_labels(predictions)
    original_refs={ (r['window_id'],r['reference_id']):r for r in R.rows(B.OUT/'references.jsonl') }
    references,records=[],[]
    for pred in predictions:
        w=windows[pred['window']['id']]; actual=pools(pred,decoded_windows[w['id']])
        targets=build_targets(w,catalogs[w['species']],pred['sequence'])
        rs=[]
        for ident in w['positive_ids']:
            chain=B.chain_key(oriented_parts(catalogs[w['species']][ident],w))
            oldref=original_refs[(w['id'],ident)]
            if chain!=B.chain_key(oldref['chain']):raise ValueError('Original reference geometry changed')
            rs.append(dict(window_id=w['id'],species=w['species'],reference_id=ident,chain=chain,
                           exons=len(chain),exon_class=oldref['exon_class'],hits={m:chain in actual[m] for m in MODES}))
        if len(rs)!=pred['expected_positive_chains'] or sum(r['exons'] for r in rs)!=pred['expected_positive_exons']:
            raise ValueError('Keep all original eligible reference/window pairs')
        references.extend(rs)
        records.append(dict(window_id=w['id'],**aggregate(rs),unique_candidates={m:len(actual[m]) for m in MODES},
                            reference_relative_labels={m:B.label_counts(actual[m],targets) for m in MODES}))
    if {(r['window_id'],r['reference_id']) for r in references}!=set(original_refs) or len(references)!=27:
        raise ValueError('Do not drop or duplicate original 27 identities')
    with (OUT/'references.jsonl').open('x') as log:
        for r in references:log.write(json.dumps(r,allow_nan=False)+'\n')
    pooled=aggregate(references)
    counts={m:sum(r['unique_candidates'][m] for r in records) for m in MODES}
    labels={m:{str(y):sum(r['reference_relative_labels'][m][str(y)] for r in records) for y in (1,0,-1)} for m in MODES}
    R.write_summary(OUT,dict(experiment='M28-A8-QUALITY-COMPLETE-TRAIN',completed=True,windows=8,limits=LIMITS,
        generator_checkpoint=C.GENERATOR,quality_checkpoint=fitted['final_checkpoint'],records=records,pooled=pooled,
        per_species={s:aggregate([r for r in references if r['species']==s]) for s in R.ALLOW},
        by_exon_class={c:aggregate([r for r in references if r['exon_class']==c]) for c in ('one','two_to_four','five_or_more')},
        unique_candidates_window_sum=counts,reference_relative_labels_window_sum=labels,
        path_control_differences=control_differences,all_returned_chains_DNA_legal=True,
        mechanism_preparation_supported=not control_differences and pooled['expanded_quality_new']>0 and
                                       pooled['quality_vs_continued_path_new']>0 and pooled['expanded_quality_lost']==0,
        reference_used_in_decoding=False,TRAIN_reference_used_after_prediction=True,optimizer_steps_this_run=0,
        checkpoints_written=False,trainable_parameters_this_run=0,new_backbone_forward=False,
        DEV_forward_used=False,test_setaria_used=False,window_seconds=timing,seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitations='Metadata-selected TRAIN cohort, not independent accuracy. Same continued histories, different complete-pool readout. Extra compute and candidate exposure; no F1/background/claim gate.'))


if __name__=='__main__':
    torch.set_num_threads(2)
    main()
