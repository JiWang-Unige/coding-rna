#!/usr/bin/env python3
"""Fixed A3 last-weight TRAIN beam diagnostic; never new fit or DEV scoring."""
import json
import math
import time
from collections import Counter
import torch
import run_full_train as R
import run_query_known_mean as A
from smoke_query_beam import check_chain, LIMITS
from src.m28.query_decode import decode_query_chains
from src.m28.query_beam import decode_query_beam
from src.m28.labels import build_targets, candidate_targets, oriented_parts

OUT = R.ROOT/'outputs/M28-A5-TRAINED-TRAIN-DIAGNOSTIC'
STATES = ('both','beam_only','greedy_only','neither')


def chain_key(chain):
    return tuple(tuple(exon) for exon in chain)


def compare_references(references, greedy, beam):
    """Every supplied reference identity is retained, even equal geometries."""
    counts = Counter({k:0 for k in STATES}); rows = []
    for ref in references:
        key = chain_key(ref['chain'])
        g,b = key in greedy,key in beam
        state = 'both' if g and b else 'greedy_only' if g else 'beam_only' if b else 'neither'
        counts[state] += 1
        rows.append({**ref,'greedy_free_hit':g,'beam_free_hit':b,'state':state})
    return {'reference_window_pairs':len(references),'states':dict(counts),
            'greedy_free_hits':counts['both']+counts['greedy_only'],
            'beam_free_hits':counts['both']+counts['beam_only']}, rows


def label_counts(chains, targets):
    labels = candidate_targets(sorted(chains),targets)
    return {str(k):labels.count(k) for k in (1,0,-1)}


def aggregate(rows):
    counts = Counter({k:0 for k in STATES})
    for r in rows:
        counts[r['state']] += 1
    return {'reference_window_pairs':len(rows),'states':dict(counts),
            'greedy_free_hits':counts['both']+counts['greedy_only'],
            'beam_free_hits':counts['both']+counts['beam_only']}


def load_selected_labels(predictions):
    selected = {r['window']['id']:r for r in predictions}
    windows,catalogs = {},{}
    for species in R.ALLOW:
        ws = [w for w in R.train_records(R.LABELS/species/'windows_train.jsonl')
              if w['id'] in selected]
        needed = {tid for w in ws for tid in w['overlapping_transcript_ids']}
        records = [t for t in R.train_records(R.LABELS/species/'transcripts.jsonl') if t['id'] in needed]
        catalog = {t['id']:t for t in records}
        if len(records)!=len(catalog) or set(catalog)!=needed:
            raise ValueError('Selected TRAIN reference join incomplete or duplicated')
        catalogs[species]=catalog
        for w in ws:
            expected=selected[w['id']]['window']
            if (w['id'] in windows or w['species']!=species or w['split']!='train' or
                w['seqid'] not in R.ALLOW[species] or
                any(w[k]!=expected[k] for k in R.GEOMETRY) or
                w['label_policy']!='CDS_assessable_primary_R4'):
                raise ValueError('Unchanged selected TRAIN R4 geometry required')
            for tid in w['overlapping_transcript_ids']:
                tx=catalog[tid]
                if tx['split']!='train' or (tx['seqid'],tx['strand'])!=(w['seqid'],w['strand']):
                    raise ValueError('Selected TRAIN transcript identity mismatch')
                if tid in w['positive_ids'] and not (tx['eligible_positive'] and tx['training_primary']):
                    raise ValueError('Existing positive eligibility changed')
            windows[w['id']]=w
    if set(windows)!=set(selected):
        raise ValueError('Do not drop selected TRAIN windows')
    return windows,catalogs


@torch.inference_mode()
def main():
    R.require_readiness()
    R.validate_feature_summary(json.loads((R.FEATURES/'summary.json').read_text()))
    witnesses=R.rows(R.READY/'witnesses.jsonl')
    plan={r['window']['id']:r for r in R.rows(R.PLAN/'feature_plan.jsonl')}
    if (len(witnesses)!=8 or len({r['window_id'] for r in witnesses})!=8 or
        Counter(r['species'] for r in witnesses)!=Counter({s:4 for s in R.ALLOW})):
        raise ValueError('Exactly the pre-existing eight TRAIN sentinels required')
    fitted=json.loads((A.STAGE_DIRS['fit']/'summary.json').read_text())
    if (A.STAGE_DIRS['fit']/'STATUS').read_text().strip()!='COMPLETED':
        raise ValueError('A3 fit must be completed')
    A.check_last_identity(fitted)
    checkpoint=R.ROOT/fitted['final_checkpoint']
    state=torch.load(checkpoint,map_location='cpu',weights_only=True)
    A.check_last_identity(state,state=True)
    model=A.new_model().eval().requires_grad_(False)
    model.load_state_dict(state['model'],strict=True);del state
    draw_counts=Counter(r['window_id'] for r in R.rows(R.Q.DRAWS))
    torch.cuda.reset_peak_memory_stats();began=time.perf_counter()
    predictions=[]
    # No R4 window/catalog reference read until every free record is saved.
    with (OUT/'predictions.jsonl').open('x') as log:
        for witness in witnesses:
            entry=plan[witness['window_id']];w=entry['window']
            if w['split']!='train' or w['seqid'] not in R.ALLOW[w['species']] or w['valid_bases']!=24576:
                raise ValueError('Fixed full-width allowed TRAIN window required')
            cached=torch.load(R.ROOT/entry['artifact'],map_location='cpu',weights_only=True)
            features,dna=R.Q.window_inputs(cached,w,'cuda')
            hidden=model(features,dna,[w['valid_bases']])['features'][0]
            sequence=cached['sequence']
            torch.cuda.synchronize();started=time.perf_counter()
            greedy=decode_query_chains(model.query,hidden,sequence,128)
            torch.cuda.synchronize();greedy_seconds=time.perf_counter()-started
            started=time.perf_counter()
            beam=decode_query_beam(model.query,hidden,sequence,128,4,4)
            torch.cuda.synchronize();beam_seconds=time.perf_counter()-started
            if len(greedy)!=16 or len(beam)!=16:
                raise ValueError('Frozen Q16 changed')
            gc={chain_key(r['chain']) for r in greedy if r['terminated']}
            bc={chain_key(r['chain']) for q in beam for r in q['completed']}
            for chain in gc|bc:
                check_chain(chain,sequence.upper())
            if any(not math.isfinite(r['query_logit']) for r in greedy):
                raise ValueError('Nonfinite original query logit')
            record={'window':w,'category':witness['category'],
                    'expected_positive_chains':witness['positive_chains'],
                    'expected_positive_exons':witness['positive_exons'],
                    'original_A3_draw_occurrences_per_pass':draw_counts[w['id']],
                    'sequence':sequence,'greedy_queries':greedy,'beam_queries':beam,
                    'greedy_seconds':greedy_seconds,'beam_seconds':beam_seconds,
                    'reference_used_in_decoding':False}
            log.write(json.dumps(record,allow_nan=False)+'\n');log.flush()
            predictions.append(record)
            print(json.dumps({'window':w['id'],'greedy_unique':len(gc),'beam_unique':len(bc),
                              'beam_seconds':beam_seconds}),flush=True)
            del cached,features,dna,hidden,greedy,beam
    prediction_seconds=time.perf_counter()-began
    windows,catalogs=load_selected_labels(predictions)
    references=[];records=[]
    for pred in predictions:
        w=windows[pred['window']['id']]
        targets=build_targets(w,catalogs[w['species']],pred['sequence'])
        gc={chain_key(r['chain']) for r in pred['greedy_queries'] if r['terminated']}
        bc={chain_key(r['chain']) for q in pred['beam_queries'] for r in q['completed']}
        refs=[]
        for ident in w['positive_ids']:
            chain=oriented_parts(catalogs[w['species']][ident],w)
            refs.append({'window_id':w['id'],'species':w['species'],'reference_id':ident,
                         'chain':chain,'exons':len(chain),
                         'exon_class':'one' if len(chain)==1 else 'two_to_four' if len(chain)<=4 else 'five_or_more',
                         'original_A3_draw_occurrences_per_pass':pred['original_A3_draw_occurrences_per_pass']})
        if (len(refs)!=pred['expected_positive_chains'] or
            sum(r['exons'] for r in refs)!=pred['expected_positive_exons']):
            raise ValueError('TRAIN reference denominator differs from frozen sentinel selection')
        comparison,observed=compare_references(refs,gc,bc)
        references.extend(observed)
        records.append({'window':pred['window'],'category':pred['category'],**comparison,
                        'original_A3_draw_occurrences_per_pass':pred['original_A3_draw_occurrences_per_pass'],
                        'greedy_unique_candidates':len(gc),'beam_unique_candidates':len(bc),
                        'candidate_unique_both':len(gc&bc),'candidate_unique_only_beam':len(bc-gc),
                        'candidate_unique_only_greedy':len(gc-bc),
                        'greedy_candidate_reference_relative_labels':label_counts(gc,targets),
                        'beam_candidate_reference_relative_labels':label_counts(bc,targets),
                        'greedy_completed_exon_histogram':dict(Counter(str(len(c)) for c in gc)),
                        'beam_completed_exon_histogram':dict(Counter(str(len(c)) for c in bc)),
                        'greedy_seconds':pred['greedy_seconds'],'beam_seconds':pred['beam_seconds']})
    with (OUT/'references.jsonl').open('x') as log:
        for r in references:log.write(json.dumps(r,allow_nan=False)+'\n')
    report={'experiment':'M28-A5-TRAINED-TRAIN-DIAGNOSTIC','completed':True,
            'scope':'fixed TRAIN pressure sentinels; not independent accuracy, DEV result, selection or new fitting',
            'checkpoint':str(checkpoint.relative_to(R.ROOT)),'source_optimizer_steps':4608,
            'limits':LIMITS,'windows':len(records),'records':records,'pooled':aggregate(references),
            'per_species':{s:aggregate([r for r in references if r['species']==s]) for s in R.ALLOW},
            'by_exon_class':{k:aggregate([r for r in references if r['exon_class']==k])
                             for k in ('one','two_to_four','five_or_more')},
            'by_original_A3_draw_exposure':{k:aggregate([r for r in references if
                bool(r['original_A3_draw_occurrences_per_pass'])==exposed])
                for k,exposed in [('drawn',True),('not_drawn',False)]},
            'optimizer_steps_this_run':0,'trainable_head_parameters':sum(p.numel() for p in model.parameters() if p.requires_grad),
            'checkpoint_read':True,'checkpoints_written':False,'new_backbone_forward':False,
            'reference_used_during_decoding':False,'TRAIN_reference_used_after_prediction':True,
            'DEV_forward_used':False,'test_setaria_used':False,
            'feature_spec':R.FEATURE_SPEC,'prediction_seconds':prediction_seconds,
            'seconds':time.perf_counter()-began,
            'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),
            'limitations':'metadata-selected small TRAIN cohort; eligibility denominator only; no precision/F1/background/gate or biological truth claim'}
    R.write_summary(OUT,report)


if __name__=='__main__':
    torch.set_num_threads(2)
    main()
