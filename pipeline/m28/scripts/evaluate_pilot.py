#!/usr/bin/env python3
"""Post-inference assembly and frozen-reference scoring, only after both complete arms."""
import argparse
import json
from collections import Counter,defaultdict
from pathlib import Path
import common_ruler as E
import manifest as M
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript,write_gff3

def rows(path):
    with path.open() as f: return [json.loads(line) for line in f if line.strip()]

def unique_chains(chains):
    return {tuple(tuple(p) for p in chain) for chain in chains}

def stage_audit(assembled,refs):
    reference=unique_chains([[(a,b) for a,b,_ in t['CDS']] for t in refs])
    hits={k:len(unique_chains(v)&reference) for k,v in assembled['stages'].items()}
    assert hits['selected_unique']<=hits['owner_unique']<=hits['free_unique']
    return {'reference':len(reference),'exact_reference_hits':hits,
            'not_in_actual_free_candidates':len(reference)-hits['free_unique'],
            'lost_at_owner':hits['free_unique']-hits['owner_unique'],
            'lost_at_conflict_or_nonpositive_gain':hits['owner_unique']-hits['selected_unique'],
            'candidate_counts':assembled['counts']}

def aggregate(per_species):
    items=list(per_species.values())
    result={}
    for field in ('exact_chain','CDS_base_unstranded_union','historical_parent_exact_chain'):
        result[field]=E.S.prf(*(sum(r[field][k] for r in items) for k in ('tp','predicted','reference')))
    result['macro_exact_chain']={k:sum(r['exact_chain'][k] for r in items)/len(items) for k in ('precision','recall','f1')}
    bg={k:sum(r['background'][k] for r in items) for k in ('bp','predicted_span_bp','predicted_CDS_bp','wholly_background_unique_chains')}
    bg['span_FPR']=bg['predicted_span_bp']/bg['bp']
    bg['wholly_background_chains_per_Mb']=bg['wholly_background_unique_chains']*1e6/bg['bp']
    result['background']=bg
    return result

def validated_scope_records(path,scope):
    count=0
    with path.open() as f:
        for i,line in enumerate(f):
            record=json.loads(line)
            expected=scope['species']+'|'+scope['seqid']+'|'+str(scope['windows'][i][0])+'|'+scope['strand']
            if record['window_index']!=i or record['window_id']!=expected or record.get('reference_used') is not False:
                raise ValueError('Window record identity or inference provenance mismatch')
            count+=1
            yield record
    if count!=len(scope['windows']): raise ValueError('Partial inference scope')

def main(pilot,out):
    pilot=pilot.resolve();out=out.resolve();out.mkdir(parents=True,exist_ok=True)
    if (out/'summary.json').exists(): raise FileExistsError('Completed analysis exists')
    fitting={arm:json.loads((pilot/arm/'summary.json').read_text()) for arm in ('C0','B1')}
    inference={arm:json.loads((pilot/('infer_'+arm)/'summary.json').read_text()) for arm in ('C0','B1')}
    expected_scopes={(s,q,strand) for s,q,_,_ in E.S.SCOPE for strand in ('+','-')}
    for arm,r in inference.items():
        if (fitting[arm]['optimizer_steps']!=4608 or not r['inference_complete'] or r['windows']!=8690
            or r['reference_used'] is not False or r['optimizer_steps']!=4608):
            raise ValueError('Both full final-checkpoint arms are required for comparison')
        if len(r['scopes'])!=4 or {(x['species'],x['seqid'],x['strand']) for x in r['scopes']}!=expected_scopes:
            raise ValueError('Inference scope mismatch')
    references={};isoforms={};parents={};annotations={};sequences={};lengths={}
    for species,seqid,length,old_count in E.S.SCOPE:
        path=E.ROOT/'data/m1_screen'/species
        sequences[species]=M.read_allowed(path/'genome.fa',{seqid:length})
        annotations[species]=E.S.E.parse_annotation(path/'reference.gff3',{seqid:length})
        coding=E.S.E.parse_annotation(path/'reference.gff3',{seqid:length},protein_coding_only=True)
        with (path/'reference.gff3').open() as f: _,_,raw,_=E.P.raw_features(f,seqid)
        rebuilt,isoforms[species]=E.select_reference(coding,raw,sequences[species])
        references[species]=rows(E.ROOT/'outputs/M28-COMMON-RULER-R3/result'/(species+'.references.jsonl'))
        if {t['id']:E.S.key(t) for t in rebuilt}!={t['id']:E.S.key(t) for t in references[species]}:
            raise ValueError('Frozen primary reference changed')
        parents[species]=E.S.E.primary_transcripts(coding)
        assert len(parents[species])==old_count
        lengths[seqid]=length
    report={'experiment':pilot.name,'status':'complete_bounded_DEV_comparison_not_independent_generalization',
            'reference_primary_chains':7728,'historical_parent_chains':6450,
            'test_setaria_used':False,'methods':{},'fit':fitting,'inference_cost':{},
            'cached_baseline_source':'reports/M28-METHOD-RESTART/common_ruler.json',
            'candidate_semantics':'C0 native decoded chains; B1 actual pruned free proposals, no reference injection',
            'scope_limit':'B1 versus C0 is a whole-system comparison, not isolated nonadditive-GRU causation'}
    for arm,r in inference.items():
        report['inference_cost'][arm]={k:v for k,v in r.items() if k!='scopes'}
        report['inference_cost'][arm]['scope_timings']=[{k:v for k,v in s.items() if k!='windows'} for s in r['scopes']]
    for method,arm,score_field in (('C0','C0','score'),('B1','B1','score'),('B1_additive','B1','additive_only_score')):
        predictions=defaultdict(list);stages={}
        for scope in inference[arm]['scopes']:
            species,seqid,strand=scope['species'],scope['seqid'],scope['strand']
            expected=[[a,min(a+M.WIDTH,lengths[seqid])] for a in M.starts(lengths[seqid])]
            if scope['windows']!=expected: raise ValueError('Inference did not cover full natural window grid')
            source=pilot/('infer_'+arm)/scope['file']
            assembled=assemble_predictions(validated_scope_records(source,scope),scope['windows'],strand,score_field)
            for i,row in enumerate(assembled['predictions']):
                identifier='M28_'+method+'_'+seqid+'_'+('p' if strand=='+' else 'm')+'_'+str(i)
                predictions[species].append(as_transcript(row['chain'],seqid,strand,identifier,row['score']))
            selected_refs=[t for t in references[species] if t['strand']==strand]
            audit=stage_audit(assembled,selected_refs)
            audit.update(species=species,seqid=seqid,strand=strand)
            stages[species+'|'+strand]=audit
            # Compact per-scope candidate counts/hits are sufficient for the result table;
            # full free lists remain in immutable inference JSONL files.
            del assembled
        gff=out/(method+'.gff3')
        txs=[t for values in predictions.values() for t in values]
        write_gff3(gff,txs,lengths,'M28_'+method)
        parsed=E.S.E.primary_transcripts(E.S.E.parse_annotation(gff,lengths),complete_only=False)
        if {E.S.key(t) for t in parsed}!={E.S.key(t) for t in txs} or len(parsed)!=len(txs):
            raise ValueError('GFF3 export lost or altered predictions')
        per_species={}
        for species,seqid,length,_ in E.S.SCOPE:
            pred=predictions[species]
            measured=E.score(pred,references[species],isoforms[species],length,annotations[species])
            pk,hk={E.S.key(t) for t in pred},{E.S.key(t) for t in parents[species]}
            measured['historical_parent_exact_chain']=E.S.prf(len(pk&hk),len(pk),len(hk))
            failures=Counter();valid=0
            for tx in pred:
                problems,_=M.chain_checks(tx,sequences[species][seqid])
                failures.update(problems);valid+=not problems
            measured['CDS_structure_diagnostic']={'fully_supported':valid,'predictions':len(pred),
                                                  'failure_counts_nonexclusive':dict(failures)}
            per_species[species]=measured
        pooled=aggregate(per_species)
        pooled['actual_free_primary_coverage']=sum(x['exact_reference_hits']['free_unique'] for x in stages.values())/7728
        pooled['owner_primary_coverage']=sum(x['exact_reference_hits']['owner_unique'] for x in stages.values())/7728
        assert sum(x['exact_reference_hits']['selected_unique'] for x in stages.values())==pooled['exact_chain']['tp']
        report['methods'][method]={'per_species':per_species,'pooled':pooled,'stages':stages,'GFF3':str(gff.relative_to(E.ROOT))}
    baseline=json.loads((E.ROOT/'reports/M28-METHOD-RESTART/common_ruler.json').read_text())
    compare={}
    for method,values in baseline['pooled'].items():
        by_species={s:r['methods'][method]['exact_chain'] for s,r in baseline['species'].items()}
        compare[method]={'per_species':by_species,'pooled':values['exact_chain'],
                         'macro':{k:sum(x[k] for x in by_species.values())/len(by_species) for k in ('precision','recall','f1')},
                         'source':'unchanged_historical_cache'}
    for method,values in report['methods'].items():
        compare[method]={'per_species':{s:r['exact_chain'] for s,r in values['per_species'].items()},
                         'pooled':values['pooled']['exact_chain'],'macro':values['pooled']['macro_exact_chain'],
                         'source':'completed_last_checkpoint' if method!='B1_additive' else 'same_B1_candidates_additive_readout',
                         'checkpoint':fitting['C0' if method=='C0' else 'B1']['final_checkpoint']}
    report['exact_chain_comparison']=compare
    with (out/'summary.json').open('x') as f: json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({'experiment':report['experiment'],'comparison':compare},indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--pilot-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();main(a.pilot_dir,a.output_dir)
