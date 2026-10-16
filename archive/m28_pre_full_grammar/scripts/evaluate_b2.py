#!/usr/bin/env python3
"""Assemble and score complete paired B2 DEV against the unchanged frozen ruler."""
import argparse
import json
from collections import Counter,defaultdict
from pathlib import Path
import common_ruler as E
import manifest as M
from evaluate_pilot import rows,stage_audit,aggregate,validated_scope_records
from fit_b2 import CONTRACT,check_output_cap
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript,write_gff3

def main(fit,inference,out):
    fit=fit.resolve();inference=inference.resolve();out=out.resolve();out.mkdir(parents=True,exist_ok=True)
    fitted=json.loads((fit/'summary.json').read_text())
    decoded=json.loads((inference/'summary.json').read_text())
    if (not fitted['completed'] or fitted['optimizer_steps']!=4608 or not decoded['inference_complete'] or
        decoded['windows']!=8690 or decoded['optimizer_steps']!=4608 or decoded['reference_used'] is not False or
        decoded['checkpoint']!=fitted['final_checkpoint'] or decoded['control_replay']['windows']!=8690 or
        not decoded['control_replay']['discrete_exact'] or decoded['control_replay']['proposal_max_abs']>1e-5 or
        not decoded['paired_common_gain']['exact_equality'] or not decoded['paired_common_gain']['single_union_readout']):
        raise ValueError('Complete matched-checkpoint paired inference and replay are required')
    expected_scopes={(s,q,strand) for s,q,_,_ in E.S.SCOPE for strand in ('+','-')}
    if set(decoded['arms'])!={'early','late'}: raise ValueError('Expected paired early/late arms')
    for arm,r in decoded['arms'].items():
        if len(r['scopes'])!=4 or {(x['species'],x['seqid'],x['strand']) for x in r['scopes']}!=expected_scopes:
            raise ValueError('Paired inference scope mismatch')
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
        if len(parents[species])!=old_count: raise ValueError('Historical parent denominator changed')
        lengths[seqid]=length
    if sum(map(len,references.values()))!=7728: raise ValueError('Primary denominator changed')
    report={'experiment':'M28-B2','status':'complete_one_seed_DEV_mechanism_not_independent_generalization',
            'contract':CONTRACT,'reference_primary_chains':7728,'historical_parent_chains':6450,
            'test_setaria_used':False,'methods':{},'fit':fitted,
            'paired_inference_cost':{k:v for k,v in decoded.items() if k!='arms'},
            'candidate_semantics':'actual pruned free pools; no reference injection in inference',
            'scope_limit':'same trained head; early versus late scoring position; not independent fits or SOTA reruns'}
    for arm,r in decoded['arms'].items():
        method='B2_'+arm;predictions=defaultdict(list);stages={}
        for scope in r['scopes']:
            species,seqid,strand=scope['species'],scope['seqid'],scope['strand']
            expected=[[a,min(a+M.WIDTH,lengths[seqid])] for a in M.starts(lengths[seqid])]
            if scope['windows']!=expected: raise ValueError('Incomplete natural DEV grid')
            assembled=assemble_predictions(validated_scope_records(inference/scope['file'],scope),scope['windows'],strand)
            for i,row in enumerate(assembled['predictions']):
                ident=method+'_'+seqid+'_'+('p' if strand=='+' else 'm')+'_'+str(i)
                predictions[species].append(as_transcript(row['chain'],seqid,strand,ident,row['score']))
            audit=stage_audit(assembled,[t for t in references[species] if t['strand']==strand])
            audit.update(species=species,seqid=seqid,strand=strand);stages[species+'|'+strand]=audit
            del assembled
        gff=out/(method+'.gff3');txs=[t for values in predictions.values() for t in values]
        write_gff3(gff,txs,lengths,method)
        parsed=E.S.E.primary_transcripts(E.S.E.parse_annotation(gff,lengths),complete_only=False)
        if {E.S.key(t) for t in parsed}!={E.S.key(t) for t in txs} or len(parsed)!=len(txs):
            raise ValueError('B2 GFF3 export changed or lost predictions')
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
        if sum(x['exact_reference_hits']['selected_unique'] for x in stages.values())!=pooled['exact_chain']['tp']:
            raise ValueError('Selected stage and final exact-chain denominator disagree')
        report['methods'][method]={'per_species':per_species,'pooled':pooled,'stages':stages,
                                   'GFF3':str(gff.relative_to(E.ROOT))}
        check_output_cap()
    baseline=json.loads((E.ROOT/'reports/M28-METHOD-RESTART/common_ruler.json').read_text())
    compare={}
    for method,values in baseline['pooled'].items():
        by_species={s:r['methods'][method]['exact_chain'] for s,r in baseline['species'].items()}
        compare[method]={'per_species':by_species,'pooled':values['exact_chain'],
                         'macro':{k:sum(x[k] for x in by_species.values())/len(by_species) for k in ('precision','recall','f1')},
                         'source':'unchanged_historical_cache_not_new_run'}
    r5=json.loads((E.ROOT/'reports/M28-METHOD-RESTART/r5_result.json').read_text())
    for method in ('C0','B1'):
        values=r5['methods'][method]
        compare[method]={'per_species':{s:r['exact_chain'] for s,r in values['per_species'].items()},
                         'pooled':values['pooled']['exact_chain'],'macro':values['pooled']['macro_exact_chain'],
                         'source':'completed_R5_cache_not_new_run'}
    for method,values in report['methods'].items():
        compare[method]={'per_species':{s:r['exact_chain'] for s,r in values['per_species'].items()},
                         'pooled':values['pooled']['exact_chain'],'macro':values['pooled']['macro_exact_chain'],
                         'source':'completed_same_new_prefix_checkpoint','checkpoint':fitted['final_checkpoint']}
    report['exact_chain_comparison']=compare
    early,late=report['methods']['B2_early'],report['methods']['B2_late']
    gate={}
    for species in references:
        counts=lambda method,stage:sum(r['exact_reference_hits'][stage] for r in method['stages'].values()
                                       if r['species']==species)
        a,b=early['per_species'][species],late['per_species'][species]
        gate[species]={'free_complete_chains_improved':counts(early,'free_unique')>counts(late,'free_unique'),
                       'owner_complete_chains_improved':counts(early,'owner_unique')>counts(late,'owner_unique'),
                       'exact_chain_F1_improved':a['exact_chain']['f1']>b['exact_chain']['f1'],
                       'background_span_not_increased':a['background']['span_FPR']<=b['background']['span_FPR']}
    report['predeclared_mechanism_gate']={'per_species':gate,'passed':all(all(v.values()) for v in gate.values()),
                                        'not_paper_success_or_independent_generalization':True}
    with (out/'summary.json').open('x') as f: json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    check_output_cap();print(json.dumps({'comparison':compare,'gate':report['predeclared_mechanism_gate']},indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--fit-dir',type=Path,required=True)
    p.add_argument('--inference-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();main(a.fit_dir,a.inference_dir,a.output_dir)
