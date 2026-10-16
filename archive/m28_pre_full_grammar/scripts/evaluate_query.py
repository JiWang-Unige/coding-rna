#!/usr/bin/env python3
"""A completed-DEV evaluation against the unchanged primary CDS ruler.

Only explicit CLI execution reads DEV references. Synthetic helper tests do
not. The predeclared gate is a pilot continuation gate, never paper success.
"""
import json
from collections import Counter, defaultdict
import common_ruler as E
from evaluate_pilot import rows, aggregate, stage_audit, validated_scope_records
import run_query as Q
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript, write_gff3

REF_COUNTS = {'arabidopsis_thaliana':5437,'oryza_sativa':2291}


def validate_complete(fitted, decoded):
    Q.check_last_identity(fitted); Q.check_last_identity(decoded)
    if (not decoded['inference_complete'] or decoded['windows'] != 8690 or
        decoded['reference_used'] is not False or decoded['test_setaria_used'] is not False or
        decoded['checkpoint'] != fitted['final_checkpoint']):
        raise ValueError('Complete matching A DEV inference is required')
    expected = {(s,q,strand):[[a,min(a+Q.M.WIDTH,length)] for a in Q.M.starts(length)]
                for s,q,length,_ in E.S.SCOPE for strand in ('+','-')}
    actual = {(s['species'],s['seqid'],s['strand']):s['windows'] for s in decoded['scopes']}
    if len(decoded['scopes']) != 4 or actual != expected:
        raise ValueError('A full natural DEV scope grid mismatch')


def query_scope_records(path, scope):
    for row in validated_scope_records(path,scope):
        if (row.get('score_kind') != 'query_existence_logit_minus_threshold' or
            row['existence_logit_threshold'] != Q.CONFIG['existence_logit_threshold'] or
            row['counts']['queries'] != Q.CONFIG['num_queries'] or
            len(row['raw_queries']) != Q.CONFIG['num_queries']):
            raise ValueError('A query operating point/count/provenance mismatch')
        yield row


def validate_gff_roundtrip(parsed, expected):
    if (len(parsed) != len(expected) or
        {E.S.key(t) for t in parsed} != {E.S.key(t) for t in expected} or
        {t['id']:t['CDS'] for t in parsed} != {t['id']:t['CDS'] for t in expected}):
        raise ValueError('A GFF3 export lost or changed identity/coordinates/phase')


def continuation_gate(current, control):
    if set(current) != set(REF_COUNTS) or set(control) != set(REF_COUNTS):
        raise ValueError('Two matching frozen species are required')
    gate = {}
    for species,count in REF_COUNTS.items():
        a,b = current[species],control[species]
        x,y = a['exact_chain'],b['exact_chain']
        ab,bb = a['background'],b['background']
        if x['reference'] != count or y['reference'] != count or ab['bp'] != bb['bp'] or ab['bp'] <= 0:
            raise ValueError('Reference/background denominator mismatch')
        # Exact integer counts, not rounded F1 or cached rounded percentages.
        improved = x['tp']*(y['predicted']+y['reference']) > y['tp']*(x['predicted']+x['reference'])
        gate[species] = {'exact_chain_F1_improved':improved,
                         'background_span_not_increased':ab['predicted_span_bp'] <= bb['predicted_span_bp']}
    return {'per_species':gate,'passed':all(all(v.values()) for v in gate.values()),
            'not_paper_success_or_independent_generalization':True}


def reference_context():
    refs,isoforms,parents,annotations,sequences,lengths = {},{},{},{},{},{}
    for species,seqid,length,old_count in E.S.SCOPE:
        path = E.ROOT/'data/m1_screen'/species
        sequences[species] = Q.M.read_allowed(path/'genome.fa',{seqid:length})
        annotations[species] = E.S.E.parse_annotation(path/'reference.gff3',{seqid:length})
        coding = E.S.E.parse_annotation(path/'reference.gff3',{seqid:length},protein_coding_only=True)
        with (path/'reference.gff3').open() as handle: _,_,raw,_ = E.P.raw_features(handle,seqid)
        rebuilt,isoforms[species] = E.select_reference(coding,raw,sequences[species])
        refs[species] = rows(E.ROOT/'outputs/M28-COMMON-RULER-R3/result'/(species+'.references.jsonl'))
        if ({t['id']:E.S.key(t) for t in rebuilt} != {t['id']:E.S.key(t) for t in refs[species]} or
            len({E.S.key(t) for t in refs[species]}) != REF_COUNTS[species]):
            raise ValueError('Frozen primary reference changed; do not filter unsupported references')
        parents[species] = E.S.E.primary_transcripts(coding)
        if len(parents[species]) != old_count: raise ValueError('Historical parent denominator changed')
        lengths[seqid] = length
    return refs,isoforms,parents,annotations,sequences,lengths


def main():
    fit,inference,out = (Q.STAGE_DIRS[k] for k in ('fit','infer','evaluate'))
    fitted = json.loads((fit/'summary.json').read_text())
    decoded = json.loads((inference/'summary.json').read_text())
    validate_complete(fitted,decoded)
    out.mkdir(parents=True,exist_ok=True)
    refs,isoforms,parents,annotations,sequences,lengths = reference_context()
    predictions,stages = defaultdict(list),{}
    for scope in decoded['scopes']:
        species,seqid,strand = scope['species'],scope['seqid'],scope['strand']
        assembled = assemble_predictions(query_scope_records(inference/scope['file'],scope),scope['windows'],strand)
        for i,row in enumerate(assembled['predictions']):
            ident = 'A_'+seqid+'_'+('p' if strand == '+' else 'm')+'_'+str(i)
            predictions[species].append(as_transcript(row['chain'],seqid,strand,ident,row['score']))
        audit = stage_audit(assembled,[t for t in refs[species] if t['strand'] == strand])
        audit.update(species=species,seqid=seqid,strand=strand); stages[species+'|'+strand] = audit
        del assembled
    transcripts = [t for values in predictions.values() for t in values]
    gff = out/'A.gff3'; write_gff3(gff,transcripts,lengths,'M28_A')
    parsed = E.S.E.primary_transcripts(E.S.E.parse_annotation(gff,lengths),complete_only=False)
    validate_gff_roundtrip(parsed,transcripts)
    per_species = {}
    for species,seqid,length,_ in E.S.SCOPE:
        pred = predictions[species]
        measured = E.score(pred,refs[species],isoforms[species],length,annotations[species])
        pk,hk = {E.S.key(t) for t in pred},{E.S.key(t) for t in parents[species]}
        measured['historical_parent_exact_chain'] = E.S.prf(len(pk&hk),len(pk),len(hk))
        failures = Counter(); valid = 0
        for tx in pred:
            problems,_ = Q.M.chain_checks(tx,sequences[species][seqid])
            failures.update(problems); valid += not problems
        measured['CDS_structure_diagnostic'] = {'fully_supported':valid,'predictions':len(pred),
                                               'failure_counts_nonexclusive':dict(failures)}
        per_species[species] = measured
    pooled = aggregate(per_species)
    pooled['actual_free_primary_coverage'] = sum(r['exact_reference_hits']['free_unique'] for r in stages.values())/7728
    pooled['owner_primary_coverage'] = sum(r['exact_reference_hits']['owner_unique'] for r in stages.values())/7728
    if sum(r['exact_reference_hits']['selected_unique'] for r in stages.values()) != pooled['exact_chain']['tp']:
        raise ValueError('Selected stage and final exact-chain counts disagree')
    historical = json.loads((E.ROOT/'reports/M28-METHOD-RESTART/common_ruler.json').read_text())
    comparison = {}
    for method,values in historical['pooled'].items():
        by_species = {s:r['methods'][method]['exact_chain'] for s,r in historical['species'].items()}
        comparison[method] = {'per_species':by_species,'pooled':values['exact_chain'],
                              'macro':{k:sum(x[k] for x in by_species.values())/2 for k in ('precision','recall','f1')},
                              'source':'unchanged_historical_cache_not_native_rerun'}
    controls = json.loads((E.ROOT/'reports/M28-METHOD-RESTART/r5_result.json').read_text())['methods']
    for method in ('C0','B1'):
        values = controls[method]
        comparison[method] = {'per_species':{s:r['exact_chain'] for s,r in values['per_species'].items()},
                              'pooled':values['pooled']['exact_chain'],'macro':values['pooled']['macro_exact_chain'],
                              'source':'completed_R5_same_exposure_cache_not_new_fit'}
    comparison['A'] = {'per_species':{s:r['exact_chain'] for s,r in per_species.items()},
                       'pooled':pooled['exact_chain'],'macro':pooled['macro_exact_chain'],
                       'source':'completed_last_A_checkpoint','checkpoint':fitted['final_checkpoint']}
    report = {'experiment':'M28-A','status':'complete_one_seed_bounded_DEV_not_independent_generalization',
              'reference_primary_chains':7728,'historical_parent_chains':6450,'test_setaria_used':False,
              'proposal':Q.PROPOSAL,'fit':fitted,'inference_cost':{k:v for k,v in decoded.items() if k != 'scopes'},
              'methods':{'A':{'per_species':per_species,'pooled':pooled,'stages':stages,
                              'query_counts':[s['query_counts'] for s in decoded['scopes']],
                              'GFF3':str(gff.relative_to(E.ROOT))}},
              'exact_chain_comparison':comparison,
              'predeclared_pilot_gate':continuation_gate(per_species,controls['C0']['per_species']),
              'scope_limit':'whole-system A versus C0/B1; query/auxiliary/shared-gradient/negative-policy changes not isolated'}
    report['artifact_bytes_before_summary'] = Q.artifact_bytes()
    with (out/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    Q.artifact_bytes()
    print(json.dumps({'comparison':comparison,'gate':report['predeclared_pilot_gate']},allow_nan=False),flush=True)


if __name__ == '__main__':
    main()
