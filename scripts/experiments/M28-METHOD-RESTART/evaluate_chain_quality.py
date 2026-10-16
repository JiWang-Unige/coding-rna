#!/usr/bin/env python3
"""A6 full DEV ruler after completed raw-FASTA inference, never a test claim."""
from collections import Counter, defaultdict
import json
import run_chain_quality as Q
import evaluate_query as A
from evaluate_query_collaboration import continuation_gate
from evaluate_pilot import aggregate, stage_audit, validated_scope_records
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript

OUT = Q.ROOT/'outputs/M28-A6-CHAIN-QUALITY-EVAL'


def validate_complete(fitted, decoded):
    Q.validate_head(fitted)
    expected = {(s,q,t):[[a,min(a+Q.I.M.WIDTH,n)] for a in Q.I.M.starts(n)]
                for s,qs in Q.I.DEV_LENGTHS.items() for q,n in qs.items() for t in ('+','-')}
    actual = {(s['species'],s['seqid'],s['strand']):s['windows'] for s in decoded['scopes']}
    if (not decoded['completed'] or not decoded['inference_complete'] or decoded['arm']!='A6' or
        decoded['configuration']!=Q.CONFIG or decoded['protocol']!=Q.PROTOCOL or
        decoded['checkpoint']!=fitted['final_checkpoint'] or decoded['optimizer_steps']!=4608 or
        decoded['feature_revision']!=Q.I.F.REVISION or decoded['feature_spec']!=Q.R.FEATURE_SPEC or
        decoded['input_mode']!='fixed_DEV_raw_FASTA' or decoded['reference_used'] or
        decoded['labels_or_feature_cache_used'] or decoded['optimizer_steps_this_run']!=0 or
        decoded['DEV_used_for_selection'] or decoded['test_setaria_used'] or
        decoded['windows']!=8690 or len(decoded['scopes'])!=4 or actual!=expected):
        raise ValueError('Complete fixed raw-FASTA A6 inference required before reference access')


def records(path, scope):
    for row in validated_scope_records(path,scope):
        if (row['score_kind']!='A6_complete_chain_logit_minus_window_null' or
            row['quality_logit_threshold']!=0. or row['reference_used'] or
            len(row['chains'])!=len(row['score']) or len(row['chains'])!=len(row['decoder_sources']) or
            any(not sources or not set(sources)<={'greedy','beam4'} for sources in row['decoder_sources']) or
            len({tuple(map(tuple,c)) for c in row['chains']})!=len(row['chains'])):
            raise ValueError('A6 score/pool/source boundary differs')
        yield row


def score_species(predictions, refs, isoforms, length, annotation, parents):
    measured = A.E.score(predictions, refs, isoforms, length, annotation)
    predicted, historical = {A.E.S.key(t) for t in predictions}, {A.E.S.key(t) for t in parents}
    measured['historical_parent_exact_chain'] = A.E.S.prf(
        len(predicted & historical), len(predicted), len(historical))
    return measured


def main():
    fitted = json.loads((Q.FIT/'summary.json').read_text())
    decoded = json.loads((Q.INFER/'summary.json').read_text())
    validate_complete(fitted,decoded)
    if any((p/'STATUS').read_text().strip()!='COMPLETED' for p in (Q.FIT,Q.INFER)):
        raise ValueError('Terminal fit/inference required')
    context = A.reference_context()
    refs,isoforms,parents,annotations,sequences,lengths = context
    predictions,stages = defaultdict(list),{}
    for scope in decoded['scopes']:
        species,seqid,strand = scope['species'],scope['seqid'],scope['strand']
        assembled = assemble_predictions(records(Q.INFER/scope['file'],scope),scope['windows'],strand)
        if assembled['counts']!=scope['assembly_counts']:
            raise ValueError('Saved inference and reconstructed A6 assembly differ')
        for index,row in enumerate(assembled['predictions']):
            predictions[species].append(as_transcript(row['chain'],seqid,strand,
                'M28_A6_'+str(scope['ordinal'])+'_'+strand+'_'+str(index),row['score']))
        audit = stage_audit(assembled,[r for r in refs[species] if r['strand']==strand])
        audit.update(species=species,seqid=seqid,strand=strand)
        stages[species+'|'+strand] = audit
    txs = [t for values in predictions.values() for t in values]
    saved = A.E.S.E.primary_transcripts(A.E.S.E.parse_annotation(Q.INFER/'A6.gff3',lengths),complete_only=False)
    A.validate_gff_roundtrip(saved,txs)
    if len(txs)!=decoded['primary_CDS_chain_predictions']:
        raise ValueError('A6 GFF prediction count changed')
    per_species = {}
    for species,seqid,length,_ in A.E.S.SCOPE:
        measured = score_species(predictions[species],refs[species],isoforms[species],
                                 length,annotations[species],parents[species])
        measured['prediction_nonexact_primary'] = measured.pop('reference_nonexact_primary')
        failures = Counter()
        for tx in predictions[species]: failures.update(Q.I.M.chain_checks(tx,sequences[species][seqid])[0])
        measured['CDS_structure_failure_counts_nonexclusive'] = dict(failures)
        per_species[species] = measured
    pooled = aggregate(per_species)
    if sum(r['exact_reference_hits']['selected_unique'] for r in stages.values())!=pooled['exact_chain']['tp']:
        raise ValueError('A6 stage hits and exact-chain TP disagree')
    previous = json.loads((Q.ROOT/'outputs/M28-A-EVALUATION/summary.json').read_text())
    controls = json.loads((Q.ROOT/'reports/M28-METHOD-RESTART/r5_result.json').read_text())
    if any(s['reference_primary_chains']!=7728 or s['test_setaria_used'] for s in (previous,controls)):
        raise ValueError('Historical gate denominator differs')
    report = {'experiment':'M28-A6-CHAIN-QUALITY','completed':True,
        'status':'complete_one_seed_bounded_DEV_not_independent_generalization',
        'protocol':Q.PROTOCOL,'reference_primary_chains':7728,'test_setaria_used':False,
        'fit':fitted,'inference_receipt':{k:v for k,v in decoded.items() if k!='scopes'},
        'methods':{'A6':{'per_species':per_species,'pooled':pooled,'stages':stages}},
        'predeclared_pilot_gate':continuation_gate(per_species,previous['methods']['A']['per_species'],
                                                  controls['methods']['C0']['per_species']),
        'limitations':'extra head exposure and expanded pool; not a matched causal contrast or paper/generalization PASS'}
    Q.R.write_summary(OUT,report)


if __name__=='__main__': main()
