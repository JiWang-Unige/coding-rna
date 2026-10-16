#!/usr/bin/env python3
"""Paired full-TRAIN last-checkpoint DEV evaluation; never an independent-test claim."""
import json
from collections import Counter, defaultdict
import evaluate_query as A
import infer_full_train as I
from evaluate_pilot import aggregate, stage_audit, validated_scope_records
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript, write_gff3

E, R, M = A.E, I.R, I.M
OUT = I.ROOT/'outputs/M28-FULL-TRAIN-EVAL'


def validate_complete(fitted, decoded, arm):
    I.validate_fit(fitted,arm)
    expected = {(s,q,t):[[a,min(a+M.WIDTH,length)] for a in M.starts(length)]
                for s,qs in I.DEV_LENGTHS.items() for q,length in qs.items() for t in ('+','-')}
    actual = {(s['species'],s['seqid'],s['strand']):s['windows'] for s in decoded['scopes']}
    if (decoded['inference_complete'] is not True or decoded['arm'] != arm or
        decoded['configuration'] != R.CONFIGS[arm] or decoded['optimizer_steps'] != R.STEPS or
        decoded['checkpoint'] != fitted['final_checkpoint'] or
        decoded['encoder_revision'] != I.F.REVISION or decoded['feature_spec'] != R.FEATURE_SPEC or
        decoded['input_mode'] != 'fixed_DEV_raw_FASTA' or
        decoded['reference_used'] is not False or decoded['labels_or_feature_cache_used'] is not False or
        decoded['optimizer_steps_this_run'] != 0 or decoded['test_setaria_used'] is not False or
        decoded['DEV_used_for_selection'] is not False or decoded['windows'] != 8690 or
        len(decoded['scopes']) != 4 or actual != expected or
        decoded['lengths'] != {q:n for qs in I.DEV_LENGTHS.values() for q,n in qs.items()}):
        raise ValueError('Both complete matching raw-FASTA final DEV arms required')


def scope_records(path, scope, arm, counts):
    for row in validated_scope_records(path,scope):
        if len(row['chains']) != len(row['score']):
            raise ValueError('Every candidate must have one score')
        if arm=='Q':
            c = R.CONFIGS['Q']
            if (row['score_kind'] != 'query_existence_logit_minus_threshold' or
                row['existence_logit_threshold'] != c['existence_logit_threshold'] or
                row['counts']['queries'] != c['num_queries'] or len(row['raw_queries']) != c['num_queries']):
                raise ValueError('Fixed Q query count/threshold/provenance mismatch')
        elif row['score_kind'] != 'native_HMM_chain_removal_gain':
            raise ValueError('Unchanged C0 native gain required')
        chains = [tuple(map(tuple,x)) for x in row['chains']]
        positive = [x for x,score in zip(chains,row['score']) if score>0]
        counts.update({'windows':1,'completed_candidates':len(chains),
                       'positive_completed_candidates':len(positive),
                       'unique_positive_chains_within_window':len(set(positive)),
                       'within_window_duplicate_positive_candidates':len(positive)-len(set(positive))})
        yield row


def measure(arm, decoded, context):
    refs,isoforms,parents,annotations,sequences,lengths = context
    predictions,stages,diversity = defaultdict(list),{},{}
    for scope in decoded['scopes']:
        species,seqid,strand = scope['species'],scope['seqid'],scope['strand']
        counts = Counter()
        assembled = assemble_predictions(scope_records(I.INFER[arm]/scope['file'],scope,arm,counts),
                                         scope['windows'],strand)
        if arm=='Q':
            q = scope['query_counts']
            if (counts['completed_candidates'] != q['completed_queries'] or
                counts['positive_completed_candidates'] != q['completed_queries']-q['completed_nonpositive_gain']):
                raise ValueError('Q diversity does not agree with saved inference')
        if assembled['counts'] != scope['assembly_counts']:
            raise ValueError('Reassembled candidates differ from inference assembly')
        for i,row in enumerate(assembled['predictions']):
            ident = I.prediction_id(arm,scope['ordinal'],strand,i)
            predictions[species].append(as_transcript(row['chain'],seqid,strand,ident,row['score']))
        audit = stage_audit(assembled,[t for t in refs[species] if t['strand']==strand])
        audit.update(species=species,seqid=seqid,strand=strand)
        key = species+'|'+strand; stages[key] = audit; diversity[key] = dict(counts)
    transcripts = [t for values in predictions.values() for t in values]
    saved = E.S.E.primary_transcripts(E.S.E.parse_annotation(I.INFER[arm]/(arm+'.gff3'),lengths),complete_only=False)
    A.validate_gff_roundtrip(saved,transcripts)
    if len(transcripts) != decoded['primary_CDS_chain_predictions']:
        raise ValueError('Final prediction count differs from inference')
    gff = OUT/(arm+'.gff3'); write_gff3(gff,transcripts,lengths,'M28_FULL_'+arm)
    parsed = E.S.E.primary_transcripts(E.S.E.parse_annotation(gff,lengths),complete_only=False)
    A.validate_gff_roundtrip(parsed,transcripts)
    per_species = {}
    for species,seqid,length,_ in E.S.SCOPE:
        pred = predictions[species]
        measured = E.score(pred,refs[species],isoforms[species],length,annotations[species])
        # The historical field contains nonexact PREDICTIONS, not missed references.
        measured['prediction_nonexact_primary'] = measured.pop('reference_nonexact_primary')
        pk,hk = {E.S.key(t) for t in pred},{E.S.key(t) for t in parents[species]}
        measured['historical_parent_exact_chain'] = E.S.prf(len(pk&hk),len(pk),len(hk))
        failures = Counter(); valid = 0
        for tx in pred:
            problems,_ = M.chain_checks(tx,sequences[species][seqid])
            failures.update(problems); valid += not problems
        measured['CDS_structure_diagnostic'] = {'fully_supported':valid,'predictions':len(pred),
                                               'failure_counts_nonexclusive':dict(failures)}
        per_species[species] = measured
    pooled = aggregate(per_species)
    pooled['actual_free_primary_coverage'] = sum(x['exact_reference_hits']['free_unique'] for x in stages.values())/7728
    pooled['owner_primary_coverage'] = sum(x['exact_reference_hits']['owner_unique'] for x in stages.values())/7728
    if sum(x['exact_reference_hits']['selected_unique'] for x in stages.values()) != pooled['exact_chain']['tp']:
        raise ValueError('Selected candidate-stage TP disagrees with exact-chain score')
    return {'per_species':per_species,'pooled':pooled,'stages':stages,
            'candidate_diversity':diversity,'GFF3':str(gff.relative_to(I.ROOT))}


def historical_exposure_comparison():
    result = {}
    sources = [('R5_C0','C0',I.ROOT/'reports/M28-METHOD-RESTART/r5_result.json')]
    sources += [(arm,arm,I.ROOT/'outputs'/('M28-'+arm+'-EVALUATION')/'summary.json') for arm in ('A','A2','A3')]
    for name,method,path in sources:
        if not path.exists():
            result[name] = {'source':str(path.relative_to(I.ROOT)),'status':'not_available_not_a_result'}
            continue
        s = json.loads(path.read_text())
        if s['reference_primary_chains'] != 7728 or s['test_setaria_used'] is not False:
            raise ValueError('Historical result does not use the same DEV ruler')
        values = s['methods'][method]
        if any(values['per_species'][species]['exact_chain']['reference'] != n for species,n in A.REF_COUNTS.items()):
            raise ValueError('Historical per-species reference denominator differs')
        result[name] = {'per_species':{species:r['exact_chain'] for species,r in values['per_species'].items()},
                        'pooled':values['pooled']['exact_chain'],'source':str(path.relative_to(I.ROOT)),
                        'status':'completed_historical_different_exposure_not_matched_control',
                        'training_draws_per_pass':1536,'unique_TRAIN_windows':1491,'optimizer_steps':4608}
    return result


def main():
    fits = {a:json.loads((I.FIT[a]/'summary.json').read_text()) for a in R.CONFIGS}
    inferences = {a:json.loads((I.INFER[a]/'summary.json').read_text()) for a in R.CONFIGS}
    for arm in R.CONFIGS:
        validate_complete(fits[arm],inferences[arm],arm)
        if any((p/'STATUS').read_text().strip()!='COMPLETED' for p in (I.FIT[arm],I.INFER[arm])):
            raise ValueError('Terminal fits and inferences required before reference access')
    OUT.mkdir(parents=True,exist_ok=True)
    if (OUT/'summary.json').exists(): raise FileExistsError('Completed comparison exists')
    context = A.reference_context()
    methods = {arm:measure(arm,inferences[arm],context) for arm in R.CONFIGS}
    historical = json.loads((I.ROOT/'reports/M28-METHOD-RESTART/common_ruler.json').read_text())
    comparison = {method:{'per_species':{s:r['methods'][method]['exact_chain'] for s,r in historical['species'].items()},
                          'pooled':values['exact_chain'],'source':'unchanged_historical_cache_not_native_rerun'}
                  for method,values in historical['pooled'].items()}
    for arm,values in methods.items():
        comparison['FULL_'+arm] = {'per_species':{s:r['exact_chain'] for s,r in values['per_species'].items()},
                                   'pooled':values['pooled']['exact_chain'],
                                   'macro':values['pooled']['macro_exact_chain'],
                                   'source':'new_matched_full_TRAIN_last_checkpoint_raw_FASTA',
                                   'checkpoint':fits[arm]['final_checkpoint']}
    gate = A.continuation_gate(methods['Q']['per_species'],methods['C0']['per_species'])
    report = {'experiment':'M28-FULL-TRAIN-MATCHED','status':'complete_one_seed_bounded_DEV_not_independent_generalization',
              'protocol':R.PROTOCOL,'reference_primary_chains':7728,'historical_parent_chains':6450,
              'test_setaria_used':False,'fit':fits,
              'inference_cost':{a:{k:v for k,v in s.items() if k!='scopes'} for a,s in inferences.items()},
              'methods':methods,'exact_chain_comparison':comparison,'predeclared_gate':gate,
              'different_exposure_historical_comparison':historical_exposure_comparison(),
              'scope_limit':'matched exposure/input/order, not matched objective/parameters/runtime; '
                            'historical A/A2/A3/R5 controls have different exposure; no paper/generalization PASS'}
    with (OUT/'summary.json').open('x') as f:
        json.dump(report,f,indent=2,allow_nan=False); f.write('\n')
    print(json.dumps({'comparison':comparison,'gate':gate},allow_nan=False),flush=True)


if __name__ == '__main__': main()
