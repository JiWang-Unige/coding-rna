#!/usr/bin/env python3
"""A10 fixed DEV ruler, original continuation gate and completed native comparators."""
import json
from collections import Counter, defaultdict
import infer_query_context as I
import evaluate_query_known_mean as A3
from evaluate_pilot import aggregate, stage_audit
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript

E = A3.E
OUT = I.ROOT/'outputs/M28-A10-CONTEXT-EVAL'


def validate_complete(fitted, decoded):
    I.validate_fit(fitted)
    expected = {(s,q,t):[[a,min(a+I.M.WIDTH,n)] for a in I.M.starts(n)]
        for s,qs in I.DEV_LENGTHS.items() for q,n in qs.items() for t in ('+','-')}
    actual = {(s['species'],s['seqid'],s['strand']):s['windows'] for s in decoded['scopes']}
    if (decoded['inference_complete'] is not True or decoded['arm'] != I.ARM or
        decoded['configuration'] != I.R.CONFIG or decoded['optimizer_steps'] != 4608 or
        decoded['checkpoint'] != fitted['final_checkpoint'] or
        decoded['encoder_model'] != I.F.MODEL or decoded['encoder_revision'] != I.F.REVISION or
        decoded['feature_spec'] != I.C.FULL_SPEC or decoded['input_mode'] != 'fixed_DEV_raw_FASTA' or
        decoded['reference_used'] is not False or decoded['labels_or_feature_cache_used'] is not False or
        decoded['optimizer_steps_this_run'] != 0 or decoded['test_setaria_used'] is not False or
        decoded['DEV_used_for_selection'] is not False or decoded['windows'] != 8690 or
        len(decoded['scopes']) != 4 or actual != expected or
        decoded['lengths'] != {q:n for qs in I.DEV_LENGTHS.values() for q,n in qs.items()} or
        any(s['reference_used'] is not False or s['labels_or_feature_cache_used'] is not False
            for s in decoded['scopes'])):
        raise ValueError('Complete matching 8690-window A10 raw-FASTA DEV grid required')


def reconstruct(decoded):
    """Reconcile every free record, assembly and GFF before opening any reference."""
    predictions,assemblies,diversity = defaultdict(list),{},{}
    for scope in decoded['scopes']:
        species,seqid,strand = scope['species'],scope['seqid'],scope['strand']
        counts = Counter()
        assembled = assemble_predictions(A3.query_scope_records(I.OUT/scope['file'],scope,counts),
            scope['windows'],strand)
        A3.validate_diversity_counts(counts,scope)
        if assembled['counts'] != scope['assembly_counts']:
            raise ValueError('Reassembled candidates differ from A10 inference')
        for index,row in enumerate(assembled['predictions']):
            predictions[species].append(as_transcript(row['chain'],seqid,strand,
                I.prediction_id(scope['ordinal'],strand,index),row['score']))
        key = species+'|'+strand; assemblies[key] = assembled; diversity[key] = dict(counts)
    transcripts = [t for values in predictions.values() for t in values]
    parsed = E.S.E.primary_transcripts(E.S.E.parse_annotation(I.OUT/(I.ARM+'.gff3'),decoded['lengths']),
        complete_only=False)
    A3.validate_gff_roundtrip(parsed,transcripts)
    if len(transcripts) != decoded['primary_CDS_chain_predictions']:
        raise ValueError('Final A10 prediction count differs from inference')
    return predictions,assemblies,diversity


def measure(predictions, assemblies, diversity, context):
    refs,isoforms,parents,annotations,sequences,_ = context
    stages,per_species = {},{}
    for key,assembled in assemblies.items():
        species,strand = key.split('|')
        audit = stage_audit(assembled,[t for t in refs[species] if t['strand'] == strand])
        seqid = next(iter(I.DEV_LENGTHS[species]))
        audit.update(species=species,seqid=seqid,strand=strand); stages[key] = audit
    for species,seqid,length,_ in E.S.SCOPE:
        pred = predictions[species]
        measured = E.score(pred,refs[species],isoforms[species],length,annotations[species])
        measured['prediction_nonexact_primary'] = measured.pop('reference_nonexact_primary')
        pk,hk = {E.S.key(t) for t in pred},{E.S.key(t) for t in parents[species]}
        measured['historical_parent_exact_chain'] = E.S.prf(len(pk&hk),len(pk),len(hk))
        failures = Counter(); valid = 0
        for tx in pred:
            problems,_ = I.M.chain_checks(tx,sequences[species][seqid])
            failures.update(problems); valid += not problems
        measured['CDS_structure_diagnostic'] = dict(fully_supported=valid,predictions=len(pred),
            failure_counts_nonexclusive=dict(failures))
        per_species[species] = measured
    pooled = aggregate(per_species)
    pooled['actual_free_primary_coverage'] = sum(x['exact_reference_hits']['free_unique'] for x in stages.values())/7728
    pooled['owner_primary_coverage'] = sum(x['exact_reference_hits']['owner_unique'] for x in stages.values())/7728
    if sum(x['exact_reference_hits']['selected_unique'] for x in stages.values()) != pooled['exact_chain']['tp']:
        raise ValueError('A10 selected-stage TP disagrees with exact-chain score')
    return dict(per_species=per_species,pooled=pooled,stages=stages,candidate_diversity=diversity,
        GFF3=str((I.OUT/(I.ARM+'.gff3')).relative_to(I.ROOT)))


def read_report(path):
    report = json.loads(path.read_text())
    if report['reference_primary_chains'] != 7728 or report['test_setaria_used'] is not False:
        raise ValueError('Comparator must retain the frozen complete DEV ruler')
    return report


def check_ruler(values, current):
    if set(values['per_species']) != set(A3.REF_COUNTS):
        raise ValueError('Comparator species differ')
    for species,count in A3.REF_COUNTS.items():
        old,new = values['per_species'][species],current['per_species'][species]
        if old['exact_chain']['reference'] != count or old['background']['bp'] != new['background']['bp']:
            raise ValueError('Comparator primary/background denominators differ')
    if values['pooled']['exact_chain']['reference'] != 7728:
        raise ValueError('Comparator pooled primary denominator differs')


def contrast(current, previous):
    def one(now, old):
        x,y = now['exact_chain'],old['exact_chain']
        return dict(exact_chain_F1_difference_percentage_points=100*(x['f1']-y['f1']),
            exact_chain_F1_strictly_above=x['tp']*(y['predicted']+y['reference']) >
                y['tp']*(x['predicted']+x['reference']),
            background_span_difference_bp=now['background']['predicted_span_bp']-old['background']['predicted_span_bp'],
            CDS_base_F1_difference_percentage_points=100*(now['CDS_base_unstranded_union']['f1']-
                old['CDS_base_unstranded_union']['f1']))
    return dict(per_species={s:one(current['per_species'][s],previous['per_species'][s]) for s in A3.REF_COUNTS},
        pooled=one(current['pooled'],previous['pooled']),not_independent_generalization=True)


def main():
    if any((p/'STATUS').read_text().strip() != 'COMPLETED' for p in (I.R.FIT,I.OUT)):
        raise ValueError('Terminal A10 fit and raw inference required before reference access')
    fitted = json.loads((I.R.FIT/'summary.json').read_text())
    decoded = json.loads((I.OUT/'summary.json').read_text()); validate_complete(fitted,decoded)
    predictions,assemblies,diversity = reconstruct(decoded)
    current = measure(predictions,assemblies,diversity,A3.reference_context())
    matched = read_report(I.ROOT/'outputs/M28-A3-EVALUATION/summary.json')
    old_fit = matched['fit']; A3.Q.check_last_identity(old_fit)
    if ({k:v for k,v in fitted['configuration'].items() if k != 'encoder_context'} != old_fit['configuration'] or
        fitted['trainable_parameters'] != old_fit['trainable_parameters'] or
        [e['steps'] for e in old_fit['epochs']] != [1536]*3):
        raise ValueError('A3 head configuration and fixed TRAIN exposure must remain matched')
    previous = read_report(I.ROOT/'outputs/M28-A-EVALUATION/summary.json')['methods']['A']
    control = read_report(I.ROOT/'reports/M28-METHOD-RESTART/r5_result.json')['methods']['C0']
    native = read_report(I.ROOT/'outputs/M28-NATIVE-BASELINE-EVAL/summary.json')
    if native['completed'] is not True or native['unsupported_reference_strata_retained'] is not True:
        raise ValueError('Completed native comparison with unsupported references retained required')
    comparisons = {'A3':matched['methods']['A3'],'A':previous,'R5_C0':control}
    for method in ('ANNEVO','Tiberius','Helixer'):
        comparisons[method] = dict(per_species={s:r['methods'][method] for s,r in native['species'].items()},
            pooled=native['pooled'][method])
    for values in comparisons.values():
        check_ruler(values,current)
    contrasts = {name:contrast(current,values) for name,values in comparisons.items()}
    gate = A3.continuation_gate(current['per_species'],previous['per_species'],control['per_species'])
    report = dict(experiment='M28-A10-CONTEXT',
        status='complete_one_seed_bounded_DEV_not_independent_generalization',proposal=I.R.PROPOSAL,
        reference_primary_chains=7728,historical_parent_chains=6450,test_setaria_used=False,
        fit=fitted,inference_cost={k:v for k,v in decoded.items() if k != 'scopes'},
        methods={I.ARM:current},predeclared_gate=gate,contrasts=contrasts,
        exact_chain_comparison={name:dict(per_species={s:r['exact_chain'] for s,r in v['per_species'].items()},
            pooled=v['pooled']['exact_chain'],macro_F1=sum(r['exact_chain']['f1'] for r in v['per_species'].values())/2)
            for name,v in {I.ARM:current,**comparisons}.items()},
        comparison_sources=dict(A3='outputs/M28-A3-EVALUATION/summary.json',
            A='outputs/M28-A-EVALUATION/summary.json',R5_C0='reports/M28-METHOD-RESTART/r5_result.json',
            native='outputs/M28-NATIVE-BASELINE-EVAL/summary.json'),
        matched_scope='A3 head, original draws/order/weights/objective/4608 updates matched; '
            'attention context and positions both changed; A3 cached inference timing is not a fair native runtime comparator',
        native_exposure_limit=native['source_limitations'],
        scope_limit='Repeated-use two-species DEV; relative A3 improvement does not replace original gate, '
            'native rankings, independent tests or novelty evidence. No paper/generalization PASS.')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps(dict(comparison=report['exact_chain_comparison'],gate=gate),allow_nan=False),flush=True)


if __name__ == '__main__':
    main()
