#!/usr/bin/env python3
"""A3 completed-DEV evaluation; exact F1 versus A, background span versus C0."""
import json
from collections import Counter, defaultdict
import common_ruler as E
from evaluate_pilot import rows, aggregate, stage_audit, validated_scope_records
import evaluate_query as A
from evaluate_query import REF_COUNTS, reference_context, validate_gff_roundtrip
import run_query_known_mean as Q
import evaluate_query_collaboration as A2
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript, write_gff3


def validate_complete(fitted, decoded):
    Q.check_last_identity(fitted); Q.check_last_identity(decoded)
    if (not decoded['inference_complete'] or decoded['windows'] != 8690 or
        decoded['reference_used'] is not False or decoded['test_setaria_used'] is not False or
        decoded['checkpoint'] != fitted['final_checkpoint']):
        raise ValueError('Complete matching A3 DEV inference is required')
    expected = {(s,q,strand):[[a,min(a+A.Q.M.WIDTH,length)] for a in A.Q.M.starts(length)]
                for s,q,length,_ in E.S.SCOPE for strand in ('+','-')}
    actual = {(s['species'],s['seqid'],s['strand']):s['windows'] for s in decoded['scopes']}
    if len(decoded['scopes']) != 4 or actual != expected:
        raise ValueError('A3 full natural DEV scope grid mismatch')


def query_scope_records(path, scope, counts):
    for row in validated_scope_records(path,scope):
        if (row.get('score_kind') != 'query_existence_logit_minus_threshold' or
            row['existence_logit_threshold'] != Q.CONFIG['existence_logit_threshold'] or
            row['counts']['queries'] != Q.CONFIG['num_queries'] or
            len(row['raw_queries']) != Q.CONFIG['num_queries']):
            raise ValueError('A3 query operating point/count/provenance mismatch')
        chains = [tuple(map(tuple, chain)) for chain in row['chains']]
        positive = [chain for chain, score in zip(chains,row['score']) if score > 0]
        counts.update({'windows':1, 'completed_queries':len(chains),
                       'unique_completed_chains_within_window':len(set(chains)),
                       'within_window_duplicate_completed_queries':len(chains)-len(set(chains)),
                       'positive_completed_queries':len(positive),
                       'unique_positive_chains_within_window':len(set(positive)),
                       'within_window_duplicate_positive_queries':len(positive)-len(set(positive))})
        yield row


def validate_diversity_counts(counts, scope):
    q = scope['query_counts']
    if (counts['windows'] != len(scope['windows']) or
        counts['completed_queries'] != q['completed_queries'] or
        counts['positive_completed_queries'] != q['completed_queries']-q['completed_nonpositive_gain']):
        raise ValueError('A3 diversity counts differ from completed inference')


def continuation_gate(current, previous, control):
    if any(set(values) != set(REF_COUNTS) for values in (current,previous,control)):
        raise ValueError('Two matching frozen species are required')
    gate = {}
    for species,count in REF_COUNTS.items():
        now,old,c0 = current[species],previous[species],control[species]
        x,y = now['exact_chain'],old['exact_chain']
        ab,bb,cb = now['background'],old['background'],c0['background']
        if (any(r['exact_chain']['reference'] != count for r in (now,old,c0)) or
            not ab['bp'] == bb['bp'] == cb['bp'] or ab['bp'] <= 0):
            raise ValueError('Reference/background denominator mismatch')
        improved = x['tp']*(y['predicted']+y['reference']) > y['tp']*(x['predicted']+x['reference'])
        gate[species] = {'exact_chain_F1_strictly_above_A':improved,
                         'background_span_not_above_C0':ab['predicted_span_bp'] <= cb['predicted_span_bp']}
    return {'per_species':gate,'passed':all(all(v.values()) for v in gate.values()),
            'not_paper_success_or_independent_generalization':True}


def mechanism_diagnostics(stages, diversity, previous, old_fit):
    saved = json.loads((E.ROOT/'outputs/M28-A-QUERY-DIVERSITY/summary.json').read_text())
    if (saved['checkpoint'] != old_fit['final_checkpoint'] or
        saved['reference_used'] is not False or saved['scores_or_candidates_changed'] is not False):
        raise ValueError('Completed A diversity source identity mismatch')
    diagnostics = {}
    for species,count in REF_COUNTS.items():
        free = sum(r['exact_reference_hits']['free_unique'] for r in stages.values() if r['species'] == species)
        old_free = sum(r['exact_reference_hits']['free_unique'] for r in previous['stages'].values() if r['species'] == species)
        current = Counter()
        old = Counter()
        for key,value in diversity.items():
            if key.startswith(species+'|'): current.update(value)
        for value in saved['scope_reports']:
            if value['species'] == species: old.update(value['counts'])
        n,d = current['within_window_duplicate_positive_queries'],current['positive_completed_queries']
        on,od = old['within_window_duplicate_positive_queries'],old['positive_completed_queries']
        diagnostics[species] = {
            'primary_reference':count,'A_free_exact_hits':old_free,'A3_free_exact_hits':free,
            'free_exact_coverage_above_A':free > old_free,
            'A3_query_diversity':dict(current),'A_query_diversity':dict(old),
            'positive_duplicate_fraction_A3':n/d if d else None,
            'positive_duplicate_fraction_A':on/od if od else None,
            'positive_duplicate_fraction_reduced':n*od < on*d if d and od else None,
            'positive_gain_is_not_biological_truth':True,'not_a_replacement_for_primary_gate':True}
    return diagnostics


def main():
    fit,inference,out = (Q.STAGE_DIRS[k] for k in ('fit','infer','evaluate'))
    fitted = json.loads((fit/'summary.json').read_text())
    decoded = json.loads((inference/'summary.json').read_text())
    validate_complete(fitted,decoded)
    out.mkdir(parents=True,exist_ok=True)
    old_fit = json.loads((A.Q.STAGE_DIRS['fit']/'summary.json').read_text())
    old_infer = json.loads((A.Q.STAGE_DIRS['infer']/'summary.json').read_text())
    A.validate_complete(old_fit,old_infer)
    old_report = json.loads((A.Q.STAGE_DIRS['evaluate']/'summary.json').read_text())
    if old_report['fit']['final_checkpoint'] != old_fit['final_checkpoint']:
        raise ValueError('Completed A comparison checkpoint mismatch')
    previous = old_report['methods']['A']
    a2_report = json.loads((E.ROOT/'outputs/M28-A2-EVALUATION/summary.json').read_text())
    a2_fit = a2_report['fit']
    a2_decoded = json.loads((E.ROOT/'outputs/M28-A2-DEV-INFER-RECOVERY/summary.json').read_text())
    A2.validate_complete(a2_fit, a2_decoded)
    if ({k:v for k,v in fitted['configuration'].items() if k != 'existence_reduction'} != a2_fit['configuration'] or
        fitted['trainable_parameters'] != a2_fit['trainable_parameters'] or
        fitted['optimizer_steps'] != a2_fit['optimizer_steps']):
        raise ValueError('A3 must preserve A2 architecture and exposure')
    for current_epoch, previous_epoch in zip(fitted['epochs'], a2_fit['epochs']):
        if current_epoch['supervision_by_species'] != previous_epoch['supervision_by_species']:
            raise ValueError('A3 changed known/unknown TRAIN supervision counts')
    a2_method = a2_report['methods']['A2']
    refs,isoforms,parents,annotations,sequences,lengths = reference_context()
    predictions,stages = defaultdict(list),{}
    diversity = {}
    for scope in decoded['scopes']:
        species,seqid,strand = scope['species'],scope['seqid'],scope['strand']
        counts = Counter()
        assembled = assemble_predictions(query_scope_records(inference/scope['file'],scope,counts),scope['windows'],strand)
        validate_diversity_counts(counts,scope)
        diversity[species+'|'+strand] = dict(counts)
        for i,row in enumerate(assembled['predictions']):
            ident = 'A3_'+seqid+'_'+('p' if strand == '+' else 'm')+'_'+str(i)
            predictions[species].append(as_transcript(row['chain'],seqid,strand,ident,row['score']))
        audit = stage_audit(assembled,[t for t in refs[species] if t['strand'] == strand])
        audit.update(species=species,seqid=seqid,strand=strand); stages[species+'|'+strand] = audit
        del assembled
    transcripts = [t for values in predictions.values() for t in values]
    gff = out/'A3.gff3'; write_gff3(gff,transcripts,lengths,'M28_A3')
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
    comparison['A'] = {'per_species':{s:r['exact_chain'] for s,r in previous['per_species'].items()},
                       'pooled':previous['pooled']['exact_chain'],'macro':previous['pooled']['macro_exact_chain'],
                       'source':'completed_fixed_A_cache_not_new_fit','checkpoint':old_fit['final_checkpoint']}
    comparison['A3'] = {'per_species':{s:r['exact_chain'] for s,r in per_species.items()},
                       'pooled':pooled['exact_chain'],'macro':pooled['macro_exact_chain'],
                       'source':'completed_last_A3_checkpoint','checkpoint':fitted['final_checkpoint']}
    comparison['A2'] = {'per_species':{s:r['exact_chain'] for s,r in a2_method['per_species'].items()},
                        'pooled':a2_method['pooled']['exact_chain'],
                        'macro':a2_method['pooled']['macro_exact_chain'],
                        'source':'completed_fixed_A2_cache_not_new_fit','checkpoint':a2_fit['final_checkpoint']}
    report = {'experiment':'M28-A3','status':'complete_one_seed_bounded_DEV_not_independent_generalization',
              'reference_primary_chains':7728,'historical_parent_chains':6450,'test_setaria_used':False,
              'proposal':Q.PROPOSAL,'fit':fitted,'inference_cost':{k:v for k,v in decoded.items() if k != 'scopes'},
              'methods':{'A3':{'per_species':per_species,'pooled':pooled,'stages':stages,
                              'query_counts':[s['query_counts'] for s in decoded['scopes']],
                              'GFF3':str(gff.relative_to(E.ROOT))}},
              'exact_chain_comparison':comparison,
              'predeclared_pilot_gate':continuation_gate(per_species,previous['per_species'],controls['C0']['per_species']),
              'scope_limit':'A3 versus A2 changes only existence reduction and its exact matching cost; DEV mechanism, not independent generalization'}
    report['query_diversity'] = diversity
    report['mechanism_diagnostics'] = mechanism_diagnostics(stages,diversity,previous,old_fit)
    report['objective_comparison_A2'] = {'previous_checkpoint':a2_fit['final_checkpoint'],
        'architecture_parameters_initialization_and_exposure_unchanged':True,
        'existence_reduction':'mean_over_known_queries',
        'per_species':{}}
    for species in REF_COUNTS:
        now, old = per_species[species], a2_method['per_species'][species]
        flow = lambda stage, rows:sum(v['exact_reference_hits'][stage] for v in rows.values() if v['species'] == species)
        report['objective_comparison_A2']['per_species'][species] = {
            'exact_chain_A2':old['exact_chain'], 'exact_chain_A3':now['exact_chain'],
            'background_A2':old['background'], 'background_A3':now['background'],
            'free_exact_A2':flow('free_unique',a2_method['stages']),
            'free_exact_A3':flow('free_unique',stages),
            'owner_exact_A2':flow('owner_unique',a2_method['stages']),
            'owner_exact_A3':flow('owner_unique',stages),
            'known_query_mean_is_not_chain_correctness_probability':True}
    report['artifact_bytes_before_summary'] = Q.artifact_bytes()
    with (out/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    Q.artifact_bytes()
    print(json.dumps({'comparison':comparison,'gate':report['predeclared_pilot_gate']},allow_nan=False),flush=True)


if __name__ == '__main__':
    main()
