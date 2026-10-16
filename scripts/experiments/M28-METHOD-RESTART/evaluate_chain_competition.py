#!/usr/bin/env python3
"""A13 complete DEV ruler, matched A6 pool and original gate/native comparators."""
import json
import math
from collections import Counter, defaultdict
from itertools import zip_longest
import infer_chain_competition as I
import evaluate_chain_quality as A6
import evaluate_query_context as H
from evaluate_pilot import aggregate, stage_audit, validated_scope_records
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript

A = A6.A
E = A.E
OUT = I.ROOT/'outputs/M28-A13-CHAIN-COMPETITION-EVAL'


def validate_complete(fitted, decoded):
    I.validate_fit(fitted)
    expected = {(s,q,t):[[a,min(a+I.M.WIDTH,n)] for a in I.M.starts(n)]
        for s,qs in I.C.I.DEV_LENGTHS.items() for q,n in qs.items() for t in ('+','-')}
    actual = {(s['species'],s['seqid'],s['strand']):s['windows'] for s in decoded['scopes']}
    if (decoded['arm'] != 'A13' or decoded['completed'] is not True or decoded['inference_complete'] is not True or
        decoded['configuration'] != I.R.CONFIG or decoded['protocol'] != I.DEV_PROTOCOL or
        decoded['TRAIN_protocol'] != I.R.PROTOCOL or decoded['TRAIN_preparation_gate_passed'] is not True or
        decoded['checkpoint'] != fitted['final_checkpoint'] or decoded['optimizer_steps'] != 4608 or
        decoded['generator_checkpoint'] != I.C.GENERATOR or decoded['feature_revision'] != I.C.I.F.REVISION or
        decoded['feature_spec'] != I.C.R.FEATURE_SPEC or decoded['input_mode'] != 'fixed_DEV_raw_FASTA' or
        decoded['reference_used'] is not False or decoded['labels_or_feature_cache_used'] is not False or
        decoded['optimizer_steps_this_run'] != 0 or decoded['DEV_used_for_selection'] is not False or
        decoded['test_setaria_used'] is not False or decoded['windows'] != 8690 or
        len(decoded['scopes']) != 4 or actual != expected or
        decoded['lengths'] != {q:n for qs in I.C.I.DEV_LENGTHS.values() for q,n in qs.items()}):
        raise ValueError('Complete fixed 8690-window A13 raw-FASTA DEV grid required')


def matched_records(path, scope, old_path, old_scope, counts):
    """Direct raw geometry/order/origin comparison; neither reference nor hashes."""
    current = validated_scope_records(path,scope)
    previous = A6.records(old_path,old_scope)
    for row,old in zip_longest(current,previous):
        if row is None or old is None:
            raise ValueError('Partial paired A6/A13 inference scope')
        chains,scores,sources = row['chains'],row['score'],row['decoder_sources']
        actual = dict(free_candidates=len(chains),nonpositive_candidates=sum(v <= 0 for v in scores))
        if (row['score_kind'] != I.SCORE_KIND or row['quality_logit_threshold'] != 0. or
            len(chains) != len(scores) or len(chains) != len(sources) or
            any(not math.isfinite(v) for v in scores) or row['counts'] != actual or
            any(not ss or not set(ss) <= {'greedy','beam4'} for ss in sources) or
            len({tuple(map(tuple,c)) for c in chains}) != len(chains)):
            raise ValueError('A13 score/pool/source/count boundary differs')
        if any(row[k] != old[k] for k in ('window_index','window_id','chains','decoder_sources')):
            raise ValueError('Matched frozen-generator A6/A13 candidate geometry/order/source differs')
        counts.update(actual)
        yield row


def reconstruct(decoded, baseline):
    """Close every raw record, matched natural pool and GFF before reference access."""
    old_scopes = {(s['species'],s['seqid'],s['strand']):s for s in baseline['scopes']}
    predictions,assemblies,counts_by_scope = defaultdict(list),{},{}
    for scope in decoded['scopes']:
        species,seqid,strand = scope['species'],scope['seqid'],scope['strand']
        old = old_scopes[(species,seqid,strand)]; counts = Counter()
        assembled = assemble_predictions(matched_records(I.OUT/scope['file'],scope,
            I.C.INFER/old['file'],old,counts),scope['windows'],strand)
        if dict(counts) != scope['quality_counts'] or assembled['counts'] != scope['assembly_counts']:
            raise ValueError('A13 actual record counts/assembly differ from inference')
        for index,row in enumerate(assembled['predictions']):
            predictions[species].append(as_transcript(row['chain'],seqid,strand,
                I.prediction_id(scope['ordinal'],strand,index),row['score']))
        key = species+'|'+strand; assemblies[key] = assembled; counts_by_scope[key] = dict(counts)
    transcripts = [t for values in predictions.values() for t in values]
    parsed = E.S.E.primary_transcripts(E.S.E.parse_annotation(I.OUT/'A13.gff3',decoded['lengths']),complete_only=False)
    A.validate_gff_roundtrip(parsed,transcripts)
    if len(transcripts) != decoded['primary_CDS_chain_predictions']:
        raise ValueError('Final A13 actual GFF prediction count differs')
    return predictions,assemblies,counts_by_scope


def measure(predictions, assemblies, counts, context):
    refs,isoforms,parents,annotations,sequences,_ = context
    per_species,stages = {},{}
    for key,assembled in assemblies.items():
        species,strand = key.split('|')
        audit = stage_audit(assembled,[r for r in refs[species] if r['strand'] == strand])
        audit.update(species=species,seqid=next(iter(I.C.I.DEV_LENGTHS[species])),strand=strand)
        stages[key] = audit
    for species,seqid,length,_ in E.S.SCOPE:
        measured = A6.score_species(predictions[species],refs[species],isoforms[species],
            length,annotations[species],parents[species])
        measured['prediction_nonexact_primary'] = measured.pop('reference_nonexact_primary')
        failures = Counter(); valid = 0
        for tx in predictions[species]:
            problems,_ = I.M.chain_checks(tx,sequences[species][seqid])
            failures.update(problems); valid += not problems
        measured['CDS_structure_diagnostic'] = dict(fully_supported=valid,
            predictions=len(predictions[species]),failure_counts_nonexclusive=dict(failures))
        per_species[species] = measured
    pooled = aggregate(per_species)
    pooled['actual_free_primary_coverage'] = sum(x['exact_reference_hits']['free_unique'] for x in stages.values())/7728
    pooled['owner_primary_coverage'] = sum(x['exact_reference_hits']['owner_unique'] for x in stages.values())/7728
    if sum(x['exact_reference_hits']['selected_unique'] for x in stages.values()) != pooled['exact_chain']['tp']:
        raise ValueError('A13 selected-stage TP disagrees with exact-chain score')
    return dict(per_species=per_species,pooled=pooled,stages=stages,quality_counts=counts,
        GFF3=str((I.OUT/'A13.gff3').relative_to(I.ROOT)))


def main():
    if any((p/'STATUS').read_text().strip() != 'COMPLETED'
        for p in (I.R.FIT,I.R.READOUT,I.OUT,I.C.FIT,I.C.INFER)):
        raise ValueError('Terminal A13 TRAIN/inference and closed A6 sources required')
    fitted = json.loads((I.R.FIT/'summary.json').read_text())
    prepared = json.loads((I.R.READOUT/'summary.json').read_text()); I.validate_train_preparation(fitted,prepared)
    decoded = json.loads((I.OUT/'summary.json').read_text()); validate_complete(fitted,decoded)
    old_fit = json.loads((I.C.FIT/'summary.json').read_text())
    baseline = json.loads((I.C.INFER/'summary.json').read_text()); A6.validate_complete(old_fit,baseline)
    if ({k:v for k,v in fitted['configuration'].items() if k not in
        ('loss','pool','positive_negative_group_mass','conflicts','temperature','margin','initialization','main_evaluation')} !=
        {k:v for k,v in old_fit['configuration'].items() if k not in ('loss','pool','main_evaluation')}):
        raise ValueError('Matched A6 generator/architecture/optimizer/exposure settings differ')
    predictions,assemblies,counts = reconstruct(decoded,baseline)
    current = measure(predictions,assemblies,counts,A.reference_context())
    previous = H.read_report(I.ROOT/'outputs/M28-A-EVALUATION/summary.json')['methods']['A']
    control = H.read_report(I.ROOT/'reports/M28-METHOD-RESTART/r5_result.json')['methods']['C0']
    matched = H.read_report(I.ROOT/'outputs/M28-A6-CHAIN-QUALITY-EVAL/summary.json')
    if matched['completed'] is not True or matched['fit']['final_checkpoint'] != old_fit['final_checkpoint']:
        raise ValueError('Completed same A6 comparator required')
    native = H.read_report(I.ROOT/'outputs/M28-NATIVE-BASELINE-EVAL/summary.json')
    if native['completed'] is not True or native['unsupported_reference_strata_retained'] is not True:
        raise ValueError('Completed native comparison retaining unsupported reference strata required')
    comparisons = dict(A6=matched['methods']['A6'],A=previous,R5_C0=control)
    for method in ('ANNEVO','Tiberius','Helixer'):
        comparisons[method] = dict(per_species={s:r['methods'][method] for s,r in native['species'].items()},
            pooled=native['pooled'][method])
    for values in comparisons.values():
        H.check_ruler(values,current)
    gate = A6.continuation_gate(current['per_species'],previous['per_species'],control['per_species'])
    report = dict(experiment='M28-A13-CHAIN-COMPETITION',completed=True,protocol=I.DEV_PROTOCOL,
        status='complete_one_seed_bounded_DEV_not_independent_generalization',
        reference_primary_chains=7728,test_setaria_used=False,fit=fitted,TRAIN_preparation=prepared,
        inference_receipt={k:v for k,v in decoded.items() if k != 'scopes'},methods={'A13':current},
        matched_A6_natural_candidate_geometry_order_sources_unchanged=True,
        predeclared_pilot_gate=gate,contrasts={name:H.contrast(current,v) for name,v in comparisons.items()},
        exact_chain_comparison={name:dict(per_species={s:r['exact_chain'] for s,r in v['per_species'].items()},
            pooled=v['pooled']['exact_chain'],macro_F1=sum(r['exact_chain']['f1'] for r in v['per_species'].values())/2)
            for name,v in {'A13':current,**comparisons}.items()},
        comparison_sources=dict(A6='outputs/M28-A6-CHAIN-QUALITY-EVAL/summary.json',
            A='outputs/M28-A-EVALUATION/summary.json',R5_C0='reports/M28-METHOD-RESTART/r5_result.json',
            native='outputs/M28-NATIVE-BASELINE-EVAL/summary.json'),
        matched_scope='Same A6 head architecture/seed0 recipe/1536 original draws x3/optimizer and '
            'frozen generator/natural pool; objective changes local competition and group mass together.',
        native_exposure_limit=native['source_limitations'],
        limitations='Repeated-use DEV, not independent generalization. Relative A6 improvement does not '
            'replace original F1/background gate, native performance, novelty evidence or paper success.')
    I.C.R.write_summary(OUT,report)
    print(json.dumps(dict(comparison=report['exact_chain_comparison'],gate=gate),allow_nan=False),flush=True)


if __name__ == '__main__':
    main()
