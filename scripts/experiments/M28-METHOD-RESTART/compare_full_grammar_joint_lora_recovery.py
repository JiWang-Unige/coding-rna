#!/usr/bin/env python3
"""Post-outcome TRAIN recovery accounting for fixed final LoRA versus fresh-joint."""
from collections import Counter
import json
import os
import time
import diagnose_native_train_reference as D
from compare_carbon_a_joint24_train import completed, checked_predictions

ROOT = D.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-JOINT-LORA-RECOVERY'
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_joint_lora_result.md'
SOURCES = {
    'baseline': (ROOT/'outputs/M28-FULL-GRAMMAR-FRESH-JOINT24-TRAIN-DIAGNOSTIC', '13494898'),
    'adapted': (ROOT/'outputs/M28-FULL-GRAMMAR-JOINT-LORA-TRAIN-DIAGNOSTIC', '13499823')}
FIELDS = ('reference_chains', 'retained', 'gained', 'lost', 'still_missed',
          'baseline_hits', 'adapted_hits', 'net_reference_hits')


def pair(rows):
    counts = Counter()
    for row in rows:
        a, b = row['baseline_hit'], row['adapted_hit']
        counts['retained' if a and b else 'gained' if b else 'lost' if a else 'still_missed'] += 1
    baseline, adapted = counts['retained']+counts['lost'], counts['retained']+counts['gained']
    n = len(rows)
    return dict(reference_chains=n, **{k:counts[k] for k in ('retained','gained','lost','still_missed')},
        baseline_hits=baseline, adapted_hits=adapted, net_reference_hits=adapted-baseline,
        baseline_reference_recall=baseline/n if n else None,
        adapted_reference_recall=adapted/n if n else None,
        adapted_minus_baseline_recall_pp=100*(adapted-baseline)/n if n else None)


def describe(rows):
    value = pair(rows)
    groups = {
        'CDS_segments': {'one':[r for r in rows if r['CDS_segments']==1],
            'two_to_five':[r for r in rows if 2<=r['CDS_segments']<=5],
            'six_or_more':[r for r in rows if r['CDS_segments']>=6]},
        'CDS_genomic_span_bp': {'at_most_24576':[r for r in rows if r['span_bp']<=24576],
            'over_24576':[r for r in rows if r['span_bp']>24576]}}
    value['reference_defined_strata'] = {
        name:{key:pair(part) for key,part in bins.items()} for name,bins in groups.items()}
    for bins in value['reference_defined_strata'].values():
        if any(sum(part[k] for part in bins.values()) != value[k] for k in FIELDS):
            raise ValueError('Every reference-defined partition must retain all primary chains')
    return value


def checks():
    rows = [dict(baseline_hit=a, adapted_hit=b, CDS_segments=segments, span_bp=span)
        for a,b,segments,span in ((True,True,1,24575), (False,True,2,24576),
                                 (True,False,5,24577), (False,False,6,60000))]
    result = describe(rows)
    expected = dict(reference_chains=4,retained=1,gained=1,lost=1,still_missed=1,
        baseline_hits=2,adapted_hits=2,net_reference_hits=0)
    if any(result[k]!=v for k,v in expected.items()):
        raise ValueError('Artificial four-way hit accounting failed')
    strata = result['reference_defined_strata']
    if (strata['CDS_segments']['two_to_five']['reference_chains']!=2
            or strata['CDS_genomic_span_bp']['at_most_24576']['reference_chains']!=2
            or pair([])['adapted_reference_recall'] is not None):
        raise ValueError('Artificial segment/span boundary or empty-bin check failed')


def main():
    began = time.perf_counter()
    checks()
    reports = {arm:completed(folder,job) for arm,(folder,job) in SOURCES.items()}
    reference = completed(D.OUT,'13461637')
    baseline = reports['baseline']
    if (reference['reference_policy']!=baseline['reference_policy']
            or reference['reference_primary_chains']!=22688 or reference['reference_filters_added']
            or not reference['exact_hit_counts_reproduce_completed_evaluations']
            or reference['DEV_used'] or reference['sealed_data_used'] or reference['optimizer_steps_this_run']):
        raise ValueError('Completed original primary-reference population required')
    for arm, report in reports.items():
        if (report['reference_policy']!=baseline['reference_policy']
                or report['pooled']['exact_chain']['reference']!=22688
                or set(report['chromosomes'])!=set(D.EXPECTED)
                or report['DEV_used'] or report['sealed_data_used'] or report['optimizer_steps_this_run']
                or report['reference_used_in_inference'] or not report['unsupported_reference_strata_retained']
                or not report['training_eligibility_does_not_filter_reference']):
            raise ValueError('Two completed original-TRAIN final9216 sources required')
    inputs = json.loads((SOURCES['adapted'][0]/'paired_input_match.json').read_text())
    if (inputs['job_id']!='13499823' or inputs['adapted_fit_job']!='13499817'
            or inputs['baseline_fit_job']!='13494895'
            or not inputs['original_six_TRAIN_inputs_directly_equal']
            or not inputs['comparison_before_reference_access'] or inputs['DEV_used'] or inputs['sealed_data_used']):
        raise ValueError('Completed direct DNA-equality receipt required')
    rows = D.read_rows(D.OUT/'references.jsonl')
    scope = {(r['species'],seqid):D.EXPECTED[seqid] for seqid,r in baseline['chromosomes'].items()}
    if (len(rows)!=22688 or len({D.reference_key(r) for r in rows})!=22688
            or Counter((r['species'],r['seqid']) for r in rows)!=Counter(scope)
            or any(not r['CDS'] or r['CDS_segments']!=len(r['CDS'])
                   or r['span_bp']!=r['CDS'][-1][1]-r['CDS'][0][0] for r in rows)):
        raise ValueError('Saved unique reference geometry and six-contig population differ')
    for seqid, official in baseline['chromosomes'].items():
        species,length = official['species'],official['length']
        part = [r for r in rows if (r['species'],r['seqid'])==(species,seqid)]
        refs = {D.reference_key(r) for r in part}
        for arm, report in reports.items():
            source = report['chromosomes'][seqid]
            if (source['species']!=species or source['seqid']!=seqid or source['length']!=length
                    or any(source[k]!=official[k] for k in ('reference_primary_records',
                        'reference_assessable_isoform_records','retained_unsupported_reference_strata'))):
                raise ValueError('Paired chromosome or reference strata changed')
            path = ROOT/source['raw_prediction_path']
            prefix = 'FRESH-JOINT24' if arm=='baseline' else 'JOINT-LORA'
            expected = ROOT/('outputs/M28-FULL-GRAMMAR-'+prefix+'-TRAIN-RAW')/seqid/'prediction/full_grammar.gff3'
            if path.resolve()!=expected.resolve():
                raise ValueError('Saved prediction path differs')
            keys = checked_predictions(path,species,seqid,length,refs,source)
            for row in part:
                row[arm+'_hit'] = D.reference_key(row) in keys
    views = {'pooled':rows, **{s:[r for r in rows if r['species']==s] for s in baseline['per_species']},
             **{q:[r for r in rows if r['seqid']==q] for q in D.EXPECTED}}
    scopes = {name:describe(part) for name,part in views.items()}
    for name, value in scopes.items():
        for arm,report in reports.items():
            source = report['pooled'] if name=='pooled' else report['per_species'][name] if name in report['per_species'] else report['chromosomes'][name]
            if (value[arm+'_hits']!=source['exact_chain']['tp']
                    or value['reference_chains']!=source['exact_chain']['reference']):
                raise ValueError('Replay must reconstruct chromosome/species/pooled primary hits')
    summary = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],protocol=PROTOCOL,
        source_evaluation_jobs={arm:job for arm,(_,job) in SOURCES.items()},
        source_reference_job='13461637',reference_policy=baseline['reference_policy'],
        reference_primary_chains=22688,scopes=scopes,
        artificial_four_way_and_bin_boundary_checks_passed=True,
        exact_hits_unique_predictions_and_duplicates_reproduce_completed_evaluations=True,
        every_reference_partition_reconciled=True,reference_filters_added=False,
        source_DNA_equality_receipt_verified=True,DNA_read=False,real_reference_GFF_read=False,
        saved_primary_reference_profile_read=True,prediction_GFF_read=True,checkpoint_read=False,
        native_forward=False,predictions_modified=False,optimizer_steps_this_run=0,
        DEV_used=False,sealed_data_used=False,frozen_DEV_family_changed=False,
        post_outcome_exploratory_diagnostic=True,seconds=time.perf_counter()-began,
        limitations=['Post-E-outcome descriptive accounting, not pre-DEV preregistration or independent confirmation.',
            'Reference intersections and fixed bins measure recall, not subset prediction precision or F1.',
            'Segment-count and span partitions overlap; do not add their contributions.',
            'Original TRAIN and one fit do not isolate biological causality, generalization or novel genes.'])
    with (OUT/'summary.json').open('x') as handle:
        json.dump(summary,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps(dict(completed=True,reference_primary_chains=22688)),flush=True)


if __name__=='__main__':
    main()

