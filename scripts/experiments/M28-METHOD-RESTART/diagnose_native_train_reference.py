#!/usr/bin/env python3
"""Pair native and original4608 exact hits using completed original TRAIN profiles."""
import itertools
import json
import os
from collections import Counter
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT/'outputs/M28-NATIVE-TRAIN-REFERENCE-DIAGNOSTIC'
CHECKS = ROOT/'outputs/M28-NATIVE-TRAIN-REFERENCE-CHECKS'
PROFILE = ROOT/'outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-CHAIN-PROFILE'
ERRORS = ROOT/'outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-ERROR-PROFILE-R1'
METHODS = ('ANNEVO','Tiberius','Helixer')
MODEL = 'original4608'
EXPECTED = {'NC_003070.9':7135,'NC_003071.7':4293,'NC_003076.8':6311,
            'NC_089043.1':1638,'NC_089044.1':1601,'NC_089046.1':1710}
CATEGORIES = ('exact_primary_chain','no_same_strand_CDS_overlap',
    'same_intron_chain_CDS_start_only_differs','same_intron_chain_CDS_stop_only_differs',
    'same_intron_chain_both_CDS_ends_differ','same_CDS_ends_internal_structure_differs',
    'CDS_ends_and_internal_structure_differ')


def reference_key(row):
    return (row['species'],row['seqid'],row['strand'],
            tuple((c[0],c[1]) for c in row['CDS']))


def record_id(row):
    return row['species'],row['seqid'],row['strand'],row['id']


def join_profiles(references, errors, scope):
    """Do not multiply, silently drop, or cross-join primary references by transcript ID."""
    indexed = {record_id(r):r for r in errors}
    ref_ids = {record_id(r) for r in references}
    if (len(indexed) != len(errors) or len(ref_ids) != len(references)
            or set(indexed) != ref_ids
            or len({reference_key(r) for r in references}) != len(references)
            or Counter((r['species'],r['seqid']) for r in references) != Counter(scope)):
        raise ValueError('Complete unique TRAIN reference identities and denominators must agree')
    records = []
    for ref in references:
        err = indexed[record_id(ref)]
        category = err['category']
        if (err['CDS_segments'] != ref['CDS_segments'] or ref['CDS_segments'] != len(ref['CDS'])
                or category not in CATEGORIES
                or type(ref['exact_primary_hit']) is not bool
                or (category == 'exact_primary_chain') != ref['exact_primary_hit']):
            raise ValueError('Original exact-hit flag and saved error profile disagree')
        records.append({**ref,'original_model_error_category':category,'native_exact_hits':{}})
    return records


def exact_counts(predicted, references):
    return {'tp':len(predicted & references),'predicted':len(predicted),'reference':len(references)}


def summarize(records):
    size = len(records)
    patterns = Counter((r['exact_primary_hit'],)+tuple(r['native_exact_hits'][m] for m in METHODS)
                       for r in records)
    hits = {MODEL:sum(r['exact_primary_hit'] for r in records)}
    hits.update({m:sum(r['native_exact_hits'][m] for r in records) for m in METHODS})
    paired = {}
    for method in METHODS:
        counts, missing = Counter(), Counter()
        for r in records:
            a,n = r['exact_primary_hit'],r['native_exact_hits'][method]
            counts['both' if a and n else 'native_only' if n else 'model_only' if a else 'neither'] += 1
            if n and not a: missing[r['original_model_error_category']] += 1
        paired[method] = {'membership':{k:counts[k] for k in ('both','native_only','model_only','neither')},
            'native_only_by_original_model_error':{k:missing[k] for k in CATEGORIES}}
    any_native = Counter()
    for r in records:
        a,n = r['exact_primary_hit'],any(r['native_exact_hits'].values())
        any_native['both' if a and n else 'native_only' if n else 'model_only' if a else 'neither'] += 1
    return {'reference_chains':size,'exact_primary_hits':hits,
        'reference_recall':{m:n/size if size else None for m,n in hits.items()},
        'native_vs_model':paired,
        'any_native_vs_model':{k:any_native[k] for k in ('both','native_only','model_only','neither')},
        'all_three_native_exact_and_model_missed':sum(
            not r['exact_primary_hit'] and all(r['native_exact_hits'].values()) for r in records),
        'membership_patterns':[{'hits':dict(zip((MODEL,)+METHODS,p)),'reference_chains':patterns[p]}
                               for p in itertools.product((False,True),repeat=4)]}


def describe(records):
    report = summarize(records)
    report['strata'] = {
        'CDS_segments':{k:summarize([r for r in records if (r['CDS_segments']==1)==single])
                       for k,single in (('single',True),('multiple',False))},
        'span_bp':{k:summarize([r for r in records if (r['span_bp']>6144)==long])
                   for k,long in (('at_most_6144',False),('over_6144',True))},
        'original_draw_R4_positive_chain':{
            k:summarize([r for r in records if r['original_draw_R4_positive_chain']==present])
            for k,present in (('present',True),('absent',False))},
        'original_model_error_category':{
            k:summarize([r for r in records if r['original_model_error_category']==k]) for k in CATEGORIES}}
    return report


def completed(folder):
    if (folder/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Required source is not complete: '+str(folder))
    row = json.loads((folder/'summary.json').read_text())
    if not row['completed']: raise ValueError('Required source summary is incomplete')
    return row


def read_rows(path):
    with path.open() as handle: return [json.loads(line) for line in handle]


def reconcile_hits(records, official):
    observed = summarize(records)
    if (observed['reference_chains'] != official['reference']
            or observed['exact_primary_hits'][MODEL] != official['tp']):
        raise ValueError('Saved original4608 reference hits differ from completed scoring')


def main():
    began = time.perf_counter()
    check = completed(CHECKS)
    if not check['fixture_only']: raise ValueError('Synthetic check source differs')
    import evaluate_native_train_baselines as E
    N,C = E.N,E.C
    native,control,profile,error = (completed(p) for p in (E.OUT,E.CONTROL,PROFILE,ERRORS))
    if (native['job_id'] != '13461307' or set(native['methods']) != set(METHODS)
            or native['reference_primary_chains'] != sum(EXPECTED.values())
            or native['status'] != 'native_original_TRAIN_comparison_not_DEV_or_generalization'
            or native['input_mode'] != N.INPUT_MODE
            or profile['job_id'] != '13449803' or error['job_id'] != '13450002'
            or error['source_profile_job'] != profile['job_id']
            or any(s['source_diagnostic_job'] != control['job_id'] for s in (profile,error))
            or any(s['checkpoint'] != control['checkpoint'] for s in (profile,error))
            or any(s['reference_policy'] != control['reference_policy'] for s in (native,profile,error))
            or any(s['DEV_used'] or s['sealed_data_used'] or s['optimizer_steps_this_run'] != 0
                   for s in (native,control,profile,error))):
        raise ValueError('Completed TRAIN source scope or original4608 checkpoint differs')
    N.completed_runs()
    scope = {(species,seqid):EXPECTED[seqid] for species,seqid,_ in N.SCOPE}
    records = join_profiles(read_rows(PROFILE/'reference_chain_profile.jsonl'),
                            read_rows(ERRORS/'reference_error_profile.jsonl'),scope)
    for species,seqid,length in N.SCOPE:
        part = [r for r in records if (r['species'],r['seqid'])==(species,seqid)]
        keys = {reference_key(r) for r in part}
        for method in METHODS:
            official = native['methods'][method]['chromosomes'][seqid]
            prediction = ROOT/official['raw_prediction_path']
            if prediction.resolve() != (N.OUT/method/seqid/'prediction.gff3').resolve():
                raise ValueError('Native prediction path differs from completed evaluation')
            txs = C.S.E.primary_transcripts(C.S.E.parse_annotation(prediction,{seqid:length}),complete_only=False)
            predicted = {(species,)+C.S.key(t) for t in txs}
            if (exact_counts(predicted,keys) != {k:official['exact_chain'][k] for k in ('tp','predicted','reference')}
                    or len(txs) != official['prediction_records']
                    or len(txs)-len(predicted) != official['duplicate_chain_records']):
                raise ValueError('Native unique-chain replay differs from completed common ruler')
            for row in part: row['native_exact_hits'][method] = reference_key(row) in predicted
        reconcile_hits(part,control['chromosomes'][seqid]['exact_chain'])
        if summarize(part)['exact_primary_hits'][MODEL] != profile['chromosomes'][seqid]['exact_chain']['tp']:
            raise ValueError('Original chain profile per-contig hit count differs')
    per_species = {s:describe([r for r in records if r['species']==s]) for s in sorted({x[0] for x in N.SCOPE})}
    pooled = describe(records)
    reconcile_hits(records,control['pooled']['exact_chain'])
    for method in METHODS:
        if pooled['exact_primary_hits'][method] != native['methods'][method]['pooled']['exact_chain']['tp']:
            raise ValueError('Pooled native hits do not reconcile')
        for species,rows in per_species.items():
            if rows['exact_primary_hits'][method] != native['methods'][method]['per_species'][species]['exact_chain']['tp']:
                raise ValueError('Species native hits do not reconcile')
    report = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],
        plan='reports/M28-METHOD-RESTART/native_train_reference_diagnostic_plan.md',
        status='descriptive_original_TRAIN_paired_reference_diagnostic_not_DEV_or_generalization',
        source_evaluation_job=native['job_id'],source_profile_job=profile['job_id'],
        source_error_profile_job=error['job_id'],source_original_evaluation_job=control['job_id'],
        checkpoint=control['checkpoint'],reference_policy=control['reference_policy'],
        reference_primary_chains=len(records),reference_filters_added=False,
        key_definition='species/seqid/strand/ordered 0-based half-open CDS coordinates; phase excluded',
        exact_hit_counts_reproduce_completed_evaluations=True,
        DNA_read=False,real_reference_GFF_read=False,checkpoint_read=False,new_model_inference=False,
        predictions_modified=False,DEV_used=False,sealed_data_used=False,optimizer_steps_this_run=0,
        chromosomes={seqid:summarize([r for r in records if r['seqid']==seqid]) for seqid in EXPECTED},
        per_species=per_species,pooled=pooled,
        limitations=[
            'Reference-level recall intersections only; no prediction precision or F1 is defined for these subsets.',
            'Native-only and model-only hits are annotated primary references, not newly discovered genes.',
            'Error categories describe saved model geometry associations, not causal failure mechanisms.',
            'Absent original R4 positive-chain membership does not mean untrained or independently held out.',
            'All reference strata remain in the full ruler; TRAIN and comparator-dependent pretraining cannot establish generalization.'])
    report['seconds'] = time.perf_counter()-began
    with (OUT/'references.jsonl').open('x') as handle:
        for row in records: handle.write(json.dumps(row,allow_nan=False)+'\n')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps({'completed':True,'reference_primary_chains':len(records),
                     'hits':pooled['exact_primary_hits']}),flush=True)


if __name__ == '__main__': main()
