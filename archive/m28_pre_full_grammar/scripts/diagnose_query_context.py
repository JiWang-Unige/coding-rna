#!/usr/bin/env python3
"""Completed A10 versus A3/native primary-reference turnover; no new model inference."""
import csv
import json
import time
from collections import Counter
import diagnose_native_reference as N
from diagnose_b2_recovery import chain_key, recovery_state
import evaluate_query_context as E
from src.m28.assembly import assemble_predictions, nearest_owner
from src.m28.c0_decode import genomic_chain

ROOT = E.I.ROOT
OUT = ROOT/'outputs/M28-A10-CONTEXT-DIAGNOSTIC'
BASELINE = N.OUT
STAGES, METHODS = N.STAGES, N.METHODS


def parse_baseline(raw):
    row = N.parse_record(raw)
    for method in METHODS:
        value = raw[method+'_exact_primary']
        if value not in ('True','False'):
            raise ValueError('Native membership must use explicit CSV booleans')
        row[method+'_exact_primary'] = value == 'True'
    return row


def replay_scope(path, scope, reference):
    windows = scope['windows']; centers = [a+b for a,b in windows]
    gains = {}; counts = Counter()
    def observed():
        for row in E.A3.query_scope_records(path,scope,counts):
            index = row['window_index']; a,b = windows[index]
            for chain,score in zip(row['chains'],row['score']):
                key = tuple(genomic_chain(chain,a,b,scope['strand']))
                if key in reference and nearest_owner(key,windows,centers) == index:
                    gains[key] = max(float(score),gains.get(key,float('-inf')))
            yield row
    assembled = assemble_predictions(observed(),windows,scope['strand'])
    E.A3.validate_diversity_counts(counts,scope)
    if assembled['counts'] != scope['assembly_counts']:
        raise ValueError('A10 reference replay differs from original assembly')
    return assembled,gains


def add_context(row, stage_sets, gains):
    key = chain_key(row['CDS']); flags = [key in stage_sets[s] for s in STAGES]
    if (flags[1] and not flags[0]) or (flags[2] and not flags[1]):
        raise ValueError('A10 stage membership must be nested')
    for stage,value in zip(STAGES,flags):
        row['A10_'+stage] = value
    gain = gains.get(key)
    row['A10_owner_gain'] = gain
    row['A10_state'] = recovery_state(*flags,gain)


def exon_stratum(count):
    return '1' if count == 1 else '2-4' if count <= 4 else '5-9' if count <= 9 else '10+'


def span_stratum(length):
    return '<=6144' if length <= 6144 else '6145-24576' if length <= 24576 else '>24576'


def summarize(records):
    n = len(records); result = dict(reference_chains=n,stages={},states={},native_vs_A10={})
    for stage in STAGES:
        counts = Counter()
        for row in records:
            old,new = row['A3_'+stage],row['A10_'+stage]
            counts['both' if old and new else 'A10_only' if new else 'A3_only' if old else 'neither'] += 1
        result['stages'][stage] = {k:counts[k] for k in ('both','A10_only','A3_only','neither')}
    for arm in ('A3','A10'):
        result['states'][arm] = dict(Counter(r[arm+'_state'] for r in records))
        result.setdefault('reference_hit_rates',{})[arm] = {
            s:sum(r[arm+'_'+s] for r in records)/n if n else None for s in STAGES}
    states = Counter((r['A3_state'],r['A10_state']) for r in records)
    result['A3_to_A10_states'] = [dict(A3=a,A10=b,reference_chains=n) for (a,b),n in sorted(states.items())]
    for method in METHODS:
        counts,missing = Counter(),Counter()
        for row in records:
            a,b = row['A10_selected_unique'],row[method+'_exact_primary']
            counts['both' if a and b else 'native_only' if b else 'A10_only' if a else 'neither'] += 1
            if b and not a:
                missing[row['A10_state']] += 1
        result['native_vs_A10'][method] = dict(
            membership={k:counts[k] for k in ('both','native_only','A10_only','neither')},
            native_only_by_A10_state=dict(missing))
    result['A10_selected_and_no_native_exact_hit'] = sum(r['A10_selected_unique'] and
        not any(r[m+'_exact_primary'] for m in METHODS) for r in records)
    result['not_subset_precision_or_F1'] = True
    return result


def main():
    began = time.perf_counter()
    if OUT.exists():
        raise FileExistsError('Preserve original or failed diagnostic output')
    if any((p/'STATUS').read_text().strip() != 'COMPLETED' for p in (E.I.R.FIT,E.I.OUT,E.OUT,BASELINE)):
        raise ValueError('A10 fit/infer/eval and accepted baseline diagnosis must all be terminal')
    fit = json.loads((E.I.R.FIT/'summary.json').read_text())
    decoded = json.loads((E.I.OUT/'summary.json').read_text())
    evaluated = E.read_report(E.OUT/'summary.json'); E.validate_complete(fit,decoded)
    baseline = json.loads((BASELINE/'summary.json').read_text())
    before = E.read_report(ROOT/'outputs/M28-A3-EVALUATION/summary.json')
    if (evaluated['experiment'] != 'M28-A10-CONTEXT' or evaluated['fit'] != fit or
        set(evaluated['methods']) != {E.I.ARM} or baseline['completed'] is not True or
        baseline['reference_primary_chains'] != 7728 or baseline['reference_filters_added'] is not False or
        baseline['test_setaria_used'] is not False or
        baseline['source']['A3_checkpoint'] != before['fit']['final_checkpoint'] or
        evaluated['comparison_sources']['native'] != 'outputs/M28-NATIVE-BASELINE-EVAL/summary.json'):
        raise ValueError('Completed matched A10/A3/native source identities required')
    with (BASELINE/'references.tsv').open() as handle:
        records = [parse_baseline(row) for row in csv.DictReader(handle,delimiter='\t')]
    if len(records) != 7728 or len({N.reference_key(r) for r in records}) != 7728:
        raise ValueError('All unique frozen primary references required')
    for species,count in N.COUNTS.items():
        rows = [r for r in records if r['species'] == species]
        official = baseline['per_species'][species]
        if (len(rows) != count or any(sum(r['A3_'+s] for r in rows) != official['A3_stage_hits'][s]
            for s in STAGES) or any(sum(r[m+'_exact_primary'] for r in rows) !=
            official['native_vs_A3'][m]['exact_primary_hits'] for m in METHODS)):
            raise ValueError('Accepted A3/native membership source is incomplete or differs')
    for scope in decoded['scopes']:
        species,seqid,strand = scope['species'],scope['seqid'],scope['strand']
        rows = [r for r in records if (r['species'],r['seqid'],r['strand']) == (species,seqid,strand)]
        reference = {chain_key(r['CDS']) for r in rows}
        assembled,gains = replay_scope(E.I.OUT/scope['file'],scope,reference)
        sets = {s:{chain_key(c) for c in assembled['stages'][s]} & reference for s in STAGES}
        official = evaluated['methods'][E.I.ARM]['stages'][species+'|'+strand]
        if ({s:len(sets[s]) for s in STAGES} != official['exact_reference_hits'] or
            assembled['counts'] != official['candidate_counts']):
            raise ValueError('Per-reference replay does not match completed A10 stage accounting')
        for row in rows:
            add_context(row,sets,gains)
    pooled = summarize(records)
    if sum(r['A10_selected_unique'] for r in records) != evaluated['methods'][E.I.ARM]['pooled']['exact_chain']['tp']:
        raise ValueError('Final A10 reference membership differs from completed exact-chain TP')
    report = dict(completed=True,experiment=OUT.name,reference_primary_chains=7728,
        reference_filters_added=False,new_training=False,new_model_inference=False,
        optimizer_steps_this_run=0,test_setaria_used=False,pooled=pooled,
        per_species={s:summarize([r for r in records if r['species'] == s]) for s in N.COUNTS},
        exon_strata={k:summarize([r for r in records if exon_stratum(r['exon_count']) == k])
            for k in ('1','2-4','5-9','10+')},
        span_strata={k:summarize([r for r in records if span_stratum(r['span_bp']) == k])
            for k in ('<=6144','6145-24576','>24576')},
        unsupported_strata={k:summarize([r for r in records if bool(r['unsupported_strata']) == flag])
            for k,flag in (('none_recorded',False),('any_recorded',True))},
        source=dict(A10_checkpoint=fit['final_checkpoint'],A10_inference=str(E.I.OUT.relative_to(ROOT)),
            A10_evaluation=str(E.OUT.relative_to(ROOT)),accepted_A3_native_membership=str(BASELINE.relative_to(ROOT)),
            A3_checkpoint=before['fit']['final_checkpoint']),
        predeclared_gate_unchanged=evaluated['predeclared_gate'],seconds=time.perf_counter()-began,
        limitations=[
            'Descriptive repeated-use DEV reference sensitivity/turnover only, not stratum precision/F1 or a new gate.',
            'Reference CDS segment count is not a count of UTR exons. Unsupported, long and overlapping references remain.',
            'Span relative to the old 6144bp encoding block is descriptive; context and position handling both changed.',
            'Stage attrition is not a unique causal mechanism; matched full-ruler F1/background and native rankings stay primary.',
            'A10-only exact hits are already annotated primary references, not new gene discoveries or independent generalization.'])
    OUT.mkdir()
    for row in records:
        row.pop('CDS')
    with (OUT/'references.tsv').open('x',newline='') as handle:
        writer = csv.DictWriter(handle,fieldnames=list(records[0]),delimiter='\t')
        writer.writeheader(); writer.writerows(records)
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps(dict(completed=True,reference_primary_chains=7728,seconds=report['seconds'])),flush=True)


if __name__ == '__main__':
    main()
