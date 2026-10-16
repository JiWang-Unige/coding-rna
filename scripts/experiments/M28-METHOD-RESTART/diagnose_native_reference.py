#!/usr/bin/env python3
"""Join fixed native exact hits to completed A3 reference-stage turnover; DEV diagnosis only."""
import csv
import json
from collections import Counter
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from diagnose_b2_recovery import chain_key, recovery_state
from diagnose_query_turnover import require_source

OUT = ROOT/'outputs/M28-NATIVE-A3-REFERENCE-DIAGNOSTIC'
TURN = ROOT/'outputs/M28-A2-A3-REFERENCE-TURNOVER'
NATIVE = ROOT/'outputs/M28-NATIVE-BASELINE-EVAL'
METHODS = ('ANNEVO', 'Tiberius', 'Helixer')
STAGES = ('free_unique', 'owner_unique', 'selected_unique')
COUNTS = {'arabidopsis_thaliana':5437, 'oryza_sativa':2291}


def reference_key(row):
    return row['species'], row['seqid'], row['strand'], chain_key(row['CDS'])


def parse_record(raw):
    row = {k:raw[k] for k in ('species','id','gene_id','seqid','strand','exon_count',
                              'span_bp','coding_bp','CDS_intervals','unsupported_strata')}
    flags = []
    for stage in STAGES:
        value = raw['A3_'+stage]
        if value not in ('True','False'):
            raise ValueError('A3 stage flags must not use CSV string truthiness')
        row['A3_'+stage] = value == 'True'
        flags.append(value == 'True')
    if (flags[1] and not flags[0]) or (flags[2] and not flags[1]):
        raise ValueError('A3 stages are not nested')
    gain = float(raw['A3_owner_gain']) if raw['A3_owner_gain'] else None
    if recovery_state(*flags,gain) != raw['A3_state']:
        raise ValueError('A3 state does not match flags and recorded owner gain')
    row['A3_owner_gain'], row['A3_state'] = gain, raw['A3_state']
    row['exon_count'], row['span_bp'], row['coding_bp'] = (
        int(row[k]) for k in ('exon_count','span_bp','coding_bp'))
    row['CDS'] = json.loads(row['CDS_intervals'])
    return row


def summarize(records):
    report = {'reference_chains':len(records),
              'A3_stage_hits':{s:sum(r['A3_'+s] for r in records) for s in STAGES},
              'A3_states':dict(Counter(r['A3_state'] for r in records)),
              'native_vs_A3':{}}
    for method in METHODS:
        counts = Counter()
        missing = Counter()
        for row in records:
            a, n = row['A3_selected_unique'], row[method+'_exact_primary']
            counts['both' if a and n else 'native_only' if n else 'A3_only' if a else 'neither'] += 1
            if n and not a:
                missing[row['A3_state']] += 1
        report['native_vs_A3'][method] = {
            'exact_primary_hits':sum(r[method+'_exact_primary'] for r in records),
            'membership':{k:counts[k] for k in ('both','native_only','A3_only','neither')},
            'native_only_by_A3_state':dict(missing)}
    patterns = Counter(tuple(r[k] for k in
                       ('A3_selected_unique',)+tuple(m+'_exact_primary' for m in METHODS))
                       for r in records)
    report['selected_membership_patterns'] = [
        {'A3':a,'ANNEVO':b,'Tiberius':c,'Helixer':d,'reference_chains':n}
        for (a,b,c,d),n in sorted(patterns.items())]
    report['A3_selected_and_no_native_exact_hit'] = sum(
        r['A3_selected_unique'] and not any(r[m+'_exact_primary'] for m in METHODS)
        for r in records)
    report['any_native_exact_hit_absent_from_A3_free'] = sum(
        not r['A3_free_unique'] and any(r[m+'_exact_primary'] for m in METHODS)
        for r in records)
    report['no_native_exact_hit'] = sum(
        not any(r[m+'_exact_primary'] for m in METHODS) for r in records)
    return report


def main():
    began = time.perf_counter()
    if OUT.exists():
        raise FileExistsError('Do not overwrite completed or failed diagnostic artifacts')
    paths = [ROOT/'outputs'/x for x in
             ('M28-A3-FIRST-FIT','M28-A3-DEV-INFER','M28-A3-EVALUATION')]
    if any((p/'STATUS').read_text().strip() != 'COMPLETED' for p in (*paths,TURN,NATIVE)):
        raise ValueError('All original inference/evaluation and turnover sources must be complete')
    fit, decoded, evaluated = [json.loads((p/'summary.json').read_text()) for p in paths]
    require_source('A3',fit,decoded,evaluated)
    turnover = json.loads((TURN/'summary.json').read_text())
    native = json.loads((NATIVE/'summary.json').read_text())
    if (turnover['source']['A3']['checkpoint'] != fit['final_checkpoint'] or
        turnover['reference_primary_chains'] != 7728 or turnover['test_setaria_used'] or
        not native['completed'] or native['reference_primary_chains'] != 7728 or
        native['test_setaria_used'] or set(native['pooled']) != set(METHODS)):
        raise ValueError('Source scope or checkpoint differs')
    with (TURN/'references.tsv').open() as handle:
        records = [parse_record(r) for r in csv.DictReader(handle,delimiter='\t')]
    indexed = {reference_key(r):r for r in records}
    if len(indexed) != len(records) or len(records) != 7728:
        raise ValueError('Frozen unique reference grain required')

    # Real reference/annotation processing occurs only inside the CPU Slurm job.
    import common_ruler as C
    from run_native_baseline import SCOPE, OUT as PRED
    for species, seqid, length in SCOPE:
        current = [r for r in records if r['species'] == species]
        ledger_path = ROOT/'outputs/M28-COMMON-RULER-R3/result'/(species+'.references.jsonl')
        with ledger_path.open() as handle:
            refs = [json.loads(line) for line in handle]
        keys = {reference_key(r) for r in refs}
        if (len(refs) != COUNTS[species] or len(keys) != COUNTS[species] or
            len(current) != COUNTS[species] or keys != {reference_key(r) for r in current}):
            raise ValueError('Complete primary ledger and turnover identities must agree')
        for ref in refs:
            row = indexed[reference_key(ref)]
            if ((row['id'],row['gene_id']) != (ref['id'],ref['gene_id']) or
                row['unsupported_strata'] != ';'.join(ref['M28_unsupported_strata_not_excluded'])):
                raise ValueError('Primary identity/unsupported stratum differs')
        for strand in ('+','-'):
            part = [r for r in current if r['strand'] == strand]
            stage_counts = {s:sum(r['A3_'+s] for r in part) for s in STAGES}
            official = evaluated['methods']['A3']['stages'][species+'|'+strand]['exact_reference_hits']
            if stage_counts != official:
                raise ValueError('A3 table must reproduce completed official stage accounting')
        for method in METHODS:
            folder = PRED/method/species
            if (folder/'STATUS').read_text().strip() != 'COMPLETED':
                raise ValueError('Native prediction is not terminal complete')
            prediction = ROOT/native['species'][species]['methods'][method]['raw_prediction_path']
            txs = C.S.E.primary_transcripts(
                C.S.E.parse_annotation(prediction,{seqid:length}),complete_only=False)
            predicted = {(species,)+C.S.key(t) for t in txs}
            observed = {'tp':len(predicted & keys),'predicted':len(predicted),'reference':len(keys)}
            official = native['species'][species]['methods'][method]['exact_chain']
            if observed != {k:official[k] for k in observed}:
                raise ValueError('Native exact keys differ from completed common-ruler metrics')
            for row in current:
                row[method+'_exact_primary'] = reference_key(row) in predicted

    report = {'completed':True,'experiment':OUT.name,
              'scope':'descriptive two-DEV reference membership, not a causal experiment or independent test',
              'reference_primary_chains':7728,'reference_filters_added':False,
              'test_setaria_used':False,'new_training':False,'new_model_inference':False,
              'key_definition':'species/seqid/strand/ordered 0-based half-open CDS coordinates; phase not part of exact key',
              'source':{'native_evaluation':str(NATIVE.relative_to(ROOT)),
                        'A3_turnover':str(TURN.relative_to(ROOT)),
                        'A3_evaluation':str(paths[2].relative_to(ROOT)),
                        'A3_checkpoint':fit['final_checkpoint']},
              'A3_predeclared_gate_unchanged':evaluated['predeclared_pilot_gate'],
              'pooled':summarize(records),
              'per_species':{s:summarize([r for r in records if r['species']==s]) for s in COUNTS},
              'exon_strata':{k:summarize([r for r in records if (r['exon_count']==1)==single])
                             for k,single in (('single_CDS_segment',True),('multiple_CDS_segments',False))},
              'limitations':[
                  'Only reference sensitivity/membership is decomposed; unmatched predictions and biological false genes are not classified.',
                  'A3-only hits are already annotated primary references, not newly discovered genes.',
                  'Pretrained exposure is comparator-dependent or unknown; no clean universal zero-shot or causal claim.',
                  'Stage states describe observed attrition; they do not prove unique loss causes.',
                  'No subset metric replaces the full 7728-chain ruler or frozen gates.']}
    if report['pooled']['A3_stage_hits']['selected_unique'] != evaluated['methods']['A3']['pooled']['exact_chain']['tp']:
        raise ValueError('A3 selected count differs')
    report['seconds'] = time.perf_counter()-began
    OUT.mkdir()
    for row in records:
        row.pop('CDS')
    with (OUT/'references.tsv').open('x',newline='') as handle:
        writer = csv.DictWriter(handle,fieldnames=list(records[0]),delimiter='\t')
        writer.writeheader()
        writer.writerows(records)
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False)
        handle.write('\n')
    print(json.dumps(report['pooled'],indent=2),flush=True)


if __name__ == '__main__':
    main()
