#!/usr/bin/env python3
"""Fixed A/A2 or A2/A3 reference turnover; no training, tuning or new inference."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from diagnose_b2_recovery import chain_key, recovery_state

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT/'outputs/M28-A-A2-REFERENCE-TURNOVER'
ARMS = ('A','A2')
STAGES = ('free_unique','owner_unique','selected_unique')
SOURCE_DIRS = {'A':('M28-A-FIRST-FIT','M28-A-DEV-INFER','M28-A-EVALUATION'),
               'A2':('M28-A2-FIRST-FIT','M28-A2-DEV-INFER-RECOVERY','M28-A2-EVALUATION'),
               'A3':('M28-A3-FIRST-FIT','M28-A3-DEV-INFER','M28-A3-EVALUATION')}
SOURCES = {arm:SOURCE_DIRS[arm] for arm in ARMS}


def require_source(arm, fit, decoded, evaluated):
    if (fit.get('arm') != arm or not fit.get('completed') or fit.get('optimizer_steps') != 4608 or
        decoded.get('arm') != arm or not decoded.get('inference_complete') or
        decoded.get('windows') != 8690 or decoded.get('reference_used') is not False or
        decoded.get('test_setaria_used') is not False or
        decoded.get('checkpoint') != fit.get('final_checkpoint') or
        decoded.get('configuration') != fit.get('configuration') or
        evaluated.get('fit',{}).get('final_checkpoint') != fit.get('final_checkpoint') or
        evaluated.get('reference_primary_chains') != 7728 or
        evaluated.get('test_setaria_used') is not False or arm not in evaluated.get('methods',{})):
        raise ValueError('Complete matching fixed '+arm+' source is required')


def reference_row(ref, stage_sets, owner_gains):
    chain = chain_key(ref['CDS'])
    row = {k:ref[k] for k in ('species','id','gene_id','seqid','strand')}
    row.update(exon_count=len(chain),span_bp=chain[-1][1]-chain[0][0],
               coding_bp=sum(b-a for a,b in chain),
               CDS_intervals=json.dumps(chain,separators=(',',':')),
               unsupported_strata=';'.join(ref.get('M28_unsupported_strata_not_excluded',[])))
    for arm in ARMS:
        flags = [chain in stage_sets[arm][stage] for stage in STAGES]
        if (flags[1] and not flags[0]) or (flags[2] and not flags[1]):
            raise ValueError('Reference stages must be nested')
        gain = owner_gains[arm].get(chain)
        for stage,flag in zip(STAGES,flags):
            row[arm+'_'+stage] = flag
        row[arm+'_owner_gain'] = gain
        row[arm+'_state'] = recovery_state(*flags,gain)
    return row


def transition_summary(records):
    before,after = ARMS
    report = {'reference_chains':len(records),'stages':{},'states':{}}
    for stage in STAGES:
        c = Counter()
        for row in records:
            old,new = row[before+'_'+stage],row[after+'_'+stage]
            c['both' if old and new else after+'_only' if new else before+'_only' if old else 'neither'] += 1
        report['stages'][stage] = {k:c[k] for k in ('both',after+'_only',before+'_only','neither')}
    for arm in ARMS:
        report['states'][arm] = dict(Counter(row[arm+'_state'] for row in records))
        # An optimistic reference-aware bound, not an executable predictor.
        hits = sum(row[arm+'_owner_unique'] for row in records)
        free = sum(row[arm+'_free_unique'] for row in records)
        report.setdefault('reference_aware_F1_upper_bound',{})[arm] = {
            'owner_hits':hits,'free_hits':free,
            'owner_F1_max_if_all_available_TP_selected_and_all_other_predictions_rejected':
                2*hits/(hits+len(records)) if records else 0.,
            'free_F1_max_ignoring_owner':2*free/(free+len(records)) if records else 0.,
            'ignores_reference_overlap_conflicts_and_is_not_model_performance':True}
    pairs = Counter((row[before+'_state'],row[after+'_state']) for row in records)
    report[before+'_to_'+after+'_states'] = [{before:a,after:b,'reference_chains':count}
                                          for (a,b),count in sorted(pairs.items())]
    return report


def main(contrast='A_A2'):
    global ARMS,SOURCES,OUT
    ARMS = {'A_A2':('A','A2'),'A2_A3':('A2','A3')}[contrast]
    SOURCES = {arm:SOURCE_DIRS[arm] for arm in ARMS}
    OUT = ROOT/('outputs/M28-'+'-'.join(ARMS)+'-REFERENCE-TURNOVER')
    if OUT.exists():
        raise FileExistsError('Use the single new turnover output directory')
    source = {}
    for arm,names in SOURCES.items():
        paths = [ROOT/'outputs'/name for name in names]
        if any((p/'STATUS').read_text().strip() != 'COMPLETED' for p in paths):
            raise ValueError('All fixed fit/inference/evaluation sources must be terminal complete')
        source[arm] = tuple(json.loads((p/'summary.json').read_text()) for p in paths)
        require_source(arm,*source[arm])

    # Heavy/project imports and real reference processing occur only in Slurm.
    import evaluate_query as A
    import evaluate_query_collaboration as A2
    import evaluate_query_known_mean as A3
    from evaluate_pilot import rows
    from src.m28.assembly import assemble_predictions, nearest_owner
    from src.m28.c0_decode import genomic_chain
    modules = {'A':A,'A2':A2,'A3':A3}
    for arm in ARMS:
        modules[arm].validate_complete(*source[arm][:2])
    scopes = {arm:{(r['species'],r['seqid'],r['strand']):r for r in source[arm][1]['scopes']}
              for arm in ARMS}
    if scopes[ARMS[0]].keys() != scopes[ARMS[1]].keys():
        raise ValueError('Fixed scopes differ')
    references = {s:rows(ROOT/'outputs/M28-COMMON-RULER-R3/result'/(s+'.references.jsonl'))
                  for s in A.REF_COUNTS}
    for species,count in A.REF_COUNTS.items():
        refs = references[species]
        if len(refs) != count or len({r['id'] for r in refs}) != count or len({(r['seqid'],r['strand'],chain_key(r['CDS'])) for r in refs}) != count:
            raise ValueError('Full frozen unique primary-reference ledger required')
    records = []
    for species,seqid,strand in sorted(scopes[ARMS[0]]):
        refs = [r for r in references[species] if (r['seqid'],r['strand']) == (seqid,strand)]
        reference = {chain_key(r['CDS']) for r in refs}
        sets,gains = {},{}
        for arm in ARMS:
            scope = scopes[arm][species,seqid,strand]
            windows,centers = scope['windows'],[a+b for a,b in scope['windows']]
            gains[arm] = {}
            inference = ROOT/'outputs'/SOURCES[arm][1]
            def observed():
                for row in A.query_scope_records(inference/scope['file'],scope):
                    i = row['window_index']; a,b = windows[i]
                    for chain,score in zip(row['chains'],row['score']):
                        key = tuple(genomic_chain(chain,a,b,strand))
                        if key in reference and nearest_owner(key,windows,centers) == i:
                            gains[arm][key] = max(float(score),gains[arm].get(key,float('-inf')))
                    yield row
            assembled = assemble_predictions(observed(),windows,strand)
            sets[arm] = {k:{chain_key(c) for c in assembled['stages'][k]}&reference for k in STAGES}
            official = source[arm][2]['methods'][arm]['stages'][species+'|'+strand]
            if ({k:len(sets[arm][k]) for k in STAGES} != official['exact_reference_hits'] or
                assembled['counts'] != official['candidate_counts']):
                raise ValueError('Turnover replay differs from completed official accounting')
            del assembled
        records.extend(reference_row(ref,sets,gains) for ref in refs)
    if len(records) != 7728:
        raise ValueError('Do not drop primary references')
    report = {'experiment':OUT.name,
              'scope':'descriptive completed DEV turnover; not a causal proof, new gate or independent test',
              'source':{arm:{'checkpoint':source[arm][0]['final_checkpoint'],
                             'inference':SOURCES[arm][1],'evaluation':SOURCES[arm][2],
                             'predeclared_gate_unchanged':source[arm][2]['predeclared_pilot_gate']}
                        for arm in ARMS},
              'reference_primary_chains':7728,'reference_filters_added':False,
              'new_model_inference':False,'test_setaria_used':False,
              'pooled':transition_summary(records),
              'per_species':{s:transition_summary([r for r in records if r['species'] == s])
                             for s in A.REF_COUNTS}}
    for arm in ARMS:
        final = report['pooled']['states'][arm].get('selected',0)
        if final != source[arm][2]['methods'][arm]['pooled']['exact_chain']['tp']:
            raise ValueError('Selected references and official final TP disagree')
    OUT.mkdir()
    with (OUT/'references.tsv').open('x',newline='') as f:
        writer = csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t')
        writer.writeheader(); writer.writerows(records)
    with (OUT/'summary.json').open('x') as f:
        json.dump(report,f,indent=2,allow_nan=False); f.write('\n')
    print(json.dumps(report['pooled'],indent=2,allow_nan=False),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contrast',choices=('A_A2','A2_A3'),default='A_A2')
    main(parser.parse_args().contrast)
