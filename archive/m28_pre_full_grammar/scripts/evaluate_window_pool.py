#!/usr/bin/env python3
"""Fixed A4 CPU readout of A/A2 immutable natural candidates; no new fit/inference."""
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
import common_ruler as E
import evaluate_query as A
import evaluate_query_collaboration as A2
from evaluate_pilot import aggregate
from diagnose_query_turnover import ARMS, SOURCES, require_source
from src.m28.window_pool import assemble_window_pool
from src.m28.export import as_transcript, write_gff3

OUT = E.ROOT/'outputs/M28-A4-WINDOW-POOL'
PROPOSAL = 'reports/M28-METHOD-RESTART/a4_window_pool_protocol.md'


def pool_audit(assembled, refs):
    reference = {tuple((a,b) for a,b,_ in ref['CDS']) for ref in refs}
    hits = {stage:len(set(chains)&reference) for stage,chains in assembled['stages'].items()}
    if not hits['selected_unique'] <= hits['window_pool_unique'] == hits['free_unique']:
        raise ValueError('Nested full window-pool reference accounting required')
    return {'reference':len(reference),'exact_reference_hits':hits,
            'not_in_actual_free_candidates':len(reference)-hits['free_unique'],
            'lost_at_window_pool':hits['free_unique']-hits['window_pool_unique'],
            'lost_at_conflict_or_nonpositive_gain':hits['window_pool_unique']-hits['selected_unique'],
            'candidate_counts':assembled['counts']}


def main():
    if OUT.exists():
        raise FileExistsError('Use the single new A4 output directory')
    source = {}
    for arm,names in SOURCES.items():
        paths = [E.ROOT/'outputs'/name for name in names]
        if any((p/'STATUS').read_text().strip() != 'COMPLETED' for p in paths):
            raise ValueError('Completed matching immutable A/A2 sources required')
        source[arm] = tuple(json.loads((p/'summary.json').read_text()) for p in paths)
        require_source(arm,*source[arm])
    A.validate_complete(*source['A'][:2]); A2.validate_complete(*source['A2'][:2])
    OUT.mkdir()
    lengths = {seqid:length for _,seqid,length,_ in E.S.SCOPE}
    assembled_scopes = {}; predictions = {}; gffs = {}
    # Both source pools are assembled/exported before any reference access.
    for arm in ARMS:
        predictions[arm] = defaultdict(list); assembled_scopes[arm] = {}
        inference = E.ROOT/'outputs'/SOURCES[arm][1]
        for scope in source[arm][1]['scopes']:
            species,seqid,strand = scope['species'],scope['seqid'],scope['strand']
            key = species+'|'+strand
            assembled = assemble_window_pool(A.query_scope_records(inference/scope['file'],scope),
                                              scope['windows'],strand)
            original = source[arm][2]['methods'][arm]['stages'][key]['candidate_counts']
            if any(assembled['counts'][k] != original[k] for k in ('free_records','free_unique')):
                raise ValueError('A4 source pool changed')
            assembled_scopes[arm][key] = (scope,assembled)
            for i,row in enumerate(assembled['predictions']):
                identifier = 'A4_'+arm+'_'+seqid+'_'+('p' if strand == '+' else 'm')+'_'+str(i)
                predictions[arm][species].append(as_transcript(row['chain'],seqid,strand,identifier,row['score']))
        gff = OUT/(arm+'_window_pool.gff3'); gffs[arm] = str(gff.relative_to(E.ROOT))
        transcripts = [tx for values in predictions[arm].values() for tx in values]
        write_gff3(gff,transcripts,lengths,'M28_A4_'+arm)
        parsed = E.S.E.primary_transcripts(E.S.E.parse_annotation(gff,lengths),complete_only=False)
        A.validate_gff_roundtrip(parsed,transcripts)

    refs,isoforms,parents,annotations,sequences,_ = A.reference_context()
    controls = json.loads((E.ROOT/'reports/M28-METHOD-RESTART/r5_result.json').read_text())['methods']
    reference_turnover = E.ROOT/'outputs/M28-A-A2-REFERENCE-TURNOVER'
    if (reference_turnover/'STATUS').read_text().strip() != 'COMPLETED':
        raise ValueError('Completed descriptive reference turnover required for attrition accounting')
    with (reference_turnover/'references.tsv').open() as f:
        turnover = list(csv.DictReader(f,delimiter='\t'))
    if len(turnover) != 7728:
        raise ValueError('Full primary-reference turnover required')
    report = {'experiment':'M28-A4-WINDOW-POOL','proposal':PROPOSAL,
              'status':'complete_fixed_DEV_readout_not_independent_generalization',
              'reference_primary_chains':7728,'historical_parent_chains':6450,
              'reference_used_in_candidate_pooling':False,'reference_filters_added':False,
              'test_setaria_used':False,'new_model_inference':False,'optimizer_steps_this_run':0,
              'pool_policy':'max_existing_gain_per_global_chain_across_all_generating_windows',
              'tau':0.,'selector':'original_positive_gain_same_strand_nonoverlapping',
              'methods':{},'predeclared_gate':{},
              'scope_limit':'ownership-policy readout; not a new model, latest native comparator or paper result'}
    for arm in ARMS:
        stages = {}
        for key,(scope,assembled) in assembled_scopes[arm].items():
            species,strand = scope['species'],scope['strand']
            audit = pool_audit(assembled,[t for t in refs[species] if t['strand'] == strand])
            original = source[arm][2]['methods'][arm]['stages'][key]
            if audit['exact_reference_hits']['free_unique'] != original['exact_reference_hits']['free_unique']:
                raise ValueError('A4 free primary-reference pool changed')
            audit.update(species=species,seqid=scope['seqid'],strand=strand)
            stages[key] = audit
        per_species = {}
        for species,seqid,length,_ in E.S.SCOPE:
            pred = predictions[arm][species]
            measured = E.score(pred,refs[species],isoforms[species],length,annotations[species])
            pk,hk = {E.S.key(t) for t in pred},{E.S.key(t) for t in parents[species]}
            measured['historical_parent_exact_chain'] = E.S.prf(len(pk&hk),len(pk),len(hk))
            failures = Counter(); valid = 0
            for tx in pred:
                problems,_ = A.Q.M.chain_checks(tx,sequences[species][seqid])
                failures.update(problems); valid += not problems
            measured['CDS_structure_diagnostic'] = {'fully_supported':valid,'predictions':len(pred),
                                                   'failure_counts_nonexclusive':dict(failures)}
            per_species[species] = measured
        pooled = aggregate(per_species)
        if sum(s['exact_reference_hits']['selected_unique'] for s in stages.values()) != pooled['exact_chain']['tp']:
            raise ValueError('Selected primary hits and final TP disagree')
        selected = {(t['seqid'],t['strand'],tuple((a,b) for a,b,_ in t['CDS']))
                    for values in predictions[arm].values() for t in values}
        transitions = Counter()
        for row in turnover:
            chain = tuple(tuple(x) for x in json.loads(row['CDS_intervals']))
            hit = (row['seqid'],row['strand'],chain) in selected
            transitions[(row[arm+'_state'],'selected' if hit else 'not_selected')] += 1
        if sum(n for (old,new),n in transitions.items() if new == 'selected') != pooled['exact_chain']['tp']:
            raise ValueError('Reference turnover and final TP disagree')
        report['methods'][arm+'_window_pool'] = {
            'checkpoint':source[arm][0]['final_checkpoint'],'source_inference':SOURCES[arm][1],
            'per_species':per_species,'pooled':pooled,'stages':stages,'GFF3':gffs[arm],
            'fixed_owner_control':source[arm][2]['methods'][arm],
            'old_state_to_A4_selection':[{'old_state':old,'A4_state':new,'reference_chains':n}
                                         for (old,new),n in sorted(transitions.items())]}
        report['predeclared_gate'][arm+'_window_pool'] = A2.continuation_gate(
            per_species,source['A'][2]['methods']['A']['per_species'],controls['C0']['per_species'])
    with (OUT/'summary.json').open('x') as f:
        json.dump(report,f,indent=2,allow_nan=False); f.write('\n')
    print(json.dumps({'gate':report['predeclared_gate'],
                      'exact_chain':{k:v['pooled']['exact_chain'] for k,v in report['methods'].items()}},
                     allow_nan=False),flush=True)


if __name__ == '__main__':
    main()
