#!/usr/bin/env python3
"""Terminal FULL raw/GFF/reference/numeric receipt plus completed native comparison."""
import csv
import json
from collections import Counter, defaultdict
import evaluate_full_train as F
import evaluate_query_context as H
from src.m28.assembly import assemble_predictions
from src.m28.export import as_transcript

ROOT, I, E = F.I.ROOT, F.I, F.E
OUT = ROOT/'outputs/M28-FULL-TRAIN-RESULT-CHECK'
JOBS = {'Q-FIT':'13358083','C0-FIT':'13358084','Q-INFER':'13358085',
        'C0-INFER':'13358086','EVAL':'13358090'}


def main():
    for suffix,job in JOBS.items():
        folder = ROOT/'outputs'/('M28-FULL-TRAIN-'+suffix)
        assert (folder/'STATUS').read_text().strip() == 'COMPLETED'
        assert (folder/'JOBID').read_text().strip() == job
    report = H.read_report(F.OUT/'summary.json')
    assert report['protocol'] == F.R.PROTOCOL and set(report['methods']) == {'Q','C0'}
    lengths = {q:n for _,q,n,_ in E.S.SCOPE}
    predicted, assemblies, actual_windows = {},{},{}
    # Reconcile all actual raw records and both exported GFFs before opening labels.
    for arm in ('Q','C0'):
        fitted = json.loads((I.FIT[arm]/'summary.json').read_text())
        decoded = json.loads((I.INFER[arm]/'summary.json').read_text())
        F.validate_complete(fitted,decoded,arm)
        assert fitted == report['fit'][arm]
        predictions = defaultdict(list); scope_assemblies = {}; windows = 0
        for scope in decoded['scopes']:
            counts = Counter(); windows += len(scope['windows'])
            assembled = assemble_predictions(F.scope_records(I.INFER[arm]/scope['file'],scope,arm,counts),
                scope['windows'],scope['strand'])
            key = scope['species']+'|'+scope['strand']
            assert dict(counts) == report['methods'][arm]['candidate_diversity'][key]
            assert assembled['counts'] == scope['assembly_counts']
            for index,row in enumerate(assembled['predictions']):
                predictions[scope['species']].append(as_transcript(row['chain'],scope['seqid'],scope['strand'],
                    I.prediction_id(arm,scope['ordinal'],scope['strand'],index),row['score']))
            scope_assemblies[key] = assembled
        txs = [t for values in predictions.values() for t in values]
        assert windows == 8690 and len(txs) == decoded['primary_CDS_chain_predictions']
        assert len({E.S.key(t) for t in txs}) == len(txs)
        assert all(E.R7.phase_chain_consistent(t) for t in txs)
        for folder in (I.INFER[arm],F.OUT):
            parsed = E.S.E.primary_transcripts(E.S.E.parse_annotation(folder/(arm+'.gff3'),lengths),complete_only=False)
            F.A.validate_gff_roundtrip(parsed,txs)
        predicted[arm],assemblies[arm],actual_windows[arm] = predictions,scope_assemblies,windows
    refs,isoforms,parents,annotations,sequences,_ = F.A.reference_context()
    with (ROOT/'outputs/M28-NATIVE-A3-REFERENCE-DIAGNOSTIC/references.tsv').open() as handle:
        native_rows = list(csv.DictReader(handle,delimiter='\t'))
    key = lambda species,ident,seqid,strand,cds:(species,ident,seqid,strand,tuple(tuple(c[:2]) for c in cds))
    current_keys = {key(s,t['id'],t['seqid'],t['strand'],t['CDS']) for s,ts in refs.items() for t in ts}
    native_keys = {key(r['species'],r['id'],r['seqid'],r['strand'],json.loads(r['CDS_intervals'])) for r in native_rows}
    assert len(native_rows) == len(native_keys) == len(current_keys) == 7728 and current_keys == native_keys
    native = H.read_report(ROOT/'outputs/M28-NATIVE-BASELINE-EVAL/summary.json')
    assert native['completed'] is True and native['unsupported_reference_strata_retained'] is True
    methods = {}; stage_counts = {}
    for arm in ('Q','C0'):
        per_species = {}; stages = {}
        for scope_key,assembled in assemblies[arm].items():
            species,strand = scope_key.split('|')
            audit = F.stage_audit(assembled,[t for t in refs[species] if t['strand'] == strand])
            assert audit['exact_reference_hits'] == report['methods'][arm]['stages'][scope_key]['exact_reference_hits']
            stages[scope_key] = audit['exact_reference_hits']
        for species,seqid,length,_ in E.S.SCOPE:
            txs = predicted[arm][species]
            measured = E.score(txs,refs[species],isoforms[species],length,annotations[species])
            measured['prediction_nonexact_primary'] = measured.pop('reference_nonexact_primary')
            old = report['methods'][arm]['per_species'][species]
            assert all(measured[k] == old[k] for k in measured)
            pk,hk = {E.S.key(t) for t in txs},{E.S.key(t) for t in parents[species]}
            measured['historical_parent_exact_chain'] = E.S.prf(len(pk&hk),len(pk),len(hk))
            failures = Counter(); valid = 0
            for tx in txs:
                problems,_ = I.M.chain_checks(tx,sequences[species][seqid])
                failures.update(problems); valid += not problems
            structure = dict(fully_supported=valid,predictions=len(txs),failure_counts_nonexclusive=dict(failures))
            assert structure == old['CDS_structure_diagnostic']
            measured['CDS_structure_diagnostic'] = structure
            per_species[species] = measured
        pooled = F.aggregate(per_species)
        pooled['actual_free_primary_coverage'] = sum(v['free_unique'] for v in stages.values())/7728
        pooled['owner_primary_coverage'] = sum(v['owner_unique'] for v in stages.values())/7728
        assert pooled == report['methods'][arm]['pooled']
        methods[arm] = dict(per_species=per_species,pooled=pooled)
        stage_counts[arm] = {s:{stage:sum(v[stage] for k,v in stages.items() if k.startswith(s+'|'))
            for stage in ('free_unique','owner_unique','selected_unique')} for s in F.A.REF_COUNTS}
    gate = F.A.continuation_gate(methods['Q']['per_species'],methods['C0']['per_species'])
    assert gate == report['predeclared_gate']
    comparators = {m:dict(per_species={s:r['methods'][m] for s,r in native['species'].items()},
        pooled=native['pooled'][m]) for m in ('ANNEVO','Tiberius','Helixer')}
    for values in comparators.values():
        H.check_ruler(values,methods['Q'])
    result = dict(completed=True,reference_primary_chains=7728,actual_raw_windows=actual_windows,
        raw_assembly_inference_and_evaluation_GFF_identity_phase_verified=True,
        actual_reference_identity_matches_all_7728_native_rows=True,
        actual_GFF_metrics_background_and_stages_recomputed=True,methods=methods,stages=stage_counts,
        original_frozen_gate=gate,
        native_comparison={name:dict(per_species={s:v['exact_chain'] for s,v in values['per_species'].items()},
            pooled=values['pooled']['exact_chain']) for name,values in comparators.items()},
        contrasts={arm:{name:H.contrast(methods[arm],values) for name,values in comparators.items()}
            for arm in methods},native_exposure_limit=native['source_limitations'],
        optimizer_steps_this_run=0,new_model_inference=False,test_setaria_used=False,
        limitation='Repeated-use bounded DEV. Failed frozen gate remains failed; no independent generalization or paper PASS.')
    assert not (OUT/'summary.json').exists()
    OUT.mkdir(parents=True,exist_ok=True)
    with (OUT/'summary.json').open('x') as handle:
        json.dump(result,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps(result,allow_nan=False),flush=True)


if __name__ == '__main__':
    main()
