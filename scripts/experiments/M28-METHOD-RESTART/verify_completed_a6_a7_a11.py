#!/usr/bin/env python3
"""One terminal numerical receipt, not a new experiment or model selection."""
import csv
import json
import math
from collections import Counter
from pathlib import Path
import evaluate_chain_quality as F
from smoke_query_beam import check_chain

ROOT = F.Q.ROOT
OUT = ROOT/'outputs/M28-A6-A7-A11-RESULT-CHECK'


def read(folder):
    return json.loads((ROOT/'outputs'/folder/'summary.json').read_text())


def rows(path):
    with path.open() as handle:
        return [json.loads(line) for line in handle]


def chain(parts):
    return tuple(tuple(p) for p in parts)


def ref_key(row):
    return row['species'], row['id'], row['seqid'], row['strand'], chain(json.loads(row['CDS_intervals'])), row['unsupported_strata']


def main():
    if OUT.exists():
        raise FileExistsError('Keep any prior receipt intact')
    OUT.mkdir()
    a6 = read('M28-A6-CHAIN-QUALITY-EVAL')
    a7 = read('M28-A7-QUALITY-WINDOW-POOL')
    native_path = ROOT/'outputs/M28-NATIVE-A3-REFERENCE-DIAGNOSTIC/references.tsv'
    with native_path.open() as handle:
        native = list(csv.DictReader(handle, delimiter='\t'))
    with (ROOT/'outputs/M28-A7-QUALITY-WINDOW-POOL/references.tsv').open() as handle:
        paired = list(csv.DictReader(handle, delimiter='\t'))
    assert len(paired) == len(native) == 7728
    assert len({ref_key(r) for r in paired}) == 7728
    assert {ref_key(r) for r in paired} == {ref_key(r) for r in native}
    lengths = {q:n for _,q,n,_ in F.A.E.S.SCOPE}
    txs, geometry = {}, {}
    gffs = {'A6':'M28-A6-CHAIN-QUALITY-INFER/A6.gff3', 'A7':'M28-A7-QUALITY-WINDOW-POOL/A7.gff3'}
    for arm,path in gffs.items():
        txs[arm] = F.A.E.S.E.primary_transcripts(F.A.E.S.E.parse_annotation(ROOT/'outputs'/path,lengths),complete_only=False)
        geometry[arm] = {F.A.E.S.key(t) for t in txs[arm]}
        assert len(geometry[arm]) == len(txs[arm])
        assert all(F.A.E.R7.phase_chain_consistent(t) for t in txs[arm])
    membership, checked_metrics = {}, {}
    for species,count in F.A.REF_COUNTS.items():
        part = [r for r in paired if r['species']==species]
        assert len(part) == count
        state = Counter()
        measured = {}
        for row in part:
            key = row['seqid'],row['strand'],chain(json.loads(row['CDS_intervals']))
            a,b = key in geometry['A6'], key in geometry['A7']
            assert row['A6_selected'] == str(a) and row['A7_selected'] == str(b)
            state['both' if a and b else 'A7_only' if b else 'A6_only' if a else 'neither'] += 1
        for arm,result in [('A6',a6),('A7',a7)]:
            value = result['methods'][arm]['per_species'][species]['exact_chain']
            tp = state['both'] + state[arm+'_only']
            pred = sum(t['seqid']==part[0]['seqid'] for t in txs[arm])
            assert (value['tp'],value['predicted'],value['reference']) == (tp,pred,count)
            assert math.isclose(value['f1'],2*tp/(pred+count),rel_tol=0,abs_tol=1e-14)
            measured[arm] = value
        membership[species] = {k:state[k] for k in ('both','A7_only','A6_only','neither')}
        assert membership[species] == a7['paired_reference_membership']['per_species'][species]
        checked_metrics[species] = measured
    a = read('M28-A-EVALUATION')['methods']['A']['per_species']
    c0 = json.loads((ROOT/'reports/M28-METHOD-RESTART/r5_result.json').read_text())['methods']['C0']['per_species']
    a6_gate = F.continuation_gate(a6['methods']['A6']['per_species'],a,c0)
    a7_gate = F.A.continuation_gate(a7['methods']['A7']['per_species'],a6['methods']['A6']['per_species'])
    assert a6_gate == a6['predeclared_pilot_gate']
    assert a7_gate == a7['predeclared_matched_owner_gate']
    s = read('M28-A11-COMPONENT-MEAN-TRAIN-PROBE')
    old = {r['window']['id']:r for r in rows(ROOT/'outputs/M28-A5-TRAINED-TRAIN-DIAGNOSTIC/predictions.jsonl')}
    predictions = rows(ROOT/'outputs/M28-A11-COMPONENT-MEAN-TRAIN-PROBE/predictions.jsonl')
    refs = rows(ROOT/'outputs/M28-A11-COMPONENT-MEAN-TRAIN-PROBE/references.jsonl')
    original_refs = rows(ROOT/'outputs/M28-A5-TRAINED-TRAIN-DIAGNOSTIC/references.jsonl')
    identity = lambda r:(r['window_id'],r['reference_id'],chain(r['chain']))
    assert len(predictions)==8 and {r['window']['id'] for r in predictions}==set(old)
    assert len(refs)==27 and {identity(r) for r in refs}=={identity(r) for r in original_refs}
    pools, delta, completed_lengths = {}, 0., {arm:Counter() for arm in s['arms'] if arm!='ORIGINAL'}
    for pred in predictions:
        ident = pred['window']['id']; source = old[ident]
        greedy = {chain(r['chain']) for r in source['greedy_queries'] if r['terminated']}
        pools[ident] = {'ORIGINAL':greedy|{chain(r['chain']) for q in source['beam_queries'] for r in q['completed']}}
        for raw,saved in zip(pred['arms']['RAW_SUM'],source['beam_queries']):
            x = {chain(r['chain']):r['path_log_score'] for r in raw['completed']}
            y = {chain(r['chain']):r['path_log_score'] for r in saved['completed']}
            assert x.keys()==y.keys()
            delta = max(delta,max((abs(x[k]-y[k]) for k in x),default=0.))
        assert delta<=1e-3
        for left,right in zip(pred['arms']['COMPONENT_SUM'],pred['arms']['COMPONENT_MEAN']):
            assert left['counts']==right['counts']
        for arm,queries in pred['arms'].items():
            assert [q['query_index'] for q in queries]==list(range(16))
            pool = set()
            for q in queries:
                assert len(q['completed'])<=4 and q['counts']['score_bound_stopped_live']==0
                for r in q['completed']:
                    check_chain(r['chain'],source['sequence'].upper())
                    score = r['path_log_score'] if arm=='RAW_SUM' else r['boundary_log_sum']+2*r['termination_log_sum']
                    if arm=='COMPONENT_MEAN': score/=2*len(r['chain'])
                    assert math.isclose(score,r['completed_score'],rel_tol=0,abs_tol=1e-10)
                    pool.add(chain(r['chain'])); completed_lengths[arm][len(r['chain'])]+=1
            pools[ident][arm] = pool|greedy
    for r in refs:
        assert r['hits']=={arm:chain(r['chain']) in pools[r['window_id']][arm] for arm in s['arms']}
    counts = {arm:sum(r['hits'][arm] for r in refs) for arm in s['arms']}
    assert counts==s['counts']
    assert {species:{arm:sum(r['hits'][arm] for r in refs if r['species']==species) for arm in s['arms']} for species in F.A.REF_COUNTS}==s['per_species']
    new = [r for r in refs if r['hits']['COMPONENT_MEAN'] and not r['hits']['COMPONENT_SUM']]
    lost = [r for r in refs if r['hits']['COMPONENT_SUM'] and not r['hits']['COMPONENT_MEAN']]
    assert (len(new),len(lost),sum(r['exons']>=5 for r in new))==(s['mean_vs_component_sum_new'],s['mean_vs_component_sum_lost'],s['new_five_or_more_exon_pairs'])
    for r in s['candidate_records']:
        assert all(len(pools[r['window_id']][arm])==r['arms'][arm]['unique_geometry'] for arm in s['arms'])
    receipt = dict(completed=True,validation='terminal_GFF_reference_membership_and_A11_actual_pool_reconciliation',
        reference_primary_chains=7728,reference_geometry_and_unsupported_equal_accepted_native=True,
        GFF_counts={arm:len(value) for arm,value in txs.items()},GFF_all_phase_consistent=True,
        per_species_metrics=checked_metrics,paired_reference_membership=membership,
        A6_frozen_gate=a6_gate,A7_frozen_gate=a7_gate,A11_counts=counts,
        A11_original27_identity_preserved=True,A11_actual_pool_hits_equal_receipt=True,
        A11_common_raw_score_max_abs_difference=delta,A11_weighted_search_shared=True,
        A11_completed_segment_histogram=completed_lengths,
        A11_pool_unique_geometry_totals={arm:sum(len(p[arm]) for p in pools.values()) for arm in s['arms']},
        optimizer_steps_this_run=0,new_model_inference=False,test_setaria_used=False)
    with (OUT/'summary.json').open('x') as handle: json.dump(receipt,handle,indent=2,allow_nan=False)
    print(json.dumps(receipt,allow_nan=False),flush=True)


if __name__=='__main__': main()
