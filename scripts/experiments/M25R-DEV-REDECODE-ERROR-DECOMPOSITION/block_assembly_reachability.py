#!/usr/bin/env python3
"""Oracle membership in a run-preserving candidate graph; no predictions."""
import csv
import json
from bisect import bisect_right
from collections import Counter, defaultdict
from pathlib import Path
import posthoc_r4_case_correction as C

D = C.diag
ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT/'outputs/M25R-E1-BLOCK-ASSEMBLY-REACHABILITY'

def choose_ordered(options):
    """Smallest increasing run-ID witness; None if no such path exists."""
    witness, previous = [], -1
    for choices in options:
        i = bisect_right(choices, previous)
        if i == len(choices):
            return None
        previous = choices[i]
        witness.append(previous)
    return witness

def run_index(states, sequence, traces):
    index = defaultdict(set)
    carriers = []
    for trace in sorted(traces, key=lambda t: t['block_span']):
        for a,b in trace['runs']:
            rid = len(carriers)
            carriers.append((trace['lineage_id'], a, b))
            for position in (a,b):
                if not 0 < position < len(states):
                    continue
                pair = (int(states[position-1]), int(states[position]))
                for event in D.m25.TRANSITIONS.get(pair, ()):
                    anchor = D.event_anchor(position,event)
                    for motif in D.motif_positions(sequence,event,anchor):
                        index[(event,motif)].add(rid)
    return index,carriers

def inspect(ref, index, carriers):
    options=[]
    for i,(a,b) in enumerate(ref['cds']):
        left=('start',a) if i==0 else ('acceptor',a-2)
        right=('stop',b-3) if i==len(ref['cds'])-1 else ('donor',b)
        options.append(sorted(index.get(left,set()) & index.get(right,set())))
    edges_present=all(index.get((event,p)) for event,ps in ref['events'].items() for p in ps)
    witness=choose_ordered(options) if all(options) else None
    cds=ref['cds']
    geometry=all(a<b for a,b in cds) and all(b+4<=c for (_,b),(c,_) in zip(cds,cds[1:]))
    if not edges_present:
        category='missing_boundary_candidate'
    elif not all(options):
        category='boundaries_present_no_same_run_fragment'
    elif witness is None or not geometry:
        category='fragments_present_chain_incompatible'
    else:
        category='full_chain_in_graph'
    return {'category':category,'all_event_positions_present':bool(edges_present),
            'fragment_carrier_options':options,'ordered_witness_run_ids':witness,
            'witness_distinct_blocks':len({carriers[r][0] for r in witness}) if witness else None,
            'skipped_run_carriers_in_witness_span':witness[-1]-witness[0]+1-len(witness) if witness else None,
            'canonical_ORF':ref['canonical'],
            'structurally_valid_chain_in_graph':category=='full_chain_in_graph' and ref['canonical']}

def main():
    # Focused order/geometry algorithm checks precede real data within allocation.
    assert choose_ordered([[1,4],[2,5],[3,6]])==[1,2,3]
    assert choose_ordered([[2],[1,2]]) is None
    assert choose_ordered([[1,2],[2],[3]])==[1,2,3]
    if (OUT/'summary.json').exists() or (OUT/'references.jsonl').exists():
        raise FileExistsError('Do not overwrite diagnostic payload')
    config=C.yaml.safe_load((ROOT/'configs/M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0.yaml').read_text())
    species=D.load_species(ROOT,config)
    refs,lengths=D.validation_truth(species)
    assert set(refs)==C.VALIDATION
    records=D.build_reference_records(species,refs)
    assert len(records)==6450
    original={}
    for line in (ROOT/'outputs/M25R-R4-GROUP-PHASE-RETROSPECTIVE/epoch_1/reference_decomposition.jsonl').open():
        r=json.loads(line)
        original[tuple(r[k] for k in ('species','seqid','strand','transcript_id'))]=r
    corrected={}
    with (ROOT/'outputs/M25R-DEV-REDECODE-ERROR-DECOMPOSITION-R4-CASE-CORRECTION/epoch_1/reference_case_correction.tsv').open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            corrected[tuple(r[k] for k in ('species','seqid','strand','transcript_id'))]=r
    assert set(original)==set(corrected)=={r['key'] for r in records}
    traces=defaultdict(list)
    for line in (ROOT/'outputs/M25R-DEV-REDECODE-ERROR-DECOMPOSITION-R4/epoch_1/candidate_lineages.jsonl').open():
        r=json.loads(line)
        traces[(r['species'],r['seqid'],r['strand'])].append(r)
    grouped=defaultdict(list)
    for r in records:
        grouped[r['key'][:3]].append(r)
    result=[]
    for key,group in grouped.items():
        seq=species[key[0]]['seqs'][key[1]]
        if key[2]=='-': seq=D.m25.reverse_complement(seq)
        states=C.reconstruct_states(len(seq),traces[key])
        index,carriers=run_index(states,seq,traces[key])
        for ref in group:
            r=inspect(ref,index,carriers)
            old=original[ref['key']]
            if old['category']=='recovered_exact' and not r['structurally_valid_chain_in_graph']:
                raise AssertionError(f"existing exact chain unreachable: {ref['key']}")
            r.update(species=key[0],seqid=key[1],strand=key[2],transcript_id=ref['transcript_id'],
                     original_category=old['category'],gene_span_bp=ref['span'][1]-ref['span'][0],
                     CDS_count=len(ref['cds']),old_motif_reachable=bool(int(corrected[ref['key']]['motif_reachable'])))
            # Event membership may be looser than legacy one-to-one anchor matching.
            # Keep both, never silently rename the legacy metric.
            result.append(r)
        del states,index,carriers
        print(f'completed {key}',flush=True)
    assert sum(r['original_category']=='recovered_exact' for r in result)==1389
    assert sum(r['original_category']=='CDS_support_fragmented_across_blocks' for r in result)==1634
    assert sum(r['old_motif_reachable'] for r in result)==3616
    summary={}
    for name,rows in [('all',result),('fragmented',[r for r in result if r['original_category']=='CDS_support_fragmented_across_blocks'])]:
        valid=sum(r['structurally_valid_chain_in_graph'] for r in rows)
        summary[name]={'n':len(rows),'categories':dict(Counter(r['category'] for r in rows)),
                       'canonical_chain_in_graph':valid,
                       'cross_block_witness':sum((r['witness_distinct_blocks'] or 0)>1 for r in rows),
                       'skipping_witness':sum((r['skipped_run_carriers_in_witness_span'] or 0)>0 for r in rows)}
    # Even perfect precision cannot exceed this F1 for the defined graph.
    maximum=summary['all']['canonical_chain_in_graph']
    summary.update(exact_positive_checks=1389,reference_total=6450,
                   optimistic_chain_F1=2*maximum/(6450+maximum),chain_gate=0.55,
                   scope='epoch1 run-preserving cross-block graph, score and phase filters relaxed',
                   actual_predictions=False,model_loaded=False,setaria_access=False)
    with (OUT/'references.jsonl').open('x') as f:
        for r in result: f.write(json.dumps(r)+'\n')
    with (OUT/'summary.json').open('x') as f: json.dump(summary,f,indent=2)
    (OUT/'STATUS').write_text('COMPLETED\n')
    print(json.dumps(summary),flush=True)

if __name__=='__main__': main()
