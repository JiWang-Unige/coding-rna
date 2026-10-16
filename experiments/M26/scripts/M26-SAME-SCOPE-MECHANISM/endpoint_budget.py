#!/usr/bin/env python3
"""R8 fixed endpoint-budget diagnostic; no prediction, fit, or reference-guided candidates."""
import argparse
import csv
import gzip
import json
import math
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path
import numpy as np
import support_partition as P

S, G = P.S, P.G
EVENTS = S.M.BOUNDARY_NAMES
CUTOFF = {e: math.log(t/(1-t)) for e,t in
          dict(start=.5,stop=.5,donor=.1,acceptor=.1).items()}
ARMS = ('original', 'unlocked', 'ranked', 'permuted')
WINDOW = 6144
SEED = 20260926


def motifs(sequence, event):
    encoded = np.frombuffer(sequence.encode('ascii'), dtype=np.uint8)
    width = len(next(iter(S.M.MOTIFS[event])))
    hits = np.zeros(len(sequence)-width+1, dtype=bool)
    for motif in S.M.MOTIFS[event]:
        hit = np.ones(len(hits), dtype=bool)
        for i,c in enumerate(motif):
            hit &= encoded[i:i+len(hits)] == ord(c)
        hits |= hit
    return np.flatnonzero(hits)


def permute_windows(positions, values, rng, window=WINDOW):
    shuffled = values.copy()
    cuts = np.r_[0, np.flatnonzero(np.diff(positions//window))+1, len(positions)]
    for lo,hi in zip(cuts,cuts[1:]):
        shuffled[lo:hi] = rng.permutation(values[lo:hi])
    return shuffled


def select(positions, values, count, cutoff):
    eligible = np.flatnonzero(values >= cutoff)
    assert len(eligible) >= count, 'Original endpoint budget exceeds eligible global motif pool'
    order = np.lexsort((positions[eligible], -values[eligible].astype(np.float64)))
    return np.sort(positions[eligible[order[:count]]])


def endpoints(exons):
    found = {e:set() for e in EVENTS}
    for x in exons:
        found[x.left].add(x.a-(2 if x.left=='acceptor' else 0))
        found[x.right].add(x.b-(3 if x.right=='stop' else 0))
    return {e:np.array(sorted(v),dtype=np.int64) for e,v in found.items()}


def exists(values, p):
    i = np.searchsorted(values,p)
    return i < len(values) and values[i] == p


def fragments(chain):
    return [(a,b,'start' if i==0 else 'acceptor',
             'stop' if i==len(chain)-1 else 'donor') for i,(a,b) in enumerate(chain)]


def events_for(chain):
    return [(event,p) for a,b,left,right in fragments(chain)
            for event,p in ((left,a-(2 if left=='acceptor' else 0)),
                            (right,b-(3 if right=='stop' else 0)))]


def grammar(chain, orf):
    if sum(b-a for a,b in chain)<6 or any(b+4>c for (_,b),(c,_) in zip(chain,chain[1:])):
        return False
    state = ''
    for i,(a,b,left,right) in enumerate(fragments(chain)):
        state = orf.advance(state,G.Exon(i,a,b,left,right,0.))
        if state is None:
            return False
    return state == ''


def viable(exon, orf):
    return any(orf.advance(state,exon) is not None
               for state in (('',) if exon.left=='start' else G.AUTOMATON))


def upper_end(a, left, orf):
    # Any accepted fragment must stop before completing the first internal stop.
    # Three possible incoming phases cover all acceptor suffix/INIT states.
    limits = []
    for offset in ((0,) if left=='start' else (0,1,2)):
        p = a+offset
        positions = orf.stops[p%3]
        i = bisect_left(positions,p)
        limits.append(int(positions[i])+3 if i<len(positions) else len(orf.sequence))
    return max(limits)


def flat_fragments(points, orf):
    for left in ('start','acceptor'):
        for p in points[left]:
            a = int(p)+(2 if left=='acceptor' else 0)
            upper = upper_end(a,left,orf)
            for right in ('donor','stop'):
                offset = 3 if right=='stop' else 0
                ends = points[right]
                lo = np.searchsorted(ends,a-offset,side='right')
                hi = np.searchsorted(ends,upper-offset,side='right')
                for q in ends[lo:hi]:
                    e = G.Exon(0,a,int(q)+offset,left,right,0.)
                    if viable(e,orf):
                        yield e


def opportunity_counts(points, orf, original=None):
    # Fragment counts are exact local-ORF-admissible vertices, NOT complete paths.
    # Intron counts are geometric vertex-pair upper bounds, NOT phase-compatible edges.
    geometric = len(original) if original is not None else sum(
        int(np.searchsorted(points[l]+(2 if l=='acceptor' else 0),
                            points[r]+(3 if r=='stop' else 0),side='left').sum())
        for l in ('start','acceptor') for r in ('donor','stop'))
    vertices = (e for e in original if viable(e,orf)) if original is not None else flat_fragments(points,orf)
    donors, acceptors, n = [], [], 0
    for e in vertices:
        n += 1
        if e.right=='donor':
            donors.append((e.b,e.run))
        if e.left=='acceptor':
            acceptors.append((e.a,e.run))
    if original is None:
        ends = np.array(sorted(b for b,_r in donors),dtype=np.int64)
        joins = sum(int(np.searchsorted(ends,a-4,side='right')) for a,_r in acceptors)
    else:
        # Enforce the old increasing-run constraint when counting original opportunities.
        donors.sort()
        acceptors.sort()
        size = max((e.run for e in original),default=-1)+2
        tree = [0]*(size+1)
        cursor, joins = 0,0
        for a,r in acceptors:
            while cursor<len(donors) and donors[cursor][0]+4<=a:
                i = donors[cursor][1]+1
                while i<=size:
                    tree[i] += 1
                    i += i & -i
                cursor += 1
            i = r
            while i:
                joins += tree[i]
                i -= i & -i
    return dict(geometric_fragments=geometric,local_ORF_admissible_fragments=n,
                geometric_intron_vertex_pair_upper_bound=joins)


def summarize(rows):
    result = {'references':len(rows),'outside_O0':sum(not r['O0'] for r in rows),'arms':{}}
    for arm in ARMS:
        result['arms'][arm] = {
            'reachable':sum(r[arm] for r in rows),
            'joint_region_supported':sum(r[arm] and r['region_supported'] for r in rows),
            'outside_O0_reachable':sum(r[arm] and not r['O0'] for r in rows),
            'gained_vs_original':sum(r[arm] and not r['original'] for r in rows),
            'lost_vs_original':sum(not r[arm] and r['original'] for r in rows)}
        result['arms'][arm]['net_vs_original'] = result['arms'][arm]['gained_vs_original']-result['arms'][arm]['lost_vs_original']
    result['ranked_vs_permuted'] = {
        'ranked_only':sum(r['ranked'] and not r['permuted'] for r in rows),
        'permuted_only':sum(r['permuted'] and not r['ranked'] for r in rows)}
    result['reasons_may_overlap'] = {
        'ORF_or_geometry_incompatible':sum(not r['grammar'] for r in rows),
        'noncanonical_event_motif':sum(not r['all_motifs'] for r in rows),
        'some_boundary_score_insufficient':sum(not r['all_raw_boundary_pass'] for r in rows),
        'some_rank_budget_pruned':sum(r['rank_budget_pruned']>0 for r in rows),
        'some_original_endpoint_absent':sum(r['original_endpoint_absent']>0 for r in rows),
        'some_original_fragment_absent':sum(r['original_fragment_absent']>0 for r in rows),
        'all_original_fragments_but_no_run_order':sum(r['run_order_loss'] for r in rows),
        'some_CDS_region_unsupported':sum(not r['region_supported'] for r in rows)}
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir',required=True)
    args = parser.parse_args()
    out = Path(args.output_dir)
    assert not (out/'endpoint_budget.json').exists()
    out.mkdir(parents=True,exist_ok=True)
    ledger, counts = [], []
    # Reference membership is consulted only after all candidates for that orientation are fixed/saved.
    with (S.ROOT/'outputs/M26-ACTUAL-GRAPH-SUPPORT-R3/support_ledger.tsv').open() as f:
        old = {(r['species'],r['strand'],r['transcript']):r for r in csv.DictReader(f,delimiter='\t')}
    assert len(old)==6450
    for si,(species,seqid,length,expected) in enumerate(S.SCOPE):
        sequence = S.V.read_fasta(S.ROOT/'data/m1_screen'/species/'genome.fa')[seqid].upper()
        assert len(sequence)==length
        for ti,strand in enumerate(('+','-')):
            seq = sequence if strand=='+' else S.M.reverse_complement(sequence)
            stem = f'{seqid}_{"plus" if strand=="+" else "minus"}'
            cache = S.ROOT/'outputs/M25R-E1-RUN-GRAPH-R1/scores'
            region = np.load(cache/f'{stem}_region.npy',mmap_mode='r')
            boundary = np.load(cache/f'{stem}_boundary.npy',mmap_mode='r')
            assert region.shape==(length,3) and boundary.shape==(length,4)
            assert region.dtype==boundary.dtype==np.float16
            exons,intron_prefix,_ = G.candidates(seq,region,boundary)
            del intron_prefix
            original = endpoints(exons)
            points = {'original':original,'unlocked':original,'ranked':{},'permuted':{}}
            event_counts = {}
            for ei,event in enumerate(EVENTS):
                positions = motifs(seq,event)
                scores = boundary[positions,ei].astype(np.float32)
                assert np.isfinite(scores).all()
                shuffled = permute_windows(positions,scores,np.random.default_rng(np.random.SeedSequence([SEED,si,ti,ei])))
                m = len(original[event])
                points['ranked'][event] = select(positions,scores,m,CUTOFF[event])
                points['permuted'][event] = select(positions,shuffled,m,CUTOFF[event])
                assert all(exists(positions,p) for p in original[event])
                assert all(len(points[a][event])==m for a in ARMS)
                event_counts[event] = dict(motif_pool=len(positions),above_cutoff=int((scores>=CUTOFF[event]).sum()),budget=m)
            np.savez_compressed(out/f'{stem}_endpoints.npz',**{a+'_'+e:points[a][e] for a in ARMS for e in EVENTS})
            del positions,scores,shuffled
            orf = G.ORF(seq)
            opportunities = {}
            for arm in ARMS:
                opportunities[arm] = opportunity_counts(points[arm],orf,exons if arm=='original' else None)
                print(stem,arm,opportunities[arm],flush=True)
            counts.append(dict(species=species,strand=strand,events=event_counts,opportunities=opportunities))
            index = defaultdict(list)
            for e in exons:
                index[e.a,e.b,e.left,e.right].append(e.run)
            # C-I and C-G interval margins are diagnostic only, never candidate filters.
            margins = [np.r_[0.,np.cumsum(region[:,1].astype(np.float64)-region[:,c].astype(np.float64))]
                       for c in (0,2)]
            ann = S.E.parse_annotation(S.ROOT/'data/m1_screen'/species/'reference.gff3',{seqid:length},protein_coding_only=True)
            refs = S.E.primary_transcripts(ann)
            assert len(refs)==expected
            for t in refs:
                if t['strand']!=strand:
                    continue
                chain = [(a,b) for a,b,_ in t['CDS']]
                if strand=='-':
                    chain = [(length-b,length-a) for a,b in reversed(chain)]
                prior = old[species,strand,t['id']]
                o1 = P.membership(chain,index,orf)
                assert o1==(prior['O1']=='True'), 'Original membership identity mismatch'
                gr = grammar(chain,orf)
                event_rows = []
                for event,p in events_for(chain):
                    width = len(next(iter(S.M.MOTIFS[event])))
                    valid = 0<=p and p+width<=length and seq[p:p+width] in S.M.MOTIFS[event]
                    raw = float(boundary[p,EVENTS.index(event)]) if 0<=p<length else None
                    event_rows.append(dict(event=event,position=p,motif=valid,raw_logit=raw,
                        raw_pass=raw is not None and raw>=CUTOFF[event],
                        **{a:bool(exists(points[a][event],p)) for a in ARMS}))
                cds_rows = []
                for a,b,left,right in fragments(chain):
                    delta = [float(prefix[b]-prefix[a]) for prefix in margins]
                    cds_rows.append(dict(a=a,b=b,C_minus_I=delta[0],C_minus_G=delta[1],
                        original_fragment=bool(index.get((a,b,left,right))),region_support=min(delta)>0))
                allfragments = all(r['original_fragment'] for r in cds_rows)
                row = dict(species=species,seqid=seqid,strand=strand,transcript=t['id'],
                    O0=prior['O0']=='True',original=o1,
                    **{a:gr and all(r[a] for r in event_rows) for a in ARMS if a!='original'},
                    grammar=gr,region_supported=all(r['region_support'] for r in cds_rows),
                    all_raw_boundary_pass=all(r['raw_pass'] for r in event_rows),
                    all_motifs=all(r['motif'] for r in event_rows),
                    rank_budget_pruned=sum(r['motif'] and r['raw_pass'] and not r['ranked'] for r in event_rows),
                    original_endpoint_absent=sum(not r['original'] for r in event_rows),
                    original_fragment_absent=sum(not r['original_fragment'] for r in cds_rows),
                    run_order_loss=allfragments and P.witness([index[x] for x in fragments(chain)]) is None,
                    events=event_rows,CDS=cds_rows)
                assert not o1 or row['unlocked']
                ledger.append(row)
            del exons,index,orf,margins,region,boundary,points
    assert len(ledger)==6450 and sum(r['original'] for r in ledger)==3561
    assert sum(r['O0'] for r in ledger)==3585
    summary = dict(status='fixed_budget_reference_assisted_DEV_diagnostic',setaria_access=False,
                   new_forward=False,training=False,predictions_changed=False,seed=SEED,
                   pooled=summarize(ledger),per_species={s:summarize([r for r in ledger if r['species']==s])
                       for s,_q,_n,_c in S.SCOPE},candidate_opportunities=counts)
    pools = summary['pooled']['arms']
    summary['investment_gate'] = dict(
        minimum_net_323=pools['ranked']['net_vs_original']>=323,
        positive_net_both_species=all(x['arms']['ranked']['net_vs_original']>0 for x in summary['per_species'].values()),
        exceeds_permuted=pools['ranked']['reachable']>pools['permuted']['reachable'],
        joint_support_exceeds_original_and_permuted=pools['ranked']['joint_region_supported']>
            max(pools[a]['joint_region_supported'] for a in ('original','permuted')))
    summary['investment_gate']['pass'] = all(summary['investment_gate'].values())
    with gzip.open(out/'chain_ledger.jsonl.gz','wt') as f:
        for row in ledger:
            f.write(json.dumps(row)+'\n')
    with (out/'endpoint_budget.json').open('x') as f:
        json.dump(summary,f,indent=2)
        f.write('\n')
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    main()
