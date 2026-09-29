"""Reference-free bounded exon DAG for B1; approximate additive proposal search.

Only DNA and neural scores enter generate(). Supervision/injection lives in
training.py. All coordinates refer to oriented 0-based half-open CDS intervals.
"""
from bisect import bisect_left
from collections import defaultdict
from dataclasses import dataclass,asdict
import math
import numpy as np

STOP={"TAA","TAG","TGA"}

@dataclass(frozen=True)
class Budget:
    bin_bp:int=1024
    sites_per_bin:int=8
    ends_per_start_kind:int=8
    intron_links_per_donor:int=8
    beam_per_carry:int=2
    chains_per_kb:int=8
    max_exons:int=128
    min_intron_bp:int=5

class ORFIndex:
    def __init__(self,sequence):
        self.sequence=sequence
        self.bad=np.concatenate([[0],np.cumsum([c not in "ACGT" for c in sequence])])
        self.stops={f:[] for f in range(3)}
        for i in range(len(sequence)-2):
            if sequence[i:i+3] in STOP: self.stops[i%3].append(i)

    def extend(self,a,b,carry,terminal=False):
        """Carry is the actual unfinished codon, not phase alone."""
        if self.bad[b]!=self.bad[a]: return None
        size=b-a
        if carry:
            need=3-len(carry)
            if size<need: return None if terminal else carry+self.sequence[a:b]
            if carry+self.sequence[a:a+need] in STOP: return None
            first=a+need
        else: first=a
        leftover=(len(carry)+size)%3
        positions=self.stops[first%3]
        lo=bisect_left(positions,first);hi=bisect_left(positions,b-2)
        if terminal:
            if leftover or hi-lo!=1 or positions[lo]!=b-3: return None
        elif hi>lo: return None
        return self.sequence[b-leftover:b] if leftover else ""

def _top_sites(coordinates,scores,budget,offset=0):
    bins=defaultdict(list)
    for x in coordinates: bins[(x+offset)//budget.bin_bp].append(x)
    kept=[]
    for group in bins.values():
        kept.extend(sorted(group,key=lambda x:(-float(scores[x]),x))[:budget.sites_per_bin])
    return sorted(kept)

def generate(sequence,endpoint_logits,cds_probability,link_score_fn=None,budget=Budget()):
    """Return the ACTUAL pruned candidate list, scores, considered links and counts."""
    L=len(sequence)
    ep=np.asarray(endpoint_logits,dtype=float);cp=np.asarray(cds_probability,dtype=float)
    if ep.shape!=(L,4) or cp.shape!=(L,) or not np.isfinite(ep).all() or not np.isfinite(cp).all():
        raise ValueError("DNA / finite neural-score alignment mismatch")
    raw={"start":[i for i in range(L-2) if sequence[i:i+3]=="ATG"],
         "stop":[i+3 for i in range(L-2) if sequence[i:i+3] in STOP],
         "donor":[i for i in range(1,L-1) if sequence[i:i+2] in ("GT","GC")],
         "acceptor":[i+2 for i in range(L-2) if sequence[i:i+2]=="AG"]}
    boundary_scores={}
    for kind,col,offset in (("start",0,0),("stop",1,-1),("donor",2,-1),("acceptor",3,0)):
        boundary_scores[kind]={x:float(ep[x+offset,col]) for x in raw[kind]}
    sites={k:_top_sites(v,boundary_scores[k],budget,-1 if k in ("stop","donor") else 0) for k,v in raw.items()}
    counts={"raw_sites":{k:len(v) for k,v in raw.items()},
            "kept_sites":{k:len(v) for k,v in sites.items()},
            "exon_edges_considered":0,"legal_exon_edges":0,"exon_edges_retained":0,
            "intron_links_considered":0,"intron_links_retained":0,
            "beam_states_considered":0,"beam_states_pruned":0,"duplicate_states":0,
            "extensions":0,"ORF_rejected_context_edges":0,"max_exon_pruned":0,
            "complete_paths":0,"duplicate_complete_paths":0,"final_chains_pruned":0}
    pairs=[(d,a) for d in sites["donor"] for a in sites["acceptor"] if a-d>=budget.min_intron_bp]
    ls=np.zeros(len(pairs)) if link_score_fn is None else np.asarray(link_score_fn(pairs),dtype=float)
    if ls.shape!=(len(pairs),) or not np.isfinite(ls).all(): raise ValueError("Invalid link scorer output")
    outgoing=defaultdict(list)
    for pair,score in zip(pairs,ls): outgoing[pair[0]].append((pair[1],float(score)))
    counts["intron_links_considered"]=len(pairs)
    for d,choices in outgoing.items():
        choices.sort(key=lambda x:(-x[1],x[0]))
        outgoing[d]=choices[:budget.intron_links_per_donor]
        counts["intron_links_retained"]+=len(outgoing[d])
    kept_pairs=[(d,a) for d,choices in outgoing.items() for a,_ in choices]
    pc=np.clip(cp,1e-5,1-1e-5)
    coding_prefix=np.concatenate([[0.],np.cumsum(np.log(pc/(1-pc)))])
    end_events=sorted([(x,"donor") for x in sites["donor"]]+[(x,"stop") for x in sites["stop"]])
    orf=ORFIndex(sequence)
    edge_cache={}
    def edges(a,carry,first_exon):
        key=(a,carry,first_exon)
        if key in edge_cache: return edge_cache[key]
        choices=defaultdict(list)
        for b,kind in end_events:
            if b<=a or ((first_exon or kind=="stop") and b-a<3): continue
            counts["exon_edges_considered"]+=1
            remainder=orf.extend(a,b,carry,kind=="stop")
            if remainder is None:
                counts["ORF_rejected_context_edges"]+=1;continue
            counts["legal_exon_edges"]+=1
            score=boundary_scores[kind][b]+float((coding_prefix[b]-coding_prefix[a])/(b-a))
            choices[kind].append((b,kind,score,remainder))
        kept=[]
        for values in choices.values():
            values.sort(key=lambda x:(-x[2],x[0]))
            kept+=values[:budget.ends_per_start_kind]
        counts["exon_edges_retained"]+=len(kept)
        edge_cache[key]=kept
        return kept
    # state = (score, exon tuple, unfinished-codon string). No annotation identifiers.
    entries=defaultdict(list);donors=defaultdict(list);completed={}
    starts=set(sites["start"]);acceptors=set(sites["acceptor"]);donor_sites=set(sites["donor"])

    def prune(states):
        counts["beam_states_considered"]+=len(states)
        unique={}
        for state in states:
            key=(state[1],state[2])
            if key in unique:
                counts["duplicate_states"]+=1
                if state[0]>unique[key][0]: unique[key]=state
            else: unique[key]=state
        groups=defaultdict(list)
        for state in unique.values(): groups[state[2]].append(state)
        result=[]
        for values in groups.values():
            values.sort(key=lambda x:(-x[0],x[1]))
            counts["beam_states_pruned"]+=max(0,len(values)-budget.beam_per_carry)
            result+=values[:budget.beam_per_carry]
        return result

    for position in sorted(starts|acceptors|donor_sites):
        if position in donor_sites:
            for score,chain,carry in prune(donors.pop(position,[])):
                for a,link_score in outgoing[position]:
                    entries[a].append((score+link_score+boundary_scores["acceptor"][a],chain,carry))
        if position not in starts and position not in acceptors: continue
        incoming=entries.pop(position,[])
        if position in starts: incoming.append((boundary_scores["start"][position],(),""))
        for score,chain,carry in prune(incoming):
            if len(chain)>=budget.max_exons:
                counts["max_exon_pruned"]+=1;continue
            for end,kind,edge_score,remainder in edges(position,carry,not chain):
                terminal=kind=="stop"
                counts["extensions"]+=1
                new_chain=chain+((position,end),)
                new_score=score+edge_score
                if terminal:
                    counts["complete_paths"]+=1
                    if new_chain in completed: counts["duplicate_complete_paths"]+=1
                    completed[new_chain]=max(completed.get(new_chain,-math.inf),new_score)
                else: donors[end].append((new_score,new_chain,remainder))
    limit=max(1,math.ceil(L/1024)*budget.chains_per_kb)
    ranked=sorted(completed.items(),key=lambda x:(-x[1],x[0]))
    counts["final_chains_pruned"]=max(0,len(ranked)-limit)
    ranked=ranked[:limit]
    return {"chains":[list(c) for c,_ in ranked],"proposal_scores":[s for _,s in ranked],
            "links":kept_pairs,"counts":counts,"budget":asdict(budget),"reference_used":False}
