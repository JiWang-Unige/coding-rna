"""Local supervision for M28. Reference-derived objects are training-only.

15-class order follows ANNEVO's documented emission mapping:
intergenic; CDS0/1/2; intron0/1/2; donor0/1/2; acceptor0/1/2;
start; stop. CDS0 is the first coding base (not GFF phase).
"""
from collections import defaultdict
import numpy as np

IGNORE=-100
CDS_COLUMNS=(1,2,3,7,8,9,10,11,12,13,14)
REGION_COLUMNS=((0,),CDS_COLUMNS,(4,5,6))

def oriented_parts(tx,window):
    a,b=window["start"],window["end"]
    parts=[(x-a,y-a) if window["strand"]=="+" else (b-y,b-x) for x,y,_ in tx["CDS"]]
    return sorted(parts)

def _render(tx,window):
    L=window["width"];n=window["valid_bases"]
    parts=oriented_parts(tx,window)
    seg=np.zeros(L,dtype=np.int64);region=seg.copy()
    endpoint=np.zeros((L,4),dtype=np.float32);known=np.ones((L,4),dtype=bool)
    a,b=parts[0][0],parts[-1][1]
    lo,hi=max(0,a),min(n,b)
    supported=tx["CDS_complete_supported"]
    offset=0
    for j,(x,y) in enumerate(parts):
        left,right=max(0,x),min(n,y)
        if left<right:
            region[left:right]=1
            seg[left:right]=1+(offset+np.arange(left-x,right-x))%3
        offset+=y-x
        if j<len(parts)-1:
            z=parts[j+1][0]
            left,right=max(0,y),min(n,z)
            if left<right:
                region[left:right]=2;seg[left:right]=4+(offset-1)%3
            if 0<=y-1<n:
                endpoint[y-1,2]=1
            if 0<=z<n:
                endpoint[z,3]=1
    offset=0
    for (x,y),(z,_) in zip(parts,parts[1:]):
        offset+=y-x
        if 0<=y-1<n: seg[y-1]=7+(offset-1)%3
        if 0<=z<n: seg[z]=10+offset%3
    if supported:
        seg[max(0,a):max(0,min(n,a+3))]=13
        seg[max(0,b-3):max(0,min(n,b))]=14
        if 0<=a<n: endpoint[a,0]=1
        if 0<=b-1<n: endpoint[b-1,1]=1
    else:
        seg[lo:hi]=IGNORE
        known[lo:hi,:2]=False
    seg[n:]=IGNORE;region[n:]=IGNORE;known[n:]=False
    return {"seg":seg,"region":region,"endpoint":endpoint,"endpoint_known":known,
            "span":(lo,hi),"parts":parts}

def build_targets(window,catalog,sequence):
    L,n=window["width"],window["valid_bases"]
    if len(sequence)!=n: raise ValueError("Expected oriented, unpadded valid DNA")
    seg=np.zeros(L,dtype=np.int64);region=seg.copy()
    endpoint=np.zeros((L,4),dtype=np.float32);known=np.ones((L,4),dtype=bool)
    groups=defaultdict(list)
    for ident in window["overlapping_transcript_ids"]: groups[catalog[ident]["gene_id"]].append(catalog[ident])
    occupied=np.zeros(L,dtype=bool)
    for txs in groups.values():
        primary=next((t for t in txs if t["training_primary"]),None)
        if primary is None:
            # Primary can lie outside this window while an alternative extends into it.
            lo=min(max(0,oriented_parts(t,window)[0][0]) for t in txs)
            hi=max(min(n,oriented_parts(t,window)[-1][1]) for t in txs)
            seg[lo:hi]=IGNORE;region[lo:hi]=IGNORE;known[lo:hi]=False;occupied[lo:hi]=True
            continue
        current=_render(primary,window)
        lo,hi=current["span"]
        for alt in txs:
            if alt["id"]==primary["id"]: continue
            other=_render(alt,window)
            lo=min(lo,other["span"][0]);hi=max(hi,other["span"][1])
            # Mask disagreements, never combine exons into an invented union chain.
            for key in ("seg","region"):
                conflict=current[key][lo:hi]!=other[key][lo:hi]
                current[key][lo:hi][conflict]=IGNORE
            conflict=(current["endpoint"][lo:hi]!=other["endpoint"][lo:hi]) | ~other["endpoint_known"][lo:hi]
            current["endpoint_known"][lo:hi][conflict]=False
        sl=slice(lo,hi);seen=occupied[sl]
        for target,key in ((seg,"seg"),(region,"region")):
            incoming=current[key][sl]
            disagreement=seen & (target[sl]!=incoming)
            target[sl][~seen]=incoming[~seen]
            target[sl][disagreement]=IGNORE
        conflict=seen[:,None] & ((endpoint[sl]!=current["endpoint"][sl]) | ~current["endpoint_known"][sl])
        endpoint[sl][~seen]=current["endpoint"][sl][~seen]
        known[sl][~seen]=current["endpoint_known"][sl][~seen]
        known[sl][conflict]=False
        occupied[sl]=True
    for a,b in window["fine_label_unknown_intervals"]: seg[a:b]=IGNORE
    callable_bases=np.array([c in "ACGT" for c in sequence],dtype=bool)
    valid=np.zeros(L,dtype=bool);valid[:n]=callable_bases
    seg[~valid]=IGNORE;region[~valid]=IGNORE;known[~valid]=False
    positives=[oriented_parts(catalog[i],window) for i in window["positive_ids"]]
    protected=[oriented_parts(catalog[i],window) for i in window["overlapping_transcript_ids"]
               if oriented_parts(catalog[i],window)[0][0]>=0 and oriented_parts(catalog[i],window)[-1][1]<=n]
    return {"seg":seg,"region":region,"endpoint":endpoint,"endpoint_known":known,
            "positive_chains":positives,"protected_chains":protected,
            "unknown_chain_spans":window["chain_unknown_intervals"],
            "callable_bases":valid,
            "reference_empty_window":window["no_annotated_gene_either_strand"] and bool(callable_bases.all())}

def candidate_targets(chains,targets):
    positives={tuple(x) for x in targets["positive_chains"]}
    protected={tuple(x) for x in targets["protected_chains"]}
    values=[]
    for chain in chains:
        key=tuple(chain);a,b=chain[0][0],chain[-1][1]
        if key in positives: values.append(1)
        elif key in protected or any(x<b and y>a for x,y in targets["unknown_chain_spans"]) or not targets["callable_bases"][a:b].all():
            values.append(-1)
        else: values.append(0)
    return values

def local_loss(logits,targets):
    """Fine labels where resolved; grouped region likelihood on uncertain fine labels."""
    import torch
    logp=logits.log_softmax(-1)
    fine=torch.as_tensor(targets["seg"],device=logits.device)
    region=torch.as_tensor(targets["region"],device=logits.device)
    valid=fine>=0
    loss=-logp[valid,fine[valid]].sum()
    count=valid.sum()
    for label,columns in enumerate(REGION_COLUMNS):
        keep=(fine<0)&(region==label)
        loss=loss-torch.logsumexp(logp[keep][:,list(columns)],dim=-1).sum()
        count=count+keep.sum()
    return loss/count.clamp_min(1)
