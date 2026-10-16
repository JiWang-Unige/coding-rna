"""M28 training-only reference supervision; never imported by candidate search."""
import numpy as np
import torch
from .labels import local_loss,candidate_targets,oriented_parts

def balanced_binary_loss(logits,labels,known=None):
    """Equal mass to observed positive/negative classes, mean per output column.

    Each window has unit mass before the window sampling weight. If only one
    class is present it has unit mass; no known examples give exactly zero.
    """
    y=torch.as_tensor(labels,device=logits.device,dtype=logits.dtype)
    mask=y>=0 if known is None else torch.as_tensor(known,device=logits.device,dtype=torch.bool)
    if logits.ndim==1: logits=logits[:,None];y=y[:,None];mask=mask[:,None]
    terms=[]
    for k in range(logits.shape[-1]):
        parts=[]
        for cls in (0,1):
            selected=mask[:,k]&(y[:,k]==cls)
            if selected.any():
                parts.append(torch.nn.functional.binary_cross_entropy_with_logits(
                    logits[selected,k],y[selected,k],reduction="mean"))
        terms.append(torch.stack(parts).mean() if parts else logits[:,k].sum()*0)
    return torch.stack(terms).mean()

def junction_sets(window,catalog):
    primary=set();protected=set();n=window["valid_bases"]
    for ident in window["overlapping_transcript_ids"]:
        tx=catalog[ident];parts=oriented_parts(tx,window)
        pairs={(y,z) for (_,y),(z,_) in zip(parts,parts[1:]) if 0<y<z<n}
        protected |= pairs
        if tx["training_primary"]: primary |= pairs
    return primary,protected

def connection_targets(pairs,window,catalog,targets):
    primary,protected=junction_sets(window,catalog)
    labels=[]
    for d,a in pairs:
        known=targets["endpoint_known"][d-1,2] and targets["endpoint_known"][a,3]
        if (d,a) in primary and known: labels.append(1)
        elif (d,a) in protected or not known or any(x<a and y>d for x,y in targets["unknown_chain_spans"]) or not targets["callable_bases"][d-1:a+1].all():
            labels.append(-1)
        else: labels.append(0)
    return labels

def inject_training(free_chains,positive_chains):
    """Copy/deduplicate, retaining an immutable record of free vs injected origin."""
    unique={};origins=[]
    for source,chains in (("free",free_chains),("reference_injected",positive_chains)):
        for chain in chains:
            key=tuple(tuple(x) for x in chain)
            if key not in unique: unique[key]=list(key);origins.append(source)
    return list(unique.values()),origins

def joint_losses(model,outputs,free,window,catalog,targets,window_weight=1.,coefficients=(1.,1.,1.,1.)):
    """One real window. Proposals were already generated and saved without labels."""
    h=outputs["features"][0,:window["valid_bases"]]
    pool,origins=inject_training(free["chains"],targets["positive_chains"])
    chain_labels=candidate_targets(pool,targets)
    scores=model.chain(h,pool)
    primary,_=junction_sets(window,catalog)
    pairs=list(dict.fromkeys(tuple(x) for x in free["links"]))
    pairs+=sorted(primary-set(pairs))
    link_labels=connection_targets(pairs,window,catalog,targets)
    link_logits=model.chain.link_logits(h,pairs)
    losses={"local":local_loss(outputs["segmentation_logits"][0],targets),
            "endpoint":balanced_binary_loss(outputs["endpoint_logits"][0],targets["endpoint"],targets["endpoint_known"]),
            "link":balanced_binary_loss(link_logits,link_labels),
            "chain":balanced_binary_loss(scores["logits"],chain_labels)}
    total=sum(c*loss for c,loss in zip(coefficients,losses.values()))*float(window_weight)
    return {"total":total,"losses":losses,"chain_scores":scores,"chain_labels":chain_labels,
            "training_chains":pool,"origins":origins,"link_pairs":pairs,
            "link_logits":link_logits,"link_labels":link_labels}

def one_hot(sequence,width):
    codes=torch.tensor(["ACGT".find(c) if c in "ACGT" else 4 for c in sequence]+[4]*(width-len(sequence)))
    return torch.nn.functional.one_hot(codes,5).float()
