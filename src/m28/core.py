"""M28 model core. Inputs are aligned frozen-GLM features, not an encoder wrapper.

C0/B1 must receive the same feature cache and window order. This module is not
yet an end-to-end annotation caller. Discrete candidate selection is outside
autograd; score gradients reach the selected exons and shared context.
"""
import math
import torch
from torch import nn
from torch.nn.utils.rnn import pack_padded_sequence, pad_packed_sequence

class SharedContext(nn.Module):
    def __init__(self, feature_dim, hidden=128, dna_hidden=32):
        super().__init__()
        self.project=nn.Linear(feature_dim,hidden)
        self.context=nn.GRU(hidden,hidden//2,batch_first=True,bidirectional=True)
        self.dna=nn.Sequential(nn.Conv1d(5,dna_hidden,9,padding=4),nn.GELU(),
                               nn.Conv1d(dna_hidden,dna_hidden,5,padding=2),nn.GELU())
        self.fuse=nn.Sequential(nn.Linear(hidden+dna_hidden,hidden),nn.GELU(),nn.LayerNorm(hidden))

    def forward(self, token_features, dna, valid_bases):
        # token_features [B,ceil(L/6),D]; DNA one-hot A,C,G,T,N [B,L,5].
        # Token features contain exactly one vector per six bases: no CLS/SEP.
        B,L,_=dna.shape
        if token_features.shape[1] != math.ceil(L/6):
            raise ValueError("Features must be aligned six-mer vectors without special tokens")
        valid_bases=torch.as_tensor(valid_bases,device=dna.device,dtype=torch.long)
        if valid_bases.numel()!=B or (valid_bases<=0).any() or (valid_bases>L).any():
            raise ValueError("Invalid per-window valid length")
        mask=torch.arange(L,device=dna.device)[None,:] < valid_bases[:,None]
        packed=pack_padded_sequence(self.project(token_features),((valid_bases+5)//6).cpu(),
                                   batch_first=True,enforce_sorted=False)
        ctx,_=self.context(packed)
        ctx,_=pad_packed_sequence(ctx,batch_first=True,total_length=token_features.shape[1])
        ctx=ctx.repeat_interleave(6,dim=1)[:,:L]
        local=self.dna((dna*mask[...,None]).transpose(1,2)).transpose(1,2)
        return self.fuse(torch.cat([ctx,local],dim=-1))*mask[...,None],mask

class ChainHead(nn.Module):
    """Non-additive ordered-exon score with a shared null competitor.

    Multiple nonoverlapping genes may coexist: never softmax all window chains
    into a single mutually exclusive class. Use candidate-vs-null BCE; a later
    interval selector handles overlapping structures with positive score gain.
    """
    def __init__(self,hidden=128):
        super().__init__()
        self.exon=nn.Sequential(nn.Linear(hidden*3+3,hidden),nn.GELU())
        self.sequence=nn.GRU(hidden,hidden,batch_first=True)
        self.additive=nn.Linear(hidden,1)
        self.nonadditive=nn.Sequential(nn.Linear(hidden,hidden),nn.Tanh(),nn.Linear(hidden,1))
        self.null=nn.Linear(hidden,1)
        self.link=nn.Sequential(nn.Linear(hidden*2+1,hidden),nn.GELU(),nn.Linear(hidden,1))

    def link_logits(self,features,pairs):
        if not pairs: return features.new_empty((0,))
        indices=torch.as_tensor(pairs,device=features.device,dtype=torch.long)
        donor,acceptor=indices[:,0],indices[:,1]
        if ((donor<=0)|(acceptor<=donor)|(acceptor>=len(features))).any():
            raise ValueError("Expected oriented exon-end / next-exon-start coordinates")
        geometry=torch.log1p((acceptor-donor).to(features.dtype))[:,None]
        values=torch.cat([features[donor-1],features[acceptor],geometry],dim=-1)
        return self.link(values).squeeze(-1)

    def forward(self,features,chains,valid_bases=None):
        L=len(features) if valid_bases is None else int(valid_bases)
        features=features[:L]
        null=self.null(features.mean(dim=0)).squeeze(-1)
        additive,nonadditive=[],[]
        for chain in chains:
            if not chain: raise ValueError("Empty chain is represented by null, not by an exon list")
            values=[]
            previous=0
            coding=0
            for j,(a,b) in enumerate(chain):
                if not 0<=a<b<=L or (j and a<previous):
                    raise ValueError("Exons must be nonoverlapping, ordered, oriented half-open intervals")
                geo=features.new_tensor([math.log1p(b-a),math.log1p(a-previous) if j else 0.,coding%3/2.])
                values.append(torch.cat([features[a:b].mean(dim=0),features[a],features[b-1],geo]))
                previous=b
                coding+=b-a
            z=self.exon(torch.stack(values))
            _,h=self.sequence(z[None])
            additive.append(self.additive(z).sum())
            nonadditive.append(self.nonadditive(h[-1,0]).squeeze(-1))
        add=torch.stack(additive) if additive else features.new_empty((0,))
        non=torch.stack(nonadditive) if nonadditive else features.new_empty((0,))
        return {"additive":add,"nonadditive":non,"null":null,
                "logits":add+non-null,"additive_only_logits":add-null}

class M28Core(nn.Module):
    def __init__(self,arm,feature_dim,hidden=128):
        super().__init__()
        if arm not in ("C0","B1"): raise ValueError(arm)
        self.arm=arm
        self.shared=SharedContext(feature_dim,hidden)
        self.segmentation=nn.Linear(hidden,15)
        self.endpoint=nn.Linear(hidden,4) if arm=="B1" else None
        self.chain=ChainHead(hidden) if arm=="B1" else None

    def forward(self,token_features,dna,valid_bases):
        h,mask=self.shared(token_features,dna,valid_bases)
        result={"features":h,"mask":mask,"segmentation_logits":self.segmentation(h)}
        if self.endpoint is not None: result["endpoint_logits"]=self.endpoint(h)
        return result

def chain_loss(scores,targets,weights=None):
    """Unknown/alternative/partial chains use target=-1, never a false-negative label."""
    targets=torch.as_tensor(targets,device=scores["logits"].device,dtype=scores["logits"].dtype)
    keep=targets>=0
    if not keep.any():
        return scores["logits"].sum()*0+scores["null"]*0
    loss=nn.functional.binary_cross_entropy_with_logits(scores["logits"][keep],targets[keep],reduction="none")
    if weights is not None: loss=loss*torch.as_tensor(weights,device=loss.device)[keep]
    return loss.mean()

def select_nonoverlapping(chains,logits):
    """Maximum positive-gain same-strand interval set. Zero gain selects null."""
    import bisect
    if len(chains)!=len(logits): raise ValueError("One score per chain required")
    order=sorted(range(len(chains)),key=lambda i:(chains[i][-1][1],chains[i][0][0],i))
    ends=[chains[i][-1][1] for i in order]
    # Backpointers avoid copying the whole gene set at every chromosome candidate.
    value=[0.]; previous_states=[]; take=[]
    for k,i in enumerate(order):
        previous=bisect.bisect_right(ends,chains[i][0][0],0,k)
        include=value[previous]+float(logits[i])
        chosen=include>value[-1]
        previous_states.append(previous);take.append(chosen)
        value.append(include if chosen else value[-1])
    picked=[];k=len(order)
    while k:
        if take[k-1]:
            picked.append(order[k-1]);k=previous_states[k-1]
        else: k-=1
    return picked[::-1]
