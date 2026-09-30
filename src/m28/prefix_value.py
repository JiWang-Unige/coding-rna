"""B2 learned absolute prefix compatibility; no residual additive search energy."""
from collections import defaultdict
import math
import torch
from torch import nn

STAGES={'entry':0,'donor':1,'final':2}

class PrefixValue(nn.Module):
    def __init__(self,hidden=128):
        super().__init__()
        self.hidden=hidden
        self.exon=nn.Sequential(nn.Linear(hidden*3+3,hidden),nn.GELU())
        self.sequence=nn.GRUCell(hidden,hidden)
        self.stage=nn.Embedding(3,16)
        self.value=nn.Sequential(nn.Linear(hidden*2+16+4,hidden),nn.GELU(),nn.Linear(hidden,1))

    def for_window(self,features):
        return PrefixWindow(self,features)

class PrefixWindow:
    """Per-window hidden-prefix cache avoids recomputing earlier exons.

    Construct separate instances for no-grad search and gradient supervision.
    It is intentionally not persistent across windows or optimizer updates.
    """
    def __init__(self,model,features):
        self.model=model
        self.features=features.detach()
        self.sums=torch.cat([self.features.new_zeros((1,model.hidden)),self.features.cumsum(0)])
        self.hidden={():self.features.new_zeros(model.hidden)}

    def _encode(self,chains):
        missing=set()
        for chain in chains:
            prefix=chain
            while prefix not in self.hidden:
                missing.add(prefix);prefix=prefix[:-1]
        by_length=defaultdict(list)
        for prefix in missing: by_length[len(prefix)].append(prefix)
        for length in sorted(by_length):
            group=sorted(by_length[length])
            a=torch.tensor([c[-1][0] for c in group],device=self.features.device)
            b=torch.tensor([c[-1][1] for c in group],device=self.features.device)
            if ((a<0)|(b<=a)|(b>len(self.features))).any():
                raise ValueError('Invalid prefix exon coordinates')
            averages=(self.sums[b]-self.sums[a])/(b-a)[:,None]
            geometry=self.features.new_tensor([[math.log1p(y-x),
                        math.log1p(x-c[-2][1]) if len(c)>1 else 0.,
                        sum(v-u for u,v in c[:-1])%3/2.] for c in group for x,y in [c[-1]]])
            tokens=self.model.exon(torch.cat([averages,self.features[a],self.features[b-1],geometry],dim=-1))
            parents=torch.stack([self.hidden[c[:-1]] for c in group])
            states=self.model.sequence(tokens,parents)
            for chain,state in zip(group,states): self.hidden[chain]=state

    def score(self,stage,position,states):
        if not states: return self.features.new_empty(0)
        if stage not in STAGES: raise ValueError(stage)
        chains=[tuple(tuple(x) for x in s[1]) for s in states]
        if stage!='entry' and any(not c for c in chains): raise ValueError('Only entry can be empty')
        self._encode(chains)
        hidden=torch.stack([self.hidden[c] for c in chains])
        endpoints=[int(position) if stage=='entry' else
                   int(position)-1 if stage=='donor' else c[-1][1]-1 for c in chains]
        if any(not 0<=x<len(self.features) for x in endpoints): raise ValueError('Invalid prefix event')
        event_features=self.features[endpoints]
        geometry=self.features.new_tensor([[math.log1p(len(c)),
                    math.log1p(c[-1][1]-c[0][0]) if c else 0.,
                    math.log1p(sum(b-a for a,b in c)),
                    math.log1p(int(position)-c[-1][1]) if stage=='entry' and c else 0.] for c in chains])
        stage_ids=torch.full((len(states),),STAGES[stage],device=self.features.device,dtype=torch.long)
        return self.model.value(torch.cat([hidden,event_features,self.model.stage(stage_ids),geometry],dim=-1)).squeeze(-1)

    def priorities(self,stage,position,states):
        return self.score(stage,position,states).detach().cpu().numpy()
