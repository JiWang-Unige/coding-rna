"""Query-conditioned non-additive values for legal gene-prefix events.

The frozen pointer supplies recurrent states, not the ranking score. Prefix CDS
features are averaged, never summed into an exon-count reward. Values are not
monotone future bounds or biological probabilities. Labels live separately.
"""
import math
import torch
from torch import nn
from .query_pointer import STAGES

START, NEXT_START, DONOR, STOP = range(4)
EVENT_NAMES = ('start','next_start','donor','stop')


class PrefixValueHead(nn.Module):
    def __init__(self, hidden=128):
        super().__init__()
        self.event = nn.Embedding(4,8)
        self.value = nn.Sequential(nn.Linear(3*hidden+6+8,hidden),nn.GELU(),
                                   nn.Linear(hidden,hidden),nn.Tanh(),nn.Linear(hidden,1))

    def forward(self, initial, state, coding_mean, geometry, events):
        # Freeze source states/features even if a caller accidentally retains a graph.
        x=torch.cat([initial.detach(),state.detach(),coding_mean.detach(),geometry.detach(),self.event(events)],dim=-1)
        return self.value(x).squeeze(-1)


def root_row(window,query_index):
    return dict(query_index=query_index,state=window.initial[query_index],chain=(),carry='',
                previous_end=-1,coding_sum=window.features.new_zeros(window.model.hidden),coding_bases=0,
                path_log_score=0.,prefix_logit=0.)


def feature_sums(window):
    return torch.cat([window.features.new_zeros((1,window.model.hidden)),window.features.cumsum(0)])


def _inputs(window,row,states,sums,sizes,coordinates,first,events,carries):
    n=len(coordinates);L=window.length
    geometry=torch.stack([first/L,(coordinates-first)/L,sizes/L,
        sizes.new_full((n,),math.log1p(len(row['chain'])+1)/math.log1p(128)),
        sizes.new_tensor([len(c)/2 for c in carries]),(events==STOP).to(sizes.dtype)],dim=-1)
    return (window.initial[row['query_index']].expand(n,-1),states,sums/sizes[:,None],geometry,events)


@torch.no_grad()
def entry_events(window,row,coordinates):
    cs=torch.tensor(coordinates,device=row['state'].device,dtype=torch.long)
    event=START if not row['chain'] else NEXT_START
    stage='start' if event==START else 'next_start'
    tokens=window.entry_features[cs]+window.positions[cs]+window.model.event.weight[STAGES[stage]]
    states=window.model.recurrent(tokens,row['state'].expand(len(cs),-1))
    sums=row['coding_sum'][None]+window.features[cs]
    sizes=window.features.new_full((len(cs),),row['coding_bases']+1)
    first=cs.to(sizes.dtype) if not row['chain'] else sizes.new_full((len(cs),),row['chain'][0][0])
    events=torch.full((len(cs),),event,device=cs.device,dtype=torch.long)
    nodes=[dict(query_index=row['query_index'],event=event,chain=row['chain'],entry=a) for a in coordinates]
    return dict(states=states,nodes=nodes,inputs=_inputs(window,row,states,sums,sizes,cs.to(sizes.dtype),first,events,[row['carry']]*len(cs)))


@torch.no_grad()
def end_events(window,prefix_sums,row,choices):
    cs=torch.tensor([b for b,_,_ in choices],device=row['state'].device,dtype=torch.long)
    tokens=window.end_features[cs]+window.positions[cs]+window.model.event.weight[STAGES['end']]
    states=window.model.recurrent(tokens,row['state'].expand(len(cs),-1))
    sums=row['coding_sum'][None]+prefix_sums[cs]-prefix_sums[row['entry']]
    sizes=(cs-row['entry']+row['coding_bases']).to(window.features.dtype)
    first=sizes.new_full((len(cs),),row['chain'][0][0] if row['chain'] else row['entry'])
    events=torch.tensor([STOP if kind=='stop' else DONOR for _,kind,_ in choices],device=cs.device)
    nodes=[dict(query_index=row['query_index'],event=int(event),chain=row['chain']+((row['entry'],b),),entry=None)
           for (b,_,_),event in zip(choices,events.tolist())]
    return dict(states=states,coding_sums=sums,coding_bases=sizes,nodes=nodes,
                inputs=_inputs(window,row,states,sums,sizes,cs.to(sizes.dtype),first,events,[c for _,_,c in choices]))
