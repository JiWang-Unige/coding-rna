"""TRAIN-only frozen-teacher query-prefix labels and fixed counterfactual samples.

This module is never imported by the reference-free decoder. Positive injection
is supervision, not free coverage. Compatibility is relative to a matched query,
not a claim that another known gene or protected isoform is biologically false.
"""
import torch
from torch.nn import functional as F
from .candidates import Budget,ORFIndex
from .query_decode import _end_events,_legal_ends
from .query_prefix_value import START,NEXT_START,DONOR,STOP,root_row,entry_events,end_events


def _chain(chain):return tuple(tuple(e) for e in chain)


def compatible(node,chain):
    if chain is None:return False
    chain=_chain(chain);prefix=_chain(node['chain']);n=len(prefix)
    if prefix!=chain[:n]:return False
    if node['event'] in (START,NEXT_START):return n<len(chain) and node['entry']==chain[n][0]
    return n>0 and ((node['event']==DONOR and n<len(chain)) or (node['event']==STOP and n==len(chain)))


def prefix_labels(nodes,assigned,query_labels,targets):
    if set(assigned)!={q for q,y in enumerate(query_labels) if y==1} or any(y not in (-1,0,1) for y in query_labels):
        raise ValueError('Every frozen-teacher positive query needs exactly its assigned chain')
    positives={_chain(c) for c in targets['positive_chains']}
    protected=[c for c in targets['protected_chains'] if _chain(c) not in positives]
    labels=[]
    for node in nodes:
        q=node['query_index']
        if query_labels[q]==-1:labels.append(-1);continue
        if compatible(node,assigned.get(q)):labels.append(1);continue
        chain=node['chain']
        a=chain[0][0] if chain else node['entry']
        b=node['entry']+1 if node['entry'] is not None else chain[-1][1]
        unknown=(any(compatible(node,c) for c in protected) or
                 any(x<b and y>a for x,y in targets['unknown_chain_spans']) or
                 not targets['callable_bases'][a:b].all())
        labels.append(-1 if unknown else 0)
    return labels


def prefix_loss(logits,labels,events):
    """Equal event-type means, not positive/negative class balancing."""
    labels=torch.as_tensor(labels,device=logits.device)
    if labels.shape!=logits.shape or events.shape!=logits.shape or not ((labels==-1)|(labels==0)|(labels==1)).all():
        raise ValueError('One 1/0/-1 label and event per prefix value required')
    losses=[]
    for event in (START,NEXT_START,DONOR,STOP):
        known=(labels>=0)&(events==event)
        if known.any():losses.append(F.binary_cross_entropy_with_logits(logits[known],labels[known].to(logits.dtype)))
    return torch.stack(losses).mean() if losses else logits.sum()*0


class PrefixSamples:
    """One-window deduplication of the actual current sample identities."""
    def __init__(self,window):self.window,self.rows=window,{}

    def add(self,batch,origin):
        for k,node in enumerate(batch['nodes']):
            key=(node['query_index'],node['event'],_chain(node['chain']),node['entry'])
            if key not in self.rows:
                self.rows[key]=dict(node,inputs=tuple(x[k] for x in batch['inputs']),origins=set())
            self.rows[key]['origins'].add(origin)

    def records(self):
        ordered=sorted(self.rows.values(),key=lambda r:(r['query_index'],r['event'],_chain(r['chain']),-1 if r['entry'] is None else r['entry']))
        return ordered

    def public_records(self):
        return [{k:(sorted(v) if k=='origins' else v) for k,v in r.items() if k!='inputs'} for r in self.records()]

    def inputs(self):
        records=self.records()
        if records:return tuple(torch.stack([r['inputs'][k] for r in records]) for k in range(5))
        f=self.window.features;h=self.window.model.hidden
        return (f.new_empty((0,h)),f.new_empty((0,h)),f.new_empty((0,h)),f.new_empty((0,6)),torch.empty(0,device=f.device,dtype=torch.long))


def _carry(carry,sequence,a,b):
    coding=carry+sequence[a:b];n=len(coding)%3
    return coding[-n:] if n else ''


def _end_row(row,batch,k,choice):
    b,_,carry=choice
    return dict(query_index=row['query_index'],state=batch['states'][k],chain=batch['nodes'][k]['chain'],carry=carry,
                previous_end=b,coding_sum=batch['coding_sums'][k],coding_bases=int(batch['coding_bases'][k]))


@torch.no_grad()
def collect_free(window,prefix_sums,sequence,greedy,beam):
    sequence=sequence.upper();samples=PrefixSamples(window)
    starts=[a for a in range(len(sequence)-2) if sequence[a:a+3]=='ATG']
    for q in range(window.model.num_queries):
        if starts:samples.add(entry_events(window,root_row(window,q),starts),'free_initial_all_ATG')
    paths=[(r,'free_greedy') for r in greedy]+[(r,'free_beam4') for q in beam for r in q['completed']]
    for path,origin in paths:
        row=root_row(window,path['query_index'])
        chain=_chain(path['chain'])
        for j,(a,b) in enumerate(chain):
            entry=entry_events(window,row,[a]);samples.add(entry,origin)
            row=dict(row,entry=a,state=entry['states'][0])
            kind='stop' if path['terminated'] and j==len(chain)-1 else 'donor'
            choice=(b,kind,_carry(row['carry'],sequence,a,b))
            ends=end_events(window,prefix_sums,row,[choice]);samples.add(ends,origin)
            row=_end_row(row,ends,0,choice)
    return samples


@torch.no_grad()
def collect_teacher(samples,prefix_sums,sequence,assigned):
    """Inject all assigned true prefixes; fixed raw-top4 entry hard alternatives.

    All ORF-legal ends of each teacher entry are sampled. A true boundary outside
    the unchanged decoder grammar is still TRAIN supervision, explicitly counted
    as grammar-unavailable injection; no positive is silently removed.
    """
    window=samples.window;sequence=sequence.upper()
    starts=[a for a in range(len(sequence)-2) if sequence[a:a+3]=='ATG']
    acceptors=[a+2 for a in range(len(sequence)-2) if sequence[a:a+2]=='AG']
    events,orf=_end_events(sequence),ORFIndex(sequence)
    unavailable=[]
    for q,chain in sorted(assigned.items()):
        row=root_row(window,q);window.validate_chain(chain)
        for j,(a,b) in enumerate(chain):
            lower=0 if j==0 else row['previous_end']+Budget().min_intron_bp
            stage='start' if j==0 else 'next_start'
            legal=starts if j==0 else [x for x in acceptors if x>=lower]
            if a not in legal:unavailable.append(dict(query_index=q,exon=j+1,event=stage))
            raw=window.pointer_logits(row['state'],stage,lower)
            top=[legal[k] for k in torch.argsort(raw[legal],descending=True,stable=True)[:4].tolist()] if legal else []
            coordinates=sorted(set(top+[a]));entry=entry_events(window,row,coordinates)
            samples.add(entry,'TRAIN_teacher_or_counterfactual')
            row=dict(row,entry=a,state=entry['states'][coordinates.index(a)])
            kind='stop' if j==len(chain)-1 else 'donor'
            choices=_legal_ends(orf,events,a,row['carry'],first_exon=j==0)
            target=(b,kind,_carry(row['carry'],sequence,a,b))
            if not any((x,k)==(b,kind) for x,k,_ in choices):
                unavailable.append(dict(query_index=q,exon=j+1,event='end_'+kind))
                choices.append(target);choices.sort(key=lambda r:(r[0],r[1]))
            index=next(k for k,(x,event,_) in enumerate(choices) if (x,event)==(b,kind))
            ends=end_events(window,prefix_sums,row,choices);samples.add(ends,'TRAIN_teacher_or_counterfactual')
            row=_end_row(row,ends,index,choices[index])
    return unavailable
