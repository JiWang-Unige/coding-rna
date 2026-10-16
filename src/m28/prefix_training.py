"""B2 training-only prefix labels and first actual pruning supervision."""
from collections import defaultdict
import torch
from .labels import candidate_targets
from .training import balanced_binary_loss,inject_training

def prefix_targets(states,stage,position,targets):
    if stage=='final':
        return candidate_targets([list(s[1]) for s in states],targets)
    positives={tuple(tuple(e) for e in c) for c in targets['positive_chains']}
    protected={tuple(tuple(e) for e in c) for c in targets['protected_chains']}
    def matches(prefix,full):
        if stage=='donor':
            return 0<len(prefix)<len(full) and prefix==full[:len(prefix)]
        return len(prefix)<len(full) and prefix==full[:len(prefix)] and full[len(prefix)][0]==position
    labels=[]
    for _,chain,_ in states:
        prefix=tuple(tuple(x) for x in chain)
        a=prefix[0][0] if prefix else position
        b=position+1 if stage=='entry' else position
        if any(matches(prefix,c) for c in positives): labels.append(1)
        elif any(matches(prefix,c) for c in protected) or any(x<b and y>a for x,y in targets['unknown_chain_spans']) or not targets['callable_bases'][a:b].all():
            labels.append(-1)
        else: labels.append(0)
    return labels

class FirstCutObserver:
    """References select loss observations only, never alter candidate states."""
    def __init__(self,positive_chains):
        self.chains=[tuple(tuple(e) for e in c) for c in positive_chains]
        self.targets=defaultdict(list);self.blocked=set();self.cuts=[]

    def bind(self,orf,budget):
        self.budget=budget
        for ident,chain in enumerate(self.chains):
            carry=''
            for j,(a,b) in enumerate(chain):
                if carry is None: self.blocked.add(ident);break
                self.targets['entry',a].append((ident,chain[:j],carry))
                carry=orf.extend(a,b,carry,j==len(chain)-1)
                if j<len(chain)-1: self.targets['donor',b].append((ident,chain[:j+1],carry))

    def prune_event(self,stage,position,groups):
        for ident,prefix,carry in self.targets.get((stage,position),[]):
            if ident in self.blocked: continue
            values=groups.get(carry,[])
            rank=next((i for i,s in enumerate(values) if s[1]==prefix),None)
            if rank is None: self.blocked.add(ident);continue
            if rank>=self.budget.beam_per_carry:
                self.cuts.append({'stage':stage,'position':position,'positive':values[rank],
                                  'competitors':list(values[:self.budget.beam_per_carry])})
                self.blocked.add(ident)

    def final_event(self,ranked,limit):
        for ident,chain in enumerate(self.chains):
            if ident in self.blocked: continue
            rank=next((i for i,(c,_) in enumerate(ranked) if c==chain),None)
            if rank is not None and rank>=limit:
                indices=list(dict.fromkeys([0,1,limit-1])) if limit>1 else [0]
                self.cuts.append({'stage':'final','position':None,
                                  'positive':(ranked[rank][1],chain,''),
                                  'competitors':[(ranked[i][1],ranked[i][0],'') for i in indices]})

def prefix_losses(model,features,free,observer,targets,window_weight=1.):
    pool,origins=inject_training(free['chains'],targets['positive_chains'])
    scorer=model.for_window(features)
    final_states=[(0.,tuple(c),'') for c in pool]
    logits=scorer.score('final',None,final_states)
    labels=candidate_targets(pool,targets)
    final_loss=balanced_binary_loss(logits,labels) if final_states else sum(p.sum()*0 for p in model.parameters())
    ranking=[];known_pairs=0;unknown_competitors=0
    for cut in observer.cuts:
        states=[cut['positive']]+cut['competitors']
        y=prefix_targets(states,cut['stage'],cut['position'],targets)
        if y[0]!=1: raise ValueError('Observed gold prefix is not positive')
        negatives=[i for i in range(1,len(y)) if y[i]==0]
        unknown_competitors+=sum(v<0 for v in y[1:])
        if not negatives: continue
        values=scorer.score(cut['stage'],cut['position'],states)
        ranking.append(torch.nn.functional.softplus(values[negatives]-values[0]).mean())
        known_pairs+=len(negatives)
    rank_loss=torch.stack(ranking).mean() if ranking else logits.sum()*0
    return {'total':(final_loss+rank_loss)*float(window_weight),'final':final_loss,'rank':rank_loss,
            'known_rank_pairs':known_pairs,'unknown_competitors':unknown_competitors,
            'first_cut_events':len(observer.cuts),'final_labels':labels,'origins':origins,
            'training_chains':pool,'final_logits':logits}
