"""DNA-only generation ranked by learned query-prefix values, no teacher labels."""
import math
import torch
from torch.nn import functional as F
from .candidates import Budget,ORFIndex
from .query_decode import _end_events,_legal_ends
from .query_prefix_value import root_row,feature_sums,entry_events,end_events


def _best(rows,k):return sorted(rows,key=lambda r:(-r['prefix_logit'],r['chain'],r.get('entry',-1)))[:k]


@torch.no_grad()
def decode_prefix_values(model,features,sequence,head,max_exons,beam_width,completed_limit,valid_bases=None):
    if min(max_exons,beam_width,completed_limit)<1:raise ValueError('Positive explicit prefix search limits required')
    window=model.for_window(features,valid_bases);sequence=sequence.upper()
    if len(sequence)!=window.length:raise ValueError('Oriented DNA and features must align')
    sums=feature_sums(window)
    starts=[a for a in range(len(sequence)-2) if sequence[a:a+3]=='ATG']
    acceptors=[a+2 for a in range(len(sequence)-2) if sequence[a:a+2]=='AG']
    events,orf=_end_events(sequence),ORFIndex(sequence)
    results=[]
    for q in range(model.num_queries):
        live=[root_row(window,q)];done={}
        counts=dict(prefix_events_scored=0,entry_candidates=0,entry_pruned=0,continuation_pruned=0,live_at_exon_cap=0,depths=0)
        for depth in range(max_exons):
            entries=[];entry_count=0
            for row in live:
                lower=0 if depth==0 else row['previous_end']+Budget().min_intron_bp
                coordinates=starts if depth==0 else [a for a in acceptors if a>=lower]
                if not coordinates:continue
                batch=entry_events(window,row,coordinates);values=head(*batch['inputs'])
                if not bool(torch.isfinite(values).all()):raise ValueError('Nonfinite entry prefix value')
                order=torch.argsort(values,descending=True,stable=True)[:beam_width]
                stage='start' if depth==0 else 'next_start'
                raw=F.log_softmax(window.pointer_logits(row['state'],stage,lower),dim=-1)
                counts['prefix_events_scored']+=len(coordinates);entry_count+=len(coordinates)
                for k in order.tolist():
                    a=coordinates[k]
                    entries.append(dict(row,entry=a,state=batch['states'][k],prefix_logit=float(values[k]),
                                        path_log_score=row['path_log_score']+float(raw[a])))
            kept=_best(entries,beam_width)
            counts['entry_candidates']+=entry_count;counts['entry_pruned']+=entry_count-len(kept)
            continuations=[]
            for row in kept:
                a=row['entry'];choices=_legal_ends(orf,events,a,row['carry'],first_exon=depth==0)
                if not choices:continue
                batch=end_events(window,sums,row,choices);values=head(*batch['inputs'])
                if not bool(torch.isfinite(values).all()):raise ValueError('Nonfinite end prefix value')
                eos=window.model.terminate(batch['states']).squeeze(-1)
                terminal=torch.tensor([kind=='stop' for _,kind,_ in choices],device=eos.device)
                branch=torch.where(terminal,F.logsigmoid(eos),F.logsigmoid(-eos))
                raw=F.log_softmax(window.pointer_logits(row['state'],'end',a),dim=-1)
                counts['prefix_events_scored']+=len(choices)
                for k,(b,kind,carry) in enumerate(choices):
                    new=dict(query_index=q,state=batch['states'][k],chain=batch['nodes'][k]['chain'],carry=carry,
                             previous_end=b,coding_sum=batch['coding_sums'][k],coding_bases=int(batch['coding_bases'][k]),
                             prefix_logit=float(values[k]),path_log_score=row['path_log_score']+float(raw[b]+branch[k]))
                    if kind=='stop':
                        compact={k:new[k] for k in ('chain','prefix_logit','path_log_score')}
                        old=done.get(new['chain'])
                        if old is None or new['prefix_logit']>old['prefix_logit']:done[new['chain']]=compact
                    else:continuations.append(new)
            done={r['chain']:r for r in _best(list(done.values()),completed_limit)}
            live=_best(continuations,beam_width)
            counts['continuation_pruned']+=len(continuations)-len(live);counts['depths']=depth+1
            if not live:break
            if depth==max_exons-1:counts['live_at_exon_cap']=len(live)
            # Prefix/STOP logits are not monotone; no likelihood early-stop bound.
        completed=[dict(r,query_index=q,candidate_rank=i,terminated=True,stop_reason='canonical_stop')
                   for i,r in enumerate(_best(list(done.values()),completed_limit))]
        if any(not math.isfinite(r['path_log_score']) for r in completed):raise ValueError('Nonfinite path receipt')
        results.append(dict(query_index=q,query_logit=float(window.query_logits[q]),completed=completed,counts=counts))
    return results


@torch.no_grad()
def complete_query_values(model,features,sequence,head,identities,valid_bases=None):
    """One union readout for (query, complete geometry); no reference labels.

    Both arms receive these same values. Query-slot values are distinct even
    when exon geometry coincides; window geometry pooling must explicitly max
    over slots rather than mistake each slot for an independent gene.
    """
    window=model.for_window(features,valid_bases);sums=feature_sums(window)
    sequence=sequence.upper();inputs=[]
    if len(sequence)!=window.length:raise ValueError('Oriented DNA and features must align')
    identities=sorted(set((int(q),tuple(tuple(e) for e in c)) for q,c in identities))
    for q,chain in identities:
        window.validate_chain(chain);row=root_row(window,q)
        for j,(a,b) in enumerate(chain):
            entry=entry_events(window,row,[a]);row=dict(row,entry=a,state=entry['states'][0])
            coding=row['carry']+sequence[a:b];n=len(coding)%3;carry=coding[-n:] if n else ''
            end=end_events(window,sums,row,[(b,'stop' if j==len(chain)-1 else 'donor',carry)])
            row=dict(query_index=q,state=end['states'][0],chain=end['nodes'][0]['chain'],carry=carry,
                     previous_end=b,coding_sum=end['coding_sums'][0],coding_bases=int(end['coding_bases'][0]))
        inputs.append(tuple(x[0] for x in end['inputs']))
    if not inputs:return {}
    values=head(*(torch.stack([row[k] for row in inputs]) for k in range(5)))
    if not torch.isfinite(values).all():raise ValueError('Nonfinite shared complete value')
    return {identity:float(value) for identity,value in zip(identities,values)}
