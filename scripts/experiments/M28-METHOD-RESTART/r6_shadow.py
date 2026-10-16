"""R6 finite search observer; ranking remains a literal copy of the frozen generator."""
from bisect import bisect_left
from collections import defaultdict
from dataclasses import asdict
import math
import numpy as np
from src.m28.candidates import Budget,ORFIndex,STOP,_top_sites

def generate_shadow(sequence,endpoint_logits,cds_probability,link_score_fn=None,budget=Budget(),observer=None):
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
    if observer is not None:
        observer.bind(sequence,raw,sites,boundary_scores,coding_prefix,dict(zip(pairs,ls)),orf,budget)
    # state = (score, exon tuple, unfinished-codon string). No annotation identifiers.
    entries=defaultdict(list);donors=defaultdict(list);completed={}
    starts=set(sites["start"]);acceptors=set(sites["acceptor"]);donor_sites=set(sites["donor"])

    def prune(states,stage,position):
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
        if observer is not None: observer.prune_event(stage,position,groups)
        return result

    for position in sorted(starts|acceptors|donor_sites):
        if position in donor_sites:
            for score,chain,carry in prune(donors.pop(position,[]),'donor',position):
                for a,link_score in outgoing[position]:
                    entries[a].append((score+link_score+boundary_scores["acceptor"][a],chain,carry))
        if position not in starts and position not in acceptors: continue
        incoming=entries.pop(position,[])
        if position in starts: incoming.append((boundary_scores["start"][position],(),""))
        for score,chain,carry in prune(incoming,'entry',position):
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
    if observer is not None: observer.final_event(ranked,limit)
    counts["final_chains_pruned"]=max(0,len(ranked)-limit)
    ranked=ranked[:limit]
    return {"chains":[list(c) for c,_ in ranked],"proposal_scores":[s for _,s in ranked],
            "links":kept_pairs,"counts":counts,"budget":asdict(budget),"reference_used":False}


class ReferenceTrace:
    """Read-only observer. Reference geometry never controls search or scoring."""
    def __init__(self,anchors):
        self.anchors=anchors
        self.events={}
        self.final={}
        self.targets=defaultdict(list)

    def bind(self,sequence,raw,sites,boundary,coding_prefix,link_scores,orf,budget):
        self.sequence=sequence;self.raw=raw;self.sites=sites
        self.boundary=boundary;self.coding_prefix=coding_prefix
        self.link_scores=link_scores;self.orf=orf;self.budget=budget
        self.kept_links=set()
        outgoing=defaultdict(list)
        for (d,a),score in link_scores.items(): outgoing[d].append((a,float(score)))
        for d,values in outgoing.items():
            self.kept_links.update((d,a) for a,_ in sorted(values,key=lambda x:(-x[1],x[0]))[:budget.intron_links_per_donor])
        for anchor in self.anchors:
            ident=anchor['reference_id'];chain=tuple(tuple(x) for x in anchor['oriented_chain'])
            carry='';self.events[ident]={}
            for j,(a,b) in enumerate(chain):
                self.targets['entry',a].append((ident,j,chain[:j],carry))
                if carry is None: break
                carry=orf.extend(a,b,carry,j==len(chain)-1)
                if j<len(chain)-1:
                    self.targets['donor',b].append((ident,j,chain[:j+1],carry))
            self.final[ident]={'present_before_final':False,'survived':False}

    def describe(self,state,stage,position):
        score,chain,carry=state
        first=chain[0][0] if chain else position
        ep={'start':self.boundary['start'][first],'acceptor':0.,'donor':0.,'stop':0.}
        region=0.;link=0.;parts=[]
        for j,(a,b) in enumerate(chain):
            kind='stop' if stage=='final' and j==len(chain)-1 else 'donor'
            endpoint=self.boundary[kind][b]
            avg=float((self.coding_prefix[b]-self.coding_prefix[a])/(b-a))
            ep[kind]+=endpoint;region+=avg
            piece={'exon':[a,b],'end_kind':kind,'end_endpoint':endpoint,'CDS_region_mean_log_odds':avg}
            if j:
                d=chain[j-1][1];v=float(self.link_scores[d,a])
                link+=v;ep['acceptor']+=self.boundary['acceptor'][a]
                piece.update(intron=[d,a],link=v,acceptor=self.boundary['acceptor'][a])
            parts.append(piece)
        if stage=='entry' and chain:
            d=chain[-1][1];v=float(self.link_scores[d,position])
            link+=v;ep['acceptor']+=self.boundary['acceptor'][position]
            parts.append({'pending_intron':[d,position],'link':v,'acceptor':self.boundary['acceptor'][position]})
        total=sum(ep.values())+region+link
        if abs(total-score)>1e-6: raise ValueError('Additive decomposition does not reconcile')
        return {'prefix':[list(x) for x in chain],'carry':carry,'proposal_score':float(score),
                'endpoint_components':ep,'CDS_region_sum':region,'link_sum':link,
                'components_sum':total,'parts':parts}

    def prune_event(self,stage,position,groups):
        for ident,j,prefix,carry in self.targets.get((stage,position),[]):
            values=groups.get(carry,[])
            rank=next((i+1 for i,s in enumerate(values) if s[1]==prefix),None)
            record={'stage':stage,'position':position,'exon_index':j,'carry':carry,
                    'present_before_beam':rank is not None,'rank':rank,
                    'same_carry_competitors':len(values),'limit':self.budget.beam_per_carry,
                    'survived':rank is not None and rank<=self.budget.beam_per_carry}
            if rank is not None and not record['survived']:
                target=values[rank-1];kept=values[:self.budget.beam_per_carry]
                record.update(target=self.describe(target,stage,position),
                              cutoff_gap=float(kept[-1][0]-target[0]),
                              surviving_competitors=[self.describe(s,stage,position) for s in kept])
            self.events[ident][stage,j]=record

    def final_event(self,ranked,limit):
        for anchor in self.anchors:
            ident=anchor['reference_id'];chain=tuple(tuple(x) for x in anchor['oriented_chain'])
            rank=next((i+1 for i,(c,_) in enumerate(ranked) if c==chain),None)
            result={'present_before_final':rank is not None,'rank':rank,'limit':limit,
                    'complete_unique_chains':len(ranked),'survived':rank is not None and rank<=limit}
            if rank is not None and rank>limit:
                score=ranked[rank-1][1]
                result.update(target=self.describe((score,chain,''),'final',None),
                              cutoff_gap=float(ranked[limit-1][1]-score),
                              surviving_competitors=[self.describe((s,c,''),'final',None) for c,s in ranked[:limit]])
            self.final[ident]=result

    def conditional_edge(self,a,b,carry,first,terminal):
        """Hypothetical reference carry, NOT a replayed actual search state."""
        kind='stop' if terminal else 'donor'
        base={'exon':[a,b],'carry_in':carry,'end_kind':kind,
              'actual_execution_claim':False,'legal':False,'retained':False,'rank':None}
        if carry is None: return dict(base,reason='invalid_previous_reference_carry')
        if b not in self.sites[kind]: return dict(base,reason='end_endpoint_not_kept')
        choices=[]
        for end in self.sites[kind]:
            if end<=a or ((first or terminal) and end-a<3): continue
            rem=self.orf.extend(a,end,carry,terminal)
            if rem is None: continue
            score=self.boundary[kind][end]+float((self.coding_prefix[end]-self.coding_prefix[a])/(end-a))
            choices.append((end,float(score),rem))
        choices.sort(key=lambda x:(-x[1],x[0]))
        rank=next((i+1 for i,x in enumerate(choices) if x[0]==b),None)
        if rank is None: return dict(base,reason='ORF_or_size_rejected')
        end,score,rem=choices[rank-1]
        result=dict(base,legal=True,rank=rank,retained=rank<=self.budget.ends_per_start_kind,
                    carry_out=rem,edge_proposal_score=score,
                    end_endpoint=self.boundary[kind][b],
                    CDS_region_mean_log_odds=float((self.coding_prefix[b]-self.coding_prefix[a])/(b-a)),
                    limit=self.budget.ends_per_start_kind)
        if not result['retained']:
            result['cutoff_gap']=float(choices[self.budget.ends_per_start_kind-1][1]-score)
            result['surviving_endpoints']=[{'end':x,'edge_proposal_score':s,'carry_out':c}
                                          for x,s,c in choices[:self.budget.ends_per_start_kind]]
        return result

    def results(self):
        results=[]
        for anchor in self.anchors:
            ident=anchor['reference_id'];chain=[tuple(x) for x in anchor['oriented_chain']]
            endpoints=[('start',chain[0][0]),('stop',chain[-1][1])]
            for (_,d),(a,_) in zip(chain,chain[1:]): endpoints.extend([('donor',d),('acceptor',a)])
            ep=[{'kind':k,'position':p,'motif_available':p in self.raw[k],
                 'kept':p in self.sites[k]} for k,p in endpoints]
            links=[{'pair':[d,a],'retained':(d,a) in self.kept_links}
                   for (_,d),(a,_) in zip(chain,chain[1:])]
            edges=[];carry='';grammar=True
            for j,(a,b) in enumerate(chain):
                edges.append(self.conditional_edge(a,b,carry,j==0,j==len(chain)-1))
                if carry is not None: carry=self.orf.extend(a,b,carry,j==len(chain)-1)
                grammar &= carry is not None
            first_cut=None
            def cut(stage,detail):
                nonlocal first_cut
                if first_cut is None: first_cut={'stage':stage,'detail':detail}
            start=ep[0]
            if not start['motif_available']: cut('start_motif',start)
            elif not start['kept']: cut('start_endpoint',start)
            for j,(a,b) in enumerate(chain):
                if first_cut is not None: break
                entry=self.events[ident].get(('entry',j))
                if entry is None or not entry['present_before_beam']: cut('entry_not_reached',entry);break
                if not entry['survived']: cut('entry_beam',entry);break
                if j>=self.budget.max_exons: cut('max_exons',{'exon_index':j});break
                edge=edges[j]
                if not edge['legal']: cut('endpoint_or_ORF_support',edge);break
                if not edge['retained']: cut('exon_edge_cut',edge);break
                if j<len(chain)-1:
                    donor=self.events[ident].get(('donor',j))
                    if donor is None or not donor['present_before_beam']: cut('donor_not_reached',donor);break
                    if not donor['survived']: cut('donor_beam',donor);break
                    if not links[j]['retained']: cut('intron_link',links[j]);break
                    acceptor=next(q for q in ep if q['kind']=='acceptor' and q['position']==chain[j+1][0])
                    if not acceptor['kept']: cut('acceptor_endpoint',acceptor);break
            final=self.final[ident]
            if first_cut is None and not final['present_before_final']: cut('complete_not_reached',final)
            if first_cut is None and not final['survived']: cut('final_chain_cut',final)
            if first_cut and first_cut['stage'] in ('entry_not_reached','donor_not_reached','complete_not_reached'):
                raise ValueError('Unexplained actual trajectory gap for '+ident+': '+first_cut['stage'])
            legal_continuation=grammar and all(q['kept'] for q in ep) and all(q['retained'] for q in links)
            all_edges=all(e['retained'] for e in edges)
            results.append({'reference_id':ident,'stratum':anchor['stratum'],'species':anchor['species'],
                            'strand':anchor['strand'],'oriented_chain':[list(x) for x in chain],
                            'raw_and_kept_endpoints':ep,'necessary_links':links,
                            'reference_ORF_grammar_legal':bool(grammar),
                            'conditional_reference_edges':edges,
                            'actual_beam_events':[v for _,v in sorted(self.events[ident].items(),key=lambda kv:(kv[0][1],kv[0][0]!='entry'))],
                            'final':final,'first_actual_cut':first_cut,
                            'legal_continuation_with_kept_endpoints_links':bool(legal_continuation),
                            'all_conditional_reference_edges_retained':bool(all_edges),
                            'mechanism_witness':bool(legal_continuation and all_edges and first_cut and
                                first_cut['stage'] in ('entry_beam','donor_beam','final_chain_cut'))})
        return results
