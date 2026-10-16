"""DNA-only A9 early versus original greedy+beam paired window records."""
import torch
from .query_decode import decode_query_chains
from .query_beam import decode_query_beam
from .query_prefix_value_decode import decode_prefix_values,complete_query_values

ARMS=('early','late')
SCORE_KIND='A9_max_query_complete_prefix_logit'


@torch.no_grad()
def paired_query_records(model,features,sequence,head,index,ident,max_exons,beam_width,completed_limit,valid_bases=None):
    """Same input/head, different candidate policies; not equal search exposure.

    Keep all completed values, including <=0. Geometry max over query slots is
    explicit; null=0 is used later by the unchanged chromosome selector.
    """
    n=len(sequence) if valid_bases is None else valid_bases
    greedy=decode_query_chains(model,features,sequence,max_exons,n)
    beam=decode_query_beam(model,features,sequence,max_exons,beam_width,completed_limit,n)
    early=decode_prefix_values(model,features,sequence,head,max_exons,beam_width,completed_limit,n)
    paths={arm:[] for arm in ARMS}
    paths['early']=[(r,'prefix4') for q in early for r in q['completed']]
    paths['late']=[(r,'greedy') for r in greedy if r['terminated']]+[(r,'beam4') for q in beam for r in q['completed']]
    sources={arm:{} for arm in ARMS}
    for arm in ARMS:
        for row,origin in paths[arm]:
            key=(row['query_index'],tuple(tuple(e) for e in row['chain']))
            sources[arm].setdefault(key,set()).add(origin)
    union=sorted(sources['early'].keys()|sources['late'].keys())
    # One fresh union STOP readout, mapped into both arms without batch drift.
    values=complete_query_values(model,features,sequence,head,union,n)
    arm_records={}
    for arm in ARMS:
        gains={};origins={};queries={}
        for (q,chain),origin in sources[arm].items():
            value=values[(q,chain)]
            gains[chain]=max(value,gains.get(chain,float('-inf')))
            origins.setdefault(chain,set()).update(origin);queries.setdefault(chain,set()).add(q)
        ordered=sorted(gains)
        arm_records[arm]=dict(window_index=index,window_id=ident,chains=ordered,
            score=[gains[c] for c in ordered],decoder_sources=[sorted(origins[c]) for c in ordered],
            query_indices=[sorted(queries[c]) for c in ordered],score_kind=SCORE_KIND,
            quality_logit_threshold=0.,reference_used=False,
            counts=dict(completed_records=len(paths[arm]),unique_query_geometry=len(sources[arm]),
                        free_candidates=len(ordered),nonpositive_candidates=sum(gains[c]<=0 for c in ordered)))
    shared=sorted(sources['early'].keys()&sources['late'].keys())
    return dict(window_index=index,window_id=ident,arms=arm_records,reference_used=False,
        joint_query_chain_values=[dict(query_index=q,chain=c,complete_value=values[(q,c)],
            arms=[a for a in ARMS if (q,c) in sources[a]]) for q,c in union],
        common_query_chain_records=len(shared),common_query_chain_values_identical=True,
        early_search_counts=[q['counts'] for q in early],late_search_counts=[q['counts'] for q in beam])
