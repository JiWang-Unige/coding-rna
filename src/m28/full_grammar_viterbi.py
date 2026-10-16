"""Checkpointed maximum-score complete gene-set automaton.

Uses the independently validated boundary state graph, not the stopped A43
endpoint/trace-arena prototype. Scores are FP64 per-base role-minus-background
increments. Ties prefer background continuation, then state-table arrival;
no promise of the endpoint prefix implementation's rounding/tie identity.
"""
import numpy as np
from numba import njit
import torch

from .full_grammar_boundary_compiled import _state_table


@njit(cache=False,fastmath=False,parallel=False)
def _step(alpha,base,row,destinations,roles,background,following,predecessor,gain):
    following.fill(-np.inf); predecessor.fill(0)
    for role in range(15): gain[role]=np.float64(row[role])-np.float64(row[0])
    following[background]=alpha[background]
    predecessor[background]=background*3+1
    for state in range(len(alpha)):
        if alpha[state]==-np.inf: continue
        for edge in range(3):
            target=destinations[base,state,edge]
            if target<0: break
            value=alpha[state]+gain[roles[base,state,edge]]
            if value>following[target]:
                following[target]=value; predecessor[target]=state*3+edge+1


@njit(cache=False,fastmath=False,parallel=False)
def _forward(codes,logits,positions,destinations,roles,background):
    size=destinations.shape[1]
    alpha=np.full(size,-np.inf,dtype=np.float64); alpha[background]=0.
    following=np.empty(size,dtype=np.float64); predecessor=np.empty(size,dtype=np.uint32)
    gain=np.empty(15,dtype=np.float64)
    checkpoints=np.empty((len(positions),size),dtype=np.float64); checkpoints[0]=alpha
    cursor=1
    for i in range(len(codes)):
        _step(alpha,codes[i],logits[i],destinations,roles,background,following,predecessor,gain)
        alpha,following=following,alpha
        if i+1==positions[cursor]:
            checkpoints[cursor]=alpha; cursor+=1
    return alpha[background],checkpoints


@njit(cache=False,fastmath=False,parallel=False)
def _trace(codes,logits,positions,checkpoints,destinations,roles,background):
    size=destinations.shape[1]; path=np.empty(len(codes),dtype=np.uint8)
    tail=background; maximum_bytes=0
    gain=np.empty(15,dtype=np.float64)
    for block in range(len(positions)-2,-1,-1):
        a,b=positions[block],positions[block+1]
        predecessors=np.empty((b-a,size),dtype=np.uint32)
        maximum_bytes=max(maximum_bytes,predecessors.nbytes)
        alpha=checkpoints[block].copy(); following=np.empty(size,dtype=np.float64)
        for i in range(a,b):
            _step(alpha,codes[i],logits[i],destinations,roles,background,following,predecessors[i-a],gain)
            alpha,following=following,alpha
        if alpha[tail]!=checkpoints[block+1,tail]:
            raise ValueError('Selected boundary score changed under exact recomputation')
        for i in range(b-1,a-1,-1):
            token=predecessors[i-a,tail]
            if token==0: raise ValueError('Missing predecessor on selected complete path')
            state=(token-1)//3; edge=(token-1)%3
            if destinations[codes[i],state,edge]!=tail: raise ValueError('Invalid selected predecessor edge')
            path[i]=roles[codes[i],state,edge]; tail=state
    if tail!=background: raise ValueError('Selected path does not begin in background')
    return path,maximum_bytes


@njit(cache=False,fastmath=False,parallel=False)
def _selected_scores(logits,path,starts,ends):
    total=0.; baseline=0.; cursor=0; gains=np.empty(len(starts),dtype=np.float64)
    for i in range(len(path)):
        if cursor<len(starts) and i==starts[cursor]: baseline=total
        total+=np.float64(logits[i,path[i]])-np.float64(logits[i,0])
        if cursor<len(ends) and i+1==ends[cursor]:
            gains[cursor]=total-baseline; cursor+=1
    if cursor!=len(starts): raise ValueError('Incomplete selected gene score intervals')
    return total,gains


@torch.no_grad()
def decode_full_grammar_viterbi(sequence,logits,max_exons,checkpoint_bp=8192):
    sequence=sequence.upper(); n=len(sequence)
    if (not n or set(sequence)-set('ACGTN') or checkpoint_bp<1
            or logits.device.type!='cpu' or logits.dtype not in (torch.float32,torch.float64)
            or tuple(logits.shape)!=(n,15) or not bool(torch.isfinite(logits).all())):
        raise ValueError('Normalized DNA, aligned finite CPU roles and positive checkpoint span required')
    states,destinations,roles,background=_state_table(max_exons)
    if len(states)*3>=2**32: raise ValueError('State table exceeds explicit uint32 predecessor representation')
    lookup=np.full(256,4,dtype=np.int8)
    for i,base in enumerate('ACGTN'): lookup[ord(base)]=i
    codes=lookup[np.frombuffer(sequence.encode('ascii'),dtype=np.uint8)]
    positions=np.array([0]+list(range(checkpoint_bp,n,checkpoint_bp))+[n],dtype=np.int64)
    values=logits.detach().numpy()
    score,checkpoints=_forward(codes,values,positions,destinations,roles,background)
    path,trace_bytes=_trace(codes,values,positions,checkpoints,destinations,roles,background)
    starts=np.flatnonzero((path==13)&np.r_[True,path[:-1]!=13])
    ends=np.flatnonzero((path==14)&np.r_[path[1:]!=14,True])+1
    if len(starts)!=len(ends) or np.any(starts>=ends): raise ValueError('Incomplete selected gene boundary roles')
    coding=(path!=0)&((path<4)|(path>6))
    opening=np.r_[True,~coding[:-1]]; opening[starts]=True
    closing=np.r_[~coding[1:],True]; closing[ends-1]=True
    exon_starts=np.flatnonzero(coding&opening)
    exon_ends=np.flatnonzero(coding&closing)+1
    chains=[]; cursor=0
    for a,b in zip(starts,ends):
        chain=[]
        while cursor<len(exon_starts) and exon_starts[cursor]<b:
            x,y=int(exon_starts[cursor]),int(exon_ends[cursor])
            if not a<=x<y<=b: raise ValueError('Selected exon outside its complete gene')
            chain.append((x,y)); cursor+=1
        if not 1<=len(chain)<=max_exons: raise ValueError('Selected complete chain exceeds exon cap')
        chains.append(tuple(chain))
    if cursor!=len(exon_starts): raise ValueError('Selected coding rows outside complete genes')
    reconstructed,gains=_selected_scores(values,path,starts,ends)
    if reconstructed!=score: raise ValueError('Selected per-base score does not equal complete-path optimum')
    return dict(chains=tuple(chains),gain=float(score),chain_gains=tuple(float(v) for v in gains),
        counts=dict(state_count=len(states),checkpoint_rows=len(positions),checkpoint_bytes=checkpoints.nbytes,
            checkpoint_bp=checkpoint_bp,maximum_block_predecessor_bytes=int(trace_bytes),
            selected_role_path_bytes=path.nbytes,input_logit_bytes=logits.numel()*logits.element_size(),
            forward_and_recomputed_bases=2*n),reference_used=False,
        tie_policy='background_continuation_then_state_table_arrival',
        score_accumulation='FP64_per_base_forward_order')
