"""C0 bridge to the unmodified ANNEVO HMM research dependency.

ANNEVO source is not vendored here; its non-commercial research license and
copyright remain in the separately obtained dependency. This bridge neither
calls its neural model nor reuses the M25 region/run candidate graph.
"""
import importlib.util
import sys
from functools import lru_cache
from pathlib import Path
import numpy as np

# Internal M28 phases increase 0,1,2 along CDS. Native input columns use 0,2,1.
NATIVE_FROM_M28=(0,1,3,2,4,6,5,7,9,8,10,12,11,13,14)
DEFAULT_HMM=Path(__file__).resolve().parents[2]/"refs/repos/annevo-2026/src/HMM.py"

def native_probabilities(probabilities):
    p=np.asarray(probabilities,dtype=np.float32)
    if p.ndim!=2 or p.shape[1]!=15 or not np.isfinite(p).all() or (p<0).any():
        raise ValueError("Expected finite nonnegative [bases,15] probabilities")
    return p[:,NATIVE_FROM_M28].copy()

@lru_cache(maxsize=1)
def load_native(path=DEFAULT_HMM):
    spec=importlib.util.spec_from_file_location("m28_native_annevo_hmm",path)
    module=importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=module
    spec.loader.exec_module(module)
    return module

# Explicit native input columns used by viterbi_decoding, not define_columns order.
NATIVE_EMISSION_COLUMN = {
    "INTERGENIC":0, "CODING_EXON_0":1, "CODING_EXON_1":3, "CODING_EXON_2":2,
    "INTRON_0":4, "INTRON_1":6, "INTRON_2":5,
    "DSS_0":7, "DSS_1":9, "DSS_2":8,
    "ASS_0":10, "ASS_1":12, "ASS_2":11, "START":13, "END":14}

def native_path_objective_inputs(module,states,n,groups,sequence,clipped_native_p,min_intron_internal):
    """Reuse native transition constructors, with the same zero penalties as decode_c0."""
    initial=np.full((n,n),-np.inf,dtype=np.float32)
    initial=module.set_transition_matrix_common_state(initial,states,min_intron_internal,0,0,0,0)
    matrices=module.set_transition_matrix_conditional_state(initial,states,min_intron_internal,0,0,0,0)
    transitions=np.stack([matrices[x] for x in ("A","T","C","G","N")]).astype(np.float32)
    columns=np.empty(n,dtype=np.int32)
    for group,members in groups.items(): columns[members]=NATIVE_EMISSION_COLUMN[group]
    return np.log(clipped_native_p),columns,transitions,module._encode_sequence(sequence)

def native_path_score(path,inputs,lo=0,hi=None):
    """Native objective: initial state intergenic, no emission at index zero."""
    path=np.asarray(path,dtype=np.int32)
    logp,columns,transitions,codes=inputs
    hi=len(path) if hi is None else hi
    t=np.arange(max(1,lo),hi)
    if path[0]!=0 or path[-1]!=0: return -np.inf
    values=transitions[codes[t],path[t-1],path[t]]+logp[t,columns[path[t]]]
    return float(values.astype(np.float64).sum())

def native_chain_gains(path,chains,inputs):
    """Replace one full gene span by intergenic; include incoming/outgoing edges."""
    path=np.asarray(path,dtype=np.int32)
    gains=[]
    logp,columns,transitions,codes=inputs
    for chain in chains:
        a,b=chain[0][0],chain[-1][1]
        t=np.arange(max(1,a),min(len(path),b+1))
        original_from,original_to=path[t-1],path[t]
        removed_from=np.where((t-1>=a)&(t-1<b),0,original_from)
        removed_to=np.where((t>=a)&(t<b),0,original_to)
        actual=transitions[codes[t],original_from,original_to]+logp[t,columns[original_to]]
        removed=transitions[codes[t],removed_from,removed_to]+logp[t,columns[removed_to]]
        if not np.isfinite(actual).all() or not np.isfinite(removed).all():
            raise ValueError("Native chain removal is not a finite legal intergenic alternative")
        gains.append(float(actual.astype(np.float64).sum()-removed.astype(np.float64).sum()))
    return gains

def decode_c0(sequence,probabilities,min_intron_internal=1,native_path=DEFAULT_HMM):
    if len(sequence)!=len(probabilities): raise ValueError("DNA/emission length mismatch")
    if not sequence: return {"chains":[],"partial_paths":0,"gains":[]}
    module=load_native(native_path)
    states,n=module.define_state(min_intron_length=min_intron_internal)
    groups=module.define_columns(states)
    p=native_probabilities(probabilities)
    # Native implementation changes np.seterr and clips its input: isolate both.
    with np.errstate():
        path=module.viterbi_decoding(p,sequence,states,n,groups,min_intron_length=min_intron_internal)
    coding={v for key,values in groups.items() if key not in ("INTERGENIC","INTRON_0","INTRON_1","INTRON_2") for v in values}
    chains=[];current=None;partial=0
    for i,state in enumerate(path):
        if state==states["start0"]:
            if current is not None: partial+=1
            current=[]
        if current is not None and state in coding:
            if current and current[-1][1]==i:
                current[-1]=(current[-1][0],i+1)
            else:
                current.append((i,i+1))
        if state==states["end2"] and current is not None:
            chains.append(current);current=None
    if current is not None: partial+=1
    inputs=native_path_objective_inputs(module,states,n,groups,sequence,p,min_intron_internal)
    gains=native_chain_gains(path,chains,inputs)
    return {"chains":chains,"partial_paths":partial,"gains":gains}

def genomic_chain(chain,window_start,window_end,strand):
    if strand=="+": return [(window_start+a,window_start+b) for a,b in chain]
    return sorted((window_end-b,window_end-a) for a,b in chain)

def owner_window(chain,windows):
    """Nearest window center to CDS-span midpoint; lower-start window wins ties."""
    midpoint=(chain[0][0]+chain[-1][1])/2
    return min(range(len(windows)),key=lambda i:(abs((windows[i][0]+windows[i][1])/2-midpoint),windows[i][0],i))
