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

def decode_c0(sequence,probabilities,min_intron_internal=1,native_path=DEFAULT_HMM):
    if len(sequence)!=len(probabilities): raise ValueError("DNA/emission length mismatch")
    if not sequence: return {"chains":[],"partial_paths":0}
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
    return {"chains":chains,"partial_paths":partial}

def genomic_chain(chain,window_start,window_end,strand):
    if strand=="+": return [(window_start+a,window_start+b) for a,b in chain]
    return sorted((window_end-b,window_end-a) for a,b in chain)

def owner_window(chain,windows):
    """Nearest window center to CDS-span midpoint; lower-start window wins ties."""
    midpoint=(chain[0][0]+chain[-1][1])/2
    return min(range(len(windows)),key=lambda i:(abs((windows[i][0]+windows[i][1])/2-midpoint),windows[i][0],i))
