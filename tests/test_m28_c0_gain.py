import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.m28.c0_decode import (decode_c0,load_native,native_probabilities,
    native_path_objective_inputs,native_path_score,native_chain_gains)
from test_m28_c0_decode import example

def test_exported_gene_gain_matches_full_native_path_counterfactual_all_phases():
    module=load_native()
    states,n=module.define_state(min_intron_length=1);groups=module.define_columns(states)
    for motif in ('GT','GC'):
        for phase in (0,1,2):
            dna,p,expected=example(12+phase,motif)
            bg=np.full((8,15),1e-6,np.float32);bg[:,0]=1-14e-6
            sequence=dna+'C'*8+dna
            probs=np.concatenate([p,bg,p])
            shift=len(dna)+8
            chains=expected+[[(a+shift,b+shift) for a,b in expected[0]]]
            native=native_probabilities(probs)
            with np.errstate():
                path=module.viterbi_decoding(native,sequence,states,n,groups,min_intron_length=1)
            inputs=native_path_objective_inputs(module,states,n,groups,sequence,native,1)
            actual=native_path_score(path,inputs)
            gains=native_chain_gains(path,chains,inputs)
            assert np.isfinite(actual)
            for chain,gain in zip(chains,gains):
                removed=np.array(path,dtype=np.int32,copy=True);removed[chain[0][0]:chain[-1][1]]=0
                np.testing.assert_allclose(gain,actual-native_path_score(removed,inputs),rtol=1e-10,atol=1e-8)
                assert gain>0
            result=decode_c0(sequence,probs)
            assert result['chains']==chains
            np.testing.assert_allclose(result['gains'],gains)
    empty=np.full((30,15),1e-6,np.float32);empty[:,0]=1-14e-6
    assert decode_c0('C'*30,empty)['gains']==[]
