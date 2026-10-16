"""Reference-free chromosome/strand assembly for C0/B1 candidate outputs."""
from bisect import bisect_left
from .c0_decode import genomic_chain
from .core import select_nonoverlapping

def nearest_owner(chain,windows,centers=None):
    centers=centers if centers is not None else [a+b for a,b in windows]
    midpoint=chain[0][0]+chain[-1][1]
    i=bisect_left(centers,midpoint)
    choices=[k for k in (i-1,i) if 0<=k<len(windows)]
    return min(choices,key=lambda k:(abs(centers[k]-midpoint),windows[k][0],k))

def assemble_predictions(records,windows,strand,score_field='score'):
    """One complete chromosome/strand, all windows in increasing genomic start.

    Reference annotations are not inputs. Stage lists permit a separate evaluator
    to measure actual reference-chain attrition. Only owner-window candidates may
    compete; a chain seen exclusively in non-owner windows is recorded as lost.
    """
    windows=[tuple(w) for w in windows]
    if not windows or windows!=sorted(set(windows)) or strand not in ('+','-'):
        raise ValueError('Expected nonempty unique sorted windows and a single strand')
    centers=[a+b for a,b in windows]
    if centers!=sorted(centers): raise ValueError('Window centers must increase')
    before={};owned={};total=owner_records=0
    for record in records:
        i=record['window_index'];a,b=windows[i]
        if len(record['chains'])!=len(record[score_field]):
            raise ValueError('One score required per candidate')
        for chain,score in zip(record['chains'],record[score_field]):
            if not chain or not 0<=chain[0][0]<chain[-1][1]<=b-a:
                raise ValueError('Candidate lies outside its oriented window')
            global_chain=tuple(genomic_chain(chain,a,b,strand))
            total+=1
            before.setdefault(global_chain,None)
            if nearest_owner(global_chain,windows,centers)!=i: continue
            owner_records+=1
            value=float(score)
            old=owned.get(global_chain)
            if old is None or value>old['score']:
                owned[global_chain]={'chain':global_chain,'score':value,'window_index':i}
    pool=[owned[k] for k in sorted(owned)]
    selected=select_nonoverlapping([r['chain'] for r in pool],[r['score'] for r in pool])
    output=[pool[i] for i in selected]
    return {'predictions':output,
            'stages':{'free_unique':sorted(before),'owner_unique':sorted(owned),
                      'selected_unique':[r['chain'] for r in output]},
            'counts':{'free_records':total,'free_unique':len(before),'owner_records':owner_records,
                      'owner_unique':len(owned),'duplicate_owner_records':owner_records-len(owned),
                      'nonowner_only_unique_lost':len(before.keys()-owned.keys()),
                      'selected_unique':len(output),'score_field':score_field}}
