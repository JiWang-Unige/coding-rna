#!/usr/bin/env python3
"""Descriptive first-two-pass TRAIN exposure counts, not held-out model recall."""
import json,pathlib,collections
r=[json.loads(x) for x in pathlib.Path('outputs/M28-PILOT-R4/B1/training.jsonl').read_text().splitlines()][:3072]
assert len(r)==3072
a,b=r[:1536],r[1536:]
assert all(x['window_id']==y['window_id'] and x['draw_ordinal']==y['draw_ordinal'] and x['chain_positive']==y['chain_positive'] for x,y in zip(a,b))
summary=[]
for epoch,rows in ((1,a),(2,b)):
    for species in ('arabidopsis_thaliana','oryza_sativa'):
        subset=[x for x in rows if x['species']==species]
        pos=sum(x['chain_positive'] for x in subset)
        inj=sum(x['reference_injected'] for x in subset)
        assert all(0<=x['reference_injected']<=x['chain_positive'] for x in subset)
        summary.append(dict(epoch=epoch,species=species,draws=len(subset),known_positive_chain_exposures=pos,
                            freely_proposed_positive_exposures=pos-inj,injected_positive_exposures=inj,
                            free_exposure_fraction=(pos-inj)/pos,
                            draws_requiring_positive_injection=sum(x['reference_injected']>0 for x in subset),
                            mean_free_candidate_count=sum(x['free_candidates'] for x in subset)/len(subset)))
print(json.dumps(dict(source='outputs/M28-PILOT-R4/B1/training.jsonl: first 3072 completed records',pairs_match=True,denominator='Repeated TRAIN window positive-chain exposures, not unique genes or final-checkpoint DEV recall',rows=summary),indent=2))
