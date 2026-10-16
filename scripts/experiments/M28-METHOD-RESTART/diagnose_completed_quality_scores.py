#!/usr/bin/env python3
"""Exploratory sign/selection attrition from immutable A6 raw scores and A7 membership."""
import csv
import json
import math
from collections import Counter, defaultdict
import evaluate_chain_quality as F
from src.m28.assembly import nearest_owner
from src.m28.c0_decode import genomic_chain

ROOT = F.Q.ROOT
OUT = ROOT/'outputs/M28-A6-A7-SCORE-ATTRITION'
STATES = ('not_free','lost_owner','nonpositive_gain','positive_unselected','selected')


def state(key, free, pool, selected):
    if key not in free: return 'not_free'
    if key not in pool: return 'lost_owner'
    if selected:
        if pool[key] <= 0: raise ValueError('Selected reference has nonpositive gain')
        return 'selected'
    return 'nonpositive_gain' if pool[key] <= 0 else 'positive_unselected'


def main():
    if OUT.exists(): raise FileExistsError('Do not overwrite an earlier analysis')
    fit = json.loads((F.Q.FIT/'summary.json').read_text())
    infer = json.loads((F.Q.INFER/'summary.json').read_text())
    a6 = json.loads((ROOT/'outputs/M28-A6-CHAIN-QUALITY-EVAL/summary.json').read_text())
    a7 = json.loads((ROOT/'outputs/M28-A7-QUALITY-WINDOW-POOL/summary.json').read_text())
    F.validate_complete(fit,infer)
    if (not a7['completed'] or a7['reference_primary_chains']!=7728 or a7['test_setaria_used'] or
        a7['parent_checkpoint']!=fit['final_checkpoint']):
        raise ValueError('Exact completed A6/A7 sources required')
    with (ROOT/'outputs/M28-A7-QUALITY-WINDOW-POOL/references.tsv').open() as handle:
        refs = list(csv.DictReader(handle,delimiter='\t'))
    if len(refs)!=7728 or len({(r['species'],r['id']) for r in refs})!=7728:
        raise ValueError('Keep every original primary reference')
    OUT.mkdir()
    counts = defaultdict(lambda:defaultdict(Counter)); strata = defaultdict(lambda:defaultdict(Counter))
    diagnostics = []
    for scope in infer['scopes']:
        windows = scope['windows']; centers = [a+b for a,b in windows]
        free,owner = {},{}
        raw_count = owner_records = 0
        for row in F.records(F.Q.INFER/scope['file'],scope):
            index = row['window_index']; a,b = windows[index]
            for parts,value in zip(row['chains'],row['score']):
                key = tuple(genomic_chain(parts,a,b,scope['strand'])); gain = float(value)
                if not math.isfinite(gain): raise ValueError('Nonfinite original score')
                raw_count += 1
                free[key] = max(free.get(key,-math.inf),gain)
                if nearest_owner(key,windows,centers)==index:
                    owner_records += 1
                    owner[key] = max(owner.get(key,-math.inf),gain)
        expected = scope['assembly_counts']
        if (raw_count,len(free),owner_records,len(owner)) != tuple(expected[k] for k in
            ('free_records','free_unique','owner_records','owner_unique')):
            raise ValueError('Original actual candidate/owner pool changed')
        species,strand = scope['species'],scope['strand']
        local = {arm:Counter() for arm in ('A6','A7')}
        for r in refs:
            if (r['species'],r['strand'])!=(species,strand): continue
            key = tuple(tuple(p) for p in json.loads(r['CDS_intervals']))
            out = dict(r)
            out.update(A6_owner_gain=owner.get(key),A7_max_gain=free.get(key))
            for arm,pool in [('A6',owner),('A7',free)]:
                selected = r[arm+'_selected']=='True'
                label = state(key,free,pool,selected)
                out[arm+'_state'] = label
                counts[species][arm][label] += 1; local[arm][label] += 1
                group = 'one' if len(key)==1 else 'two_to_four' if len(key)<=4 else 'five_or_more'
                strata[group][arm][label] += 1
            diagnostics.append(out)
        saved = a6['methods']['A6']['stages'][species+'|'+strand]
        if (local['A6']['not_free'],local['A6']['lost_owner'],
            local['A6']['nonpositive_gain']+local['A6']['positive_unselected'],local['A6']['selected']) != (
            saved['not_in_actual_free_candidates'],saved['lost_at_owner'],
            saved['lost_at_conflict_or_nonpositive_gain'],saved['exact_reference_hits']['selected_unique']):
            raise ValueError('New classification must replay original A6 stage counts')
        saved = a7['methods']['A7']['stages'][species+'|'+strand]
        if (local['A7']['not_free'],local['A7']['lost_owner'],
            local['A7']['nonpositive_gain']+local['A7']['positive_unselected'],local['A7']['selected']) != (
            saved['not_in_actual_free_candidates'],0,
            saved['lost_at_conflict_or_nonpositive_gain'],saved['exact_reference_hits']['selected_unique']):
            raise ValueError('New classification must replay original A7 stage counts')
    if len(diagnostics)!=7728: raise ValueError('A reference was dropped or duplicated')
    normalized = {s:{arm:{k:c[k] for k in STATES} for arm,c in values.items()} for s,values in counts.items()}
    pooled = {arm:{k:sum(values[arm][k] for values in normalized.values()) for k in STATES} for arm in ('A6','A7')}
    report = dict(completed=True,experiment='M28-A6-A7-SCORE-ATTRITION',reference_primary_chains=7728,
        source_inference='outputs/M28-A6-CHAIN-QUALITY-INFER',source_evaluation='outputs/M28-A6-CHAIN-QUALITY-EVAL',
        source_A7='outputs/M28-A7-QUALITY-WINDOW-POOL',source_checkpoint=fit['final_checkpoint'],
        reference_filters_added=False,new_training=False,new_model_inference=False,new_predictions=False,
        optimizer_steps_this_run=0,tau=0.,test_setaria_used=False,
        per_species=normalized,pooled=pooled,
        by_exon_class={g:{arm:{k:c[k] for k in STATES} for arm,c in values.items()} for g,values in strata.items()},
        interpretation='exploratory DEV attrition only; positive_unselected is original positive-gain selector exclusion, not proof of one causal training defect; no threshold sweep or new F1/gate')
    with (OUT/'references.tsv').open('x',newline='') as handle:
        writer = csv.DictWriter(handle,fieldnames=list(diagnostics[0]),delimiter='\t'); writer.writeheader(); writer.writerows(diagnostics)
    with (OUT/'summary.json').open('x') as handle: json.dump(report,handle,indent=2,allow_nan=False)
    print(json.dumps(report,allow_nan=False),flush=True)


if __name__=='__main__': main()
