#!/usr/bin/env python3
"""Reference-known geometric F1 ceiling of fixed pools; not model performance."""
import csv
import json
from collections import defaultdict
from pathlib import Path
from src.m28.core import select_nonoverlapping

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT/'outputs/M28-A-A2-POOL-CEILING'
ARMS = ('A','A2')
STAGES = ('free_unique','owner_unique','selected_unique')


def ceiling(records, arm, stage):
    scopes = defaultdict(list)
    for row in records:
        flag = row[arm+'_'+stage]
        if flag not in ('True','False'):
            raise ValueError('Exact saved reference-membership flags required')
        if flag == 'True':
            chain = tuple(tuple(x) for x in json.loads(row['CDS_intervals']))
            scopes[row['species'],row['seqid'],row['strand']].append(chain)
    available = sum(len(chains) for chains in scopes.values())
    maximum = sum(len(select_nonoverlapping(chains,[1.]*len(chains))) for chains in scopes.values())
    return {'reference_chains':len(records),'available_exact_reference_chains':available,
            'maximum_nonoverlapping_TP':maximum,'overlap_conflict_loss':available-maximum,
            'F1_ceiling_with_reference_known_perfect_reranking':
                2*maximum/(maximum+len(records)) if records else 0.,
            'not_model_performance':True}


def main():
    if OUT.exists():
        raise FileExistsError('Use the single new geometric-ceiling directory')
    source = ROOT/'outputs/M28-A-A2-REFERENCE-TURNOVER'
    if any((p/'STATUS').read_text().strip() != 'COMPLETED'
           for p in (source,ROOT/'outputs/M28-A4-WINDOW-POOL')):
        raise ValueError('Complete original turnover and fixed A4 readout required')
    saved = json.loads((source/'summary.json').read_text())
    with (source/'references.tsv').open() as f:
        records = list(csv.DictReader(f,delimiter='\t'))
    if len(records) != 7728 or saved['reference_primary_chains'] != 7728:
        raise ValueError('Full unchanged primary-reference denominator required')
    report = {'experiment':'M28-A-A2-POOL-CEILING',
              'scope':'post-A4 descriptive reference-known maximum; not a predictor or independent result',
              'source':str(source.relative_to(ROOT)),'source_checkpoints':saved['source'],
              'reference_used_for_diagnosis':True,'reference_filters_added':False,
              'new_model_inference':False,'test_setaria_used':False,'optimizer_steps_this_run':0,
              'rule':'Reference hits only, unit gains; same-strand span nonoverlap; all other predictions ideally rejected',
              'pooled':{},'per_species':{}}
    for species in saved['per_species']:
        rows = [row for row in records if row['species'] == species]
        report['per_species'][species] = {
            arm:{stage:ceiling(rows,arm,stage) for stage in STAGES} for arm in ARMS}
    for arm in ARMS:
        report['pooled'][arm] = {stage:ceiling(records,arm,stage) for stage in STAGES}
        for stage in STAGES:
            value = report['pooled'][arm][stage]
            counted = saved['pooled']['stages'][stage]['both']+saved['pooled']['stages'][stage][arm+'_only']
            if value['available_exact_reference_chains'] != counted:
                raise ValueError('Geometric pool differs from saved actual reference hits')
        selected = report['pooled'][arm]['selected_unique']
        if selected['maximum_nonoverlapping_TP'] != selected['available_exact_reference_chains']:
            raise ValueError('Original selected TP should be mutually nonoverlapping')
    OUT.mkdir()
    with (OUT/'summary.json').open('x') as f:
        json.dump(report,f,indent=2,allow_nan=False); f.write('\n')
    print(json.dumps(report['pooled'],indent=2,allow_nan=False),flush=True)


if __name__ == '__main__':
    main()
