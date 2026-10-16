#!/usr/bin/env python3
"""Independent primary-pool ranking reconciliation, without DNA or model reads."""
from collections import defaultdict
import json
import math
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT/'outputs/M28-FULL-GRAMMAR-ENERGY-POOL-POSTERIOR-TRAIN'
PROFILE = ROOT/'outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-CHAIN-PROFILE'
JOBS = ('13556785','13556786','13556787','13556788','13556789','13556790')
SEQIDS = ('NC_003070.9','NC_003071.7','NC_003076.8','NC_089043.1','NC_089044.1','NC_089046.1')


def read(path):
    return json.loads(path.read_text())


def same(a, b):
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(same(v,b[k]) for k,v in a.items())
    if isinstance(a, list):
        return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    if isinstance(a, float):
        return math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12)
    return a==b


def primary_check(rows, ref_count, reported):
    n=len(rows); positives=sum(r['primary'] for r in rows); negatives=n-positives
    for column, rank in (('log_probability','posterior'),('gain','internal_gain')):
        groups=defaultdict(lambda:[0,0])
        for row in rows:
            if not math.isfinite(row[column]): raise ValueError('Nonfinite saved query score')
            groups[row[column]][int(row['primary'])]+=1
        below=0; auc=0.
        for score in sorted(groups):
            gn,gp=groups[score]; auc+=gp*(below+gn/2); below+=gn
        total=tp=0; contributions=[]; fractions={}
        for score in sorted(groups,reverse=True):
            gn,gp=groups[score]; total+=gn+gp; tp+=gp
            contributions.append(gp*tp/total)
            for f in (.25,.5,.75,.9,1.):
                if str(f) not in fractions and total>=math.ceil(n*f):
                    fractions[str(f)]=dict(score=score,retained_predictions=total,true_reference_matches=tp,
                        precision=tp/total,recall_within_selected_positive_pool=tp/positives,
                        recall_of_all_primary_references=tp/ref_count,requested_fraction=f,
                        actual_fraction=total/n,ties_kept_together=True)
        actual=dict(predictions=n,positive=positives,negative=negatives,
            prediction_pool_average_precision=math.fsum(contributions)/positives,
            prediction_pool_AUROC=auc/(positives*negatives),retained_fraction_diagnostics=fractions,
            unpredicted_reference_genes_in_AP_or_AUROC=False,threshold_selected=False)
        if not same(actual,reported[rank]): raise ValueError('Independent primary AP/AUC/tie-budget mismatch')
    probabilities=[math.exp(r['log_probability']) for r in rows]
    if any(not 0<=p<=1 for p in probabilities): raise ValueError('Invalid structural probability')
    brier=math.fsum((p-r['primary'])**2 for p,r in zip(probabilities,rows))/n
    if not same(brier,reported['probability']['reference_relative_Brier']): raise ValueError('Brier mismatch')
    for i,saved in enumerate(reported['probability']['fixed_probability_bins']):
        pairs=[(p,r['primary']) for p,r in zip(probabilities,rows) if min(9,int(p*10))==i]
        count=len(pairs); hits=sum(y for _,y in pairs)
        if (count,hits)!=(saved['count'],saved['matches']): raise ValueError('Primary fixed-bin counts differ')
        if count and not same(math.fsum(p for p,_ in pairs)/count,saved['mean_structural_probability']):
            raise ValueError('Primary fixed-bin probability mean differs')


def main():
    report=read(SOURCE/'evaluation/summary.json')
    if (report['job_id']!='13556802' or not report['completed'] or tuple(report['posterior_jobs'])!=JOBS
            or any(report[k] for k in ('DEV_used','sealed_data_used','predictions_modified','weights_modified',
                'native_forward','checkpoint_read','optimizer_steps_this_run','threshold_selected','calibration_fitted'))):
        raise ValueError('Completed immutable fixed TRAIN readout required')
    profile=read(PROFILE/'summary.json')
    if ((PROFILE/'STATUS').read_text().strip()!='COMPLETED' or not profile['completed']
            or profile['job_id']!='13449803' or profile['DEV_used'] or profile['sealed_data_used']):
        raise ValueError('Original completed TRAIN profile required')
    refs=defaultdict(set)
    with (PROFILE/'reference_chain_profile.jsonl').open() as handle:
        for line in handle:
            r=json.loads(line); refs[r['seqid']].add((r['strand'],tuple((a,b) for a,b,_ in r['CDS'])))
    all_rows=[]; by_species=defaultdict(list)
    for task,(seqid,job) in enumerate(zip(SEQIDS,JOBS)):
        folder=SOURCE/seqid/'posterior'; raw=read(folder/'summary.json')
        if ((folder/'STATUS').read_text().strip()!='COMPLETED' or (folder/'JOBID').read_text().strip()!=job
                or not raw['completed'] or raw['job_id']!=job or raw['task']!=task or raw['seqid']!=seqid
                or raw['identity']!=report['identity'] or {s['strand'] for s in raw['records']}!={'+','-'}
                or any(s['supplied_bases']!=raw['length'] for s in raw['records'])):
            raise ValueError('All six actual double-strand readouts required')
        rows=raw['genes']; keys={(r['strand'],tuple(tuple(p) for p in r['genomic_CDS'])) for r in rows}
        if len(keys)!=len(rows) or len(refs[seqid])!=report['retained_primary_reference_counts'][seqid]:
            raise ValueError('Unique query keys and complete reference denominator required')
        for row in rows: row['primary']=(row['strand'],tuple(tuple(p) for p in row['genomic_CDS'])) in refs[seqid]
        primary_check(rows,len(refs[seqid]),report['energy_pool']['chromosomes'][seqid]['exact_primary'])
        by_species[raw['species']].extend(rows); all_rows.extend(rows)
    for species,rows in by_species.items():
        count=sum(report['retained_primary_reference_counts'][q] for q in SEQIDS if
            (q.startswith('NC_003') and species=='arabidopsis_thaliana') or
            (q.startswith('NC_089') and species=='oryza_sativa'))
        primary_check(rows,count,report['energy_pool']['per_species'][species]['exact_primary'])
    if sum(len(refs[q]) for q in SEQIDS)!=22688 or len(all_rows)!=22066 or sum(r['primary'] for r in all_rows)!=16527:
        raise ValueError('Original full pool counts must replay')
    primary_check(all_rows,22688,report['energy_pool']['pooled']['exact_primary'])
    out=ROOT/'outputs/M28-FULL-GRAMMAR-ENERGY-POOL-POSTERIOR-READOUT-VERIFIED'/os.environ['SLURM_JOB_ID']
    out.mkdir(parents=True,exist_ok=False)
    compact=dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],source_evaluation_job='13556802',
        primary_AP_AUROC_fixed_fractions_Brier_bins_independently_reconciled=True,
        any_isoform_metrics_projected_from_completed_evaluation_not_independently_recomputed=True,
        energy_pool=report['energy_pool'],pool_minus_parent=report['pool_minus_parent'],
        reference_counts=report['retained_primary_reference_counts'],program_seconds=report['seconds'],
        DNA_read=False,reference_GFF_read=False,checkpoint_read=False,native_forward=False,
        predictions_modified=False,DEV_used=False,sealed_data_used=False,optimizer_steps_this_run=0)
    with (out/'compact.json').open('x') as handle: json.dump(compact,handle,indent=2,allow_nan=False)
    print(json.dumps(dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],output=str(out/'compact.json'))),flush=True)


if __name__=='__main__': main()
