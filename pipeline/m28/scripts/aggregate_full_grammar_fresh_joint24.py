#!/usr/bin/env python3
"""Completed original-TRAIN sampling x encoder 2x2; no fitting or DEV."""
import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import common_ruler as C

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-FRESH-JOINT24-FOUR-CELL'
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_fresh_joint24_plan.md'
SOURCE = (
    ('F00_original_split','OBJECTIVE-CONTRAST-TRAIN-DIAGNOSTIC/STRUCTURED_CONTINUE','13451625','OBJECTIVE-CONTRAST-TRAIN-RAW/STRUCTURED_CONTINUE'),
    ('F10_fresh_split','FRESH-DRAWS-TRAIN-DIAGNOSTIC','13453583','FRESH-DRAWS-TRAIN-RAW'),
    ('F01_original_joint','JOINT24-TRAIN-DIAGNOSTIC','13453780','JOINT24-TRAIN-RAW'),
    ('F11_fresh_joint','FRESH-JOINT24-TRAIN-DIAGNOSTIC',None,'FRESH-JOINT24-TRAIN-RAW'),
)
SPECIES = ('arabidopsis_thaliana','oryza_sativa')
SEQIDS = ('NC_003070.9','NC_003071.7','NC_003076.8','NC_089043.1','NC_089044.1','NC_089046.1')


def checked_prf(value):
    tp,predicted,reference = (value[k] for k in ('tp','predicted','reference'))
    if not (0 <= tp <= min(predicted,reference)):
        raise ValueError('Invalid complete-chain counts')
    expected = dict(precision=tp/predicted if predicted else 0.,
                    recall=tp/reference if reference else 0.,
                    f1=2*tp/(predicted+reference) if predicted+reference else 0.)
    if any(not math.isclose(value[k],v,rel_tol=1e-12,abs_tol=1e-12) for k,v in expected.items()):
        raise ValueError('Count-derived complete-chain PRF required')
    return dict(value)


def contrasts(values):
    a,b,c,d = values
    return dict(encoding_effect_original=difference(c,a),encoding_effect_fresh=difference(d,b),
                fresh_effect_split=difference(b,a),fresh_effect_joint=difference(d,c),
                interaction=100*(d-c-b+a))


def difference(a,b):
    return 100*(a-b)


def main(evaluation_job):
    reports,roots,source_jobs = {},{},{}
    for cell,suffix,job,raw in SOURCE:
        directory = ROOT/'outputs'/('M28-FULL-GRAMMAR-'+suffix)
        value = json.loads((directory/'summary.json').read_text())
        expected_job = job or evaluation_job
        if ((directory/'STATUS').read_text().strip() != 'COMPLETED'
                or (directory/'JOBID').read_text().strip() != expected_job
                or not value['completed'] or value['job_id'] != expected_job
                or value['DEV_used'] or value['sealed_data_used'] or value['optimizer_steps_this_run']
                or set(value['chromosomes']) != set(SEQIDS) or set(value['per_species']) != set(SPECIES)
                or not value['unsupported_reference_strata_retained']
                or value['reference_used_in_inference']):
            raise ValueError('All four complete original-TRAIN sources required')
        if reports and value['reference_policy'] != next(iter(reports.values()))['reference_policy']:
            raise ValueError('Primary reference policy differs')
        reports[cell],roots[cell],source_jobs[cell] = value,ROOT/'outputs'/('M28-FULL-GRAMMAR-'+raw),expected_job
    # The existing ruler checks geometry and exported DNA, but does not itself
    # compare caller FASTA across branches. Do that before claiming matched inputs.
    for seqid in SEQIDS:
        reference = reports[SOURCE[0][0]]['chromosomes'][seqid]
        species,length = reference['species'],reference['length']
        original = C.M.read_allowed(ROOT/'data/m1_screen'/species/'genome.fa',{seqid:length})[seqid]
        normalized = re.sub('[^ACGT]','N',original.upper())
        for cell in reports:
            row = reports[cell]['chromosomes'][seqid]
            if (row['species'] != species or row['length'] != length
                    or row['exact_chain']['reference'] != reference['exact_chain']['reference']):
                raise ValueError('Matched chromosome or reference denominator changed')
            if C.M.read_allowed(roots[cell]/seqid/'input.fasta',{seqid:length})[seqid] != normalized:
                raise ValueError('Four caller inputs must match each other and normalized original DNA')
        del original,normalized
    scopes = {}
    for scope in (*SEQIDS,*SPECIES,'pooled'):
        cells = {}
        for cell,report in reports.items():
            row = report['pooled'] if scope=='pooled' else report['per_species'][scope] if scope in SPECIES else report['chromosomes'][scope]
            exact = checked_prf(row['exact_chain'])
            if scope in SPECIES or scope=='pooled':
                parts = list(report['chromosomes'].values()) if scope=='pooled' else [r for r in report['chromosomes'].values() if r['species']==scope]
                if any(exact[k] != sum(r['exact_chain'][k] for r in parts) for k in ('tp','predicted','reference')):
                    raise ValueError('Species and pooled counts must reconstruct all chromosomes')
            cells[cell] = dict(exact_chain=exact,CDS_base_unstranded_union=row['CDS_base_unstranded_union'],background=row['background'])
        if len({v['exact_chain']['reference'] for v in cells.values()}) != 1:
            raise ValueError('Keep common retained reference denominator')
        scopes[scope] = dict(cells=cells,chain_PRF_contrasts_percentage_points={metric:contrasts([cells[k]['exact_chain'][metric] for k in reports]) for metric in ('precision','recall','f1')})
    for endpoint in ('macro_species','worst_species'):
        values = [sum(r['per_species'][s]['exact_chain']['f1'] for s in SPECIES)/2 if endpoint=='macro_species' else min(r['per_species'][s]['exact_chain']['f1'] for s in SPECIES) for r in reports.values()]
        scopes[endpoint] = dict(chain_F1={k:v for k,v in zip(reports,values)},contrasts_percentage_points=contrasts(values))
    if scopes['pooled']['cells'][SOURCE[0][0]]['exact_chain']['reference'] != 22688:
        raise ValueError('Retain all22688 original primary references')
    report = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],generated_UTC=datetime.now(timezone.utc).isoformat(),
        protocol=PROTOCOL,source_jobs=source_jobs,cell_order=list(reports),scopes=scopes,
        all_four_caller_inputs_and_normalized_original_DNA_directly_equal=True,
        DNA_comparison_policy='original_uppercase_nonACGT_to_N_preserving_length',
        counts_and_species_pooled_PRF_reconciled=True,reference_policy=next(iter(reports.values()))['reference_policy'],
        optimizer_steps_this_run=0,native_forward=False,DEV_used=False,sealed_data_used=False,
        limitation='single-fit descriptive TRAIN interaction on nonlinear PRF scale; sampling changes exposure and gradients; not causal, independent generalization or model selection')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False);handle.write('\n')
    lines=['# Original TRAIN sampling by encoding four cell comparison','',
           '| Scope | Original split F1 % | Fresh split F1 % | Original joint F1 % | Fresh joint F1 % | Interaction pp |',
           '| --- | ---: | ---: | ---: | ---: | ---: |']
    for scope in (*SPECIES,'pooled'):
        s=scopes[scope];values=[100*s['cells'][k]['exact_chain']['f1'] for k in reports]
        lines.append('| '+scope+' | '+' | '.join(f'{v:.4f}' for v in values)+f" | {s['chain_PRF_contrasts_percentage_points']['f1']['interaction']:+.4f} |")
    lines.extend(['',report['limitation']])
    with (OUT/'comparison.md').open('x') as handle:
        handle.write('\n'.join(lines)+'\n')
    print(json.dumps(dict(completed=True,job_id=report['job_id'],source_jobs=source_jobs)),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evaluation-job',required=True)
    main(parser.parse_args().evaluation_job)
