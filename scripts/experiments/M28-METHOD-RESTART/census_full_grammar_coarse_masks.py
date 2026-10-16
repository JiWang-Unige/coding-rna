#!/usr/bin/env python3
"""Original TRAIN coarse concordance and read-only fine-mask contrasts."""
import argparse
from collections import Counter, defaultdict
import json
import os
import time

import numpy as np
import torch
import census_full_grammar_population_supervision as S

P = S.P
ROOT = P.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-COARSE-MASKS'
CHECKS = ROOT/'outputs/M28-FULL-GRAMMAR-COARSE-MASKS-CHECKS'
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_coarse_mask_census.md'
REGIONS = np.array([0,1,1,1,2,2,2,1,1,1,1,1,1,1,1])
PATTERNS = ['neither_restores', 'explicit_only_but_not_joint',
    'alternatives_only_but_not_joint', 'either_single_but_not_joint',
    'joint_removal_required', 'explicit_removal_sufficient',
    'alternative_removal_sufficient', 'either_single_removal_sufficient']


def primary_only(window, catalog):
    groups = defaultdict(list)
    for ident in window['overlapping_transcript_ids']:
        groups[catalog[ident]['gene_id']].append(ident)
    retained = set()
    for ids in groups.values():
        primary = next((i for i in ids if catalog[i]['training_primary']), None)
        # Do not turn alternative-only coverage into invented background.
        retained.update([primary] if primary is not None else ids)
    return {**window, 'overlapping_transcript_ids':[
        i for i in window['overlapping_transcript_ids'] if i in retained]}


def contrasts(window, catalog, sequence):
    no_explicit = {**window, 'fine_label_unknown_intervals':[]}
    no_alternatives = primary_only(window, catalog)
    joint = {**no_alternatives, 'fine_label_unknown_intervals':[]}
    return [P.O.build_targets(w, catalog, sequence)
            for w in (window, no_explicit, no_alternatives, joint)]


def counts(targets, roles, window, lo, hi):
    old = S.role_counts(targets[0]['seg'], targets[0]['region'], roles, window, lo, hi)
    offset = lo-window['start'] if window['strand'] == '+' else window['end']-hi
    segs = [t['seg'][offset:offset+len(roles)] for t in targets]
    region = targets[0]['region'][offset:offset+len(roles)]
    coarse = (segs[0] < 0) & (region >= 0)
    expected = REGIONS[roles]
    matrix = np.bincount(3*expected[coarse]+region[coarse], minlength=9).reshape(3,3)
    restored = [s[coarse] == roles[coarse] for s in segs[1:]]
    patterns = restored[0].astype(np.int64)+2*restored[1]+4*restored[2]
    histogram = np.bincount(patterns, minlength=8)
    by_region = np.bincount(8*expected[coarse]+patterns, minlength=24).reshape(3,8)
    return dict(old, matching_coarse_rows=int(np.trace(matrix)),
        disagreeing_coarse_rows=int(matrix.sum()-np.trace(matrix)),
        expected_by_observed_region_confusion=matrix.tolist(),
        coarse_rows_by_primary_role=np.bincount(roles[coarse], minlength=15).tolist(),
        restoration_patterns=histogram.tolist(),
        restoration_by_expected_region=by_region.tolist(),
        joint_unmask_valid_fine_disagreement_rows=int(np.count_nonzero(
            coarse & (segs[3] >= 0) & (segs[3] != roles))))


def fixtures():
    from src.m28.labels import REGION_COLUMNS
    for region, roles in enumerate(REGION_COLUMNS):
        assert np.array_equal(REGIONS[list(roles)], np.full(len(roles), region))
    roles = np.array([13,1,4,0,14,7])
    original = dict(seg=np.array([-100,-100,-100,-100,14,-100]),
        region=np.array([1,2,2,1,2,-100]))
    for strand, window, lo, hi, off in (
        ('+', dict(start=10,end=20,strand='+'),12,18,2),
        ('-', dict(start=10,end=20,strand='-'),11,17,3)):
        target = {k:np.pad(v,(off,10-off-len(v)),constant_values=-100)
                  for k,v in original.items()}
        result = counts([target]*4, roles, window, lo, hi)
        assert (result['matching_coarse_rows'],result['disagreeing_coarse_rows'],
                result['matching_fine_rows'],result['fully_unknown_rows']) == (2,2,1,1)
        assert result['restoration_patterns'] == [4,0,0,0,0,0,0,0]
    roles8 = np.ones(8,dtype=np.int64)
    targets = [dict(seg=np.full(8,-100),region=np.ones(8,dtype=np.int64))]
    for bit in (1,2,4):
        targets.append(dict(seg=np.where(np.arange(8)&bit,1,-100),region=np.ones(8,dtype=np.int64)))
    assert counts(targets,roles8,dict(start=0,end=8,strand='+'),0,8)['restoration_patterns'] == [1]*8
    old_fixture = S.fixture()
    genome = 'CCATGCCCTAACCC'
    primary = dict(id='t',gene_id='g',strand='+',CDS=[(2,11,0)],
        training_primary=True,CDS_complete_supported=True)
    window = dict(start=0,end=14,width=14,valid_bases=14,strand='+',
        overlapping_transcript_ids=['t'],positive_ids=['t'],
        fine_label_unknown_intervals=[],chain_unknown_intervals=[],
        no_annotated_gene_either_strand=False)
    roles = np.array([13,13,13,1,2,3,14,14,14])
    alternative = {**primary,'id':'alt','training_primary':False,'CDS_complete_supported':False}
    catalog = {'t':primary,'alt':alternative}
    cases = [
        ({**window,'fine_label_unknown_intervals':[(2,4)]},catalog,{5:2}),
        ({**window,'overlapping_transcript_ids':['t','alt']},catalog,{6:9}),
        ({**window,'overlapping_transcript_ids':['t','alt'],
            'fine_label_unknown_intervals':[(2,4)]},catalog,{4:2,6:7}),
        (window,{'t':{**primary,'CDS_complete_supported':False}},{0:9})]
    for w,c,expected in cases:
        actual = counts(contrasts(w,c,genome),roles,w,2,11)
        assert actual['disagreeing_coarse_rows'] == 0
        assert actual['restoration_patterns'] == [expected.get(i,0) for i in range(8)]
    # An alternative-only gene must remain protected, not become background.
    protected = {**window,'overlapping_transcript_ids':['alt'],'positive_ids':[]}
    assert primary_only(protected,catalog)['overlapping_transcript_ids'] == ['alt']
    return dict(cases=9,existing_geometry_fixture=old_fixture,
        checks='role_mapping_fine_first_both_strands_eight_patterns_explicit_alternative_joint_unsupported_alternative_only')


def describe(rows):
    cores = [c for r in rows for c in r['cores']]
    confusion = sum((np.array(c['expected_by_observed_region_confusion'],dtype=np.int64)
                     for c in cores),np.zeros((3,3),dtype=np.int64))
    patterns = sum((np.array(c['restoration_patterns'],dtype=np.int64)
                    for c in cores),np.zeros(8,dtype=np.int64))
    by_region = sum((np.array(c['restoration_by_expected_region'],dtype=np.int64)
                     for c in cores),np.zeros((3,8),dtype=np.int64))
    role_hist = sum((np.array(c['coarse_rows_by_primary_role'],dtype=np.int64)
                     for c in cores),np.zeros(15,dtype=np.int64))
    coarse = sum(c['coarse_only_rows'] for c in cores)
    if (int(confusion.sum()) != coarse or int(patterns.sum()) != coarse
            or int(role_hist.sum()) != coarse or not np.array_equal(by_region.sum(0),patterns)
            or not np.array_equal(by_region.sum(1),confusion.sum(1))):
        raise ValueError('Every original coarse-only reference-core row must be counted once per partition')
    with_coarse = [r for r in rows if any(c['coarse_only_rows'] for c in r['cores'])]
    return dict(reference_chains=len(rows),supported=sum(r['single_chain_supported'] for r in rows),
        source_supervision_categories=dict(Counter(r['population_supervision_category'] for r in rows)),
        unassessed_categories=dict(Counter(r['population_supervision_category'] for r in rows if not r['cores'])),
        assessed_references=sum(bool(r['cores']) for r in rows),containing_reference_core_pairs=len(cores),
        references_with_coarse_rows=len(with_coarse),
        references_with_any_coarse_disagreement=sum(any(c['disagreeing_coarse_rows'] for c in r['cores']) for r in rows),
        references_with_all_cores_coarse_disagreeing=sum(bool(r['cores']) and all(
            c['disagreeing_coarse_rows']>0 for c in r['cores']) for r in rows),
        coarse_row_exposures=coarse,matching_coarse_rows=int(np.trace(confusion)),
        disagreeing_coarse_rows=int(confusion.sum()-np.trace(confusion)),
        expected_by_observed_region_confusion=confusion.tolist(),
        coarse_rows_by_primary_role=role_hist.tolist(),restoration_patterns=patterns.tolist(),
        restoration_by_expected_region=by_region.tolist(),
        nonmonotone_restoration_rows=int(patterns[1:4].sum()),
        joint_unmask_valid_fine_disagreement_rows=sum(c['joint_unmask_valid_fine_disagreement_rows'] for c in cores))


def main(fixture_only,checks_job):
    began = time.perf_counter()
    checked = fixtures()
    if fixture_only:
        P.write(CHECKS/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],
            protocol=PROTOCOL,fixture=checked,fixture_only=True,research_data_read=False,
            optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False))
        return
    ready = S.G.D.completed(CHECKS)
    if ready['job_id'] != checks_job or not ready['fixture_only'] or ready['research_data_read']:
        raise ValueError('Completed matching data-free mask checks required')
    source = S.G.D.completed(S.OUT)
    if (source['job_id'] != '13462343' or source['reference_primary_chains'] != 22688
            or source['full_existing_TRAIN_windows'] != 25512 or source['DEV_used'] or source['sealed_data_used']):
        raise ValueError('Completed unchanged original population census required')
    source_rows = list(P.rows(S.OUT/'references.jsonl'))
    if len(S.G.indexed(source_rows)) != 22688 or sum(r['single_chain_supported'] for r in source_rows) != 22479:
        raise ValueError('All original references, including unsupported, required')
    windows,_,catalogs,_ = P.load_population()
    sequences,descriptors = P.F.P.read_train_DNA()
    rows,tasks,roles_by_id = [],defaultdict(list),{}
    identity = ('species','seqid','id','gene_id','strand','CDS','single_chain_supported','population_supervision_category')
    for old in source_rows:
        row = {k:old[k] for k in identity}
        row['cores'] = []
        rows.append(row)
        saved = old['population_containing_core_row_counts']
        if len(saved) != old['population_containing_core_count']:
            raise ValueError('Saved containing-core count differs')
        if saved:
            roles_by_id[id(row)] = P.E.role_path(P.E.oriented_parts(row,len(sequences[row['species']][row['seqid']])))
        for counts_old in saved:
            tasks[counts_old['window_id']].append((row,counts_old))
    if len(tasks) != source['containing_windows_evaluated']:
        raise ValueError('All saved containing windows required')
    P.write(OUT/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],protocol=PROTOCOL,
        source_population_job=source['job_id'],checks_job=checks_job,
        full_existing_TRAIN_windows=len(windows),containing_windows=len(tasks),reference_primary_chains=len(rows)))
    for ordinal,(ident,refs) in enumerate(tasks.items()):
        w = windows[ident]
        raw = sequences[w['species']][w['seqid']][w['start']:w['end']]
        sequence = raw if w['strand'] == '+' else P.F.Q.M.rc(raw)
        targets = contrasts(w,catalogs[w['species']],sequence)
        for row,old in refs:
            if (row['species'],row['seqid'],row['strand']) != (w['species'],w['seqid'],w['strand']):
                raise ValueError('Original reference/core identity differs')
            lo,hi = row['CDS'][0][0],row['CDS'][-1][1]
            result = counts(targets,roles_by_id[id(row)],w,lo,hi)
            if any(result[k] != old[k] for k in ('matching_fine_rows','disagreeing_fine_rows',
                'coarse_only_rows','fully_unknown_rows','complete_matching_fine')):
                raise ValueError('Original fine/coarse/unknown counts must replay before interpreting mask contrasts')
            row['cores'].append(dict(window_id=ident,**result))
        if (ordinal+1)%2048 == 0:
            print(json.dumps(dict(windows=ordinal+1,total=len(tasks),seconds=time.perf_counter()-began)),flush=True)
    pooled = describe(rows)
    per_species = {s:describe([r for r in rows if r['species']==s]) for s in P.Q.ALLOW}
    per_chromosome = {q:describe([r for r in rows if r['seqid']==q]) for _,q,_ in P.E.I.SCOPE}
    by_source_category = {name:describe([r for r in rows if r['population_supervision_category']==name])
                          for name in S.CATEGORIES}
    if pooled['source_supervision_categories'] != source['pooled']['all_references']['categories']:
        raise ValueError('Original five-way reference partition must remain unchanged')
    with (OUT/'references.jsonl').open('x') as handle:
        for row in rows:
            handle.write(json.dumps(row,allow_nan=False)+'\n')
    report = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],protocol=PROTOCOL,
        source_population_job=source['job_id'],checks_job=checks_job,fixture=checked,
        reference_policy=source['reference_policy'],full_existing_TRAIN_windows=len(windows),
        containing_windows_evaluated=len(tasks),reference_primary_chains=len(rows),
        restoration_pattern_order=PATTERNS,restoration_bits='1 remove_explicit_fine_intervals; 2 remove_alternatives; 4 remove_both',
        pooled=pooled,per_species=per_species,per_chromosome=per_chromosome,by_source_category=by_source_category,
        all_original_reference_core_geometry_and_fine_coarse_unknown_counts_replayed=True,
        original_TRAIN_DNA_sources=descriptors,original_training_labels_changed=False,
        counterfactual_targets_only_used_for_diagnostic_counts=True,counterfactual_labels_used_for_training=False,
        actual_reference_GFF_read=False,checkpoint_read=False,feature_tensors_read=False,native_forward=False,
        predictions_modified=False,optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,
        seconds=time.perf_counter()-began,
        limitation='same_annotation_primary_concordance_and_mask_counterfactual_not_independent_truth_gradient_accuracy_or_causality')
    P.write(OUT/'summary.json',report)
    if sum(p.stat().st_size for p in OUT.iterdir() if p.is_file()) > 96*1024*1024:
        raise ValueError('Declared 96 MiB diagnostic artifact ceiling exceeded')
    print(json.dumps(dict(completed=True,seconds=report['seconds'],pooled=pooled),allow_nan=False),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture-only',action='store_true')
    parser.add_argument('--checks-job')
    args = parser.parse_args()
    torch.set_num_threads(2)
    main(args.fixture_only,args.checks_job)
