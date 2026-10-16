#!/usr/bin/env python3
"""Fixed original-TRAIN reference-relative chain geometry; no new predictions."""
from collections import Counter, defaultdict
import json
import os
import time
import unittest

import diagnose_full_grammar_e_m_missed_reference as D
from diagnose_full_grammar_train_energy import PROFILE
import readout_full_grammar_gene_posterior_train as P
import run_full_grammar_e_m_energy_pool_train as Q
from src.m28.full_grammar_chain_geometry import CATEGORIES, PredictionIndex, key

ROOT = P.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-E-M-GEOMETRY-TRAIN'
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_e_m_geometry_train.md'
ARMS = ('E_RANK4', 'M_DENSE', 'ENERGY_POOL', 'CONSENSUS')


def sources():
    profile, parents = D.sources()
    pool = P.read(Q.OUT/'evaluation/summary.json')
    folder = ROOT/'outputs/M28-FULL-GRAMMAR-E-M-CONSENSUS-TRAIN'
    consensus = P.read(folder/'summary.json')
    for path, report, job in ((Q.OUT/'evaluation', pool, '13551150'), (folder, consensus, '13550067')):
        if ((path/'STATUS').read_text().strip() != 'COMPLETED'
                or (path/'JOBID').read_text().strip() != job
                or not report['completed'] or report['job_id'] != job
                or report['reference_policy'] != profile['reference_policy']
                or not report['all_original_primary_references_retained']
                or report['DEV_used'] or report['sealed_data_used'] or report['optimizer_steps_this_run']
                or report['threshold_selected'] or report['frozen_DEV_family_changed']):
            raise ValueError('Completed fixed pool and consensus on the unchanged full reference required')
    if (pool['identity'] != Q.IDENTITY or pool['protocol'] != Q.PROTOCOL
            or pool['raw_jobs'] != ['13551143','13551145','13551146','13551147','13551148','13551149']
            or not pool['full_reference_metadata_matches_all_comparators']
            or pool['original_predictions_modified'] or pool['reference_used_in_inference']
            or not consensus['full_parent_metrics_replayed'] or consensus['predictions_modified']
            or not consensus['consensus_selected_without_reference_or_score']):
        raise ValueError('Use the completed immutable fixed E/M energy pool and exact intersection')
    return profile, dict(parents, ENERGY_POOL=pool['energy_pool'], CONSENSUS=consensus['consensus']), pool


def counts(rows):
    result = {arm:dict(Counter(r['arms'][arm]['category'] for r in rows)) for arm in ARMS}
    return dict(reference=len(rows), categories={arm:{c:result[arm].get(c, 0) for c in CATEGORIES} for arm in ARMS})


def transitions(rows, parent):
    matrix = Counter((r['arms'][parent]['category'], r['arms']['ENERGY_POOL']['category']) for r in rows)
    return dict(reference=len(rows), from_arm=parent, to_arm='ENERGY_POOL',
        category_counts={a:{b:matrix[a, b] for b in CATEGORIES} for a in CATEGORIES},
        exact_gained=sum(n for (a, b), n in matrix.items() if a != CATEGORIES[0] and b == CATEGORIES[0]),
        exact_lost=sum(n for (a, b), n in matrix.items() if a == CATEGORIES[0] and b != CATEGORIES[0]))


def main():
    began = time.perf_counter()
    P.write(OUT/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL,
        arms=ARMS, source_profile_job='13449803', source_evaluation_jobs=dict(D.EVALUATIONS,
            ENERGY_POOL='13551150', CONSENSUS='13550067'),
        research_DNA_read=False, reference_GFF_read=False, checkpoint_read=False,
        native_forward=False, optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False))
    import test_m28_full_grammar_chain_geometry as T
    checked = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(T))
    if not checked.wasSuccessful():
        raise RuntimeError('Strand, interval, internal-CDS or retained-reference geometry failure')
    profile, reports, pool = sources()
    references = defaultdict(list)
    scope = {seqid:(species, length) for species, seqid, length in P.E.I.SCOPE}
    with (PROFILE/'reference_chain_profile.jsonl').open() as handle:
        for line in handle:
            ref = json.loads(line)
            if ref['seqid'] not in scope or ref['species'] != scope[ref['seqid']][0]:
                raise ValueError('Only original-six retained primary reference rows are allowed')
            references[ref['seqid']].append(ref)
    rows, chromosomes = [], {}
    for task, (species, seqid, length) in enumerate(P.E.I.SCOPE):
        refs = references[seqid]
        if len(refs) != len({key(r) for r in refs}) or len(refs) != profile['chromosomes'][seqid]['exact_chain']['reference']:
            raise ValueError('Every unique original primary reference must be retained')
        predictions = {}
        for name, (module, _) in P.MODELS.items():
            predictions[name] = P.C.S.E.primary_transcripts(P.C.S.E.parse_annotation(
                module.OUT/seqid/'prediction/full_grammar.gff3', {seqid:length}), complete_only=False)
        predictions['ENERGY_POOL'] = P.C.S.E.primary_transcripts(P.C.S.E.parse_annotation(
            Q.OUT/seqid/'full_grammar.gff3', {seqid:length}), complete_only=False)
        raw = P.read(Q.OUT/seqid/'summary.json')
        if (raw['job_id'] != pool['raw_jobs'][task] or not raw['completed']
                or (Q.OUT/seqid/'STATUS').read_text().strip() != 'COMPLETED'
                or raw['identity'] != Q.IDENTITY or (raw['seqid'], raw['length']) != (seqid, length)
                or len(predictions['ENERGY_POOL']) != raw['primary_CDS_chain_predictions']):
            raise ValueError('Completed pooled GFF identity and exported count must replay')
        common = {key(t) for t in predictions['E_RANK4']} & {key(t) for t in predictions['M_DENSE']}
        predictions['CONSENSUS'] = [t for t in predictions['E_RANK4'] if key(t) in common]
        indices = {name:PredictionIndex(predictions[name]) for name in ARMS}
        block = [dict(species=species, seqid=seqid, id=r['id'], strand=r['strand'],
            CDS_segments=len(r['CDS']), span_bp=r['CDS'][-1][1]-r['CDS'][0][0],
            arms={name:indices[name].classify(r) for name in ARMS}) for r in refs]
        result = counts(block)
        for name in ARMS:
            replay = P.C.S.prf(result['categories'][name]['exact_primary'], len(predictions[name]), len(refs))
            if replay != reports[name]['chromosomes'][seqid]['exact_chain']:
                raise ValueError('All four existing chromosome exact-chain PRF rows must replay')
        result.update(species=species, seqid=seqid,
            predicted={name:len(predictions[name]) for name in ARMS},
            transitions={name:transitions(block, name) for name in ('E_RANK4','M_DENSE')})
        chromosomes[seqid] = result
        rows.extend(block)
        print(json.dumps(dict(completed_geometry_TRAIN=seqid, reference=len(refs))), flush=True)
    pooled = counts(rows)
    per_species = {species:counts([r for r in rows if r['species'] == species]) for species in
                   ('arabidopsis_thaliana','oryza_sativa')}
    if len(rows) != 22688:
        raise ValueError('Retain all22688 primary TRAIN references including unsupported strata')
    for name in ARMS:
        if (pooled['categories'][name]['exact_primary'] != reports[name]['pooled']['exact_chain']['tp']
                or any(per_species[s]['categories'][name]['exact_primary'] != reports[name]['per_species'][s]['exact_chain']['tp']
                       for s in per_species)):
            raise ValueError('Species and pooled exact hits must match completed original evaluation')
    with (OUT/'reference_geometry.jsonl').open('x') as handle:
        for row in rows:
            handle.write(json.dumps(row, allow_nan=False)+'\n')
    P.write(OUT/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL,
        source_profile_job=profile['job_id'], source_evaluation_jobs=dict(D.EVALUATIONS,
            ENERGY_POOL='13551150', CONSENSUS='13550067'), reference_policy=profile['reference_policy'],
        categories=CATEGORIES, chromosomes=chromosomes, per_species=per_species, pooled=pooled,
        transitions={name:transitions(rows, name) for name in ('E_RANK4','M_DENSE')},
        strata={name:counts(block) for name, block in (
            ('single_CDS_segment', [r for r in rows if r['CDS_segments'] == 1]),
            ('multiple_CDS_segments', [r for r in rows if r['CDS_segments'] > 1]),
            ('span_longer_than_native_window_not_excluded', [r for r in rows if r['span_bp'] > 24576]))},
        tests_run=checked.testsRun, failures=len(checked.failures), errors=len(checked.errors),
        all_four_existing_exact_chain_PRF_rows_replayed=True, all_original_primary_references_retained=True,
        same_strand_prediction_nonoverlap_verified=True, reference_GFF_read=False, reference_profile_rows_read=True,
        research_DNA_read=False, role_arrays_read=False, checkpoint_read=False, native_forward=False,
        optimizer_steps_this_run=0, predictions_modified=False, weights_modified=False,
        DEV_used=False, sealed_data_used=False, threshold_selected=False, frozen_DEV_family_changed=False,
        interpretation='post_outcome_TRAIN_reference_relative_geometry_not_causal_or_independent_validation',
        seconds=time.perf_counter()-began))


if __name__ == '__main__':
    main()
