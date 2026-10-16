#!/usr/bin/env python3
"""Post-readout TRAIN score/grammar diagnosis using existing E/M role arrays."""
import argparse
from collections import Counter
import json
import math
import os
from pathlib import Path
import time

import numpy as np
import torch

import readout_full_grammar_gene_posterior_train as P
import diagnose_full_grammar_train_energy as D
import diagnose_full_grammar_matched_energy as R
import diagnose_full_grammar_energy_terms as T

ROOT = P.ROOT
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_e_m_missed_reference_plan.md'
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-E-M-MISSED-REFERENCE'
CHECKS = ROOT/'outputs/M28-FULL-GRAMMAR-E-M-MISSED-REFERENCE-CHECKS'
EVALUATIONS = {'E_RANK4':'13499823', 'M_DENSE':'13510351'}


def counts(rows):
    categories = Counter(r['category'] for r in rows)
    return dict(reference=len(rows), exact_primary_hits=categories['exact_primary_chain'],
        individually_supported_reference=sum(r['single_chain_supported'] for r in rows),
        categories=dict(categories),
        replacement_improvement_counterexamples=categories['positive_reference_improves_saved_selected_set'],
        strata={name:dict(reference=len(group), categories=dict(Counter(r['category'] for r in group)))
            for name, group in (
                ('single_CDS_segment', [r for r in rows if r['CDS_segments'] == 1]),
                ('multiple_CDS_segments', [r for r in rows if r['CDS_segments'] > 1]),
                ('span_longer_than_native_window_not_excluded', [r for r in rows if r['span_bp'] > 24576]))})


def paired_rows(e, m):
    """Identical retained reference geometry; no cross-model score-scale inference."""
    a, b = ({R.key(r):r for r in rows} for rows in (e, m))
    if len(a) != len(e) or len(b) != len(m) or set(a) != set(b):
        raise ValueError('Both models must retain every identical unique primary reference')
    groups, transitions = Counter(), Counter()
    for key, left in a.items():
        right = b[key]
        if (left['single_chain_supported'] != right['single_chain_supported']
                or left['rows_by_role'] != right['rows_by_role']):
            raise ValueError('Shared DNA grammar and reference role rows must remain identical')
        group = ('both_hit' if right['exact_primary_hit'] else 'E_only') if left['exact_primary_hit'] else (
            'M_only' if right['exact_primary_hit'] else 'both_missed')
        groups[group] += 1
        transitions[left['category']+' -> '+right['category']] += 1
    names = ('both_hit', 'E_only', 'M_only', 'both_missed')
    return dict(reference=len(a), hit_groups={n:groups[n] for n in names},
        M_minus_E_exact_hits=groups['M_only']-groups['E_only'], category_transitions=dict(transitions))


def sources():
    profile = P.read(D.PROFILE/'summary.json')
    if ((D.PROFILE/'STATUS').read_text().strip() != 'COMPLETED' or not profile['completed']
            or profile['job_id'] != '13449803' or not profile['primary_reference_unchanged']
            or profile['reference_policy'] != 'M26_CDS_assessable_longest_eligible_CDS_tie_transcript_ID'
            or profile['DEV_used'] or profile['sealed_data_used']):
        raise ValueError('Completed unchanged original-six primary reference profile required')
    evaluations = {}
    for name in P.MODELS:
        module, _, _ = P.source(name)
        report = P.read(module.EVAL/'summary.json')
        if ((module.EVAL/'STATUS').read_text().strip() != 'COMPLETED' or not report['completed']
                or report['job_id'] != EVALUATIONS[name] or report['checkpoint'] != str(module.CHECKPOINT)
                or report['reference_policy'] != profile['reference_policy']
                or report['DEV_used'] or report['sealed_data_used'] or report['optimizer_steps_this_run']
                or report['reference_used_in_inference'] or not report['unsupported_reference_strata_retained']):
            raise ValueError('Fixed E/M full retained TRAIN evaluation required')
        for _, seqid, _ in P.E.I.SCOPE:
            if report['chromosomes'][seqid]['exact_chain']['reference'] != profile['chromosomes'][seqid]['exact_chain']['reference']:
                raise ValueError('Reference profile denominator must match both fixed models')
        evaluations[name] = report
    if any(e['pooled']['exact_chain']['reference'] != 22688 for e in evaluations.values()):
        raise ValueError('Keep all22688 primary TRAIN references')
    return profile, evaluations


def checks():
    import unittest
    import test_m28_e_m_missed_reference as F
    began = time.perf_counter()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(F))
    if not result.wasSuccessful():
        raise RuntimeError('New E/M paired diagnosis arithmetic failed')
    sources()
    P.write(CHECKS/os.environ['SLURM_JOB_ID']/'summary.json', dict(completed=True,
        job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        original_E_M_evaluation_and_profile_metadata_match=True, research_DNA_read=False,
        real_reference_GFF_read=False, reference_profile_rows_read=False, checkpoint_read=False,
        native_forward=False, optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
        seconds=time.perf_counter()-began))


def diagnose_chromosome(species, seqid, length, reference_count, evaluations, out):
    refs = []
    with (D.PROFILE/'reference_chain_profile.jsonl').open() as handle:
        for line in handle:
            row = json.loads(line)
            if row['seqid'] == seqid:
                refs.append(row)
    rk = {P.C.S.key(r) for r in refs}
    if len(refs) != len(rk) or len(refs) != reference_count:
        raise ValueError('Keep original unique primary profile, including unsupported references')
    dna_by_model = {}
    for name, (module, _) in P.MODELS.items():
        dna_by_model[name] = P.C.M.read_allowed(module.OUT/seqid/'input.fasta', {seqid:length})[seqid]
    if dna_by_model['E_RANK4'] != dna_by_model['M_DENSE']:
        raise ValueError('Shared literal original TRAIN DNA required for grammar comparison')
    sequence = dna_by_model['E_RANK4']; del dna_by_model
    blocks = {}
    for name, (module, _) in P.MODELS.items():
        raw = P.read(module.OUT/seqid/'prediction/summary.json')
        predicted = P.predictions(module, seqid, length, raw)
        pk = {P.C.S.key(t) for t in predicted}
        if P.C.S.prf(len(pk & rk), len(pk), len(rk)) != evaluations[name]['chromosomes'][seqid]['exact_chain']:
            raise ValueError('Fixed full primary endpoint must replay before score diagnosis')
        rows, scopes = [], []
        folder = out/seqid/name; folder.mkdir(parents=True)
        with (folder/'reference_energy.jsonl').open('x') as handle:
            for strand, filename in (('+', 'plus.npy'), ('-', 'minus.npy')):
                dna = sequence if strand == '+' else P.C.M.rc(sequence)
                values = np.load(P.OUT/name/seqid/'emissions'/filename, mmap_mode='r', allow_pickle=False)
                if values.shape != (length, 15) or values.dtype != np.float32:
                    raise ValueError('Exact verified full oriented-contig role array required')
                selected = []
                for tx in predicted:
                    if tx['strand'] != strand:
                        continue
                    parts = D.oriented_parts(tx, length); path = D.role_path(parts)
                    if not D.supported(dna, parts, path):
                        raise ValueError('Original selected chain must be grammar-supported')
                    gain, terms, role_rows = D.chain_score(values, parts, path)
                    if not math.isclose(gain, tx['score'], rel_tol=6e-8, abs_tol=D.EPSILON):
                        raise ValueError('Array score must preserve original eight-digit GFF gain')
                    selected.append(dict(id=tx['id'], parts=parts, gain=gain,
                        gain_by_role=terms, rows_by_role=role_rows))
                selected.sort(key=lambda r:r['parts'][0][0])
                if any(a['parts'][-1][1] > b['parts'][0][0] for a,b in zip(selected, selected[1:])):
                    raise ValueError('Single-reference replacement needs nonoverlapping selected spans')
                index = selected, [r['parts'][-1][1] for r in selected]
                total = math.fsum(r['gain'] for r in selected)
                original_scope = next(s for s in raw['scopes'] if s['strand'] == strand)
                if not math.isclose(total, original_scope['gain'], rel_tol=1e-8, abs_tol=D.EPSILON):
                    raise ValueError('Replayed complete selected-set gain differs')
                for ref in refs:
                    if ref['strand'] != strand:
                        continue
                    parts = D.oriented_parts(ref, length)
                    legal = (parts[0][1]-parts[0][0] >= 3 and parts[-1][1]-parts[-1][0] >= 3
                        and D.supported(dna, parts, D.role_path(parts)))
                    gain = terms = role_rows = competing = cterms = crows = margin = mterms = None
                    conflicts = []
                    if legal:
                        gain, terms, role_rows = D.chain_score(values, parts, D.role_path(parts))
                        conflicts = D.conflicting(parts, index)
                        competing, cterms, crows, margin, mterms = R.margin_terms(gain, terms, conflicts)
                    hit = P.C.S.key(ref) in pk
                    if hit and not legal:
                        raise ValueError('An exact selected reference cannot be outside the literal grammar')
                    row = dict(species=species, seqid=seqid, id=ref['id'], strand=strand,
                        reference_CDS=[(a,b) for a,b,_ in ref['CDS']], CDS_segments=len(parts),
                        span_bp=parts[-1][1]-parts[0][0], exact_primary_hit=hit,
                        single_chain_supported=legal, category=D.classify(hit, legal, gain, competing),
                        reference_gain_vs_background=gain, gain_by_role=terms, rows_by_role=role_rows,
                        conflicting_selected_ids=[r['id'] for r in conflicts], conflicting_selected_gain=competing,
                        conflicting_selected_gain_by_role=cterms, conflicting_selected_rows_by_role=crows,
                        replacement_gain_difference=margin, replacement_margin_by_role=mterms)
                    rows.append(row); handle.write(json.dumps(row, allow_nan=False)+'\n')
                handle.flush()
                scopes.append(dict(strand=strand, source_array=str(P.OUT/name/seqid/'emissions'/filename),
                    selected_set_gain_replayed=True))
                del values, selected, index, dna
        summary = dict(counts=counts(rows), role_terms=T.summarize(rows), scopes=scopes)
        P.write(folder/'summary.json', summary); blocks[name] = rows
    chromosome = dict(species=species, length=length,
        models={name:counts(rows) for name,rows in blocks.items()},
        paired=paired_rows(blocks['E_RANK4'], blocks['M_DENSE']))
    print(json.dumps(dict(completed_original_TRAIN_mechanism=seqid)), flush=True)
    return chromosome, blocks



def main(checks_job, posterior_evaluation_job, posterior_evaluation_dir=None, output_dir=None):
    began = time.perf_counter()
    parent_dir = P.OUT/'evaluation' if posterior_evaluation_dir is None else posterior_evaluation_dir
    out = OUT if output_dir is None else output_dir
    checked = P.read(CHECKS/checks_job/'summary.json')
    if ((CHECKS/checks_job/'STATUS').read_text().strip() != 'COMPLETED' or not checked['completed']
            or checked['job_id'] != checks_job or checked['protocol'] != PROTOCOL
            or checked['tests_run'] != 4 or checked['failures'] or checked['errors']
            or not checked['original_E_M_evaluation_and_profile_metadata_match']):
        raise ValueError('Completed new paired/reference-denominator qualification required')
    parent = P.read(parent_dir/'summary.json')
    if ((parent_dir/'STATUS').read_text().strip() != 'COMPLETED' or not parent['completed']
            or parent['job_id'] != posterior_evaluation_job
            or parent['protocol'] != P.PROTOCOL
            or not parent['all_twelve_posteriors_verified_before_reference_access']
            or not parent['original_exact_chain_counts_replayed'] or parent['predictions_modified']
            or parent['DEV_used'] or parent['sealed_data_used']):
        raise ValueError('All twelve original posteriors and fixed reference ranking must finish first')
    profile, evaluations = sources()
    P.write(out/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL,
        checks_job=checks_job, source_posterior_evaluation_job=posterior_evaluation_job,
        source_profile_job=profile['job_id'], source_evaluation_jobs=EVALUATIONS,
        checkpoint_read=False, native_forward=False, predictions_modified=False,
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False))
    all_rows = {name:[] for name in P.MODELS}; chromosomes = {}
    for species, seqid, length in P.E.I.SCOPE:
        chromosome, blocks = diagnose_chromosome(species, seqid, length,
            parent['retained_primary_reference_counts'][seqid], evaluations, out)
        chromosomes[seqid] = chromosome
        for name, rows in blocks.items():
            all_rows[name].extend(rows)
    scopes = {}
    for scope in ('arabidopsis_thaliana', 'oryza_sativa', 'pooled'):
        blocks = {name:[r for r in rows if scope == 'pooled' or r['species'] == scope] for name, rows in all_rows.items()}
        scopes[scope] = dict(models={name:dict(counts=counts(rows), role_terms=T.summarize(rows)) for name,rows in blocks.items()},
            paired=paired_rows(blocks['E_RANK4'], blocks['M_DENSE']))
        for name, report in evaluations.items():
            expected = (report['pooled'] if scope == 'pooled' else report['per_species'][scope])['exact_chain']
            measured = scopes[scope]['models'][name]['counts']
            if (measured['reference'], measured['exact_primary_hits']) != (expected['reference'], expected['tp']):
                raise ValueError('Full species/pooled reference and exact TP changed')
    counterexamples = sum(counts(rows)['replacement_improvement_counterexamples'] for rows in all_rows.values())
    artifact_bytes = sum(p.stat().st_size for p in out.rglob('*') if p.is_file())
    if artifact_bytes > 256*1024*1024:
        raise ValueError('Declared compact reference-diagnostic output budget exceeded')
    P.write(out/'summary.json', dict(completed=not counterexamples, job_id=os.environ['SLURM_JOB_ID'],
        protocol=PROTOCOL, checks_job=checks_job, source_posterior_evaluation_job=posterior_evaluation_job,
        source_profile_job=profile['job_id'], source_evaluation_jobs=EVALUATIONS,
        chromosomes=chromosomes, scopes=scopes, all_original_primary_references_retained=True,
        fixed_grammar_and_reference_role_rows_shared=True, original_exact_chain_counts_replayed=True,
        unexpected_improving_single_reference_replacements=counterexamples,
        checkpoint_read=False, native_forward=False, real_reference_GFF_read=False,
        reference_profile_rows_read=True, reference_used_for_posthoc_scoring_only=True,
        predictions_modified=False, decoder_rerun=False, optimizer_steps_this_run=0,
        DEV_used=False, sealed_data_used=False, scientific_threshold_changed=False,
        numerical_near_zero_tolerance=D.EPSILON, artifact_bytes_before_summary=artifact_bytes,
        seconds=time.perf_counter()-began,
        interpretation='descriptive_TRAIN_score_and_grammar_miss_not_capacity_LR_causal_or_generalization_claim'))
    if counterexamples:
        raise ValueError('Legal improving replacement: inspect saved MAP or replay before scientific interpretation')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('checks', 'diagnose'))
    parser.add_argument('--checks-job')
    parser.add_argument('--posterior-evaluation-job')
    parser.add_argument('--posterior-evaluation-dir', type=Path)
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args(); torch.set_num_threads(2)
    if args.mode == 'checks':
        checks()
    elif args.checks_job is None or args.posterior_evaluation_job is None:
        parser.error('Both --checks-job and --posterior-evaluation-job required')
    else:
        main(args.checks_job, args.posterior_evaluation_job, args.posterior_evaluation_dir, args.output_dir)
