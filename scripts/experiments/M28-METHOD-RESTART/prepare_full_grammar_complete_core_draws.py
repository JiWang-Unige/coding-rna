#!/usr/bin/env python3
"""Freeze TRAIN complete-fine-chain draws with uniform-objective importance correction."""
import argparse
from collections import Counter
import json
import os
import time

import numpy as np
import torch
import prepare_full_grammar_fresh_draws as P
import census_full_grammar_population_supervision as C

ROOT = P.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-COMPLETE-CORE-DRAWS-PREPARATION'
CHECKS = ROOT/'outputs/M28-FULL-GRAMMAR-COMPLETE-CORE-DRAWS-PREPARATION-CHECKS'
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_complete_core_draws_plan.md'
CONFIG = dict(P.CONFIG, sampling='half_uniform_half_complete_fine_chain_incidence_with_replacement',
    sampling_seed=1, uniform_mixture_probability=.5, census_job='13462343',
    importance_weight_uniform_within_stratum='1/(stratum_size*proposal_probability)',
    objective_multiplier='importance_weight/known_observation_denominator',
    baseline_or_model_outcomes_used_for_sampling=False,
    existing_DEV_candidate_families_changed=False)


def proposal(pool, eligibility):
    incidence = np.array([len(eligibility[w['id']]['complete_reference_ids']) for w in pool])
    if incidence.sum() <= 0:
        raise ValueError('This fixed proposal requires a complete fine chain in every original stratum')
    q = .5/len(pool) + .5*incidence/incidence.sum()
    if not np.isclose(q.sum(), 1.) or np.any(q <= 0):
        raise ValueError('Proposal must be normalized and retain every natural window')
    return q


def schedule(original, strata, probabilities):
    rng = np.random.default_rng(CONFIG['sampling_seed'])
    result = []
    for epoch in range(3):
        for ordinal, old in enumerate(original):
            key = old['species'], *old['stratum']
            pool, q = strata[key], probabilities[key]
            index = int(rng.choice(len(pool), p=q))
            result.append(dict(additional_step=len(result)+1, pass_index=epoch+1,
                draw_ordinal=ordinal, window_id=pool[index]['id'], species=key[0],
                stratum=list(key[1:]), stratum_size=len(pool),
                stratum_probability=float(q[index]),
                importance_weight=float(1./(len(pool)*q[index]))))
    return result


def fixture():
    pool = [dict(id=str(i)) for i in range(3)]
    eligible = {str(i):dict(complete_reference_ids=list(range(k))) for i,k in enumerate((0,1,3))}
    q = proposal(pool, eligible)
    assert np.allclose(q, [1/6, 7/24, 13/24])
    weight = 1/(len(pool)*q)
    values = np.array([2., 11., -4.])
    assert np.isclose(np.sum(q*weight*values), np.mean(values))
    assert not np.isclose(np.sum(q*values), np.mean(values))
    assert weight[0] == 2 and np.all(weight <= 2)
    old = [dict(species='a', stratum=['q','+']), dict(species='b', stratum=['r','-'])]*4
    strata = {('a','q','+'):pool, ('b','r','-'):pool[:2]}
    probabilities = {key:proposal(group, eligible) for key,group in strata.items()}
    draws = schedule(old, strata, probabilities)
    assert draws == schedule(old, strata, probabilities) and len(draws) == 24
    assert [r['additional_step'] for r in draws] == list(range(1,25))
    for epoch in range(3):
        assert [(r['species'],r['stratum']) for r in draws[epoch*8:(epoch+1)*8]] == [
            (r['species'],r['stratum']) for r in old]
    return 'unequal_strata_positive_background_probability_exact_importance_expectation_fixed_schedule'


def describe(refs):
    return dict(reference_chains=len(refs), supported=sum(r['single_chain_supported'] for r in refs),
        population_complete_core_available=sum(bool(r['population_matching_core_window_ids']) for r in refs),
        original_single_full_core=sum(r['original_single_full_core'] for r in refs),
        fresh_single_full_core=sum(r['fresh_single_full_core'] for r in refs),
        complete_core_proposal_single_full_core=sum(r['proposal_full_core_unique_windows'] > 0 for r in refs),
        expected_proposal_single_full_core=sum(r['expected_proposal_single_full_core_probability'] for r in refs),
        newly_covered_relative_to_old_and_fresh=sum(r['proposal_full_core_unique_windows'] > 0
            and not r['original_single_full_core'] and not r['fresh_single_full_core'] for r in refs),
        old_or_fresh_or_proposal_single_full_core=sum(r['original_single_full_core']
            or r['fresh_single_full_core'] or r['proposal_full_core_unique_windows'] > 0 for r in refs))


def main(fixture_only):
    began = time.perf_counter()
    checked = fixture()
    if fixture_only:
        P.write(CHECKS/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'],
            fixture=checked, configuration=CONFIG, research_data_read=False,
            optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False))
        return
    census = json.loads((C.OUT/'summary.json').read_text())
    if ((C.OUT/'STATUS').read_text().strip() != 'COMPLETED'
            or not census['completed'] or census['job_id'] != CONFIG['census_job']
            or census['reference_primary_chains'] != 22688 or census['full_existing_TRAIN_windows'] != 25512
            or not census['all_saved_original_and_fresh_unique_core_counts_reproduced']
            or census['optimizer_steps_this_run'] != 0 or census['DEV_used'] or census['sealed_data_used']):
        raise ValueError('Completed whole original TRAIN census required')
    windows, strata, catalogs, original = P.load_population()
    eligibility_rows = list(P.rows(C.OUT/'window_eligibility.jsonl'))
    eligibility = {r['window_id']:r for r in eligibility_rows}
    if len(eligibility_rows) != len(windows) or set(eligibility) != set(windows):
        raise ValueError('All original natural windows must remain in the proposal')
    for ident,w in windows.items():
        if any(eligibility[ident][k] != w[k] for k in ('species','seqid','strand','start','end')):
            raise ValueError('Census window identity differs')
    probabilities = {key:proposal(pool, eligibility) for key,pool in strata.items()}
    draws = schedule(original, strata, probabilities)
    counts = Counter(r['window_id'] for r in draws)
    updates = Counter((r['species'],*r['stratum']) for r in draws)
    if len(draws) != 4608 or updates != Counter({key:3*n for key,n in Counter(
            (r['species'],*r['stratum']) for r in original).items()}):
        raise ValueError('Original 4608 update and stratum budgets must remain')
    by_window = {w['id']:float(q) for key,pool in strata.items() for w,q in zip(pool,probabilities[key])}
    P.write(OUT/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        protocol=PROTOCOL, prospective_draws=len(draws), selected_unique_windows=len(counts)))
    source_refs = list(P.rows(C.OUT/'references.jsonl'))
    references = {tuple(C.G.D.record_id(r)):r for r in source_refs}
    if len(references) != 22688 or sum(r['single_chain_supported'] for r in source_refs) != 22479:
        raise ValueError('All primary references including unsupported must remain')
    sequences, descriptors = P.F.P.read_train_DNA()
    stats = {s:Counter() for s in P.Q.ALLOW}
    fine_hist = {s:np.zeros(15,dtype=np.int64) for s in P.Q.ALLOW}
    coarse_hist = {s:np.zeros(3,dtype=np.int64) for s in P.Q.ALLOW}
    weighted_fine = {s:np.zeros(15) for s in P.Q.ALLOW}
    weighted_coarse = {s:np.zeros(3) for s in P.Q.ALLOW}
    views = {}
    for ordinal,(ident,copies) in enumerate(counts.items()):
        w = windows[ident]
        species = w['species']
        chromosome = sequences[species][w['seqid']]
        raw = chromosome[w['start']:w['end']]
        core = raw if w['strand'] == '+' else P.F.Q.M.rc(raw)
        target = P.O.build_targets(w, catalogs[species], core)
        n = w['valid_bases']
        fine,region = target['seg'][:n],target['region'][:n]
        exact = fine >= 0
        coarse = (~exact) & (region >= 0)
        unknown = (~exact) & (~coarse)
        for key in eligibility[ident]['complete_reference_ids']:
            ref = references[tuple(key)]
            lo,hi = ref['CDS'][0][0],ref['CDS'][-1][1]
            roles = P.E.role_path(P.E.oriented_parts(ref,len(chromosome)))
            if not C.role_counts(fine,region,roles,w,lo,hi)['complete_matching_fine']:
                raise ValueError('Selected complete fine event differs from census; do not fit')
        q = by_window[ident]
        iw = 1/(len(strata[P.stratum(w)])*q)
        denominator = max(1,int(exact.sum()+coarse.sum()))
        geometry = P.F.G.context_geometry(w['start'],w['end'],len(chromosome))
        a,b = geometry['genomic_context']
        offset = w['start']-a if w['strand'] == '+' else b-w['end']
        views[ident] = dict(fine_positions=int(exact.sum()),coarse_positions=int(coarse.sum()),
            fully_unknown_positions=int(unknown.sum()),known_observation_denominator=denominator,
            objective_multiplier=iw/denominator,genomic_context=[a,b],core_oriented_interval=[offset,offset+n],
            complete_fine_chain_incidence=len(eligibility[ident]['complete_reference_ids']),
            no_annotated_gene_either_strand=w['no_annotated_gene_either_strand'])
        f = np.bincount(fine[exact],minlength=15)
        c = np.bincount(region[coarse],minlength=3)
        fine_hist[species] += copies*f
        coarse_hist[species] += copies*c
        weighted_fine[species] += copies*iw*f
        weighted_coarse[species] += copies*iw*c
        stats[species].update(unique_windows=1,valid_base_exposures=copies*n,
            unknown_exposures=copies*int(unknown.sum()),importance_weight_sum=copies*iw,
            importance_weight_square_sum=copies*iw*iw,
            importance_weighted_valid_base_exposures=copies*iw*n,
            importance_weighted_unknown_exposures=copies*iw*int(unknown.sum()),
            no_annotated_gene_either_strand_draws=copies*int(w['no_annotated_gene_either_strand']),
            no_complete_fine_chain_draws=copies*int(not eligibility[ident]['complete_reference_ids']))
        if (ordinal+1) % 512 == 0:
            print(json.dumps(dict(checked_unique_windows=ordinal+1,total=len(counts))),flush=True)
    coverage = []
    for ref in source_refs:
        ids = ref['population_matching_core_window_ids']
        key = ref['species'],ref['seqid'],ref['strand']
        p = min(1.,sum(by_window[ident] for ident in ids))
        row = {k:ref[k] for k in ('species','seqid','strand','id','single_chain_supported',
            'population_matching_core_window_ids','population_supervision_category',
            'original_single_full_core','fresh_single_full_core')}
        row.update(proposal_full_core_unique_windows=sum(ident in counts for ident in ids),
            proposal_full_core_draw_exposures=sum(counts[ident] for ident in ids),
            expected_proposal_single_full_core_probability=float(-np.expm1(updates[key]*np.log1p(-p))) if p < 1 else 1.)
        coverage.append(row)
    with (OUT/'reference_coverage.jsonl').open('x') as handle:
        for row in coverage:
            handle.write(json.dumps(row,allow_nan=False)+'\n')
    with (OUT/'training_plan.jsonl').open('x') as handle:
        for draw in draws:
            handle.write(json.dumps(dict(draw,**views[draw['window_id']]),allow_nan=False)+'\n')
    stratum_reports = []
    for key,pool in strata.items():
        q = probabilities[key]
        incidence = [len(eligibility[w['id']]['complete_reference_ids']) for w in pool]
        background = np.array([w['no_annotated_gene_either_strand'] for w in pool],dtype=bool)
        stratum_reports.append(dict(species=key[0],seqid=key[1],strand=key[2],windows=len(pool),
            complete_fine_chain_incidence=sum(incidence),updates=updates[key],
            proposal_probability_sum=float(q.sum()),importance_expectation=float(np.sum(q/(len(pool)*q))),
            importance_weight_min=float(np.min(1/(len(pool)*q))),importance_weight_max=float(np.max(1/(len(pool)*q))),
            background_windows=int(background.sum()),uniform_background_probability=float(background.mean()),
            proposal_background_probability=float(q[background].sum())))
    exposure = {}
    for s in P.Q.ALLOW:
        v = stats[s]
        if int(fine_hist[s].sum()+coarse_hist[s].sum())+v['unknown_exposures'] != v['valid_base_exposures']:
            raise ValueError('Every selected row must remain fine, coarse, or unknown')
        exposure[s] = dict(v,fine_exposures_by_role=fine_hist[s].tolist(),
            coarse_exposures_by_region=coarse_hist[s].tolist(),
            importance_weighted_fine_exposures_by_role=weighted_fine[s].tolist(),
            importance_weighted_coarse_exposures_by_region=weighted_coarse[s].tolist(),
            realized_importance_effective_sample_size=v['importance_weight_sum']**2/v['importance_weight_square_sum'])
    summary = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],configuration=CONFIG,protocol=PROTOCOL,
        fixture=checked,source_census_job=census['job_id'],full_existing_TRAIN_windows=len(windows),
        prospective_draws=len(draws),selected_unique_windows=len(counts),strata=stratum_reports,exposure=exposure,
        reference_coverage={'pooled':describe(coverage),**{s:describe([r for r in coverage if r['species'] == s]) for s in P.Q.ALLOW}},
        by_pass=[dict(pass_index=e+1,unique_windows=len({r['window_id'] for r in draws if r['pass_index'] == e+1})) for e in range(3)],
        all_selected_declared_complete_fine_events_reconstructed=True,
        all_new_exact_fine_events_feasible='NOT_YET_VERIFIED',original_TRAIN_DNA_sources=descriptors,
        labels_changed=False,reference_filters_added=False,prospective_fit_submitted=False,
        native_forward=False,checkpoint_read=False,feature_tensors_read=False,optimizer_steps_this_run=0,
        DEV_used=False,sealed_data_used=False,
        limitation='Importance correction restores expected pre-clipping uniform objective, not finite gradients, AdamW trajectories, model accuracy, or causal isolation.',
        seconds=time.perf_counter()-began)
    P.write(OUT/'summary.json',summary)
    print(json.dumps(summary,allow_nan=False),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture-only',action='store_true')
    args = parser.parse_args()
    torch.set_num_threads(2)
    main(args.fixture_only)
