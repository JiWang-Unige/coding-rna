#!/usr/bin/env python3
"""Original TRAIN: fixed complete-core proposal versus matched FRESH and STRUCT9216."""
import argparse
import json
import math
import os
import time

import compare_full_grammar_supervision_recovery as S
import prepare_full_grammar_complete_core_draws as P
import infer_full_grammar_complete_core_draws_train as I

R, G = S.R, S.G
CORE = 'M28_COMPLETE_CORE_DRAWS_final9216'
CONTROLS = (S.FRESH, R.CONTROL)
STRATA = ('both_complete', 'proposal_only_complete', 'control_only_complete', 'neither_complete')
CHECKS = R.ROOT/'outputs/M28-FULL-GRAMMAR-COMPLETE-CORE-RECOVERY-CHECKS'
OUT = R.ROOT/'outputs/M28-FULL-GRAMMAR-COMPLETE-CORE-RECOVERY'
FIELDS = ('reference_chains', 'control_hits', 'candidate_hits', 'retained', 'gained',
          'lost', 'still_missed', 'net_reference_hits')


def stratum(row, control):
    other = row['fresh_single_full_core'] if control == S.FRESH else row['original_single_full_core']
    proposal = row['proposal_single_full_core']
    return ('both_complete' if other and proposal else 'proposal_only_complete' if proposal
            else 'control_only_complete' if other else 'neither_complete')


def describe(rows, candidate, control):
    def grouped(part):
        return dict(all_references=R.pair(part, candidate, control),
            supervision_strata={k:R.pair([r for r in part if stratum(r, control) == k], candidate, control)
                                for k in STRATA})
    value = grouped(rows)
    value['fixed_native_groups'] = {k:grouped([r for r in rows if R.fixed_group(r) == k]) for k in R.GROUPS}
    value['all_three_native_exact_anchor_gap'] = grouped([
        r for r in rows if not r['exact_primary_hit'] and all(r['native_exact_hits'].values())])
    for name, field in (('grammar_support', 'single_chain_supported'),
                        ('population_complete_core', 'population_complete_fine_core_available')):
        value[name] = {k:R.pair([r for r in rows if r[field] == flag], candidate, control)
                       for k,flag in (('available', True), ('unavailable', False))}
    S.verify_partition(value)
    for name in ('grammar_support', 'population_complete_core'):
        for field in FIELDS:
            if sum(part[field] for part in value[name].values()) != value['all_references'][field]:
                raise ValueError('Eligibility strata must retain the whole reference denominator')
    return value


def fixture():
    rows = []
    for coverage, proposal in ((True, True), (False, True), (True, False), (False, False)):
        for a,b in ((True, True), (False, True), (True, False), (False, False)):
            rows.append(dict(original_single_full_core=coverage, fresh_single_full_core=not coverage,
                proposal_single_full_core=proposal, single_chain_supported=True,
                population_complete_fine_core_available=coverage or proposal,
                exact_primary_hit=a, native_exact_hits={m:b for m in R.D.METHODS},
                model_hits={R.CONTROL:a, S.FRESH:a, CORE:b}))
    for control in CONTROLS:
        result = describe(rows, CORE, control)
        assert result['all_references']['reference_chains'] == 16
        assert [result['all_references'][k] for k in ('retained','gained','lost','still_missed')] == [4]*4
        assert all([v[k] for k in ('retained','gained','lost','still_missed')] == [1]*4
                   for v in result['supervision_strata'].values())
    probe = dict(original_single_full_core=False, fresh_single_full_core=True, proposal_single_full_core=True)
    assert stratum(probe, S.FRESH) == 'both_complete'
    assert stratum(probe, R.CONTROL) == 'proposal_only_complete'
    unsupported = dict(rows[-1], single_chain_supported=False,
        original_single_full_core=False, fresh_single_full_core=False,
        proposal_single_full_core=False, population_complete_fine_core_available=False)
    result = describe(rows+[unsupported], CORE, R.CONTROL)
    assert result['all_references']['reference_chains'] == 17
    assert result['grammar_support']['unavailable']['reference_chains'] == 1
    assert result['all_references']['still_missed'] == 5
    assert R.pair([], CORE, R.CONTROL)['reference_recall'] is None
    return 'two_control_coverage_bases_four_hit_states_unsupported_retention_empty_groups'


def join_proposal(rows):
    saved = R.D.completed(P.OUT)
    if ((P.OUT/'JOBID').read_text().strip() != '13462399' or saved['job_id'] != '13462399'
            or saved['source_census_job'] != '13462343' or saved['protocol'] != P.PROTOCOL
            or saved['selected_unique_windows'] != 4107 or saved['prospective_draws'] != 4608
            or saved['DEV_used'] or saved['sealed_data_used'] or saved['optimizer_steps_this_run'] != 0):
        raise ValueError('Fixed completed complete-core preparation required')
    coverage = R.D.read_rows(P.OUT/'reference_coverage.jsonl')
    indexed, original = G.indexed(coverage), G.indexed(rows)
    if len(indexed) != 22688 or set(indexed) != set(original):
        raise ValueError('Proposal coverage must join every primary reference exactly once')
    expected = {'pooled':P.describe(coverage), **{s:P.describe([r for r in coverage if r['species'] == s])
                for s in P.P.Q.ALLOW}}
    if expected != saved['reference_coverage']:
        raise ValueError('Fixed preparation coverage counts must reproduce')
    for key,row in original.items():
        other = indexed[key]
        if any(row[k] != other[k] for k in ('single_chain_supported','original_single_full_core','fresh_single_full_core')):
            raise ValueError('Original and FRESH supervision identities differ from preparation')
        row.update(proposal_single_full_core=other['proposal_full_core_unique_windows'] > 0,
            proposal_full_core_unique_windows=other['proposal_full_core_unique_windows'],
            proposal_full_core_draw_exposures=other['proposal_full_core_draw_exposures'],
            population_complete_fine_core_available=bool(other['population_matching_core_window_ids']),
            population_supervision_category=other['population_supervision_category'])
    return saved


def core_report(anchor, C):
    # Read completion summaries only, never weights. All evaluations close before reference/GFF replay.
    I.source('13462537')
    report = R.D.completed(I.EVAL)
    if ((I.EVAL/'JOBID').read_text().strip() != '13462628' or report['job_id'] != '13462628'
            or report['status'] != 'original_TRAIN_end_to_end_diagnostic_not_DEV_or_generalization'
            or report['checkpoint'] != str(I.CHECKPOINT) or report['reference_policy'] != anchor['reference_policy']
            or report['reference_used_in_inference'] or report['DEV_used'] or report['sealed_data_used']
            or not report['training_eligibility_does_not_filter_reference']
            or not report['unsupported_reference_strata_retained'] or report['optimizer_steps_this_run'] != 0
            or set(report['chromosomes']) != set(C.REFERENCE_COUNTS) or set(report['per_species']) != set(C.SPECIES)):
        raise ValueError('Completed original-six complete-core final9216 evaluation required')
    for species,seqid,length in C.E.I.SCOPE:
        row = report['chromosomes'][seqid]
        if ((row['species'],row['seqid'],row['length'],row['reference_primary_records'],row['exact_chain']['reference'])
                != (species,seqid,length,C.REFERENCE_COUNTS[seqid],C.REFERENCE_COUNTS[seqid])
                or C.reference_metadata(row) != C.reference_metadata(anchor['chromosomes'][seqid])):
            raise ValueError('Fixed reference denominator, background or retained strata differ')
        for field in ('exact_chain','CDS_base_unstranded_union'):
            counts = row[field]
            if counts != C.E.C.S.prf(*(counts[k] for k in ('tp','predicted','reference'))):
                raise ValueError('Completed candidate count arithmetic differs')
        if (row['phase_checked_prediction_records'] != row['prediction_records']
                or row['prediction_records'] != row['exact_chain']['predicted']+row['duplicate_chain_records']):
            raise ValueError('Candidate duplicate or phase-check counts differ')
    for species in C.SPECIES:
        if report['per_species'][species] != C.E.aggregate([
                report['chromosomes'][q] for s,q,_ in C.E.I.SCOPE if s == species]):
            raise ValueError('Candidate species aggregation differs')
    f1s = [report['per_species'][s]['exact_chain']['f1'] for s in C.SPECIES]
    if (report['pooled'] != C.E.aggregate(list(report['chromosomes'].values()))
            or not math.isclose(report['macro_species_exact_chain_F1'],sum(f1s)/len(f1s),rel_tol=1e-12,abs_tol=1e-12)
            or report['worst_species_exact_chain_F1'] != min(f1s)):
        raise ValueError('Candidate pooled, macro or worst-species aggregation differs')
    return report


def replay_core(rows, report, C):
    for species,seqid,length in C.E.I.SCOPE:
        official = report['chromosomes'][seqid]
        prediction = R.ROOT/official['raw_prediction_path']
        if prediction.resolve() != (I.OUT/seqid/'prediction/full_grammar.gff3').resolve():
            raise ValueError('Completed complete-core prediction path differs')
        txs = C.E.C.S.E.primary_transcripts(C.E.C.S.E.parse_annotation(prediction,{seqid:length}),complete_only=False)
        predicted = {(species,)+C.E.C.S.key(t) for t in txs}
        R.assign_hits([r for r in rows if (r['species'],r['seqid']) == (species,seqid)],
                      CORE, predicted, len(txs), official)
    for species in C.SPECIES:
        if R.pair([r for r in rows if r['species'] == species],CORE)['candidate_hits'] != report['per_species'][species]['exact_chain']['tp']:
            raise ValueError('Candidate per-species reference hit replay differs')
    if R.pair(rows,CORE)['candidate_hits'] != report['pooled']['exact_chain']['tp']:
        raise ValueError('Candidate pooled reference hit replay differs')


def main(check_local):
    began = time.perf_counter()
    checked = fixture()
    labels = (R.ANCHOR,R.LOCAL) if check_local else (R.ANCHOR,*CONTROLS)
    C,specs,reports = R.load_reports(labels)
    if not check_local:
        reports[CORE] = core_report(reports[R.ANCHOR],C)
    rows,paired = R.source_records(reports[R.ANCHOR],C)
    supervision = S.join_supervision(rows,reports[R.ANCHOR])
    preparation = join_proposal(rows)
    for label in labels[1:]:
        R.replay(rows,label,reports[label],C)
    candidate,controls = (R.LOCAL,(R.ANCHOR,)) if check_local else (CORE,CONTROLS)
    if not check_local:
        replay_core(rows,reports[CORE],C)
    comparisons = {}
    for control in controls:
        pooled = describe(rows,candidate,control)
        per_species = {s:describe([r for r in rows if r['species'] == s],candidate,control) for s in C.SPECIES}
        for scope,value in [('pooled',pooled)]+list(per_species.items()):
            a = reports[candidate]['pooled'] if scope == 'pooled' else reports[candidate]['per_species'][scope]
            b = reports[control]['pooled'] if scope == 'pooled' else reports[control]['per_species'][scope]
            if value['all_references']['net_reference_hits'] != a['exact_chain']['tp']-b['exact_chain']['tp']:
                raise ValueError('Paired recovery/loss must match whole-evaluation TP differences')
        comparisons[control] = dict(pooled=pooled,per_species=per_species,
            cumulative_optimizer_steps_matched=not check_local,
            control_coverage_basis='fresh_single_full_core' if control == S.FRESH else 'original_single_full_core')
    report = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],plan=P.PROTOCOL,fixture=checked,
        status='completed_LOCAL_routing_rehearsal_not_complete_core_result' if check_local
               else 'completed_matched_complete_core_vs_FRESH_and_STRUCT_original_TRAIN_recovery',
        check_local=check_local,fixture_only=False,formal_comparison_completed=not check_local,
        candidate=candidate,controls=list(controls),source_evaluation_jobs={m:v['job_id'] for m,v in reports.items()},
        source_preparation_job=preparation['job_id'],source_supervision_job=supervision['job_id'],
        source_paired_reference_job=paired['job_id'],reference_primary_chains=len(rows),
        reference_policy=reports[R.ANCHOR]['reference_policy'],fixed_coverage_counts_reproduce_preparation=True,
        exact_hits_and_unique_predictions_reproduce_completed_evaluations=True,
        fixed_strata_unchanged_by_new_predictions=True,comparisons=comparisons,
        measurements={m:dict(cumulative_optimizer_steps=9216 if m == CORE else specs[m][3],
            **{k:reports[m][k] for k in ('checkpoint','chromosomes','per_species','pooled',
                'macro_species_exact_chain_F1','worst_species_exact_chain_F1')}) for m in reports},
        DNA_read=False,real_reference_GFF_read=False,checkpoint_read=False,prediction_GFF_read=True,
        saved_paired_reference_profile_read=True,new_model_inference=False,predictions_modified=False,
        optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,
        scientific_winner_or_DEV_candidate_selected=False,statistical_significance_or_mechanism_causality_claimed=False,
        limitations=[
            'Coverage strata use additional-fit schedules, excluding common starting-fit exposure.',
            'Matched9216 updates and expected uniform loss do not isolate coverage from importance, variance or AdamW trajectories.',
            'Reference subsets measure recall accounting, not subset precision or F1; full-output background and losses remain relevant.',
            'LOCAL routing rehearsal did not train with the proposal; its proposal flags are prospective grouping fixtures only.',
            'Native-gap groups are fixed original outcome-defined diagnostics, not causal strata or independent replicates.',
            'Single-seed original TRAIN results do not establish independent generalization or alter any frozen DEV family.'],
        seconds=time.perf_counter()-began)
    out = CHECKS if check_local else OUT
    with (out/'references.jsonl').open('x') as handle:
        for row in rows:
            handle.write(json.dumps(row,allow_nan=False)+'\n')
    with (out/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(completed=True,status=report['status'],reference_primary_chains=len(rows),
                          formal_comparison_completed=not check_local)),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-local',action='store_true',help='Completed LOCAL routing rehearsal only')
    main(parser.parse_args().check_local)
