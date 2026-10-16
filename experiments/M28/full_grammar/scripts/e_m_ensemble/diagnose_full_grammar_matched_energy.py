#!/usr/bin/env python3
"""Fixed final9216 STRUCTURED/JOINT original-TRAIN score replay, not new prediction."""
import argparse
from collections import Counter
import json
import math
import os
import time

import torch

import annotate_full_grammar_objective_contrast as S
import annotate_full_grammar_joint24 as J
import common_ruler as C
import diagnose_full_grammar_train_energy as D
import diagnose_full_grammar_energy_terms as T

ROOT = D.I.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-MATCHED-ENERGY-DIAGNOSTIC'
CHECKS = ROOT/'outputs/M28-FULL-GRAMMAR-MATCHED-ENERGY-CHECKS'
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_matched_energy_diagnostic.md'
ARMS = ('STRUCTURED_CONTINUE', 'JOINT24576')
FIT = (ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST/STRUCTURED_CONTINUE',
       ROOT/'outputs/M28-FULL-GRAMMAR-JOINT24-FIT')
RAW = (ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-TRAIN-RAW/STRUCTURED_CONTINUE',
       ROOT/'outputs/M28-FULL-GRAMMAR-JOINT24-TRAIN-RAW')
EVAL = (ROOT/'outputs/M28-FULL-GRAMMAR-OBJECTIVE-CONTRAST-TRAIN-DIAGNOSTIC/STRUCTURED_CONTINUE',
        ROOT/'outputs/M28-FULL-GRAMMAR-JOINT24-TRAIN-DIAGNOSTIC')
EVAL_JOBS = ('13451625', '13453780')
FIT_JOBS = ('13451136', '13453777')


def competing_terms(conflicts):
    return ([math.fsum(r['gain_by_role'][i] for r in conflicts) for i in range(15)],
            [sum(r['rows_by_role'][i] for r in conflicts) for i in range(15)])


def margin_terms(reference_gain, reference_terms, conflicts):
    gains, counts = competing_terms(conflicts)
    total = math.fsum(r['gain'] for r in conflicts)
    margin = reference_gain-total
    difference = [a-b for a, b in zip(reference_terms, gains)]
    if (len(reference_terms) != 15 or not math.isclose(math.fsum(gains), total, rel_tol=1e-9, abs_tol=D.EPSILON)
            or not math.isclose(math.fsum(reference_terms), reference_gain, rel_tol=1e-9, abs_tol=D.EPSILON)
            or not math.isclose(math.fsum(difference), margin, rel_tol=1e-9, abs_tol=D.EPSILON)):
        raise ValueError('Reference, competing roles and margin must reconstruct actual scores')
    return total, gains, counts, margin, difference


def checks():
    began = time.perf_counter()
    a, b, reference = [0.]*15, [0.]*15, [0.]*15
    a[1], a[4], b[1], b[13] = 4., -1., 2., -1.
    reference[1], reference[4], reference[13] = 5., -2., -1.
    counts = [0]+[1]*14
    conflicts = [dict(gain=3., gain_by_role=a, rows_by_role=counts),
                 dict(gain=1., gain_by_role=b, rows_by_role=counts)]
    total, gains, rows, margin, difference = margin_terms(2., reference, conflicts)
    assert total == 4. and gains[1] == 6. and rows[1] == 2 and rows[0] == 0
    assert margin == -2. and difference[1] == -1. and difference[4] == -1.
    empty = margin_terms(2., reference, [])
    assert empty[0] == 0. and empty[3] == 2. and empty[4] == reference
    try:
        margin_terms(3., reference, conflicts)
    except ValueError:
        pass
    else:
        raise AssertionError('Wrong reference role total accepted')
    index = ([dict(parts=[(0, 9)], gain=1.), dict(parts=[(12, 21)], gain=2.)], [9, 21])
    assert len(D.conflicting([(3, 18)], index)) == 2
    assert D.conflicting([(9, 12)], index) == []
    control, candidate = {}, {}
    for n in range(4):
        for mapping, hit in ((control, n in (0,1)), (candidate, n in (0,2))):
            legal = n != 3
            mapping[n] = dict(single_chain_supported=legal, exact_primary_hit=hit,
                category='exact_primary_chain' if hit else ('scoring_miss' if legal else 'unsupported'),
                rows_by_role=counts if legal else None, gain_by_role=reference if legal else None,
                conflicting_selected_gain_by_role=gains if legal else None,
                conflicting_selected_rows_by_role=rows if legal else None,
                replacement_margin_by_role=difference if legal else None)
    paired = paired_summary(control, candidate)
    assert (paired['reference'], paired['gained'], paired['lost'], paired['net_hits']) == (4,1,1,0)
    assert all(group['reference'] == 1 for group in paired['hit_groups'].values())
    assert paired['hit_groups']['both_missed']['supported_reference'] == 0
    try:
        paired_summary(control, {k:v for k,v in candidate.items() if k != 3})
    except ValueError:
        pass
    else:
        raise AssertionError('Dropped unsupported reference accepted')
    D.write_json(CHECKS/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'],
        protocol=PROTOCOL, checks='two_conflicts_role_sums_margin_empty_set_wrong_sum_touching_spans_four_hit_groups_and_unsupported_retention',
        existing_production_score_support_and_minus_geometry_helpers_unchanged=True,
        research_data_read=False, checkpoint_read=False, native_forward=False,
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
        seconds=time.perf_counter()-began))


def completed_sources():
    checked = json.loads((CHECKS/'summary.json').read_text())
    profile = json.loads((D.PROFILE/'summary.json').read_text())
    if ((CHECKS/'STATUS').read_text().strip() != 'COMPLETED' or not checked['completed']
            or checked['protocol'] != PROTOCOL or checked['research_data_read']
            or (D.PROFILE/'STATUS').read_text().strip() != 'COMPLETED' or not profile['completed']
            or profile['DEV_used']):
        raise ValueError('New arithmetic checks and original TRAIN reference profile required')
    evaluations = []
    for i in range(2):
        value = json.loads((EVAL[i]/'summary.json').read_text())
        fitted = json.loads((FIT[i]/'summary.json').read_text())
        checkpoint = FIT[i]/'step_009216.pt'
        if ((EVAL[i]/'STATUS').read_text().strip() != 'COMPLETED' or not value['completed']
                or value['job_id'] != EVAL_JOBS[i] or value['checkpoint'] != str(checkpoint)
                or value['DEV_used'] or value['sealed_data_used'] or value['optimizer_steps_this_run']
                or (FIT[i]/'STATUS').read_text().strip() != 'COMPLETED' or not fitted['completed']
                or fitted['job_id'] != FIT_JOBS[i] or fitted['cumulative_optimizer_steps'] != 9216
                or fitted['DEV_used'] or fitted['sealed_data_used']):
            raise ValueError('The two fixed completed final9216 TRAIN sources are required')
        evaluations.append(value)
    if (evaluations[0]['reference_policy'] != evaluations[1]['reference_policy']
            or any(e['pooled']['exact_chain']['reference'] != 22688 for e in evaluations)):
        raise ValueError('Keep the original retained primary-reference population')
    return profile, evaluations


@torch.inference_mode()
def task(number):
    began = time.perf_counter()
    species, seqid, length = D.I.SCOPE[number]
    out = OUT/seqid
    profile, evaluations = completed_sources()
    refs = []
    with (D.PROFILE/'reference_chain_profile.jsonl').open() as handle:
        for line in handle:
            row = json.loads(line)
            if row['seqid'] == seqid:
                refs.append(row)
    rk = {C.S.key(r) for r in refs}
    if len(refs) != len(rk):
        raise ValueError('Original reference keys must remain unique')
    sequences = [C.M.read_allowed(root/seqid/'input.fasta', {seqid:length})[seqid] for root in RAW]
    if sequences[0] != sequences[1]:
        raise ValueError('Both frozen raw models must use identical original TRAIN DNA')
    sequence = sequences[0]
    del sequences
    locations = [(a, min(a+C.M.WIDTH, length)) for a in C.M.starts(length)]
    D.write_json(out/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], task=number,
        species=species, seqid=seqid, length=length, arms=list(ARMS),
        checkpoints=[str(p/'step_009216.pt') for p in FIT], source_evaluation_jobs=list(EVAL_JOBS),
        protocol=PROTOCOL, native_forward=True, reference_used_by_encoder_or_head=False,
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False))
    backbone, tokenizer = D.I.R.load_native()
    torch.cuda.reset_peak_memory_stats()
    runs = {}
    for i, arm in enumerate(ARMS):
        folder = out/arm
        folder.mkdir()
        raw_root = RAW[i]/seqid
        raw = json.loads((raw_root/'prediction/summary.json').read_text())
        receipt = json.loads((raw_root/'summary.json').read_text())
        checkpoint = FIT[i]/'step_009216.pt'
        if ((raw_root/'STATUS').read_text().strip() != 'COMPLETED' or not raw['inference_complete']
                or not receipt['inference_complete'] or raw['checkpoint'] != str(checkpoint)
                or receipt['checkpoint'] != str(checkpoint) or receipt['source_fit_job'] != FIT_JOBS[i]
                or (receipt['species'], receipt['seqid'], receipt['length']) != (species, seqid, length)
                or receipt['DEV_used'] or receipt['sealed_data_used'] or receipt['optimizer_steps_this_run']
                or raw['reference_used'] or raw['lengths'] != {seqid:length} or len(raw['scopes']) != 2
                or {s['strand'] for s in raw['scopes']} != {'+','-'}
                or any(s['windows'] != [list(x) for x in locations] or s['context_policy'] != J.R.CONTEXT
                       for s in raw['scopes'])):
            raise ValueError('Matching fixed raw receipts and production grid required')
        state = torch.load(checkpoint, map_location='cpu', weights_only=True)
        configuration = S.validate_checkpoint(state, arm) if i == 0 else J.validate_checkpoint(state)
        if configuration != raw['configuration']:
            raise ValueError('Frozen prediction and strict-loaded weight configuration differs')
        model = J.R.M28Core('C0', 2048, 128)
        model.load_state_dict(state['model'], strict=True)
        del state
        model.eval().requires_grad_(False).to('cuda')
        emit = J.R.contig_emissions if i == 0 else J.contig_emissions
        predpath = raw_root/'prediction/full_grammar.gff3'
        predictions = C.S.E.primary_transcripts(C.S.E.parse_annotation(predpath, {seqid:length}), complete_only=False)
        pk = {C.S.key(t) for t in predictions}
        exact = C.S.prf(len(pk & rk), len(pk), len(rk))
        if len(predictions) != len(pk) or exact != evaluations[i]['chromosomes'][seqid]['exact_chain']:
            raise ValueError('Reference and prediction counts must reproduce completed evaluation')
        scores = D.gff_scores(predpath)
        rows, scopes = [], []
        arm_started = time.perf_counter()
        with (folder/'reference_energy.jsonl').open('x') as handle:
            for strand in ('+', '-'):
                start = time.perf_counter()
                joined = emit(model, backbone, tokenizer, sequence, locations, strand, 'cuda')
                encoding_seconds = time.perf_counter()-start
                values, dna = joined['local_potentials'].numpy(), joined['sequence']
                prior = next(s for s in raw['scopes'] if s['strand'] == strand)
                if [list(segment) for segment in joined['genomic_context_segments']] != prior['genomic_context_segments']:
                    raise ValueError('Replayed per-base ownership differs from frozen production')
                selected, max_error = [], 0.
                for tx in predictions:
                    if tx['strand'] != strand:
                        continue
                    parts = D.oriented_parts(tx, length)
                    path = D.role_path(parts)
                    if not D.supported(dna, parts, path):
                        raise ValueError('Frozen selected chain is outside the production graph')
                    gain, terms, counts = D.chain_score(values, parts, path)
                    error = abs(gain-scores[tx['id']])
                    max_error = max(max_error, error)
                    if error > 5.1e-8*max(1., abs(scores[tx['id']]))+D.EPSILON:
                        raise ValueError('Fixed native replay differs from eight-digit GFF gain: '+tx['id'])
                    selected.append(dict(id=tx['id'], parts=parts, gain=gain,
                                         gain_by_role=terms, rows_by_role=counts))
                selected.sort(key=lambda r: r['parts'][0][0])
                if any(a['parts'][-1][1] > b['parts'][0][0] for a,b in zip(selected, selected[1:])):
                    raise ValueError('Saved same-strand spans must be nonoverlapping')
                index = (selected, [r['parts'][-1][1] for r in selected])
                total = math.fsum(r['gain'] for r in selected)
                if not math.isclose(total, prior['gain'], rel_tol=1e-8, abs_tol=D.EPSILON):
                    raise ValueError('Replayed selected-set gain differs from frozen contig receipt')
                for ref in refs:
                    if ref['strand'] != strand:
                        continue
                    parts = D.oriented_parts(ref, length)
                    legal = (parts[0][1]-parts[0][0] >= 3 and parts[-1][1]-parts[-1][0] >= 3
                             and D.supported(dna, parts, D.role_path(parts)))
                    gain = terms = counts = competing = cterms = crows = margin = mterms = None
                    conflicts = []
                    if legal:
                        gain, terms, counts = D.chain_score(values, parts, D.role_path(parts))
                        conflicts = D.conflicting(parts, index)
                        competing, cterms, crows, margin, mterms = margin_terms(gain, terms, conflicts)
                    hit = C.S.key(ref) in pk
                    if hit and not legal:
                        raise ValueError('Exact selected primary must be grammar-supported')
                    row = dict(species=species, seqid=seqid, id=ref['id'], strand=strand,
                        reference_CDS=[[a,b] for a,b,_ in ref['CDS']], CDS_segments=len(parts),
                        category=D.classify(hit, legal, gain, competing), single_chain_supported=legal,
                        original_draw_R4_positive_chain=ref['original_draw_R4_positive_chain'],
                        exact_primary_hit=hit, reference_gain_vs_background=gain,
                        conflicting_selected_ids=[r['id'] for r in conflicts], conflicting_selected_gain=competing,
                        replacement_gain_difference=margin, gain_by_role=terms, rows_by_role=counts,
                        conflicting_selected_gain_by_role=cterms, conflicting_selected_rows_by_role=crows,
                        replacement_margin_by_role=mterms)
                    rows.append(row)
                    handle.write(json.dumps(row, allow_nan=False)+'\n')
                handle.flush()
                scopes.append(dict(strand=strand, context_policy=J.R.CONTEXT, windows=len(locations),
                    selected_chains=len(selected), replayed_selected_gain=total, original_selected_gain=prior['gain'],
                    maximum_chain_score_difference_from_8_significant_digit_GFF=max_error,
                    encoding_seconds=encoding_seconds, seconds=time.perf_counter()-start))
                print(json.dumps(dict(arm=arm, seqid=seqid, completed_strand=strand)), flush=True)
                del joined, values, dna
        if any(p.requires_grad or p.grad is not None for module in (model,backbone) for p in module.parameters()):
            raise ValueError('Diagnostic parameters must remain frozen without gradients')
        run = dict(completed=True, arm=arm, job_id=os.environ['SLURM_JOB_ID'], task=number,
            species=species, seqid=seqid, checkpoint=str(checkpoint), cumulative_optimizer_steps=9216,
            source_evaluation_job=EVAL_JOBS[i], source_raw_job=receipt['job_id'], source_fit_job=FIT_JOBS[i],
            configuration=configuration, reference_policy=evaluations[i]['reference_policy'],
            exact_chain_unchanged=exact, counts=D.summarize(rows), scopes=scopes,
            reference_used_for_posthoc_scoring=True, reference_used_by_encoder_or_head=False,
            predictions_modified=False, decoder_rerun=False, potentials_saved=False,
            optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
            seconds=time.perf_counter()-arm_started)
        if run['counts']['replacement_improvement_counterexamples']:
            run['completed'] = False
        D.write_json(folder/'summary.json', run)
        if not run['completed']:
            raise ValueError('Legal improving replacement: inspect replay before scientific interpretation')
        runs[arm] = run
        del model, predictions, rows, selected
    D.write_json(out/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], task=number,
        species=species, seqid=seqid, length=length, arms=list(ARMS), protocol=PROTOCOL,
        source_profile_job=profile['job_id'], source_evaluation_jobs=list(EVAL_JOBS),
        per_arm={arm:dict(counts=r['counts'], checkpoint=r['checkpoint'], exact_chain_unchanged=r['exact_chain_unchanged'])
                 for arm,r in runs.items()},
        predictions_modified=False, optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(), seconds=time.perf_counter()-began))


def key(row):
    return row['species'], row['seqid'], row['strand'], tuple(tuple(p) for p in row['reference_CDS'])


def paired_summary(control, candidate):
    if set(control) != set(candidate):
        raise ValueError('Both models must retain identical reference geometry')
    groups, migrations = {}, Counter()
    for k in control:
        a, b = control[k], candidate[k]
        if a['single_chain_supported'] != b['single_chain_supported'] or a['rows_by_role'] != b['rows_by_role']:
            raise ValueError('Fixed DNA grammar and reference role rows must remain identical')
        migrations[a['category']+' -> '+b['category']] += 1
        name = ('both_hit' if b['exact_primary_hit'] else 'control_only') if a['exact_primary_hit'] else ('candidate_only' if b['exact_primary_hit'] else 'both_missed')
        groups.setdefault(name, []).append((a,b))
    results = {}
    for name in ('both_hit', 'control_only', 'candidate_only', 'both_missed'):
        pairs = groups.get(name, [])
        legal = [(a,b) for a,b in pairs if a['single_chain_supported']]
        results[name] = dict(reference=len(pairs), supported_reference=len(legal),
            role_groups={g:{arm:dict(
                mean_reference_gain=math.fsum(math.fsum(p[i]['gain_by_role'][r] for r in indices) for p in legal)/len(legal) if legal else None,
                mean_conflicting_selected_gain=math.fsum(math.fsum(p[i]['conflicting_selected_gain_by_role'][r] for r in indices) for p in legal)/len(legal) if legal else None,
                mean_replacement_margin=math.fsum(math.fsum(p[i]['replacement_margin_by_role'][r] for r in indices) for p in legal)/len(legal) if legal else None,
                reference_base_rows=sum(sum(p[i]['rows_by_role'][r] for r in indices) for p in legal),
                conflicting_selected_base_rows=sum(sum(p[i]['conflicting_selected_rows_by_role'][r] for r in indices) for p in legal))
                for i,arm in enumerate(ARMS)} for g,indices in T.GROUPS.items()})
    return dict(reference=len(control), hit_groups=results, category_transitions=dict(migrations),
        gained=results['candidate_only']['reference'], lost=results['control_only']['reference'],
        net_hits=results['candidate_only']['reference']-results['control_only']['reference'])


def aggregate():
    began = time.perf_counter()
    _, evaluations = completed_sources()
    rows = {arm:[] for arm in ARMS}
    jobs = {}
    for number,(species,seqid,_) in enumerate(D.I.SCOPE):
        folder = OUT/seqid
        outer = json.loads((folder/'summary.json').read_text())
        if ((folder/'STATUS').read_text().strip() != 'COMPLETED' or not outer['completed']
                or (outer['task'],outer['species'],outer['seqid']) != (number,species,seqid)
                or outer['DEV_used'] or outer['sealed_data_used'] or outer['optimizer_steps_this_run']):
            raise ValueError('All six fixed TRAIN tasks must complete both models')
        jobs[seqid] = outer['job_id']
        for i,arm in enumerate(ARMS):
            report = json.loads((folder/arm/'summary.json').read_text())
            block = [json.loads(line) for line in (folder/arm/'reference_energy.jsonl').read_text().splitlines()]
            if (not report['completed'] or report['source_evaluation_job'] != EVAL_JOBS[i]
                    or D.summarize(block) != report['counts']
                    or report['exact_chain_unchanged'] != evaluations[i]['chromosomes'][seqid]['exact_chain']
                    or len(block) != len({key(r) for r in block})):
                raise ValueError('Completed reference rows must reproduce their frozen source')
            for r in block:
                if r['single_chain_supported']:
                    T.terms(r)
                    if (not math.isclose(math.fsum(r['conflicting_selected_gain_by_role']), r['conflicting_selected_gain'], rel_tol=1e-9, abs_tol=D.EPSILON)
                            or not math.isclose(math.fsum(r['replacement_margin_by_role']), r['replacement_gain_difference'], rel_tol=1e-9, abs_tol=D.EPSILON)):
                        raise ValueError('Saved role terms do not reconstruct selected conflicts or margin')
            rows[arm].extend(block)
    paired, by_arm = {}, {}
    for scope in ('arabidopsis_thaliana','oryza_sativa','pooled'):
        selected = {arm:[r for r in rows[arm] if scope == 'pooled' or r['species'] == scope] for arm in ARMS}
        paired[scope] = paired_summary(*({key(r):r for r in selected[arm]} for arm in ARMS))
        by_arm[scope] = {arm:dict(counts=D.summarize(selected[arm]), reference_role_terms=T.summarize(selected[arm])) for arm in ARMS}
        for i,arm in enumerate(ARMS):
            expected = evaluations[i]['pooled'] if scope == 'pooled' else evaluations[i]['per_species'][scope]
            value = by_arm[scope][arm]['counts']
            if value['reference'] != expected['exact_chain']['reference'] or value['exact_primary_hits'] != expected['exact_chain']['tp']:
                raise ValueError('Species/pool reference denominator and TP changed')
    if paired['pooled']['net_hits'] != 2230:
        raise ValueError('Paired exact-hit change must reproduce the fixed TRAIN comparison')
    report = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL,
        arms=list(ARMS), source_evaluation_jobs=list(EVAL_JOBS), source_jobs=jobs,
        cumulative_optimizer_steps=9216, reference_policy=evaluations[0]['reference_policy'],
        per_scope=by_arm, paired_reference_outcomes=paired,
        full_primary_measurements={arm:{scope:(e['pooled'] if scope == 'pooled' else e['per_species'][scope])['exact_chain']
                                   for scope in paired} for arm,e in zip(ARMS,evaluations)},
        role_groups={k:list(v) for k,v in T.GROUPS.items()}, near_zero_reporting_tolerance=D.EPSILON,
        predictions_modified=False, native_forward=False, optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
        limitation='descriptive TRAIN score arithmetic; outcome strata, uncalibrated cross-head gain scales; no causal or generalization claim',
        seconds=time.perf_counter()-began)
    D.write_json(OUT/'summary.json', report)
    lines = ['# Fixed final9216 TRAIN energy comparison', '',
             '| Scope | Control TP | Joint TP | Gained | Lost | Net |', '| --- | ---: | ---: | ---: | ---: | ---: |']
    for scope,p in paired.items():
        lines.append(f"| {scope} | {by_arm[scope][ARMS[0]]['counts']['exact_primary_hits']} | {by_arm[scope][ARMS[1]]['counts']['exact_primary_hits']} | {p['gained']} | {p['lost']} | {p['net_hits']:+d} |")
    lines.extend(['', report['limitation']])
    with (OUT/'comparison.md').open('x') as handle:
        handle.write('\n'.join(lines)+'\n')
    print(json.dumps(dict(completed=True, job_id=report['job_id'], paired_counts={s:{k:p[k] for k in ('reference','gained','lost','net_hits')} for s,p in paired.items()})), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('checks','task','aggregate'))
    parser.add_argument('--task', type=int, choices=range(6))
    args = parser.parse_args()
    torch.set_num_threads(2)
    if args.mode == 'checks':
        checks()
    elif args.mode == 'aggregate':
        aggregate()
    elif args.task is None:
        parser.error('--task is required')
    else:
        task(args.task)
