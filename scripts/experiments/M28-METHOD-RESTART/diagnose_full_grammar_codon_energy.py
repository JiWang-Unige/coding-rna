#!/usr/bin/env python3
"""Final9216 TRAIN role/codon replay on frozen predictions, not new inference."""
import argparse
from collections import Counter
import json
import math
import os
import time
import unittest

import numpy as np
import torch

import annotate_full_grammar_joint24 as J
import common_ruler as C
import diagnose_full_grammar_train_energy as D
import diagnose_full_grammar_energy_terms as T
import infer_full_grammar_fresh_joint24_train as B
import infer_full_grammar_codon_train as K
from src.m28.full_grammar_codon import CODONS, CODON_INDEX, centered_codon_scores

ROOT = B.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-CODON-ENERGY-DIAGNOSTIC'
CHECKS = ROOT/'outputs/M28-FULL-GRAMMAR-CODON-ENERGY-CHECKS'
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_codon_energy_diagnostic.md'
ARMS = ('FRESH_JOINT24576', 'CODON_ENERGY')
SOURCES = (B, K)
FIT_JOBS = ('13494895', '13497914')
EVAL_JOBS = ('13494898', '13497946')
PAIRED = ROOT/'outputs/M28-FULL-GRAMMAR-CODON-PAIRED'
PAIRED_JOB = '13497947'


def close(actual, expected):
    if not math.isclose(actual, expected, rel_tol=1e-9, abs_tol=D.EPSILON):
        raise ValueError('Role and codon terms do not reconstruct the fixed score')


def chain_terms(sequence, values, parts, coefficients):
    """Use spliced CDS counts independently of the decoder's closure-edge indices."""
    path = D.role_path(parts)
    a, b = parts[0][0], parts[-1][1]
    positions = np.concatenate([np.arange(x, y) for x, y in parts])
    coding = ''.join(sequence[x:y] for x, y in parts)
    if len(coding) < 6 or len(coding) % 3:
        raise ValueError('A supported complete CDS must have whole start/stop triplets')
    try:
        indices = np.array([CODON_INDEX[coding[i:i+3]] for i in range(0, len(coding), 3)])
    except KeyError as error:
        raise ValueError('Supported coding sequence cannot contain ambiguous bases') from error
    counts = np.bincount(indices, minlength=64)
    codon_terms = counts * coefficients
    bonus = np.zeros(b-a, dtype=np.float64)
    bonus[positions[2::3]-a] = coefficients[indices]
    local = values[a:b]
    delta = local[np.arange(b-a), path].astype(np.float64)-local[:, 0].astype(np.float64)
    role_terms = np.bincount(path, weights=delta, minlength=15)
    # A codon is added at its third coding base, in production forward order.
    total = float(np.cumsum(delta + bonus)[-1])
    close(math.fsum(role_terms)+math.fsum(codon_terms), total)
    return dict(gain=total, gain_by_role=role_terms.tolist(),
        rows_by_role=np.bincount(path, minlength=15).tolist(),
        codon_counts=counts.tolist(), codon_gain_by_type=codon_terms.tolist(),
        codon_gain=math.fsum(codon_terms))


def margin_terms(reference, conflicts):
    role = [math.fsum(r['gain_by_role'][i] for r in conflicts) for i in range(15)]
    rows = [sum(r['rows_by_role'][i] for r in conflicts) for i in range(15)]
    codon = [math.fsum(r['codon_gain_by_type'][i] for r in conflicts) for i in range(64)]
    counts = [sum(r['codon_counts'][i] for r in conflicts) for i in range(64)]
    total = math.fsum(r['gain'] for r in conflicts)
    role_margin = [a-b for a, b in zip(reference['gain_by_role'], role)]
    codon_margin = [a-b for a, b in zip(reference['codon_gain_by_type'], codon)]
    margin = reference['gain']-total
    close(math.fsum(reference['gain_by_role'])+reference['codon_gain'], reference['gain'])
    close(math.fsum(role)+math.fsum(codon), total)
    close(math.fsum(role_margin)+math.fsum(codon_margin), margin)
    return dict(conflicting_selected_gain=total, conflicting_selected_gain_by_role=role,
        conflicting_selected_rows_by_role=rows, conflicting_selected_codon_counts=counts,
        conflicting_selected_codon_gain_by_type=codon,
        replacement_gain_difference=margin, replacement_margin_by_role=role_margin,
        replacement_margin_by_codon=codon_margin,
        role_only_fixed_set_margin=math.fsum(role_margin),
        codon_fixed_set_margin=math.fsum(codon_margin))


def check_row(row, coefficients):
    if not row['single_chain_supported']:
        if row['reference_gain_vs_background'] is not None or row['exact_primary_hit']:
            raise ValueError('Unsupported references must remain unscored and not exact hits')
        return
    if (len(row['gain_by_role']) != 15 or row['gain_by_role'][0] != 0
            or row['rows_by_role'][0] != 0 or len(row['codon_counts']) != 64
            or sum(row['codon_counts'])*3 != sum(row['rows_by_role'][i] for i in (1,2,3,7,8,9,10,11,12,13,14))):
        raise ValueError('Complete codons and coding role rows must retain the same CDS geometry')
    for counts, terms in ((row['codon_counts'], row['codon_gain_by_type']),
                          (row['conflicting_selected_codon_counts'], row['conflicting_selected_codon_gain_by_type'])):
        if len(counts) != 64 or len(terms) != 64 or any(n < 0 or int(n) != n for n in counts):
            raise ValueError('All codon counts must be nonnegative integer exposures')
        for n, score, coefficient in zip(counts, terms, coefficients):
            close(n*coefficient, score)
    close(math.fsum(row['gain_by_role'])+math.fsum(row['codon_gain_by_type']), row['reference_gain_vs_background'])
    close(math.fsum(row['conflicting_selected_gain_by_role'])+math.fsum(row['conflicting_selected_codon_gain_by_type']), row['conflicting_selected_gain'])
    for a, b, difference in zip(row['gain_by_role'], row['conflicting_selected_gain_by_role'], row['replacement_margin_by_role']):
        close(a-b, difference)
    for a, b, difference in zip(row['codon_gain_by_type'], row['conflicting_selected_codon_gain_by_type'], row['replacement_margin_by_codon']):
        close(a-b, difference)
    close(row['reference_gain_vs_background']-row['conflicting_selected_gain'], row['replacement_gain_difference'])
    close(math.fsum(row['replacement_margin_by_role']), row['role_only_fixed_set_margin'])
    close(math.fsum(row['replacement_margin_by_codon']), row['codon_fixed_set_margin'])
    close(row['role_only_fixed_set_margin']+row['codon_fixed_set_margin'], row['replacement_gain_difference'])


def key(row):
    return row['species'], row['seqid'], row['strand'], tuple(tuple(p) for p in row['reference_CDS'])


def paired_summary(control, candidate):
    if set(control) != set(candidate):
        raise ValueError('Keep identical reference geometry including unsupported chains')
    groups, transitions = {}, Counter()
    for identity, a in control.items():
        b = candidate[identity]
        if any(a[field] != b[field] for field in ('single_chain_supported', 'rows_by_role', 'codon_counts')):
            raise ValueError('Fixed DNA, grammar, codons and reference role rows must match')
        transitions[a['category']+' -> '+b['category']] += 1
        group = ('both_hit' if b['exact_primary_hit'] else 'control_only') if a['exact_primary_hit'] else ('candidate_only' if b['exact_primary_hit'] else 'both_missed')
        groups.setdefault(group, []).append((a, b))
    result = {}
    for name in ('both_hit', 'control_only', 'candidate_only', 'both_missed'):
        pairs = groups.get(name, [])
        legal = [p for p in pairs if p[0]['single_chain_supported']]
        arms = {}
        for i, arm in enumerate(ARMS):
            rows = [p[i] for p in legal]
            means = {field:math.fsum(r[field] for r in rows)/len(rows) if rows else None
                     for field in ('reference_gain_vs_background', 'conflicting_selected_gain',
                                   'replacement_gain_difference', 'role_only_fixed_set_margin', 'codon_fixed_set_margin')}
            arms[arm] = dict(mean_scores=means,
                role_groups={g:math.fsum(math.fsum(r['replacement_margin_by_role'][k] for k in indices) for r in rows)/len(rows) if rows else None
                             for g, indices in T.GROUPS.items()},
                reference_codon_exposures=[sum(r['codon_counts'][k] for r in rows) for k in range(64)],
                conflicting_codon_exposures=[sum(r['conflicting_selected_codon_counts'][k] for r in rows) for k in range(64)],
                full_margin_nonpositive_but_role_only_margin_positive=sum(
                    r['replacement_gain_difference'] <= D.EPSILON and r['role_only_fixed_set_margin'] > D.EPSILON for r in rows))
        result[name] = dict(reference=len(pairs), supported_reference=len(legal), per_arm=arms)
    return dict(reference=len(control), hit_groups=result, category_transitions=dict(transitions),
        gained=result['candidate_only']['reference'], lost=result['control_only']['reference'],
        net_hits=result['candidate_only']['reference']-result['control_only']['reference'])


def completed_sources():
    checked = json.loads((CHECKS/'summary.json').read_text())
    if ((CHECKS/'STATUS').read_text().strip() != 'COMPLETED' or not checked['completed']
            or checked['protocol'] != PROTOCOL or checked['research_data_read']):
        raise ValueError('Completed artificial scoring checks required before research access')
    paired = json.loads((PAIRED/'summary.json').read_text())
    if ((PAIRED/'STATUS').read_text().strip() != 'COMPLETED' or not paired['completed']
            or paired['job_id'] != PAIRED_JOB or paired['DEV_used'] or paired['sealed_data_used']
            or paired['optimizer_steps_this_run']
            or paired['source_jobs'] != dict(baseline=EVAL_JOBS[0], codon=EVAL_JOBS[1])
            or paired['source_fit_jobs'] != dict(baseline=FIT_JOBS[0], codon=FIT_JOBS[1])
            or not paired['all_paired_caller_inputs_and_normalized_original_DNA_directly_equal']):
        raise ValueError('The fixed final9216 codon/fresh-joint paired TRAIN closure is required')
    profile = json.loads((D.PROFILE/'summary.json').read_text())
    if (not profile['completed'] or (D.PROFILE/'STATUS').read_text().strip() != 'COMPLETED' or profile['DEV_used']):
        raise ValueError('Original retained TRAIN primary-reference profile required')
    evaluations = []
    for source, fit_job, evaluation_job in zip(SOURCES, FIT_JOBS, EVAL_JOBS):
        source.source(fit_job)
        value = json.loads((source.EVAL/'summary.json').read_text())
        if ((source.EVAL/'STATUS').read_text().strip() != 'COMPLETED' or not value['completed']
                or value['job_id'] != evaluation_job or value['checkpoint'] != str(source.CHECKPOINT)
                or value['DEV_used'] or value['sealed_data_used'] or value['optimizer_steps_this_run']
                or value['reference_used_in_inference'] or not value['unsupported_reference_strata_retained']
                or value['pooled']['exact_chain']['reference'] != 22688):
            raise ValueError('Completed original six-TRAIN evaluation and retained denominator required')
        evaluations.append(value)
    if evaluations[0]['reference_policy'] != evaluations[1]['reference_policy']:
        raise ValueError('Keep the matched primary reference policy')
    return profile, evaluations


@torch.inference_mode()
def task(number):
    began = time.perf_counter()
    profile, evaluations = completed_sources()
    species, seqid, length = B.I.SCOPE[number]
    out = OUT/seqid
    sequences = [C.M.read_allowed(source.OUT/seqid/'input.fasta', {seqid:length})[seqid] for source in SOURCES]
    if sequences[0] != sequences[1]:
        raise ValueError('Both saved raw inputs must directly match')
    sequence = sequences[0]
    del sequences
    refs = []
    with (D.PROFILE/'reference_chain_profile.jsonl').open() as handle:
        for line in handle:
            row = json.loads(line)
            if row['seqid'] == seqid:
                refs.append(row)
    rk = {C.S.key(r) for r in refs}
    if len(rk) != len(refs):
        raise ValueError('Retain unique reference geometry')
    locations = [(a, min(a+C.M.WIDTH, length)) for a in C.M.starts(length)]
    D.write_json(out/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], task=number,
        species=species, seqid=seqid, length=length, protocol=PROTOCOL, arms=list(ARMS),
        source_evaluation_jobs=list(EVAL_JOBS), source_paired_job=PAIRED_JOB,
        checkpoints=[str(s.CHECKPOINT) for s in SOURCES], reference_used_by_encoder_or_head=False,
        predictions_modified=False, decoder_rerun=False, optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False))
    backbone, tokenizer = B.I.R.load_native()
    torch.cuda.reset_peak_memory_stats()
    for ordinal, (arm, source) in enumerate(zip(ARMS, SOURCES)):
        folder = out/arm
        folder.mkdir()
        raw_root = source.OUT/seqid
        receipt = json.loads((raw_root/'summary.json').read_text())
        raw = json.loads((raw_root/'prediction/summary.json').read_text())
        if ((raw_root/'STATUS').read_text().strip() != 'COMPLETED'
                or not receipt['inference_complete'] or receipt['source_fit_job'] != FIT_JOBS[ordinal]
                or (receipt['species'], receipt['seqid'], receipt['length']) != (species, seqid, length)
                or receipt['checkpoint'] != str(source.CHECKPOINT) or receipt['DEV_used'] or receipt['sealed_data_used']
                or receipt['optimizer_steps_this_run'] or not raw['inference_complete']
                or raw['checkpoint'] != str(source.CHECKPOINT) or raw['reference_used'] or raw['lengths'] != {seqid:length}
                or len(raw['scopes']) != 2 or {s['strand'] for s in raw['scopes']} != {'+', '-'}
                or any(s['windows'] != [list(x) for x in locations] or s['context_policy'] != J.R.CONTEXT for s in raw['scopes'])):
            raise ValueError('Completed frozen raw identity and production grid required')
        state = torch.load(source.CHECKPOINT, map_location='cpu', weights_only=True)
        configuration = source.A.validate_checkpoint(state)
        if configuration != raw['configuration']:
            raise ValueError('Frozen output and replayed checkpoint identity differ')
        if ordinal == 0:
            model, coefficients = source.A.model_from_state(state, 'cuda'), np.zeros(64)
        else:
            model, raw_scores = source.A.model_from_state(state, 'cuda')
            K.validate_raw_report(raw, seqid, length, raw_scores)
            coefficients = centered_codon_scores(raw_scores).numpy()
        del state
        predpath = raw_root/'prediction/full_grammar.gff3'
        predictions = C.S.E.primary_transcripts(C.S.E.parse_annotation(predpath, {seqid:length}), complete_only=False)
        pk = {C.S.key(t) for t in predictions}
        exact = C.S.prf(len(pk & rk), len(pk), len(rk))
        if len(pk) != len(predictions) or exact != evaluations[ordinal]['chromosomes'][seqid]['exact_chain']:
            raise ValueError('Prediction and retained reference counts must reproduce evaluation')
        saved_scores = D.gff_scores(predpath)
        rows, scopes = [], []
        with (folder/'reference_energy.jsonl').open('x') as handle:
            for strand in ('+', '-'):
                start = time.perf_counter()
                joined = J.contig_emissions(model, backbone, tokenizer, sequence, locations, strand, 'cuda')
                encoding_seconds = time.perf_counter()-start
                values, dna = joined['local_potentials'].numpy(), joined['sequence']
                prior = next(s for s in raw['scopes'] if s['strand'] == strand)
                if [list(p) for p in joined['genomic_context_segments']] != prior['genomic_context_segments']:
                    raise ValueError('Per-base ownership must equal the production raw call')
                selected, maximum_error = [], 0.
                for tx in predictions:
                    if tx['strand'] != strand:
                        continue
                    parts = D.oriented_parts(tx, length)
                    if not D.supported(dna, parts, D.role_path(parts)):
                        raise ValueError('Saved prediction lies outside the unchanged grammar')
                    terms = chain_terms(dna, values, parts, coefficients)
                    error = abs(terms['gain']-saved_scores[tx['id']])
                    maximum_error = max(maximum_error, error)
                    if error > 5.1e-8*max(1., abs(saved_scores[tx['id']]))+D.EPSILON:
                        raise ValueError('Role-plus-codon replay differs from the saved eight-digit GFF score')
                    selected.append(dict(id=tx['id'], parts=parts, **terms))
                selected.sort(key=lambda r:r['parts'][0][0])
                if any(a['parts'][-1][1] > b['parts'][0][0] for a,b in zip(selected, selected[1:])):
                    raise ValueError('Saved same-strand selected spans must not overlap')
                index = selected, [r['parts'][-1][1] for r in selected]
                total = math.fsum(r['gain'] for r in selected)
                if not math.isclose(total, prior['gain'], rel_tol=1e-8, abs_tol=D.EPSILON):
                    raise ValueError('Selected-set score must reproduce the frozen contig receipt')
                for ref in refs:
                    if ref['strand'] != strand:
                        continue
                    parts = D.oriented_parts(ref, length)
                    legal = (parts[0][1]-parts[0][0] >= 3 and parts[-1][1]-parts[-1][0] >= 3
                             and D.supported(dna, parts, D.role_path(parts)))
                    hit = C.S.key(ref) in pk
                    if hit and not legal:
                        raise ValueError('Exact selected reference must be grammar-supported')
                    row = dict(species=species, seqid=seqid, id=ref['id'], strand=strand,
                        reference_CDS=[[a,b] for a,b,_ in ref['CDS']], CDS_segments=len(parts),
                        single_chain_supported=legal, exact_primary_hit=hit,
                        original_draw_R4_positive_chain=ref['original_draw_R4_positive_chain'],
                        reference_gain_vs_background=None, rows_by_role=None, codon_counts=None)
                    if legal:
                        terms = chain_terms(dna, values, parts, coefficients)
                        conflicts = D.conflicting(parts, index)
                        row.update({k:v for k,v in terms.items() if k != 'gain'})
                        row.update(reference_gain_vs_background=terms['gain'],
                                   conflicting_selected_ids=[r['id'] for r in conflicts],
                                   **margin_terms(terms, conflicts))
                    row['category'] = D.classify(hit, legal, row['reference_gain_vs_background'], row.get('conflicting_selected_gain'))
                    check_row(row, coefficients)
                    if row['category'] == 'positive_reference_improves_saved_selected_set':
                        raise ValueError('Legal improving replacement: inspect replay before scientific interpretation')
                    rows.append(row)
                    handle.write(json.dumps(row, allow_nan=False)+'\n')
                scopes.append(dict(strand=strand, selected_chains=len(selected),
                    replayed_selected_gain=total, original_selected_gain=prior['gain'],
                    maximum_chain_score_difference_from_8_significant_digit_GFF=maximum_error,
                    encoding_seconds=encoding_seconds, seconds=time.perf_counter()-start))
                handle.flush()
                del joined, values, dna
                print(json.dumps(dict(arm=arm, seqid=seqid, completed_strand=strand)), flush=True)
        if any(p.requires_grad or p.grad is not None for module in (model, backbone) for p in module.parameters()):
            raise ValueError('Scoring replay must not create parameter gradients')
        D.write_json(folder/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'],
            arm=arm, task=number, species=species, seqid=seqid, checkpoint=str(source.CHECKPOINT),
            source_fit_job=FIT_JOBS[ordinal], source_raw_job=receipt['job_id'], source_evaluation_job=EVAL_JOBS[ordinal],
            exact_chain_unchanged=exact, counts=D.summarize(rows), scopes=scopes,
            centered_codon_scores=coefficients.tolist(), codon_order=list(CODONS),
            predictions_modified=False, decoder_rerun=False, optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False))
        del model, predictions, rows, selected
    D.write_json(out/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], task=number,
        species=species, seqid=seqid, length=length, protocol=PROTOCOL, source_profile_job=profile['job_id'],
        source_evaluation_jobs=list(EVAL_JOBS), source_paired_job=PAIRED_JOB,
        predictions_modified=False, decoder_rerun=False, potentials_saved=False, optimizer_steps_this_run=0,
        DEV_used=False, sealed_data_used=False, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        seconds=time.perf_counter()-began))


def aggregate():
    _, evaluations = completed_sources()
    rows = {arm:[] for arm in ARMS}
    jobs, coefficients = {}, {}
    for number, (species, seqid, length) in enumerate(B.I.SCOPE):
        folder = OUT/seqid
        outer = json.loads((folder/'summary.json').read_text())
        if ((folder/'STATUS').read_text().strip() != 'COMPLETED' or not outer['completed']
                or (outer['task'], outer['species'], outer['seqid'], outer['length']) != (number, species, seqid, length)
                or outer['protocol'] != PROTOCOL or outer['DEV_used'] or outer['sealed_data_used']
                or outer['optimizer_steps_this_run'] or outer['predictions_modified'] or outer['decoder_rerun']):
            raise ValueError('All six unchanged original-TRAIN replay tasks required')
        jobs[seqid] = outer['job_id']
        for i, arm in enumerate(ARMS):
            report = json.loads((folder/arm/'summary.json').read_text())
            block = [json.loads(line) for line in (folder/arm/'reference_energy.jsonl').open()]
            if (not report['completed'] or report['source_evaluation_job'] != EVAL_JOBS[i]
                    or report['checkpoint'] != str(SOURCES[i].CHECKPOINT) or report['codon_order'] != list(CODONS)
                    or report['counts'] != D.summarize(block)
                    or report['exact_chain_unchanged'] != evaluations[i]['chromosomes'][seqid]['exact_chain']
                    or len(block) != len({key(r) for r in block})):
                raise ValueError('Saved replay rows must reconstruct frozen sources')
            if arm in coefficients and coefficients[arm] != report['centered_codon_scores']:
                raise ValueError('All chromosomes must use identical final codon coefficients')
            coefficients[arm] = report['centered_codon_scores']
            for row in block:
                if (row['species'], row['seqid']) != (species, seqid):
                    raise ValueError('Reference rows cannot move between chromosomes')
                check_row(row, coefficients[arm])
            rows[arm].extend(block)
    scopes = {}
    for scope in (*[x[1] for x in B.I.SCOPE], 'arabidopsis_thaliana', 'oryza_sativa', 'pooled'):
        selected = {arm:[r for r in rows[arm] if scope == 'pooled' or r['seqid'] == scope or r['species'] == scope] for arm in ARMS}
        paired = paired_summary(*({key(r):r for r in selected[arm]} for arm in ARMS))
        counts = {arm:D.summarize(selected[arm]) for arm in ARMS}
        for i, arm in enumerate(ARMS):
            source = evaluations[i]
            expected = source['pooled'] if scope == 'pooled' else source['per_species'][scope] if scope in source['per_species'] else source['chromosomes'][scope]
            if counts[arm]['reference'] != expected['exact_chain']['reference'] or counts[arm]['exact_primary_hits'] != expected['exact_chain']['tp']:
                raise ValueError('TP and reference denominator must reproduce frozen evaluations')
        if paired['net_hits'] != counts[ARMS[1]]['exact_primary_hits']-counts[ARMS[0]]['exact_primary_hits']:
            raise ValueError('Paired hit accounting must reconstruct the whole-evaluation TP change')
        scopes[scope] = dict(counts=counts, paired=paired)
    if scopes['pooled']['paired']['reference'] != 22688:
        raise ValueError('Do not drop unsupported primary references')
    D.write_json(OUT/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL,
        arms=list(ARMS), source_jobs=jobs, source_evaluation_jobs=list(EVAL_JOBS), source_paired_job=PAIRED_JOB,
        reference_policy=evaluations[0]['reference_policy'], codon_order=list(CODONS), centered_codon_scores=coefficients,
        scopes=scopes, role_groups={k:list(v) for k,v in T.GROUPS.items()},
        reporting_tolerance=D.EPSILON, predictions_modified=False, decoder_rerun=False,
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
        limitation='fixed-set descriptive TRAIN score arithmetic; zero-codon scores are not counterfactual predictions; no causal, calibrated-probability or generalization claim'))


def checks():
    began = time.perf_counter()
    suite = unittest.defaultTestLoader.discover(str(ROOT/'tests'), pattern='test_m28_full_grammar_codon_energy.py')
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.testsRun != 7 or not result.wasSuccessful():
        raise ValueError('All seven new artificial scoring checks must pass')
    D.write_json(CHECKS/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], protocol=PROTOCOL,
        tests_run=result.testsRun, errors=len(result.errors), failures=len(result.failures),
        research_data_read=False, checkpoint_read=False, native_forward=False,
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False, seconds=time.perf_counter()-began))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('checks', 'task', 'aggregate'))
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
