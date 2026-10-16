#!/usr/bin/env python3
"""Fixed-final-weight TRAIN score diagnosis; no fitting or prediction changes."""
import argparse
from bisect import bisect_right
from collections import Counter, defaultdict
import json
import math
import os
import time

import numpy as np
from numba import njit
import torch
import common_ruler as C
import infer_full_grammar_original_train as I
from src.m28.full_grammar_boundary_compiled import _state_table
from src.m28.phase_evidence_chain import candidate_roles

OUT = I.ROOT/'outputs/M28-FULL-GRAMMAR-TRAIN-ENERGY-DIAGNOSTIC'
CHECKS = I.ROOT/'outputs/M28-FULL-GRAMMAR-TRAIN-ENERGY-CHECKS-R1'
PROFILE = I.ROOT/'outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-CHAIN-PROFILE'
DIAGNOSTIC = I.ROOT/'outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-DIAGNOSTIC'
EPSILON = 1e-5  # Numerical near-zero reporting only; not a model operating point.


@njit(cache=False, fastmath=False, parallel=False)
def follow_roles(codes, path, destinations, roles, background):
    """Existence of the supplied fixed role path in the actual production graph."""
    size = destinations.shape[1]
    active = np.empty(size, dtype=np.int32)
    following = np.empty(size, dtype=np.int32)
    seen = np.full(size, -1, dtype=np.int64)
    active[0] = background
    count = 1
    for position in range(len(codes)):
        kept = 0
        for j in range(count):
            state = active[j]
            for edge in range(3):
                target = destinations[codes[position], state, edge]
                if target < 0:
                    break
                if roles[codes[position], state, edge] == path[position] and seen[target] != position:
                    seen[target] = position
                    following[kept] = target
                    kept += 1
        if not kept:
            return False
        active, following = following, active
        count = kept
    for j in range(count):
        if active[j] == background:
            return True
    return False


def oriented_parts(tx, length):
    parts = [(a, b) for a, b, _ in tx['CDS']]
    return parts if tx['strand'] == '+' else [(length-b, length-a) for a, b in parts[::-1]]


def role_path(parts):
    a, b = parts[0][0], parts[-1][1]
    shifted = [(x-a, y-a) for x, y in parts]
    return candidate_roles(shifted, b-a).numpy()


def supported(sequence, parts, path, cap=128):
    _, destinations, roles, background = _state_table(cap)
    local = sequence[parts[0][0]:parts[-1][1]]
    lookup = np.full(256, 4, dtype=np.int8)
    for i, base in enumerate('ACGTN'):
        lookup[ord(base)] = i
    codes = lookup[np.frombuffer(local.encode('ascii'), dtype=np.uint8)]
    return bool(follow_roles(codes, path, destinations, roles, background))


def chain_score(values, parts, path):
    a, b = parts[0][0], parts[-1][1]
    selected = values[a:b]
    delta = selected[np.arange(b-a), path].astype(np.float64)-selected[:, 0].astype(np.float64)
    # Forward accumulation matches production, without FP32 subtraction first.
    gain = float(np.cumsum(delta)[-1])
    by_role = np.bincount(path, weights=delta, minlength=15)
    rows = np.bincount(path, minlength=15)
    return gain, by_role.tolist(), rows.tolist()


def conflicting(parts, predictions):
    rows, ends = predictions
    a, b = parts[0][0], parts[-1][1]
    result = []
    for j in range(bisect_right(ends, a), len(rows)):
        item = rows[j]
        if item['parts'][0][0] >= b:
            break
        result.append(item)
    return result


def classify(hit, legal, gain=None, competing=None):
    if hit:
        return 'exact_primary_chain'
    if not legal:
        return 'outside_fixed_single_chain_grammar'
    if gain < -EPSILON:
        return 'supported_reference_negative_gain_vs_background'
    if abs(gain) <= EPSILON:
        return 'supported_reference_near_zero_gain_vs_background'
    if competing-gain > EPSILON:
        return 'supported_positive_reference_outscored_by_selected_conflicts'
    if abs(competing-gain) <= EPSILON:
        return 'supported_positive_reference_near_replacement_tie'
    return 'positive_reference_improves_saved_selected_set'


def summarize(rows):
    bins = Counter(r['category'] for r in rows)
    strata = {}
    for name, predicate in (
            ('single_CDS_segment', lambda r: r['CDS_segments'] == 1),
            ('multiple_CDS_segments', lambda r: r['CDS_segments'] > 1),
            ('R4_positive_metadata_present', lambda r: r['original_draw_R4_positive_chain']),
            ('R4_positive_metadata_absent', lambda r: not r['original_draw_R4_positive_chain'])):
        group = [r for r in rows if predicate(r)]
        strata[name] = dict(reference=len(group), categories=dict(Counter(r['category'] for r in group)))
    return dict(reference=len(rows), exact_primary_hits=bins['exact_primary_chain'],
        categories=dict(bins), strata=strata,
        individually_supported_reference=sum(r['single_chain_supported'] for r in rows),
        replacement_improvement_counterexamples=bins['positive_reference_improves_saved_selected_set'])


def gff_scores(path):
    result = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            if fields[2] == 'mRNA':
                attrs = dict(part.split('=', 1) for part in fields[8].split(';'))
                result[attrs['ID']] = float(fields[5])
    return result


def write_json(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def checks():
    from src.m28.full_grammar_skyline import decode_full_grammar_skyline
    began = time.perf_counter()
    sequence = 'ATGCGTTAGAGTTAGAGCTAGTAA'
    values = (torch.arange(len(sequence)*15).reshape(-1, 15).float() % 11)-5
    actual = decode_full_grammar_skyline(sequence, values, 4, 3)
    assert actual['chains']
    for parts, expected in zip(actual['chains'], actual['chain_gains']):
        path = role_path(parts)
        assert supported(sequence, parts, path, 4)
        gain, terms, counts = chain_score(values.numpy(), parts, path)
        assert gain == expected == sum(terms) and sum(counts) == parts[-1][1]-parts[0][0]
    # Five-base intron, actual carried codon and short-exon precedence.
    sequence = 'ATGCGTNAGTTTAA'
    parts = [(0, 4), (9, 14)]
    assert supported(sequence, parts, role_path(parts), 2)
    assert not supported(sequence, parts, role_path(parts), 1)
    assert not supported('ATGCCCNAGTTTAA', parts, role_path(parts), 2)
    assert not supported('ATGTAATAA', [(0, 9)], role_path([(0, 9)]), 1)
    tx = dict(CDS=[(2, 5, '0'), (10, 16, '0')], strand='-')
    assert oriented_parts(tx, 20) == [(4, 10), (15, 18)]
    entries = [dict(parts=[(0, 9)], gain=3.), dict(parts=[(9, 18)], gain=7.)]
    index = (entries, [9, 18])
    assert conflicting([(0, 9)], index) == entries[:1]
    assert conflicting([(8, 12)], index) == entries
    assert classify(False, False) == 'outside_fixed_single_chain_grammar'
    assert classify(False, True, -1., 0.) == 'supported_reference_negative_gain_vs_background'
    assert classify(False, True, 2., 3.) == 'supported_positive_reference_outscored_by_selected_conflicts'
    assert classify(False, True, 4., 3.) == 'positive_reference_improves_saved_selected_set'
    report = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'],
        checks='actual_skyline_score_role_support_carry_cap_motif_stop_minus_overlap_and_classification',
        implementation_repair='explicit_loop_instead_of_unsupported_numba_generator',
        failed_check_job_retained='13450149',
        actual_research_DNA_read=False, checkpoint_read=False, native_forward=False,
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
        seconds=time.perf_counter()-began)
    write_json(CHECKS/'summary.json', report)
    print(json.dumps(report), flush=True)


@torch.inference_mode()
def task(number):
    species, seqid, length = I.SCOPE[number]
    out = OUT/seqid
    began = time.perf_counter()
    checked = json.loads((CHECKS/'summary.json').read_text())
    original = json.loads((DIAGNOSTIC/'summary.json').read_text())
    profile = json.loads((PROFILE/'summary.json').read_text())
    source = I.OUT/seqid
    raw = json.loads((source/'prediction/summary.json').read_text())
    if (not checked['completed'] or (CHECKS/'STATUS').read_text().strip() != 'COMPLETED'
            or (PROFILE/'STATUS').read_text().strip() != 'COMPLETED'
            or (DIAGNOSTIC/'STATUS').read_text().strip() != 'COMPLETED'
            or (source/'STATUS').read_text().strip() != 'COMPLETED'
            or not profile['completed'] or not original['completed']
            or profile['DEV_used'] or original['DEV_used'] or raw['reference_used']
            or original['checkpoint'] != str(I.CHECKPOINT) or raw['checkpoint'] != str(I.CHECKPOINT)):
        raise ValueError('Completed final-weight original-TRAIN sources required')
    expected_windows = [(a, min(a+C.M.WIDTH, length)) for a in C.M.starts(length)]
    if (raw['lengths'] != {seqid:length} or len(raw['scopes']) != 2
            or any(s['windows'] != [list(x) for x in expected_windows]
                   or s['context_policy'] != I.R.CONTEXT for s in raw['scopes'])):
        raise ValueError('Whole original TRAIN window ownership differs')
    refs = []
    with (PROFILE/'reference_chain_profile.jsonl').open() as handle:
        for line in handle:
            row = json.loads(line)
            if row['seqid'] == seqid:
                refs.append(row)
    predpath = source/'prediction/full_grammar.gff3'
    predicted = C.S.E.primary_transcripts(C.S.E.parse_annotation(predpath, {seqid:length}), complete_only=False)
    pk, rk = ({C.S.key(t) for t in rows} for rows in (predicted, refs))
    exact = C.S.prf(len(pk & rk), len(pk), len(rk))
    if (len(refs) != len(rk) or len(predicted) != len(pk)
            or exact != original['chromosomes'][seqid]['exact_chain']):
        raise ValueError('Do not change frozen primary denominator or predictions')
    write_json(out/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], task=number,
        species=species, seqid=seqid, length=length, checkpoint=str(I.CHECKPOINT),
        DEV_used=False, sealed_data_used=False, optimizer_steps_this_run=0))
    sequence = C.M.read_allowed(source/'input.fasta', {seqid:length})[seqid]
    state = torch.load(I.CHECKPOINT, map_location='cpu', weights_only=True)
    configuration = I.R.validate_checkpoint(state)
    model = I.R.M28Core('C0', 2048, 128)
    model.load_state_dict(state['model'], strict=True)
    del state
    model.eval().requires_grad_(False).to('cuda')
    backbone, tokenizer = I.R.load_native()
    stored_scores = gff_scores(predpath)
    torch.cuda.reset_peak_memory_stats()
    rows, scopes = [], []
    with (out/'reference_energy.jsonl').open('x') as handle:
        for strand in ('+', '-'):
            t = time.perf_counter()
            joined = I.R.contig_emissions(model, backbone, tokenizer, sequence, expected_windows, strand, 'cuda')
            encoding_seconds = time.perf_counter()-t
            values = joined['local_potentials'].numpy()
            dna = joined['sequence']
            selected = []
            max_error = 0.
            for tx in predicted:
                if tx['strand'] != strand:
                    continue
                parts = oriented_parts(tx, length)
                path = role_path(parts)
                if not supported(dna, parts, path):
                    raise ValueError('Saved selected chain not supported by actual production graph')
                gain, _, _ = chain_score(values, parts, path)
                receipt = stored_scores[tx['id']]
                error = abs(gain-receipt)
                max_error = max(max_error, error)
                # GFF export has eight significant digits; allow that rounding
                # plus 1e-5 accumulation noise, not a fitted score correction.
                if error > 5.1e-8*max(1., abs(receipt))+EPSILON:
                    raise ValueError('Native fixed-weight replay disagrees with saved chain score: '+tx['id'])
                selected.append(dict(id=tx['id'], parts=parts, gain=gain))
            selected.sort(key=lambda r: r['parts'][0][0])
            if any(a['parts'][-1][1] > b['parts'][0][0] for a, b in zip(selected, selected[1:])):
                raise ValueError('Selected same-strand span overlap prevents replacement comparison')
            index = (selected, [r['parts'][-1][1] for r in selected])
            total = math.fsum(r['gain'] for r in selected)
            prior = next(s for s in raw['scopes'] if s['strand'] == strand)
            if not math.isclose(total, prior['gain'], rel_tol=1e-8, abs_tol=EPSILON):
                raise ValueError('Replayed selected-set gain differs from frozen contig receipt')
            for ref in refs:
                if ref['strand'] != strand:
                    continue
                parts = oriented_parts(ref, length)
                failures, _ = C.M.chain_checks(ref, sequence)
                if len(parts) > 128:
                    failures.append('over_128_exons')
                if any(c-b < 5 for (_, b), (c, _) in zip(parts, parts[1:])):
                    failures.append('intron_shorter_than_5')
                if parts[0][1]-parts[0][0] < 3 or parts[-1][1]-parts[-1][0] < 3:
                    path, legal = None, False
                else:
                    path = role_path(parts)
                    legal = supported(dna, parts, path)
                gain = terms = counts = competitors = margin = None
                conflicts = []
                if legal:
                    gain, terms, counts = chain_score(values, parts, path)
                    conflicts = conflicting(parts, index)
                    competitors = math.fsum(r['gain'] for r in conflicts)
                    margin = gain-competitors
                hit = C.S.key(ref) in pk
                if hit and not legal:
                    raise ValueError('Exact selected reference must be individually grammar-supported')
                category = classify(hit, legal, gain, competitors)
                row = dict(species=species, seqid=seqid, id=ref['id'], strand=strand,
                    CDS_segments=len(parts), category=category, single_chain_supported=legal,
                    reported_sequence_structure_failures=failures,
                    original_draw_R4_positive_chain=ref['original_draw_R4_positive_chain'],
                    exact_primary_hit=hit, reference_gain_vs_background=gain,
                    conflicting_selected_ids=[r['id'] for r in conflicts],
                    conflicting_selected_gain=competitors, replacement_gain_difference=margin,
                    gain_by_role=terms, rows_by_role=counts)
                rows.append(row)
                handle.write(json.dumps(row, allow_nan=False)+'\n')
            handle.flush()
            scopes.append(dict(strand=strand, windows=len(expected_windows),
                context_policy=I.R.CONTEXT, selected_chains=len(selected),
                replayed_selected_gain=total, original_selected_gain=prior['gain'],
                maximum_chain_score_difference_from_8_significant_digit_GFF=max_error,
                encoding_seconds=encoding_seconds, seconds=time.perf_counter()-t))
            print(json.dumps(dict(completed_strand=strand, seqid=seqid, selected=len(selected))), flush=True)
            del joined, values, dna
    summary = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], task=number,
        species=species, seqid=seqid, length=length, checkpoint=str(I.CHECKPOINT), configuration=configuration,
        source_diagnostic_job=original['job_id'], source_profile_job=profile['job_id'],
        exact_chain_unchanged=exact, reference_policy=original['reference_policy'],
        DEV_used=False, sealed_data_used=False, optimizer_steps_this_run=0,
        reference_used_for_posthoc_scoring=True, reference_used_by_encoder_or_head=False,
        labels_or_feature_cache_used=False, predictions_modified=False, native_forward=True,
        decoder_rerun=False, potentials_saved=False, score_kind=raw['score_kind'],
        near_zero_reporting_tolerance=EPSILON,
        replacement_definition='remove_all_selected_same_strand_spans_overlapping_reference_then_insert_one_legal_reference',
        individual_support_does_not_prove_joint_reference_set_support=True,
        diagnosis_not_encoder_head_objective_or_data_causal_attribution=True,
        counts=summarize(rows), scopes=scopes, seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved())
    if summary['counts']['replacement_improvement_counterexamples']:
        summary['completed'] = False
    write_json(out/'summary.json', summary)
    print(json.dumps(summary, allow_nan=False), flush=True)
    if not summary['completed']:
        raise ValueError('A positive legal replacement improves the saved optimum; inspect before interpreting')


def aggregate():
    original = json.loads((DIAGNOSTIC/'summary.json').read_text())
    all_rows, runs = [], {}
    for number, (species, seqid, _) in enumerate(I.SCOPE):
        folder = OUT/seqid
        run = json.loads((folder/'summary.json').read_text())
        if ((folder/'STATUS').read_text().strip() != 'COMPLETED' or not run['completed']
                or (run['task'], run['species'], run['seqid']) != (number, species, seqid)
                or run['DEV_used'] or run['sealed_data_used'] or run['optimizer_steps_this_run']):
            raise ValueError('All six fixed TRAIN energy tasks must finish')
        rows = [json.loads(line) for line in (folder/'reference_energy.jsonl').read_text().splitlines()]
        if summarize(rows) != run['counts'] or run['exact_chain_unchanged'] != original['chromosomes'][seqid]['exact_chain']:
            raise ValueError('Completed score rows disagree with summary or original ruler')
        runs[seqid] = run
        all_rows.extend(rows)
    result = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'],
        checkpoint=str(I.CHECKPOINT), reference_policy=original['reference_policy'],
        DEV_used=False, sealed_data_used=False, optimizer_steps_this_run=0,
        predictions_modified=False, checkpoint_read=False, native_forward=False,
        chromosomes={key:value['counts'] for key, value in runs.items()},
        source_jobs={key:value['job_id'] for key, value in runs.items()},
        per_species={species:summarize([r for r in all_rows if r['species'] == species])
                     for species in {r['species'] for r in all_rows}},
        pooled=summarize(all_rows), exact_chain_unchanged=original['pooled']['exact_chain'],
        interpretation='fixed_final_weight_original_TRAIN_support_and_energy_diagnosis_not_new_model_or_generalization')
    if (result['pooled']['reference'] != original['pooled']['exact_chain']['reference']
            or result['pooled']['exact_primary_hits'] != original['pooled']['exact_chain']['tp']):
        raise ValueError('Pooled frozen reference denominator or TP changed')
    write_json(OUT/'summary.json', result)
    print(json.dumps(result, allow_nan=False), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('check', 'task', 'aggregate'))
    parser.add_argument('--task', type=int, choices=range(6))
    args = parser.parse_args()
    torch.set_num_threads(2)
    if args.mode == 'check':
        checks()
    elif args.mode == 'aggregate':
        aggregate()
    elif args.task is None:
        parser.error('task mode needs --task')
    else:
        task(args.task)
