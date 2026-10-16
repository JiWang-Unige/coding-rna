#!/usr/bin/env python3
"""Original1491 TRAIN views, four fixed-context conditions, zero updates."""
import argparse
from collections import Counter
import json
import math
import os
import time
import unittest

import numpy as np
import torch
import diagnose_full_grammar_native_context as N

G, D, ROOT = N.G, N.D, N.ROOT
READY = ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CONTEXT-TRAIN-PREPARATION'
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CONTEXT-TRAIN'
READOUT = ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CONTEXT-TRAIN-READOUT'
AGGREGATE = ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CONTEXT-TRAIN-AGGREGATE'
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_native_context_train.md'
ARMS = N.ARMS
SCOPE = D.I.SCOPE
REFERENCE_POLICY = 'M26_CDS_assessable_longest_eligible_CDS_tie_transcript_ID'
CONFIG = dict(N.CONFIG, scope='all_original1491_unique_TRAIN_views_all1536_draws_retained',
    grammar_view='unchanged_original_neighbor_view_including_four_clipped_endpoints',
    neural_readout='unchanged_original_24576_windows_and_base_ownership',
    JOINT49152='joint_actual_unchanged_view_at_most49152bp',
    JOINT98304='view_plus24576_each_side_clipped_to_real_TRAIN_chromosome_at_most98304bp',
    reference_policy=REFERENCE_POLICY, reference_unit='contained_primary_chain_overlapping_original_core_window_pair',
    training_eligibility_does_not_filter_reference=True, unsupported_reference_strata_retained=True,
    artifact_limit_bytes=20*1024**3)


def canonical(value):
    return json.loads(json.dumps(value))


def rows(path):
    with path.open() as handle:
        return [json.loads(line) for line in handle if line.strip()]


def extended_bounds(a, b, length):
    return max(0, a-24576), min(length, b+24576)


class EndpointChecks(unittest.TestCase):
    def test_clipping_and_oriented_crop_then_ownership(self):
        self.assertEqual(extended_bounds(0, 36864, 100000), (0, 61440))
        self.assertEqual(extended_bounds(63136, 100000, 100000), (38560, 100000))
        self.assertEqual(extended_bounds(30000, 79152, 200000), (5424, 103728))
        # Actual-sized left/right terminal grids, including an irregular terminal anchor.
        for start, length in ((0, 120000), (73728, 100002)):
            geometry = G.context_geometry(start, start+24576, length)
            a, b = geometry['genomic_context']; c, d = extended_bounds(a, b, length)
            for strand in ('+', '-'):
                features = torch.arange((d-c)//6).reshape(1, -1, 1)
                pieces = [N.crop_features(features, l, r, c, d, strand)[0]
                    .repeat_interleave(6, 0).expand(-1, 15) for l, r in geometry['genomic_windows']]
                actual = G.owned_context_logits(pieces, geometry, strand)
                desired = torch.arange((a-c)//6, (b-c)//6) if strand == '+' else torch.arange((d-b)//6, (d-a)//6)
                self.assertTrue(torch.equal(actual, desired.repeat_interleave(6)[:, None].expand(-1, 15)))


def population():
    _, fit, _, original, windows, catalogs = D.source_inputs()
    if fit['job_id'] != '13434603':
        raise ValueError('Original final4608 fit identity required')
    expected = {r['window_id']: r for r in original}
    count = Counter(r['window_id'] for r in original)
    weight = Counter()
    for r in original:
        weight[r['window_id']] += r['importance_weight']
    plan = []
    for ordinal, w in enumerate(windows.values()):
        task = next(i for i, (s, q, _) in enumerate(SCOPE) if (s, q) == (w['species'], w['seqid']))
        length = SCOPE[task][2]
        geometry = G.context_geometry(w['start'], w['end'], length)
        a, b = geometry['genomic_context']; c, d = extended_bounds(a, b, length)
        offset = w['start']-a if w['strand'] == '+' else b-w['end']
        old = expected[w['id']]
        if ([a, b] != old['genomic_context'] or [offset, offset+w['valid_bases']] != old['core_oriented_interval']
                or w['valid_bases'] != 24576 or (d-c) % 6 or (b-a) % 6
                or any(r-l != 24576 or (l-c) % 6 or (d-r) % 6 for l, r in geometry['genomic_windows'])):
            raise ValueError('Original view, real clipped context or local six-mer registration changed: '+w['id'])
        plan.append(dict(ordinal=ordinal, task=task, window_id=w['id'], species=w['species'],
            seqid=w['seqid'], strand=w['strand'], geometry=canonical(geometry),
            extended_genomic_bounds=[c, d], core_oriented_interval=[offset, offset+24576],
            original_draw_multiplicity=count[w['id']], original_importance_mean=weight[w['id']]/count[w['id']],
            known_observation_denominator=old['known_observation_denominator'],
            fine_positions=old['fine_positions'], coarse_positions=old['coarse_positions'],
            fully_unknown_positions=old['fully_unknown_positions']))
    return plan, windows, catalogs


def prepare():
    began = time.perf_counter()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(EndpointChecks))
    if not result.wasSuccessful():
        raise SystemExit(1)
    plan, _, _ = population()
    for root, source_job in ((D.E.PROFILE, '13449803'), (D.E.OUT, '13450172')):
        source = json.loads((root/'summary.json').read_text())
        if ((root/'STATUS').read_text().strip() != 'COMPLETED' or not source['completed']
                or source['job_id'] != source_job or source['reference_policy'] != REFERENCE_POLICY
                or source['DEV_used'] or source['sealed_data_used']):
            raise ValueError('Completed frozen TRAIN primary-reference profile and support strata required')
    # Reconcile the existing baseline population before new encoding, not its metrics.
    for task, (_, seqid, _) in enumerate(SCOPE):
        root = D.OUT/seqid
        summary = json.loads((root/'summary.json').read_text())
        if ((root/'STATUS').read_text().strip() != 'COMPLETED' or not summary['completed']
                or summary['configuration'] != D.F.CONFIG or summary['source_fit_job'] != '13434603'
                or [r['window_id'] for r in rows(root/'windows.jsonl')] != [r['window_id'] for r in plan if r['task'] == task]):
            raise ValueError('Completed unchanged original finite-view baseline required: '+seqid)
    previous = json.loads((N.READOUT/'summary.json').read_text())
    if (not previous['completed'] or previous['job_id'] != '13452596'
            or previous['configuration'] != N.CONFIG or not previous['original_short_block_NLL_and_MAP_reproduced']):
        raise ValueError('Completed preceding eight-witness native readout required')
    with (READY/'geometry_plan.jsonl').open('x') as handle:
        for row in plan:
            handle.write(json.dumps(row, allow_nan=False)+'\n')
    N.write(READY/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        protocol=PROTOCOL, tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        unique_windows=len(plan), original_draws=sum(r['original_draw_multiplicity'] for r in plan),
        expected_conditions=len(plan)*4, contig_windows={q:sum(r['seqid'] == q for r in plan) for _, q, _ in SCOPE},
        view_lengths=dict(Counter(r['geometry']['genomic_context'][1]-r['geometry']['genomic_context'][0] for r in plan)),
        all_local_windows_full24576_and_six_mer_registered=True,
        raw_FASTA_or_feature_tensors_or_checkpoint_read=False, original_TRAIN_metadata_read=True,
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False, seconds=time.perf_counter()-began))


def prepared_population():
    summary = json.loads((READY/'summary.json').read_text())
    if ((READY/'STATUS').read_text().strip() != 'COMPLETED' or not summary['completed']
            or summary['configuration'] != CONFIG or summary['unique_windows'] != 1491):
        raise ValueError('Completed exact full-population endpoint preparation required')
    plan, windows, catalogs = population()
    if rows(READY/'geometry_plan.jsonl') != plan:
        raise ValueError('Original TRAIN order, observations, weights or geometry changed')
    return plan, windows, catalogs


@torch.inference_mode()
def native(task):
    began = time.perf_counter(); species, seqid, length = SCOPE[task]
    plan, windows, _ = prepared_population(); selected = [r for r in plan if r['task'] == task]
    out = OUT/seqid
    chromosome = D.I.R.M.read_allowed(ROOT/'data/m1_screen'/species/'genome.fa', {seqid:length})[seqid]
    saved = torch.load(D.I.CHECKPOINT, map_location='cpu', weights_only=True)
    D.I.R.validate_checkpoint(saved)
    model = G.M28Core('C0', 2048, 128); model.load_state_dict(saved['model'], strict=True); del saved
    model.eval().requires_grad_(False).to('cuda'); backbone, tokenizer = G.load_native()
    previous = json.loads((N.OUT/'summary.json').read_text())
    if previous['job_id'] != '13452550' or previous['configuration'] != N.CONFIG or not previous['completed']:
        raise ValueError('Actual preceding native potentials required for overlapping replay witnesses')
    witnesses = {r['window_id']:r['artifact'] for r in previous['records']}
    N.write(out/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], task=task, seqid=seqid, species=species,
        configuration=CONFIG, expected_windows=len(selected), window_ids=[r['window_id'] for r in selected]))
    records = []; replayed = []; torch.cuda.reset_peak_memory_stats()
    with (out/'windows.jsonl').open('x') as log:
        for local_ordinal, row in enumerate(selected):
            w = windows[row['window_id']]; geometry = G.context_geometry(w['start'], w['end'], length)
            a, b = geometry['genomic_context']; c, d = row['extended_genomic_bounds']; values = {}
            for arm in ARMS:
                whole = token_ids = None
                if arm in ('JOINT49152', 'JOINT98304'):
                    lo, hi = (a, b) if arm == 'JOINT49152' else (c, d)
                    whole, token_ids = N.joint_encode(backbone, tokenizer, N.oriented(chromosome, lo, hi, w['strand']))
                pieces = []
                for left, right in geometry['genomic_windows']:
                    dna_string = N.oriented(chromosome, left, right, w['strand'])
                    if arm == 'BLOCK6144':
                        features = G.encode_window(backbone, tokenizer, dna_string, 'cuda')
                    elif arm == 'JOINT24576':
                        features, _ = N.joint_encode(backbone, tokenizer, dna_string)
                    else:
                        features = N.crop_features(whole, left, right, lo, hi, w['strand'])
                        offset = left-lo if w['strand'] == '+' else hi-right
                        ids = tokenizer(dna_string, add_special_tokens=False, return_tensors='pt', truncation=False)['input_ids']
                        if not torch.equal(ids, token_ids[:, offset//6:(offset+24576)//6]):
                            raise ValueError('Actual native/local token IDs differ: '+w['id'])
                    dna = G.one_hot(dna_string, 24576)[None].to('cuda')
                    pieces.append(model(features, dna, [24576])['segmentation_logits'][0].cpu())
                joined = G.owned_context_logits(pieces, geometry, w['strand'])
                if joined.dtype != torch.float32 or joined.shape != (b-a, 15) or not bool(torch.isfinite(joined).all()):
                    raise ValueError('Finite unchanged neural readout geometry required')
                values[arm] = joined
                del whole, token_ids, pieces, features, dna
            if w['id'] in witnesses:
                old = torch.load(N.OUT/witnesses[w['id']], map_location='cpu', weights_only=True)
                if any(not torch.equal(values[arm], old['arms'][arm]) for arm in ARMS):
                    raise ValueError('Full-population implementation differs from prior native witness potentials')
                replayed.append(w['id'])
            artifact = 'window_%04d.pt' % local_ordinal
            torch.save(dict(configuration=CONFIG, window_id=w['id'], geometry=geometry, arms=values), out/artifact)
            record = dict(ordinal=row['ordinal'], window_id=w['id'], artifact=artifact,
                bytes=(out/artifact).stat().st_size, context_bp=b-a, extended_encoder_bp=d-c,
                six_mer_registration_exact=True)
            records.append(record); log.write(json.dumps(record)+'\n'); log.flush()
            if (local_ordinal+1) % 16 == 0 or local_ordinal+1 == len(selected):
                print(json.dumps(dict(stage='native', task=task, completed_windows=local_ordinal+1,
                    total_windows=len(selected), seconds=time.perf_counter()-began)), flush=True)
            del values
    if any(p.requires_grad or p.grad is not None for module in (model, backbone) for p in module.parameters()):
        raise ValueError('Frozen parameters must have no gradients')
    total_bytes = sum(r['bytes'] for r in records)
    if total_bytes > CONFIG['artifact_limit_bytes']:
        raise ValueError('Declared full-cohort artifact ceiling exceeded by one contig')
    N.write(out/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], task=task, species=species,
        seqid=seqid, configuration=CONFIG, unique_windows=len(records), native_conditions=len(records)*4,
        records=records, prior_eight_witnesses_replayed=replayed, potentials_bytes=total_bytes,
        actual_last4608_checkpoint_read=True, all_native_token_registration_exact=True,
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(), seconds=time.perf_counter()-began))


def summarize(records):
    fine = np.sum([r['fine_confusion'] for r in records], axis=0)
    coarse = np.sum([r['coarse_confusion'] for r in records], axis=0)
    known = sum(r['fine_positions']+r['coarse_positions'] for r in records)
    draws = sum(r['original_draw_multiplicity'] for r in records)
    return dict(windows=len(records), known_core_rows=known, raw_NLL=math.fsum(r['raw_NLL'] for r in records),
        NLL_per_known_core_row=math.fsum(r['raw_NLL'] for r in records)/max(1, known),
        original_draw_weighted_mean_loss=math.fsum(r['weighted_training_loss']*r['original_draw_multiplicity'] for r in records)/draws,
        fine_reference_rows_by_role=fine.sum(1).tolist(), fine_matches_by_role=np.diag(fine).tolist(),
        fine_matches=int(np.trace(fine)), coarse_matches=int(np.trace(coarse)),
        reference_window_pairs=sum(r['reference_window_pairs'] for r in records),
        exact_reference_chain_hits=sum(r['exact_reference_chain_hits'] for r in records),
        unsupported_reference_window_pairs=sum(r['unsupported_reference_window_pairs'] for r in records),
        unsupported_exact_hits=sum(r['unsupported_exact_hits'] for r in records),
        windows_NLL_worse_than_short_block=sum(r['raw_NLL'] > r['short_block_NLL']+1e-6 for r in records),
        windows_NLL_better_than_short_block=sum(r['raw_NLL'] < r['short_block_NLL']-1e-6 for r in records))


@torch.no_grad()
def readout(task):
    began = time.perf_counter(); species, seqid, length = SCOPE[task]
    plan, windows, catalogs = prepared_population(); selected = [r for r in plan if r['task'] == task]
    native_source = json.loads((OUT/seqid/'summary.json').read_text())
    if ((OUT/seqid/'STATUS').read_text().strip() != 'COMPLETED' or not native_source['completed']
            or native_source['configuration'] != CONFIG or native_source['unique_windows'] != len(selected)
            or len(native_source['records']) != len(selected)):
        raise ValueError('Completed matching native chromosome outputs required')
    chromosome = D.I.R.M.read_allowed(ROOT/'data/m1_screen'/species/'genome.fa', {seqid:length})[seqid]
    legacy = {r['window_id']:r for r in rows(D.OUT/seqid/'windows.jsonl')}
    references = [r for r in rows(D.E.PROFILE/'reference_chain_profile.jsonl') if r['seqid'] == seqid]
    flags = {(r['id'], r['strand']):r['single_chain_supported'] for r in rows(D.E.OUT/seqid/'reference_energy.jsonl')}
    if {(r['id'], r['strand']) for r in references} != set(flags):
        raise ValueError('Frozen original primary population and support strata differ')
    out = READOUT/seqid; results = []; maximum_baseline_difference = 0.
    N.write(out/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], task=task, seqid=seqid, species=species,
        configuration=CONFIG, native_source_job=native_source['job_id'], expected_conditions=len(selected)*4))
    with (out/'windows.jsonl').open('x') as log, (out/'reference_views.jsonl').open('x') as ref_log:
        for local_ordinal, (row, emitted) in enumerate(zip(selected, native_source['records'])):
            if row['window_id'] != emitted['window_id']:
                raise ValueError('Native chromosome outputs are not the original ordered population')
            w = windows[row['window_id']]; a, b = row['geometry']['genomic_context']; offset = row['core_oriented_interval'][0]
            raw = chromosome[w['start']:w['end']]; core = raw if w['strand'] == '+' else G.Q.M.rc(raw)
            target = G.P.P.build_targets(w, catalogs[species], core); context = N.oriented(chromosome, a, b, w['strand'])
            relevant = [r for r in references if r['strand'] == w['strand'] and a <= r['CDS'][0][0]
                and r['CDS'][-1][1] <= b and r['CDS'][0][0] < w['end'] and r['CDS'][-1][1] > w['start']]
            saved = torch.load(OUT/seqid/emitted['artifact'], map_location='cpu', weights_only=True)
            if (saved['configuration'] != CONFIG or saved['window_id'] != w['id']
                    or canonical(saved['geometry']) != row['geometry']):
                raise ValueError('Actual potentials must match the declared original geometry')
            old = legacy[w['id']]
            for arm in ARMS:
                logits = saved['arms'][arm]; cuts = tuple(x for x in (offset, offset+24576) if 0 < x < len(context))
                likelihood = D.coarse_channel_nll(context, logits, 128, offset, offset+24576, target, .1,
                    cuts=cuts, checkpoint_bp=256)
                raw_NLL = float(likelihood['loss'])
                if not math.isfinite(raw_NLL) or raw_NLL < -1e-8:
                    raise ValueError('Unchanged observation NLL must remain finite and nonnegative')
                decoded = D.decode_full_grammar_skyline(context, logits, 128)
                roles = D.map_roles(decoded['chains'], len(context))
                agreement = D.core_agreement(roles[offset:offset+24576], target)
                if any(agreement[k] != row[k] for k in ('fine_positions', 'coarse_positions', 'fully_unknown_positions')):
                    raise ValueError('Original observation row denominator changed')
                if arm == 'BLOCK6144':
                    maximum_baseline_difference = max(maximum_baseline_difference, abs(raw_NLL-old['raw_NLL']))
                    if (abs(raw_NLL-old['raw_NLL']) > 1e-6 or agreement['fine_confusion'] != old['fine_confusion']
                            or agreement['coarse_confusion'] != old['coarse_confusion']):
                        raise ValueError('Actual whole-population baseline NLL/MAP failed to replay')
                predictions = {tuple(D.I.R.genomic_chain(chain, a, b, w['strand'])) for chain in decoded['chains']}
                hits = unsupported_hits = unsupported = 0
                for ref in relevant:
                    hit = tuple((x, y) for x, y, _ in ref['CDS']) in predictions
                    supported = flags[(ref['id'], ref['strand'])]
                    hits += hit; unsupported += not supported; unsupported_hits += hit and not supported
                    ref_log.write(json.dumps(dict(window_id=w['id'], arm=arm, id=ref['id'], strand=w['strand'],
                        single_chain_supported=supported, exact_chain_hit=hit))+'\n')
                result = dict(window_id=w['id'], seqid=seqid, species=species, strand=w['strand'], arm=arm,
                    raw_NLL=raw_NLL, short_block_NLL=old['raw_NLL'],
                    weighted_training_loss=raw_NLL*row['original_importance_mean']/row['known_observation_denominator'],
                    original_draw_multiplicity=row['original_draw_multiplicity'], **agreement,
                    reference_window_pairs=len(relevant), exact_reference_chain_hits=hits,
                    unsupported_reference_window_pairs=unsupported, unsupported_exact_hits=unsupported_hits,
                    selected_view_chains=len(decoded['chains']))
                results.append(result); log.write(json.dumps(result, allow_nan=False)+'\n')
            log.flush(); ref_log.flush()
            if (local_ordinal+1) % 16 == 0 or local_ordinal+1 == len(selected):
                print(json.dumps(dict(stage='readout', task=task, completed_windows=local_ordinal+1,
                    total_windows=len(selected), seconds=time.perf_counter()-began)), flush=True)
    N.write(out/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], task=task, species=species,
        seqid=seqid, configuration=CONFIG, unique_windows=len(selected), readout_conditions=len(results),
        native_source_job=native_source['job_id'], maximum_short_block_NLL_abs_difference=maximum_baseline_difference,
        original_short_block_NLL_and_MAP_reproduced=True, unsupported_reference_strata_retained=True,
        aggregate={arm:summarize([r for r in results if r['arm'] == arm]) for arm in ARMS},
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False, seconds=time.perf_counter()-began))


def aggregate():
    began = time.perf_counter(); plan, _, _ = prepared_population(); all_rows = []; native_summaries = []; read_summaries = []
    for task, (_, seqid, _) in enumerate(SCOPE):
        for root, target in ((OUT/seqid, native_summaries), (READOUT/seqid, read_summaries)):
            value = json.loads((root/'summary.json').read_text())
            if ((root/'STATUS').read_text().strip() != 'COMPLETED' or not value['completed']
                    or value['configuration'] != CONFIG or value['task'] != task):
                raise ValueError('Every native/readout task must complete with unchanged identity')
            target.append(value)
        all_rows.extend(rows(READOUT/seqid/'windows.jsonl'))
    if len(all_rows) != 5964:
        raise ValueError('All1491 times four conditions required')
    for arm in ARMS:
        group = [r for r in all_rows if r['arm'] == arm]
        if set(r['window_id'] for r in group) != set(r['window_id'] for r in plan) or len(group) != 1491:
            raise ValueError('Every unchanged original view must appear once per condition')
    for ident in {r['window_id'] for r in plan}:
        paired = [r for r in all_rows if r['window_id'] == ident]
        if len({(r['fine_positions'], r['coarse_positions'], r['fully_unknown_positions'],
                     r['reference_window_pairs'], r['unsupported_reference_window_pairs']) for r in paired}) != 1:
            raise ValueError('Condition-dependent observation or primary-reference denominator')
    if set(i for s in native_summaries for i in s['prior_eight_witnesses_replayed']) != set(
            r['window_id'] for r in json.loads((N.OUT/'summary.json').read_text())['records']):
        raise ValueError('All prior native potentials must be reconciled in the full cohort')
    total_bytes = sum(s['potentials_bytes'] for s in native_summaries)
    if total_bytes > CONFIG['artifact_limit_bytes']:
        raise ValueError('Declared full-cohort artifact ceiling exceeded')
    result = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG, protocol=PROTOCOL,
        unique_windows=1491, original_draws=1536, readout_conditions=5964,
        native_source_jobs=[s['job_id'] for s in native_summaries], readout_source_jobs=[s['job_id'] for s in read_summaries],
        all_original_TRAIN_views_retained=True, original_short_block_NLL_and_MAP_reproduced=True,
        all_eight_previous_native_witness_potentials_exact=True, unsupported_reference_strata_retained=True,
        maximum_short_block_NLL_abs_difference=max(s['maximum_short_block_NLL_abs_difference'] for s in read_summaries),
        potentials_bytes=total_bytes,
        per_species={species:{arm:summarize([r for r in all_rows if r['species'] == species and r['arm'] == arm])
            for arm in ARMS} for species in sorted(G.Q.ALLOW)},
        per_contig={seqid:{arm:summarize([r for r in all_rows if r['seqid'] == seqid and r['arm'] == arm])
            for arm in ARMS} for _, seqid, _ in SCOPE},
        optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
        checkpoint_or_context_selected=False, seconds=time.perf_counter()-began,
        interpretation='complete_original_TRAIN_finite_view_context_comparison_not_full_grid_F1_or_generalization')
    N.write(AGGREGATE/'summary.json', result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('prepare', 'native', 'readout', 'aggregate'))
    parser.add_argument('--task', type=int)
    args = parser.parse_args(); torch.set_num_threads(2)
    if args.mode in ('native', 'readout'):
        if args.task is None or not 0 <= args.task < 6:
            raise ValueError('Original six TRAIN tasks only')
        {'native': native, 'readout': readout}[args.mode](args.task)
    else:
        {'prepare': prepare, 'aggregate': aggregate}[args.mode]()
