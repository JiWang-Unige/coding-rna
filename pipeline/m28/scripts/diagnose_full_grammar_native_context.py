#!/usr/bin/env python3
"""Same eight TRAIN witnesses, fixed last4608 head, native context only; no fit."""
import argparse
import json
import math
import os
import re
import time
import unittest

import numpy as np
import torch
from torch.nn.attention import SDPBackend, sdpa_kernel

import preflight_full_grammar_glm_context as G
import diagnose_full_grammar_final_training_views as D
from src.m28.local_coarse_channel import local_coarse_channel_nll

ROOT = G.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CONTEXT-DIAGNOSTIC'
READOUT = ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CONTEXT-READOUT'
CHECKS = ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CONTEXT-CHECKS'
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_native_context_diagnostic.md'
ARMS = ('BLOCK6144', 'JOINT24576', 'JOINT49152', 'JOINT98304')
CONFIG = dict(scope='same_original_eight_TRAIN_witnesses_no_resampling',
    encoder_model=G.Q.F.MODEL, encoder_revision=G.Q.F.REVISION,
    encoder_frozen=True, encoder_dtype='BF16_hidden_then_FP32_head',
    arms=list(ARMS), joint_attention='original_native_bidirectional_SDPA_flash_only_no_RoPE_change',
    head_checkpoint=str(D.I.CHECKPOINT.relative_to(ROOT)), head_step=4608,
    head='unchanged_frozen_M28Core_C0_2048_128_FP32',
    grammar_view='unchanged_raw_coordinate_49152_neighbor_view',
    neural_readout='same_three_24576_windows_same_nearest_center_base_ownership',
    JOINT24576='one_native_forward_per_existing_head_window',
    JOINT49152='one_native_forward_shared_by_three_existing_head_windows',
    JOINT98304='extend_fixed_49152_view_by_24576_each_side_crop_same_head_windows',
    encoder_bounds_chosen_without_annotation=True, reference_denominator_changed=False,
    observation='original_fine_else_coarse_epsilon_0.1_else_unknown',
    max_exons=128, checkpoint_bp=256, score_accumulation='FP64',
    optimizer_steps=0, DEV_used=False, sealed_data_used=False)


def write(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def long_bounds(a, b, chromosome_length):
    left, right = a-24576, b+24576
    if b-a != 49152 or not 0 <= left < right <= chromosome_length:
        raise ValueError('These fixed interior witnesses require real 98304bp TRAIN context')
    return left, right


def crop_features(features, left, right, a, b, strand):
    offset = left-a if strand == '+' else b-right
    if (strand not in ('+', '-') or not a <= left < right <= b
            or offset % 6 or (right-left) % 6 or features.shape[1] != (b-a)//6):
        raise ValueError('Aligned six-mer crop in the unchanged oriented context required')
    return features[:, offset//6:(offset+right-left)//6]


class CropChecks(unittest.TestCase):
    def test_plus_and_minus_asymmetric_intervals(self):
        values = torch.arange(8).reshape(1, 8, 1)
        self.assertEqual(crop_features(values, 100, 124, 100, 148, '+').flatten().tolist(), [0, 1, 2, 3])
        self.assertEqual(crop_features(values, 100, 124, 100, 148, '-').flatten().tolist(), [4, 5, 6, 7])
        self.assertEqual(crop_features(values, 112, 136, 100, 148, '-').flatten().tolist(), [2, 3, 4, 5])
        with self.assertRaises(ValueError):
            crop_features(values, 101, 125, 100, 148, '+')

    def test_crop_then_owned_stitch_matches_whole_oriented_rows(self):
        # Independent desired global token-index rows, not a second crop implementation.
        windows = ((100, 124), (112, 136), (124, 148))
        geometry = dict(local_windows=((0, 24), (12, 36), (24, 48)),
            base_owned_segments=G.context_segments(48, ((0, 24), (12, 36), (24, 48))))
        for strand in ('+', '-'):
            values = torch.arange(8).reshape(1, 8, 1)
            pieces = [crop_features(values, l, r, 100, 148, strand)[0].repeat_interleave(6, 0).expand(-1, 15)
                      for l, r in windows]
            actual = G.owned_context_logits(pieces, geometry, strand)
            self.assertTrue(torch.equal(actual, torch.arange(8).repeat_interleave(6)[:, None].expand(-1, 15)))
        self.assertEqual(long_bounds(30000, 79152, 200000), (5424, 103728))
        with self.assertRaises(ValueError):
            long_bounds(0, 49152, 200000)


def checks():
    began = time.perf_counter()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CropChecks))
    write(CHECKS/'summary.json', dict(completed=result.wasSuccessful(), job_id=os.environ['SLURM_JOB_ID'],
        configuration=CONFIG, tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        research_data_read=False, models_or_weights_read=False, optimizer_steps_this_run=0,
        seconds=time.perf_counter()-began))
    if not result.wasSuccessful():
        raise SystemExit(1)


def inputs():
    # Reuse the actual original witness identities, not annotation-derived bounds.
    source = json.loads((G.OUT/'summary.json').read_text())
    if ((G.OUT/'STATUS').read_text().strip() != 'COMPLETED' or source['job_id'] != '13429053'
            or not source['completed'] or source['configuration'] != G.CONFIG or source['DEV_used']):
        raise ValueError('Actual completed original eight-witness source required')
    draws, index, windows, catalogs = G.Q.load_train_inputs()
    if len(draws) != 1536 or len(windows) != 1491 or len(source['records']) != 8:
        raise ValueError('Unchanged original TRAIN population and all eight witnesses required')
    sequences, descriptors = G.P.read_train_DNA()
    return source, index, windows, catalogs, sequences, descriptors


def oriented(sequence, a, b, strand):
    value = re.sub('[^ACGT]', 'N', sequence[a:b])
    return value if strand == '+' else G.Q.M.rc(value)


@torch.inference_mode()
def joint_encode(backbone, tokenizer, sequence):
    ids = tokenizer(sequence, add_special_tokens=False, return_tensors='pt', truncation=False)['input_ids']
    if tuple(ids.shape) != (1, len(sequence)//6) or len(sequence) % 6:
        raise ValueError('Native long sequence must retain all six-mer tokens without specials')
    # Fail if a fused kernel is unavailable rather than silently forming an N^2 attention matrix.
    with sdpa_kernel(SDPBackend.FLASH_ATTENTION):
        hidden = backbone(input_ids=ids.to('cuda'), attention_mask=torch.ones_like(ids, device='cuda')).last_hidden_state
    if tuple(hidden.shape) != (1, len(sequence)//6, 2048) or not bool(torch.isfinite(hidden).all()):
        raise ValueError('Actual finite untruncated native long-context hidden states required')
    return hidden.to(torch.bfloat16).float(), ids


@torch.inference_mode()
def native():
    began = time.perf_counter()
    checked = json.loads((CHECKS/'summary.json').read_text())
    if (not checked['completed'] or checked['configuration'] != CONFIG
            or (CHECKS/'STATUS').read_text().strip() != 'COMPLETED'):
        raise ValueError('New coordinate/crop checks must pass first')
    source, index, windows, _, sequences, descriptors = inputs()
    saved = torch.load(D.I.CHECKPOINT, map_location='cpu', weights_only=True)
    D.I.R.validate_checkpoint(saved)
    model = G.M28Core('C0', 2048, 128)
    model.load_state_dict(saved['model'], strict=True)
    del saved
    model.eval().requires_grad_(False).to('cuda')
    backbone, tokenizer = G.load_native()
    write(OUT/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        protocol=PROTOCOL, source_job='13429053', witness_ids=[r['window_id'] for r in source['records']]))
    records, opened = [], set()
    torch.cuda.reset_peak_memory_stats()
    with (OUT/'windows.jsonl').open('x') as log:
        for ordinal, original in enumerate(source['records']):
            started = time.perf_counter()
            w = windows[original['window_id']]
            chromosome = sequences[w['species']][w['seqid']]
            geometry = G.context_geometry(w['start'], w['end'], len(chromosome))
            a, b = geometry['genomic_context']; c, d = long_bounds(a, b, len(chromosome))
            if (json.loads(json.dumps(geometry)) != original['geometry'] or w['valid_bases'] != 24576
                    or any(r-l != 24576 for l, r in geometry['genomic_windows'])):
                raise ValueError('Original fixed full-length three-window view required')
            raw_core = chromosome[w['start']:w['end']]
            core = raw_core if w['strand'] == '+' else G.Q.M.rc(raw_core)
            cached = G.W.load_cache(w, index, opened)
            if cached['sequence'] != core:
                raise ValueError('Original core DNA differs from raw TRAIN FASTA')
            arms, core_features, encoding = {}, {}, {}
            for arm in ARMS:
                t = time.perf_counter(); whole = ids = None
                if arm in ('JOINT49152', 'JOINT98304'):
                    lo, hi = (a, b) if arm == 'JOINT49152' else (c, d)
                    whole, ids = joint_encode(backbone, tokenizer, oriented(chromosome, lo, hi, w['strand']))
                pieces = []
                for left, right in geometry['genomic_windows']:
                    dna_string = oriented(chromosome, left, right, w['strand'])
                    if arm == 'BLOCK6144':
                        features = G.encode_window(backbone, tokenizer, dna_string, 'cuda')
                    elif arm == 'JOINT24576':
                        features, _ = joint_encode(backbone, tokenizer, dna_string)
                    else:
                        features = crop_features(whole, left, right, lo, hi, w['strand'])
                        offset = left-lo if w['strand'] == '+' else hi-right
                        expected = tokenizer(dna_string, add_special_tokens=False,
                            return_tensors='pt', truncation=False)['input_ids']
                        if not torch.equal(ids[:, offset//6:(offset+24576)//6], expected):
                            raise ValueError('Long/native local six-mer registration changed')
                    if (left, right) == (w['start'], w['end']):
                        core_features[arm] = features.detach().cpu()
                        if arm == 'BLOCK6144' and not torch.equal(features[0].cpu().to(torch.bfloat16), cached['features']):
                            raise ValueError('6144 baseline fails exact original BF16 cache replay')
                    dna = G.one_hot(dna_string, 24576)[None].to('cuda')
                    pieces.append(model(features, dna, [24576])['segmentation_logits'][0].cpu())
                joined = G.owned_context_logits(pieces, geometry, w['strand'])
                if joined.shape != (49152, 15) or not bool(torch.isfinite(joined).all()):
                    raise ValueError('Fixed finite FP32 fifteen-role view required for every arm')
                arms[arm] = joined
                torch.cuda.synchronize(); encoding[arm] = time.perf_counter()-t
                del whole, ids, features, dna, pieces
            feature_changes = {}
            baseline = core_features['BLOCK6144']
            for arm in ARMS[1:]:
                value = core_features[arm]
                feature_changes[arm] = dict(relative_RMS=float((value-baseline).square().mean().sqrt()/baseline.square().mean().sqrt()),
                    mean_token_cosine=float(torch.nn.functional.cosine_similarity(value, baseline, dim=-1).mean()),
                    logits_mean_abs_difference=float((arms[arm]-arms['BLOCK6144']).abs().mean()))
            path = OUT/('window_%02d.pt' % ordinal)
            torch.save(dict(configuration=CONFIG, window_id=w['id'], geometry=geometry, arms=arms), path)
            row = dict(ordinal=ordinal, window_id=w['id'], species=w['species'], seqid=w['seqid'], strand=w['strand'],
                geometry=geometry, long_encoder_genomic_bounds=[c, d], artifact=path.name,
                artifact_bytes=path.stat().st_size, encoding_and_head_seconds=encoding,
                core_features_relative_to_original_short_blocks=feature_changes,
                baseline_core_BF16_cache_replay_exact=True, six_mer_registration_exact=True,
                native_max_tokens=16384, seconds=time.perf_counter()-started)
            records.append(row); log.write(json.dumps(row, allow_nan=False)+'\n'); log.flush()
            print(json.dumps(row, allow_nan=False), flush=True)
            del cached, arms, core_features, baseline, value
    if any(p.requires_grad or p.grad is not None for module in (model, backbone) for p in module.parameters()):
        raise ValueError('All parameters must remain frozen without gradients')
    write(OUT/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        protocol=PROTOCOL, records=records, actual_witnesses=len(records), native_conditions=32,
        source_job='13429053', actual_last4608_checkpoint_read=True, all_baseline_core_replays_exact=True,
        original_TRAIN_DNA_sources=descriptors, cached_artifacts_opened=sorted(opened),
        native_model_code_or_RoPE_changed=False, optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(), peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved(),
        GPU=torch.cuda.get_device_name(), seconds=time.perf_counter()-began,
        interpretation='native_context_wiring_and_paired_potentials_not_performance_or_generalization'))


@torch.no_grad()
def readout():
    began = time.perf_counter()
    native_source = json.loads((OUT/'summary.json').read_text())
    if (not native_source['completed'] or native_source['configuration'] != CONFIG
            or (OUT/'STATUS').read_text().strip() != 'COMPLETED' or native_source['actual_witnesses'] != 8):
        raise ValueError('All eight actual native context outputs required')
    source, _, windows, catalogs, sequences, _ = inputs()
    if [r['window_id'] for r in source['records']] != [r['window_id'] for r in native_source['records']]:
        raise ValueError('Same ordered eight witnesses required')
    # A changed baseline would invalidate attribution to native context length.
    legacy = {}
    for seqid in {r['seqid'] for r in native_source['records']}:
        root = D.OUT/seqid
        if (root/'STATUS').read_text().strip() != 'COMPLETED':
            raise ValueError('Completed original short-block finite-view baseline required')
        with (root/'windows.jsonl').open() as handle:
            for line in handle:
                row = json.loads(line); legacy[row['window_id']] = row
    rows = []
    with (READOUT/'windows.jsonl').open('x') as log:
        for original, emitted in zip(source['records'], native_source['records']):
            w = windows[original['window_id']]; n = w['valid_bases']
            chromosome = sequences[w['species']][w['seqid']]
            geometry = G.context_geometry(w['start'], w['end'], len(chromosome))
            a, b = geometry['genomic_context']; context = oriented(chromosome, a, b, w['strand'])
            core = chromosome[w['start']:w['end']]
            if w['strand'] == '-': core = G.Q.M.rc(core)
            target = G.P.P.build_targets(w, catalogs[w['species']], core)
            offset = w['start']-a if w['strand'] == '+' else b-w['end']
            saved = torch.load(OUT/emitted['artifact'], map_location='cpu', weights_only=True)
            if saved['configuration'] != CONFIG or saved['window_id'] != w['id'] or saved['geometry'] != geometry:
                raise ValueError('Native potentials must match the exact fixed readout geometry')
            references = [catalogs[w['species']][i] for i in w['overlapping_transcript_ids']
                if catalogs[w['species']][i]['training_primary']
                and a <= catalogs[w['species']][i]['CDS'][0][0]
                and catalogs[w['species']][i]['CDS'][-1][1] <= b]
            for arm in ARMS:
                started = time.perf_counter(); logits = saved['arms'][arm]
                likelihood = D.coarse_channel_nll(context, logits, 128, offset, offset+n, target, .1,
                    cuts=(offset, offset+n), checkpoint_bp=256)
                raw = float(likelihood['loss'])
                if not math.isfinite(raw) or raw < -1e-8:
                    raise ValueError('Actual unchanged observation NLL must be finite and nonnegative')
                local = local_coarse_channel_nll(logits[offset:offset+n], target, .1)
                predicted = D.decode_full_grammar_skyline(context, logits, 128)
                roles = D.map_roles(predicted['chains'], len(context))
                agreement = D.core_agreement(roles[offset:offset+n], target)
                if tuple(agreement[k] for k in ('fine_positions', 'coarse_positions', 'fully_unknown_positions')) != tuple(
                        original[k] for k in ('fine_positions', 'coarse_positions', 'fully_unknown_positions')):
                    raise ValueError('Every original fine/coarse/unknown core row must remain unchanged')
                baseline_difference = None
                if arm == 'BLOCK6144':
                    old = legacy[w['id']]
                    baseline_difference = abs(raw-old['raw_NLL'])
                    if (baseline_difference > 1e-6 or agreement['fine_confusion'] != old['fine_confusion']
                            or agreement['coarse_confusion'] != old['coarse_confusion']):
                        raise ValueError('Short-block NLL/MAP no longer reproduces the completed fixed-weight baseline')
                selected = {tuple(D.I.R.genomic_chain(chain, a, b, w['strand'])) for chain in predicted['chains']}
                hits = sum(tuple((x, y) for x, y, _ in tx['CDS']) in selected for tx in references)
                row = dict(window_id=w['id'], species=w['species'], seqid=w['seqid'], strand=w['strand'], arm=arm,
                    raw_NLL=raw, NLL_per_known_core_row=raw/max(1, agreement['fine_positions']+agreement['coarse_positions']),
                    local_observation_NLL=float(local['loss']), **agreement,
                    contained_original_primary_core_overlapping_reference_pairs=len(references),
                    exact_reference_chain_hits=hits, selected_view_chains=len(predicted['chains']),
                    original_short_block_NLL_abs_difference=baseline_difference,
                    seconds=time.perf_counter()-started)
                rows.append(row); log.write(json.dumps(row, allow_nan=False)+'\n'); log.flush()
                print(json.dumps({k:row[k] for k in ('window_id', 'arm', 'NLL_per_known_core_row', 'exact_reference_chain_hits', 'seconds')}), flush=True)
    aggregate = {}
    for species in sorted(G.Q.ALLOW):
        aggregate[species] = {}
        for arm in ARMS:
            group = [r for r in rows if r['species'] == species and r['arm'] == arm]
            fine = np.sum([r['fine_confusion'] for r in group], axis=0)
            known = sum(r['fine_positions']+r['coarse_positions'] for r in group)
            aggregate[species][arm] = dict(witnesses=len(group), raw_NLL=sum(r['raw_NLL'] for r in group),
                NLL_per_known_core_row=sum(r['raw_NLL'] for r in group)/max(1, known),
                known_core_rows=known, fine_positions=int(fine.sum()), fine_matches=int(np.trace(fine)),
                fine_reference_rows_by_role=fine.sum(1).tolist(), fine_matches_by_role=np.diag(fine).tolist(),
                reference_window_pairs=sum(r['contained_original_primary_core_overlapping_reference_pairs'] for r in group),
                exact_reference_chain_hits=sum(r['exact_reference_chain_hits'] for r in group))
    write(READOUT/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        native_source_job=native_source['job_id'], readout_conditions=len(rows), aggregate=aggregate,
        records=rows, optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
        original_short_block_NLL_and_MAP_reproduced=True,
        checkpoint_or_context_selected=False, seconds=time.perf_counter()-began,
        interpretation='eight_fixed_TRAIN_witness_pairs_no_unique_gene_F1_no_full_grid_or_generalization_claim'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('checks', 'native', 'readout'))
    args = parser.parse_args(); torch.set_num_threads(2)
    {'checks': checks, 'native': native, 'readout': readout}[args.mode]()
