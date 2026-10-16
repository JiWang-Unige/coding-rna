#!/usr/bin/env python3
"""Frozen native CDS heads: all six within-token rotations on eight TRAIN cores."""
import json
import math
import os
import time
import unittest

import numpy as np
import torch
from torch import nn
from torch.nn.attention import SDPBackend, sdpa_kernel
from transformers import AutoModelForTokenClassification, AutoTokenizer
import diagnose_full_grammar_native_context as N

G, D, ROOT = N.G, N.D, N.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CDS-SLOT-DIAGNOSTIC'
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_native_cds_slots.md'
STRATA = ('fine_known', 'coarse_only', 'fine_complete_homogeneous_tokens', 'fine_complete_mixed_tokens')
CONFIG = dict(scope='same_original_eight_TRAIN_witness_cores_no_resampling',
    encoder_model=G.Q.F.MODEL, encoder_revision=G.Q.F.REVISION,
    encoder_frozen=True, native_classification_heads_frozen=True,
    encoder_dtype='BF16_hidden_and_native_linear_then_FP32_saved_logits',
    arms=list(N.ARMS), head_indices=[0, 1], head_selected_from_labels=False,
    head_semantics='report_native_indices_separately_no_assumed_strand_assignment',
    readout='native_2048_to_6x2_linear_core_only_no_M28_head_or_grammar',
    native_attention='unchanged_SDPA_flash_only_no_RoPE_change',
    rotations=list(range(6)), rotation='cyclic_within_each_oriented_six_mer_preserve_score_multiset',
    target='original_fine_first_else_coarse_else_unknown_binary_CDS_vs_non_CDS',
    fine_role_to_region='unchanged_src_m28_labels_REGION_COLUMNS',
    binary_CE='descriptive_per_base_not_epsilon_channel_or_full_grammar_NLL',
    decision='native_argmax_two_labels_no_threshold_tuning',
    all_six_fine_known_tokens_stratified_by_binary_label_homogeneity=True,
    potential_artifact_limit_bytes=16*1024**2,
    optimizer_steps=0, DEV_used=False, sealed_data_used=False)


def binary_targets(targets):
    fine = torch.as_tensor(targets['seg'], dtype=torch.long)
    region = torch.as_tensor(targets['region'], dtype=torch.long)
    if fine.ndim != 1 or region.shape != fine.shape or len(fine) % 6:
        raise ValueError('Original complete six-mer core label rows required')
    exact = fine >= 0; coarse = (~exact) & (region >= 0)
    binary = torch.full_like(fine, -100)
    binary[exact] = torch.as_tensor(D.REGION)[fine[exact]].eq(1).long()
    binary[coarse] = region[coarse].eq(1).long()
    complete = exact.reshape(-1, 6).all(1)
    same = binary.reshape(-1, 6).amin(1) == binary.reshape(-1, 6).amax(1)
    masks = dict(fine_known=exact, coarse_only=coarse,
        fine_complete_homogeneous_tokens=(complete & same).repeat_interleave(6),
        fine_complete_mixed_tokens=(complete & ~same).repeat_interleave(6))
    counts = dict(fine_positions=int(exact.sum()), coarse_positions=int(coarse.sum()),
        fully_unknown_positions=int((binary < 0).sum()),
        fine_role_rows=torch.bincount(fine[exact], minlength=15).tolist())
    return binary, masks, counts


def rotate_slots(logits, shift):
    if logits.ndim != 2 or logits.shape[1] != 2 or len(logits) % 6 or not 0 <= shift < 6:
        raise ValueError('Complete two-class six-mer prediction rows and one declared rotation required')
    return logits.reshape(-1, 6, 2).roll(shift, 1).reshape(-1, 2)


def measure(logits, binary, masks):
    logp = logits.double().log_softmax(-1); prediction = logits.argmax(-1)
    values = {}
    for name, mask in masks.items():
        labels = binary[mask]; predictions = prediction[mask]
        table = torch.bincount(2*labels+predictions, minlength=4).reshape(2, 2)
        count = int(mask.sum())
        values[name] = dict(rows=count, binary_CE_sum=float(-logp[mask, labels].sum()),
            confusion=table.tolist())
    return values


def pooled(entries):
    values = {}
    for name in STRATA:
        count = sum(r['strata'][name]['rows'] for r in entries)
        ce = math.fsum(r['strata'][name]['binary_CE_sum'] for r in entries)
        table = np.sum([r['strata'][name]['confusion'] for r in entries], axis=0)
        tn, fp, fn, tp = map(int, table.flatten())
        values[name] = dict(rows=count, binary_CE_sum=ce,
            binary_CE_per_row=ce/count if count else None,
            confusion=table.tolist(), CDS_base_F1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None,
            CDS_rows=tp+fn, non_CDS_rows=tn+fp)
    return values


class SlotChecks(unittest.TestCase):
    def test_original_fine_precedence_and_unknown(self):
        labels, masks, counts = binary_targets(dict(seg=[13, 14, 4, 0, -100, -100],
            region=[2, 0, 1, 1, 1, -100]))
        self.assertEqual(labels.tolist(), [1, 1, 0, 0, 1, -100])
        self.assertEqual((counts['fine_positions'], counts['coarse_positions'], counts['fully_unknown_positions']), (4, 1, 1))
        self.assertEqual(int(masks['fine_complete_mixed_tokens'].sum()), 0)

    def test_all_rotations_preserve_scores_and_detect_misregistration(self):
        targets = dict(seg=[0, 1, 1, 0, 0, 1]+[1]*6, region=[0]*12)
        binary, masks, _ = binary_targets(targets)
        self.assertEqual(int(masks['fine_complete_mixed_tokens'].sum()), 6)
        self.assertEqual(int(masks['fine_complete_homogeneous_tokens'].sum()), 6)
        score = torch.stack([2.-4.*binary, 4.*binary-2.], 1)
        score[6:] += torch.arange(6)[:, None]/10
        reference = measure(score, binary, masks)
        alternatives = []
        for shift in range(6):
            shifted = rotate_slots(score, shift)
            self.assertTrue(torch.equal(score.reshape(-1, 6, 2).sort(1).values,
                shifted.reshape(-1, 6, 2).sort(1).values))
            measured = measure(shifted, binary, masks)
            self.assertAlmostEqual(reference['fine_complete_homogeneous_tokens']['binary_CE_sum'],
                measured['fine_complete_homogeneous_tokens']['binary_CE_sum'], places=10)
            if shift:
                alternatives.append(measured['fine_complete_mixed_tokens']['binary_CE_sum'])
        self.assertGreater(sum(alternatives)/5, reference['fine_complete_mixed_tokens']['binary_CE_sum'])


@torch.inference_mode()
def main():
    began = time.perf_counter()
    tested = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(SlotChecks))
    if not tested.wasSuccessful():
        raise SystemExit(1)
    source, _, windows, catalogs, sequences, _ = N.inputs()
    previous = json.loads((N.OUT/'summary.json').read_text())
    if (not previous['completed'] or previous['job_id'] != '13452550' or previous['configuration'] != N.CONFIG
            or [r['window_id'] for r in previous['records']] != [r['window_id'] for r in source['records']]):
        raise ValueError('The completed original eight native identities must be retained')
    N.write(OUT/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        protocol=PROTOCOL, witness_ids=[r['window_id'] for r in source['records']]))
    tokenizer = AutoTokenizer.from_pretrained(G.Q.F.MODEL, revision=G.Q.F.REVISION,
        trust_remote_code=True, local_files_only=True)
    full = AutoModelForTokenClassification.from_pretrained(G.Q.F.MODEL, revision=G.Q.F.REVISION,
        trust_remote_code=True, local_files_only=True, attn_implementation='sdpa', torch_dtype=torch.bfloat16)
    if (full.k != 6 or full.num_labels != 2 or full.num_prediction_heads != 2 or full.feature_layer != -1
            or full.config.id2label != {0:'NON_CODING', 1:'CDS'} or not isinstance(full.score, nn.ModuleList)
            or any(not isinstance(head, nn.Linear) or head.in_features != 2048 or head.out_features != 12
                   or head.bias is not None for head in full.score)):
        raise ValueError('Exact released two six-slot linear CDS heads required')
    full.eval().requires_grad_(False).to('cuda'); backbone = full.model
    # Capture the actual hidden states used by the public forward once; do not
    # replace its forward or infer layout from parameter shapes alone.
    first = windows[source['records'][0]['window_id']]
    core = N.oriented(sequences[first['species']][first['seqid']], first['start'], first['end'], first['strand'])
    ids = tokenizer(core[:6144], add_special_tokens=False, return_tensors='pt', truncation=False)['input_ids'].to('cuda')
    captured = {}
    hook = backbone.register_forward_hook(lambda module, args, output: captured.update(hidden=output[0]))
    try:
        with sdpa_kernel(SDPBackend.FLASH_ATTENTION):
            public = full(input_ids=ids, attention_mask=torch.ones_like(ids)).logits
        rebuilt = torch.cat([head(captured['hidden']).reshape(1, 6144, 2) for head in full.score], 1)
        if public.shape != (1, 12288, 2) or not torch.equal(public, rebuilt) or not bool(torch.isfinite(public).all()):
            raise ValueError('Actual released forward layout differs from two native six-slot projections')
    finally:
        hook.remove()
    del public, rebuilt, captured, ids
    torch.cuda.reset_peak_memory_stats(); results = []; descriptors = []; artifact_bytes = 0
    with (OUT/'measurements.jsonl').open('x') as log:
        for ordinal, (original, preceding) in enumerate(zip(source['records'], previous['records'])):
            w = windows[original['window_id']]; chromosome = sequences[w['species']][w['seqid']]
            geometry = G.context_geometry(w['start'], w['end'], len(chromosome)); a, b = geometry['genomic_context']
            c, d = N.long_bounds(a, b, len(chromosome))
            if json.loads(json.dumps(geometry)) != preceding['geometry'] or w['valid_bases'] != 24576:
                raise ValueError('Unchanged original eight interior witness views required')
            raw = chromosome[w['start']:w['end']]; oriented_raw = raw if w['strand'] == '+' else G.Q.M.rc(raw)
            target = G.P.P.build_targets(w, catalogs[w['species']], oriented_raw)
            binary, masks, counts = binary_targets(target)
            if any(counts[key] != original[key] for key in ('fine_positions', 'coarse_positions', 'fully_unknown_positions')):
                raise ValueError('Original witness observation denominators changed')
            dna = N.oriented(chromosome, w['start'], w['end'], w['strand']); arms = {}
            for arm in N.ARMS:
                if arm == 'BLOCK6144':
                    features = G.encode_window(backbone, tokenizer, dna, 'cuda')
                elif arm == 'JOINT24576':
                    features, _ = N.joint_encode(backbone, tokenizer, dna)
                else:
                    lo, hi = (a, b) if arm == 'JOINT49152' else (c, d)
                    whole, token_ids = N.joint_encode(backbone, tokenizer, N.oriented(chromosome, lo, hi, w['strand']))
                    features = N.crop_features(whole, w['start'], w['end'], lo, hi, w['strand'])
                    offset = w['start']-lo if w['strand'] == '+' else hi-w['end']
                    expected = tokenizer(dna, add_special_tokens=False, return_tensors='pt', truncation=False)['input_ids']
                    if not torch.equal(expected, token_ids[:, offset//6:(offset+24576)//6]):
                        raise ValueError('Native/core six-mer token registration changed')
                scores = torch.stack([head(features.to(torch.bfloat16)).reshape(24576, 2).float().cpu()
                    for head in full.score])
                if scores.shape != (2, 24576, 2) or not bool(torch.isfinite(scores).all()):
                    raise ValueError('Actual finite six-slot native binary core logits required')
                arms[arm] = scores
                for head_index, score in enumerate(scores):
                    reference = measure(score, binary, masks)
                    for shift in range(6):
                        shifted = rotate_slots(score, shift)
                        measured = reference if shift == 0 else measure(shifted, binary, masks)
                        if (any(measured[name]['rows'] != reference[name]['rows'] for name in STRATA)
                                or abs(measured['fine_complete_homogeneous_tokens']['binary_CE_sum']-
                                       reference['fine_complete_homogeneous_tokens']['binary_CE_sum']) > 1e-8):
                            raise ValueError('Rotation changed observation counts or homogeneous-token CE')
                        result = dict(window_id=w['id'], species=w['species'], seqid=w['seqid'], strand=w['strand'],
                            arm=arm, head_index=head_index, rotation=shift, strata=measured)
                        results.append(result); log.write(json.dumps(result, allow_nan=False)+'\n')
                del features
                if arm in ('JOINT49152', 'JOINT98304'):
                    del whole, token_ids
            artifact = 'window_%02d.pt' % ordinal
            torch.save(dict(configuration=CONFIG, window_id=w['id'], arms=arms, binary_labels=binary,
                strata_masks=masks), OUT/artifact)
            size = (OUT/artifact).stat().st_size; artifact_bytes += size
            descriptor = dict(ordinal=ordinal, window_id=w['id'], species=w['species'], seqid=w['seqid'],
                strand=w['strand'], geometry=geometry, **counts, artifact=artifact, artifact_bytes=size)
            descriptors.append(descriptor); log.flush()
            print(json.dumps(dict(stage='native_CDS_slot_diagnostic', completed_witnesses=ordinal+1,
                total_witnesses=8, seconds=time.perf_counter()-began)), flush=True)
            del arms, scores
    if len(descriptors) != 8 or len(results) != 384 or artifact_bytes > CONFIG['potential_artifact_limit_bytes']:
        raise ValueError('All eight witnesses, four contexts, two heads and six rotations within declared output size required')
    if any(p.requires_grad or p.grad is not None for p in full.parameters()):
        raise ValueError('Released backbone and both native CDS heads must stay frozen without gradients')
    N.write(OUT/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        protocol=PROTOCOL, tests_run=tested.testsRun, failures=len(tested.failures), errors=len(tested.errors),
        actual_witnesses=8, native_conditions=32, head_context_conditions=64, rotation_conditions=384,
        records=descriptors, native_public_forward_projection_replay_exact=True,
        token_registration_and_observation_denominators_unchanged=True, homogeneous_token_CE_rotation_invariant=True,
        per_species={species:{arm:{str(head):{str(shift):pooled([r for r in results if r['species'] == species
            and r['arm'] == arm and r['head_index'] == head and r['rotation'] == shift])
            for shift in range(6)} for head in range(2)} for arm in N.ARMS} for species in sorted(G.Q.ALLOW)},
        potential_artifact_bytes=artifact_bytes, peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        optimizer_steps_this_run=0, M28_checkpoint_read=False, native_head_selected=False,
        DEV_used=False, sealed_data_used=False, seconds=time.perf_counter()-began,
        interpretation='purpose_selected_TRAIN_native_binary_position_evidence_not_M28_causal_failure_or_gene_F1'))


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
