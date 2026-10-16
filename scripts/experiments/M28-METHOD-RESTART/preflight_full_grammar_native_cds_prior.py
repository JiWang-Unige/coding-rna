#!/usr/bin/env python3
"""Eight original TRAIN views: zero-prior C0 replay and real full-grammar VJP."""
import json
import math
import os
import time
import unittest

import torch
from torch import nn
from transformers import AutoModelForTokenClassification, AutoTokenizer
import diagnose_full_grammar_native_context as N
import diagnose_full_grammar_native_cds_slots as S
from native_cds_role_prior import NativeCDSRolePrior, native_cds_log_odds

G, D, ROOT = N.G, N.D, N.ROOT
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-NATIVE-CDS-PRIOR-PREFLIGHT'
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_native_cds_prior_preflight.md'
CONFIG = dict(scope='same_original_eight_TRAIN_witness_views_no_resampling',
    head_checkpoint=str(D.I.CHECKPOINT.relative_to(ROOT)), starting_head_step=4608,
    encoder_model=G.Q.F.MODEL, encoder_revision=G.Q.F.REVISION,
    encoder_blocks_bp=6144, encoder_frozen=True, native_heads_frozen=True,
    native_head_indices=[0, 1], head_selected_from_labels=False,
    evidence='native_BF16_two_head_6x2_CDS_minus_NON_CODING_FP32_per_base',
    prior='zero_initialized_Linear_2_to14_no_bias_pad_unchanged_background_role0',
    additional_parameters=28, total_trainable_parameters=366219,
    original_C0_core_or_DNA_branch_changed=False,
    grammar_view='unchanged_original_three24576_nearest_center_owned49152',
    observation='unchanged_fine_first_else_coarse_epsilon_0.1_else_unknown',
    max_exons=128, checkpoint_bp=256, score_accumulation='FP64',
    real_baseline_replay='exact_eval_owned_potentials_vs_original13452550_BLOCK6144',
    gradient_mode='candidate_train_mode_for_cuDNN_GRU_backward_native_model_eval_frozen',
    gradient_check='independent_local_role_VJP_transpose_times_both_native_margins',
    VJP_rtol=1e-4, VJP_atol=1e-3, optimizer_steps=0,
    trained_or_formal_checkpoint_created=False, DEV_used=False, sealed_data_used=False)


class PriorChecks(unittest.TestCase):
    def test_zero_original_replay_relative_background_and_padding(self):
        torch.manual_seed(0)
        baseline = G.M28Core('C0', 12, 16).eval()
        candidate = NativeCDSRolePrior(12, 16).eval(); candidate.load_c0_weights(baseline.state_dict())
        features = torch.randn(2, 2, 12)
        dna = torch.nn.functional.one_hot(torch.randint(5, (2, 12)), 5).float()
        evidence = torch.ones(2, 12, 2); evidence[..., 1] = 2
        valid = [7, 12]
        with torch.no_grad():
            old = baseline(features, dna, valid)['segmentation_logits']
            new = candidate(features, dna, valid, evidence)['segmentation_logits']
            self.assertTrue(torch.equal(old, new))
            candidate.native_cds_role_prior.weight.fill_(1)
            delta = candidate(features, dna, valid, evidence)['segmentation_logits']-old
        self.assertTrue(torch.equal(delta[..., 0], torch.zeros(2, 12)))
        self.assertTrue(torch.equal(delta[0, 7:], torch.zeros(5, 15)))
        self.assertTrue(torch.allclose(delta[0, :7, 1:], torch.full((7, 14), 3.)))

    def test_two_head_and_six_slot_projection_order(self):
        heads = nn.ModuleList([nn.Linear(12, 12, bias=False).bfloat16() for _ in range(2)])
        with torch.no_grad():
            heads[0].weight.copy_(torch.eye(12)); heads[1].weight.copy_(2*torch.eye(12))
        heads.eval().requires_grad_(False)
        features = torch.zeros(1, 2, 12)
        features.reshape(1, 12, 2)[..., 1] = torch.arange(1, 13)
        actual = native_cds_log_odds(heads, features)
        expected = torch.stack([torch.arange(1, 13), 2*torch.arange(1, 13)], -1).float()[None]
        self.assertTrue(torch.equal(actual, expected))
        self.assertFalse(actual.requires_grad or torch.is_inference(actual))


def required_source(root, job_id, configuration):
    value = json.loads((root/'summary.json').read_text())
    if ((root/'STATUS').read_text().strip() != 'COMPLETED' or not value['completed']
            or value['job_id'] != job_id or value['configuration'] != configuration):
        raise ValueError('Actual completed unchanged source required: '+str(root))
    return value


def main():
    began = time.perf_counter()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(PriorChecks))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    native_source = required_source(N.OUT, '13452550', N.CONFIG)
    readout_source = required_source(N.READOUT, '13452596', N.CONFIG)
    slot_source = required_source(S.OUT, '13452867', S.CONFIG)
    if not slot_source['native_public_forward_projection_replay_exact']:
        raise ValueError('Actual public native-head layout must already be verified')
    source, _, windows, catalogs, sequences, _ = N.inputs()
    ids = [r['window_id'] for r in source['records']]
    if (ids != [r['window_id'] for r in native_source['records']]
            or ids != [r['window_id'] for r in slot_source['records']]):
        raise ValueError('Same complete ordered eight native witnesses required')
    old_NLL = {r['window_id']:r['raw_NLL'] for r in readout_source['records'] if r['arm'] == 'BLOCK6144'}
    N.write(OUT/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        protocol=PROTOCOL, witness_ids=ids, native_potential_source_job='13452550', native_slot_source_job='13452867'))
    saved = torch.load(D.I.CHECKPOINT, map_location='cpu', weights_only=True)
    D.I.R.validate_checkpoint(saved)
    initial = {k:v.clone() for k, v in saved['model'].items()}; del saved
    baseline = G.M28Core('C0', 2048, 128); baseline.load_state_dict(initial, strict=True)
    baseline.eval().requires_grad_(False).to('cuda')
    candidate = NativeCDSRolePrior(); candidate.load_c0_weights(initial); candidate.to('cuda')
    if sum(p.numel() for p in candidate.parameters()) != CONFIG['total_trainable_parameters']:
        raise ValueError('Original366191 plus exactly28 prior coefficients required')
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
    records = []; torch.cuda.reset_peak_memory_stats()
    with (OUT/'windows.jsonl').open('x') as log:
        for ordinal, (original, emitted, slot) in enumerate(zip(source['records'], native_source['records'], slot_source['records'])):
            t = time.perf_counter(); w = windows[original['window_id']]; chromosome = sequences[w['species']][w['seqid']]
            geometry = G.context_geometry(w['start'], w['end'], len(chromosome)); a, b = geometry['genomic_context']
            if json.loads(json.dumps(geometry)) != emitted['geometry'] or w['valid_bases'] != 24576:
                raise ValueError('Unchanged original full-length witness view required')
            offset = w['start']-a if w['strand'] == '+' else b-w['end']
            core = chromosome[w['start']:w['end']]
            if w['strand'] == '-':
                core = G.Q.M.rc(core)
            target = G.P.P.build_targets(w, catalogs[w['species']], core)
            _, _, counts = S.binary_targets(target)
            if any(counts[k] != original[k] for k in ('fine_positions', 'coarse_positions', 'fully_unknown_positions')):
                raise ValueError('Original witness fine/coarse/unknown denominator changed')
            context = N.oriented(chromosome, a, b, w['strand'])
            old_potentials = torch.load(N.OUT/emitted['artifact'], map_location='cpu', weights_only=True)
            old_slots = torch.load(S.OUT/slot['artifact'], map_location='cpu', weights_only=True)
            features_list = []; dna_list = []; priors = []; eval_pieces = []; candidate.eval(); core_replayed = False
            for left, right in geometry['genomic_windows']:
                dna_string = N.oriented(chromosome, left, right, w['strand'])
                # Clone outside the frozen inference encoder: trainable layers
                # must be able to save these normal constant inputs for backward.
                features = G.encode_window(backbone, tokenizer, dna_string, 'cuda').clone()
                evidence = native_cds_log_odds(full.score, features)
                if (features.requires_grad or evidence.requires_grad or torch.is_inference(features)
                        or torch.is_inference(evidence) or not bool(torch.isfinite(evidence).all())):
                    raise ValueError('Finite normal frozen encoder/prior inputs required for trainable layers')
                dna = G.one_hot(dna_string, 24576)[None].to('cuda')
                with torch.no_grad():
                    old = baseline(features, dna, [24576])['segmentation_logits'][0].cpu()
                    current = candidate(features, dna, [24576], evidence)['segmentation_logits'][0].cpu()
                if not torch.equal(old, current):
                    raise ValueError('Zero prior changed actual original local C0 predictions')
                if (left, right) == (w['start'], w['end']):
                    old_scores = old_slots['arms']['BLOCK6144']
                    expected = (old_scores[..., 1]-old_scores[..., 0]).T[None]
                    if not torch.equal(evidence.cpu(), expected):
                        raise ValueError('Fresh two-head per-base margins differ from preceding native slot outputs')
                    core_replayed = True
                features_list.append(features); dna_list.append(dna); priors.append(evidence); eval_pieces.append(current)
            eval_joined = G.owned_context_logits(eval_pieces, geometry, w['strand'])
            if not core_replayed or not torch.equal(eval_joined, old_potentials['arms']['BLOCK6144']):
                raise ValueError('Zero prior failed exact preceding last4608 owned-potential replay')
            candidate.train(); candidate.zero_grad(set_to_none=True); tensors = []
            for features, dna, evidence in zip(features_list, dna_list, priors):
                logits = candidate(features, dna, [24576], evidence)['segmentation_logits'][0]
                logits.retain_grad(); tensors.append(logits)
            joined = G.owned_context_logits(tensors, geometry, w['strand']).cpu()
            train_eval_max = float((joined.detach()-eval_joined).abs().max())
            forward_started = time.perf_counter()
            likelihood = D.coarse_channel_nll(context, joined, 128, offset, offset+24576, target, .1,
                cuts=(offset, offset+24576), checkpoint_bp=256)
            raw_NLL = float(likelihood['loss'].detach()); forward_seconds = time.perf_counter()-forward_started
            if not math.isfinite(raw_NLL) or raw_NLL < -1e-8:
                raise ValueError('Actual unchanged full-grammar observation NLL must be finite and nonnegative')
            if any(likelihood[k] != counts[k] for k in ('fine_positions', 'coarse_positions', 'fully_unknown_positions')):
                raise ValueError('Actual structured channel observation row counts changed')
            backward_started = time.perf_counter(); likelihood['loss'].backward()
            torch.cuda.synchronize(); backward_seconds = time.perf_counter()-backward_started
            analytic = torch.zeros(14, 2, dtype=torch.float64)
            for logits, evidence in zip(tensors, priors):
                if logits.grad is None or not bool(torch.isfinite(logits.grad).all()):
                    raise ValueError('Finite actual local fifteen-role VJP required')
                analytic += logits.grad[:, 1:].double().cpu().T @ evidence[0].double().cpu()
            gradient = candidate.native_cds_role_prior.weight.grad
            if gradient is None or not bool(torch.isfinite(gradient).all()):
                raise ValueError('Both native-head coefficients must receive the actual structural derivative')
            actual = gradient.detach().double().cpu()
            torch.testing.assert_close(actual, analytic, rtol=CONFIG['VJP_rtol'], atol=CONFIG['VJP_atol'])
            row = dict(ordinal=ordinal, window_id=w['id'], species=w['species'], seqid=w['seqid'], strand=w['strand'],
                geometry=geometry, **counts, actual_prior_inputs=3, original_local_and_owned_eval_replay_exact=True,
                fresh_both_head_core_margin_replay_exact=True, training_eval_logit_max_abs_difference=train_eval_max,
                raw_NLL=raw_NLL, training_NLL_abs_difference_from_prior_eval_readout=abs(raw_NLL-old_NLL[w['id']]),
                NLL_forward_seconds=forward_seconds, NLL_backward_seconds=backward_seconds,
                native_prior_gradient=actual.tolist(), analytic_VJP_max_abs_difference=float((actual-analytic).abs().max()),
                prior_head_column_gradient_l2=actual.square().sum(0).sqrt().tolist(),
                shared_parameter_gradient_l2={name:math.sqrt(sum(float(p.grad.detach().double().square().sum())
                    for p in module.parameters() if p.grad is not None)) for name, module in (
                    ('project', candidate.shared.project), ('context', candidate.shared.context),
                    ('dna', candidate.shared.dna), ('fuse', candidate.shared.fuse), ('segmentation', candidate.segmentation))},
                row_seconds=time.perf_counter()-t)
            records.append(row); log.write(json.dumps(row, allow_nan=False)+'\n'); log.flush()
            print(json.dumps({k:row[k] for k in ('ordinal', 'window_id', 'raw_NLL',
                'prior_head_column_gradient_l2', 'analytic_VJP_max_abs_difference', 'row_seconds')}), flush=True)
            del features_list, dna_list, priors, tensors, joined, likelihood, eval_joined, eval_pieces, old_potentials, old_slots
        if (len(records) != 8 or any(not torch.equal(candidate.state_dict()[k].cpu(), v) for k, v in initial.items())
                or bool(candidate.native_cds_role_prior.weight.detach().count_nonzero())):
            raise ValueError('All original weights and zero prior must remain unchanged without optimizer updates')
    if any(p.requires_grad or p.grad is not None for p in full.parameters()):
        raise ValueError('Released backbone and both native heads must remain frozen without gradients')
    N.write(OUT/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        protocol=PROTOCOL, tests_run=tests.testsRun, failures=len(tests.failures), errors=len(tests.errors),
        actual_witnesses=8, records=records, actual_trainable_parameters=366219, new_prior_parameters=28,
        all_original_eval_predictions_replayed_exact=True, both_native_core_margins_replayed_exact=True,
        all_real_prior_VJPs_match_analytic=True,
        all_prior_head_columns_nonzero=all(all(v > 0 for v in r['prior_head_column_gradient_l2']) for r in records),
        original_C0_weights_and_zero_prior_unchanged=True, pretrained_encoder_and_heads_frozen_no_grad=True,
        optimizer_steps_this_run=0, checkpoint_created=False, DEV_used=False, sealed_data_used=False,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(), seconds=time.perf_counter()-began,
        interpretation='engineering_initialization_and_gradient_preflight_not_trained_gene_annotation_result'))


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
