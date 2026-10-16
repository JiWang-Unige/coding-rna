#!/usr/bin/env python3
"""Two fixed fresh TRAIN views; disposable joint24576 low-rank gradient preflight."""
import argparse
import json
import os
import re
import time
import torch
import fit_full_grammar_fresh_joint24 as B
from src.m28 import joint_encoder_qvo_lora as L

C, J, G, Q, ROOT = B.C, B.J, B.G, B.Q, B.ROOT
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_joint_qvo_lora_plan.md'
OUT = ROOT/'outputs/M28-FULL-GRAMMAR-JOINT-QVO-LORA-PREFLIGHT'
CHECKS = ROOT/'outputs/M28-FULL-GRAMMAR-JOINT-QVO-LORA-TRAIN-CHECKS'
CONFIG = dict(B.CONFIG, encoder_frozen=False, encoder_base_weights_frozen=True,
    model='source4608_C0_plus167936_joint_encoder_adapter_parameters_disposable',
    optimizer_steps=4609, additional_optimizer_steps=1, adapter_optimizer_steps=1,
    formal_checkpoint='none_disposable_preflight_no_checkpoint_written',
    encoder_adaptation='last_four_layers22_25_q_proj_v_proj_o_proj_only', adapter_rank=4,
    adapter_alpha=4, adapter_scale=1., adapter_seed=0, adapter_dtype='FP32', native_adapter_lr=3e-4,
    adapter_initialization_streams='qv_original_seed0_unchanged_o_separate_seed0',
    adapter_A='seed0_uniform_plus_minus_inverse_sqrt_input_dim', adapter_B='zero',
    adapter_parameters=167936, adapter_dropout=0., gradient_checkpointing=False,
    adapter_output='FP32_residual_rounded_to_BF16_added_to_BF16_projection',
    encoder_materialization='differentiable_joint24576_BF16_then_FP32_no_detach',
    optimizer='restored_head_AdamW_plus_new_adapter_group_same_hyperparameters',
    scope='two_fixed_first_fresh_TRAIN_draws_disposable_updates_no_saved_checkpoint',
    formal_training_started=False)


def potentials(model, backbone, tokenizer, chromosome, geometry, strand):
    pieces = []
    for left, right in geometry['genomic_windows']:
        sequence = re.sub('[^ACGT]', 'N', chromosome[left:right])
        if strand == '-':
            sequence = Q.M.rc(sequence)
        features = L.encode_joint_window(backbone, tokenizer, sequence, 'cuda')
        dna = G.one_hot(sequence, Q.M.WIDTH)[None].to('cuda')
        pieces.append(model(features, dna, [len(sequence)])['segmentation_logits'][0, :len(sequence)])
    return G.owned_context_logits(pieces, geometry, strand).cpu()


def main(checks_job):
    began = time.perf_counter()
    checked = CHECKS/checks_job
    summary = json.loads((checked/'summary.json').read_text())
    if ((checked/'STATUS').read_text().strip() != 'COMPLETED' or not summary['completed']
            or summary['job_id'] != checks_job or summary['tests_run'] != 4
            or summary['errors'] or summary['failures']):
        raise ValueError('Completed four new q/v/o interface and checkpoint checks required')
    plan, windows, catalogs = B.plan_and_population()
    selected = [next(i for i, row in enumerate(plan[:1536]) if row['species'] == species) for species in Q.ALLOW]
    C.write(OUT/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'], configuration=CONFIG,
        protocol=PROTOCOL, checks_job=checks_job, selected_additional_steps=[i+1 for i in selected],
        selected_rule='same_first_original_fresh_draw_per_species_as_B_preflight'))
    source = C.starting_state()
    sequences, descriptors = C.P.read_train_DNA()
    backbone, tokenizer = G.load_native()
    # Record the frozen parent before wrapping. Only the same fixed views are read.
    frozen = []
    parent_model, parent_optimizer = C.load_model_optimizer(source)
    with torch.no_grad():
        for ordinal in selected:
            expected, w = plan[ordinal], windows[plan[ordinal]['window_id']]
            chromosome, geometry, context, offset, targets = C.example(w, expected, sequences, catalogs)
            frozen.append(J.potentials(parent_model, backbone, tokenizer, chromosome, geometry, w['strand']))
    del parent_model, parent_optimizer
    adapters = L.attach_adapters(backbone)
    parameters = L.adapter_parameters(adapters)
    if sum(p.numel() for p in parameters) != L.PARAMETERS:
        raise ValueError('Exactly167936 new trainable adapter parameters required')
    adapter_ids = {id(p) for p in parameters}
    if {id(p) for p in backbone.parameters() if p.requires_grad} != adapter_ids:
        raise ValueError('Only the fixed twenty-four adapter tensors may be trainable')
    initial = L.adapter_state(adapters)
    records = []
    for parent_potentials, ordinal in zip(frozen, selected):
        started = time.perf_counter()
        L.restore_adapters(adapters, initial)
        for parameter in backbone.parameters():
            parameter.grad = None
        expected, w = plan[ordinal], windows[plan[ordinal]['window_id']]
        chromosome, geometry, context, offset, targets = C.example(w, expected, sequences, catalogs)
        model, optimizer = C.load_model_optimizer(source)
        optimizer.add_param_group({'params': parameters})
        if (optimizer.param_groups[0]['lr'] != CONFIG['lr']
                or optimizer.param_groups[1]['lr'] != CONFIG['native_adapter_lr']):
            raise ValueError('Actual disposable head/native learning rates must match the fixed q/v/o setting')
        torch.cuda.reset_peak_memory_stats()
        joined = potentials(model, backbone, tokenizer, chromosome, geometry, w['strand'])
        if not torch.equal(joined.detach(), parent_potentials):
            raise ValueError('Zero adapters changed actual native frozen-parent potentials')
        raw, _ = C.objective('STRUCTURED_CONTINUE', context, joined, offset, w['valid_bases'], targets)
        parent = parent_potentials.clone().requires_grad_()
        parent_raw, _ = C.objective('STRUCTURED_CONTINUE', context, parent, offset, w['valid_bases'], targets)
        if float(raw.detach()) != float(parent_raw.detach()):
            raise ValueError('Zero adapters changed actual native observation NLL')
        gradient, = torch.autograd.grad(raw, joined, retain_graph=True)
        parent_gradient, = torch.autograd.grad(parent_raw, parent)
        if not torch.equal(gradient, parent_gradient):
            raise ValueError('Zero adapters changed the actual role VJP')
        total = raw * expected['objective_multiplier']
        if not bool(torch.isfinite(total)) or float(raw.detach()) < -1e-8:
            raise ValueError('Finite nonnegative structured observation NLL required')
        total.backward()
        trainable = list(model.parameters()) + parameters
        if any(p.grad is None or not bool(torch.isfinite(p.grad).all()) for p in trainable):
            raise ValueError('Finite C0 and all adapter gradients required')
        if any(p.grad is not None or p.requires_grad for p in backbone.parameters() if id(p) not in adapter_ids):
            raise ValueError('Released encoder base weights must remain frozen and gradient-free')
        a_norm = math_norm([module.A.grad for module in adapters.values()])
        b_norms = {name: float(module.B.grad.double().norm()) for name, module in adapters.items()}
        if a_norm != 0 or not all(value > 0 for value in b_norms.values()):
            raise ValueError('Zero-B initialization requires zero A and nonzero B gradient in each projection')
        norm = torch.nn.utils.clip_grad_norm_(trainable, 1., error_if_nonfinite=True)
        optimizer.step()
        torch.cuda.synchronize()
        if (any(int(optimizer.state[p]['step']) != 4609 for p in model.parameters())
                or any(int(optimizer.state[p]['step']) != 1 for p in parameters)
                or any(not bool(torch.isfinite(p).all()) for p in trainable)):
            raise ValueError('Finite original head4609 and new adapter1 update required')
        updated = L.adapter_state(adapters)
        if not all(torch.count_nonzero(updated[name+'.B']) > 0 for name in adapters):
            raise ValueError('Each adapter B must receive a nonzero disposable update')
        if all(torch.equal(value.cpu(), source['model'][key]) for key, value in model.state_dict().items()):
            raise ValueError('Restored original head must receive its disposable update')
        restored = G.M28Core('C0', 2048, 128).eval().requires_grad_(False).to('cuda')
        restored.load_state_dict(model.state_dict(), strict=True)
        model.eval()
        with torch.no_grad():
            trained = potentials(model, backbone, tokenizer, chromosome, geometry, w['strand'])
        # Restore only the new FP32 state, not the released base weights.
        L.restore_adapters(adapters, initial)
        L.restore_adapters(adapters, updated)
        raw_view = J.raw_view(restored, backbone, tokenizer, chromosome, geometry, w['strand'])
        if not torch.equal(trained, raw_view):
            raise ValueError('Updated differentiable and strict raw native potentials differ')
        left = J.A.R.decode_full_grammar_skyline(context, trained, 128)
        right = J.A.R.decode_full_grammar_skyline(context, raw_view, 128)
        if any(left[k] != right[k] for k in ('chains', 'gain', 'chain_gains')):
            raise ValueError('Actual updated context training/raw MAP differs')
        txs, scopes = J.A.annotate_contig(restored, backbone, tokenizer, 'toy14',
            'ATGCGTNAGTTTAA', ordinal, 'cuda')
        if len(scopes) != 2 or any(s['windows'] != [(0, 14)] for s in scopes):
            raise ValueError('Actual short raw native caller must retain both strands and all DNA rows')
        row = dict(additional_step=ordinal+1, window_id=w['id'], species=w['species'], strand=w['strand'],
            native_context_bases=len(context), native_windows=len(geometry['genomic_windows']),
            raw_NLL=float(raw.detach()), weighted_total=float(total.detach()),
            original_known_observation_denominator=expected['known_observation_denominator'],
            gradient_norm_before_clip=float(norm), normalized_A_gradient_l2=a_norm,
            normalized_B_gradient_l2_by_projection=b_norms, updated_B_l2=math_norm([
                module.B for module in adapters.values()]), adapter_parameters=L.PARAMETERS,
            initial_head_optimizer_step=4608, final_disposable_head_optimizer_step=4609,
            initial_adapter_optimizer_step=0, final_adapter_optimizer_step=1,
            head_learning_rate=optimizer.param_groups[0]['lr'],
            native_adapter_learning_rate=optimizer.param_groups[1]['lr'],
            zero_adapter_native_potentials_NLL_and_role_VJP_parent_exact=True,
            all_twelve_projection_B_gradients_finite_nonzero=True,
            all_four_output_projection_B_updates_nonzero=True,
            actual_updated_native_logits_raw_exact=True, actual_updated_context_MAP_raw_exact=True,
            short_native_raw_both_strands_complete=True, short_complete_CDS_predictions=len(txs),
            peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
            peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved(), seconds=time.perf_counter()-started)
        records.append(row)
        print(json.dumps(row, allow_nan=False), flush=True)
        del model, optimizer, joined, raw, total, gradient, parent_gradient, parent_raw, parent
        del restored, trained, raw_view, left, right, txs, scopes, updated
    C.write(OUT/'summary.json', dict(completed=True, job_id=os.environ['SLURM_JOB_ID'],
        configuration=CONFIG, protocol=PROTOCOL, checks_job=checks_job, records=records,
        source_fit_job='13434603', support_job='13453408', original_TRAIN_DNA_sources=descriptors,
        disposable_updates_only=True, optimizer_steps_this_run=2, formal_training_started=False,
        checkpoints_written=False, saved_source_weights_and_optimizer_unchanged=True,
        pretrained_encoder_base_weights_frozen_no_grad=True, cached_features_used=False,
        DEV_used=False, sealed_data_used=False, seconds=time.perf_counter()-began,
        interpretation='native_gradient_and_memory_feasibility_not_chain_accuracy_or_novelty'))


def math_norm(values):
    return float(torch.stack([value.detach().double().square().sum() for value in values]).sum().sqrt())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checks-job', required=True)
    args = parser.parse_args()
    torch.set_num_threads(2)
    main(args.checks_job)
