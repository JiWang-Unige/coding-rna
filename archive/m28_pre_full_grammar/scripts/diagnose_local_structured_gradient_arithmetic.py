#!/usr/bin/env python3
"""Artificial CUDA arithmetic only; no research tensors or model fitting."""
import json
import time
import torch

import preflight_local_structured_emission as P
from src.m28.local_structured_emission import local_structured_emission_loss

OUT = P.ROOT/'outputs/M28-A44-GRADIENT-ARITHMETIC'


def probe(features, dna, mode):
    model = P.A.A36.new_model('cuda').train()
    before = {k: v.detach().clone() for k, v in model.state_dict().items()}
    params = list(model.parameters()); n = dna.shape[1]
    logits = model(features, dna, [n])['segmentation_logits'][0]
    targets = dict(seg=torch.zeros(n, dtype=torch.long, device='cuda'),
                   region=torch.zeros(n, dtype=torch.long, device='cuda'))
    result = local_structured_emission_loss(logits,
        [[(256, 12256)], [(256, 15256)]], [1, 0], targets, n,
        local_coefficient=1., structured_coefficient=1.)
    keys = ('local_loss', 'structured_loss', 'loss')
    emission = {k: torch.autograd.grad(result[k], logits, retain_graph=True)[0] for k in keys}
    expected = emission['local_loss']+emission['structured_loss']
    emission_ok = bool(torch.isclose(emission['loss'], expected, rtol=1e-4, atol=1e-6).all())
    gradients = {k: torch.autograd.grad(result[k], params, retain_graph=k != 'loss') for k in keys}
    per_parameter = []; squared_error = 0.; squared_total = 0.
    for (name, _), a, b, total in zip(model.named_parameters(), gradients['local_loss'],
                                      gradients['structured_loss'], gradients['loss']):
        error = total-(a+b)
        if not bool(torch.isfinite(total).all()) or not bool(torch.isfinite(error).all()):
            raise ValueError('Nonfinite artificial CUDA arithmetic')
        mismatch = int((~torch.isclose(total, a+b, rtol=1e-4, atol=1e-6)).sum())
        squared_error += float(error.double().square().sum())
        squared_total += float(total.double().square().sum())
        per_parameter.append(dict(parameter=name, elements=total.numel(),
            mismatched_at_original_guard=mismatch, maximum_absolute_error=float(error.abs().max())))
    if any(not torch.equal(v, before[k]) for k, v in model.state_dict().items()):
        raise ValueError('Artificial zero-update arithmetic changed model state')
    return dict(mode=mode, cuda_matmul_allow_tf32=torch.backends.cuda.matmul.allow_tf32,
        cudnn_allow_tf32=torch.backends.cudnn.allow_tf32,
        matmul_precision=torch.get_float32_matmul_precision(),
        emission_additivity_passes_original_guard=emission_ok,
        maximum_emission_additivity_error=float((emission['loss']-expected).abs().max()),
        upstream_parameters_pass_original_guard=not any(r['mismatched_at_original_guard'] for r in per_parameter),
        upstream_relative_L2_additivity_error=(squared_error/squared_total)**.5,
        per_parameter=per_parameter, model_weights_unchanged=True,
        optimizer_steps=0, parameters_with_accumulated_grad=sum(p.grad is not None for p in params))


def main():
    torch.manual_seed(0); torch.cuda.reset_peak_memory_stats(); began = time.perf_counter()
    features = torch.randn((1, 4096, 2048), device='cuda', dtype=torch.float32)
    dna = torch.nn.functional.one_hot(torch.arange(24576, device='cuda') % 5, 5).float()[None]
    original = (torch.backends.cuda.matmul.allow_tf32, torch.backends.cudnn.allow_tf32)
    try:
        rows = [probe(features, dna, 'unchanged_imported_runtime')]
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
        rows.append(probe(features, dna, 'artificial_TF32_disabled_diagnostic_only'))
    finally:
        torch.backends.cuda.matmul.allow_tf32, torch.backends.cudnn.allow_tf32 = original
    report = dict(completed=True, experiment='M28-A44-GRADIENT-ARITHMETIC',
        protocol=P.PROTOCOL, failed_real_preflight_job='13425634',
        original_additivity_guard=dict(rtol=1e-4, atol=1e-6),
        torch_version=torch.__version__, cuda_device=torch.cuda.get_device_name(), records=rows,
        original_runtime_flags_restored=True, optimizer_steps=0,
        research_inputs_opened=False, checkpoints_written=False, DEV_used=False,
        seconds=time.perf_counter()-began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='Synthetic FP32 arithmetic diagnosis only; no real-data retry, tolerance change, '
            'research kernel change, selected fit coefficients or model performance evidence.')
    with (OUT/'summary.json').open('x') as handle:
        json.dump(report, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(report, allow_nan=False), flush=True)


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
