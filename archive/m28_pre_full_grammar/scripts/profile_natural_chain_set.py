#!/usr/bin/env python3
"""A15 first-32-draw zero-update timing after a verified scheduler timeout."""
from collections import defaultdict
import time
import torch
import run_natural_chain_set as R

OUT = R.ROOT/'outputs/M28-A15-NATURAL-CHAIN-SET-TIMING'


def sync():
    torch.cuda.synchronize()
    return time.perf_counter()


def main():
    assert (R.FIT/'JOBID').read_text().strip() == '13367306'
    assert not (R.FIT/'step_001536.pt').exists() and not (R.FIT/'summary.json').exists()
    began = time.perf_counter()
    draws, index, windows, records, _ = R.population()
    population_seconds = time.perf_counter()-began
    generator = R.C.load_generator(); head = R.C.new_quality('cuda').train()
    initial = {k: v.detach().clone() for k, v in head.state_dict().items()}
    timings = []; totals = defaultdict(float)
    for ordinal, draw in enumerate(draws[:32]):
        ident = draw['window_id']; w = windows[ident]; local = records[ident]['candidates']; row = {}
        head.zero_grad(set_to_none=True); t = sync()
        cached = torch.load(R.ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
        row['cache_read_seconds'] = sync()-t; t = sync()
        features, dna = R.C.R.Q.window_inputs(cached, w, 'cuda')
        row['input_conversion_transfer_seconds'] = sync()-t; t = sync()
        with torch.no_grad(): hidden = generator(features, dna, [w['valid_bases']])['features'][0]
        row['frozen_generator_seconds'] = sync()-t; t = sync()
        score = head(hidden, [r['chain'] for r in local], w['valid_bases'])
        row['chain_head_seconds'] = sync()-t; t = sync()
        result = R.chain_set_loss(score, [r['chain'] for r in local], [r['label'] for r in local])
        loss = result['loss']*draw['importance_weight_uniform_within_stratum']
        if not bool(torch.isfinite(loss)):
            raise ValueError('Nonfinite unchanged natural-domain timing loss')
        row['set_loss_seconds'] = sync()-t; t = sync()
        loss.backward(); norm = torch.nn.utils.clip_grad_norm_(head.parameters(), 1., error_if_nonfinite=True)
        row['backward_clip_seconds'] = sync()-t
        assert not hidden.requires_grad and not features.requires_grad and not dna.requires_grad
        assert not any(p.grad is not None for p in generator.parameters())
        totals.update({k: totals[k]+v for k, v in row.items()})
        timings.append(dict(draw_ordinal=ordinal, window_id=ident, timings=row, hidden_shape=list(hidden.shape),
                            hidden_dtype=str(hidden.dtype), hidden_bytes=hidden.numel()*hidden.element_size(),
                            gradient_norm_before_clip=float(norm)))
        print(dict(draw_ordinal=ordinal, timings=row), flush=True)
        del cached, features, dna, hidden, score, result, loss
    assert all(torch.equal(initial[k], v) for k, v in head.state_dict().items())
    R.C.R.write_summary(OUT, dict(completed=True, configuration=R.CONFIG, protocol=R.PROTOCOL,
        original_failed_fit_job='13367306', first_original_draws=32, timings=timings,
        population_seconds=population_seconds, mean_stage_seconds={k: v/32 for k, v in totals.items()},
        predicted_FP32_all_unique_hidden_bytes=sum(w['width']*128*4 for w in windows.values()),
        assigned_GPU_total_bytes=torch.cuda.get_device_properties(0).total_memory,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(), optimizer_steps_this_run=0,
        head_weights_unchanged=True, generator_weights_frozen=True, new_backbone_forward=False,
        new_candidate_search=False, DEV_forward_used=False, test_setaria_used=False,
        seconds=time.perf_counter()-began,
        limitation='Current first-32-draw timing, not proof of historical timeout cause or model performance.'))


if __name__ == '__main__':
    torch.set_num_threads(2); main()
