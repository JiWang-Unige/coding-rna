#!/usr/bin/env python3
"""R2: scoped chain TF32-off candidate vs original R4 baseline; discarded updates."""
import argparse
import copy
import json
import math
import os
import statistics
import time
from collections import Counter
from pathlib import Path
from types import MethodType
import torch
import feature_smoke as F
import manifest as M
from src.m28.core import M28Core, ChainHead, select_nonoverlapping
from src.m28.packed_chain_precision import packed_chain_forward_no_tf32
from src.m28.candidates import generate, Budget
from src.m28.labels import build_targets, CDS_COLUMNS
from src.m28.training import joint_losses, one_hot

ATOL, RTOL = 2e-6, 2e-4
CHECKPOINT = M.C.ROOT / 'outputs/M28-PILOT-R4/B1/step_001536.pt'


def rows(path):
    with path.open() as f:
        for line in f:
            if line.strip():
                yield json.loads(line)


def compare_tensors(key, new, old, errors):
    if new.shape != old.shape or new.dtype != old.dtype:
        raise ValueError('Tensor identity mismatch: '+key)
    delta = (new-old).abs()
    finite = bool(torch.isfinite(new).all() and torch.isfinite(old).all())
    close = torch.isclose(new, old, atol=ATOL, rtol=RTOL)
    item = {'shape': list(old.shape), 'dtype': str(old.dtype),
            'finite': finite, 'within_tolerance': finite and bool(close.all()),
            'outside_tolerance': int((~close).sum()),
            'max_abs': float(delta.max()) if delta.numel() else 0.,
            'max_relative_with_atol_floor': float((delta/old.abs().clamp_min(ATOL)).max()) if delta.numel() else 0.}
    if not item['within_tolerance']:
        item['bad_flat_indices_first20'] = torch.where(~close.reshape(-1))[0][:20].tolist()
        if old.numel() <= 512:
            item['baseline_values'] = old.tolist()
            item['candidate_values'] = new.tolist()
    errors[key] = item


def distribution(values):
    ordered = sorted(values)
    return {'n': len(values), 'median': statistics.median(values),
            'p90_nearest_rank': ordered[math.ceil(.9*len(ordered))-1],
            'min': ordered[0], 'max': ordered[-1]}


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    if (out/'steps.jsonl').exists():
        raise FileExistsError('Engineering output already exists')
    environment = {'experiment': 'M28-PACKED-GRU-GPU-R2', 'job_id': os.environ.get('SLURM_JOB_ID'),
                   'GPU': torch.cuda.get_device_name(), 'torch': str(torch.__version__),
                   'CUDA': torch.version.cuda, 'cudnn': torch.backends.cudnn.version(),
                   'cudnn_allow_tf32': torch.backends.cudnn.allow_tf32,
                   'matmul_allow_tf32': torch.backends.cuda.matmul.allow_tf32,
                   'cudnn_benchmark': torch.backends.cudnn.benchmark,
                   'cudnn_deterministic': torch.backends.cudnn.deterministic,
                   'autocast_cuda': torch.is_autocast_enabled(),
                   'candidate_precision_scope': 'only packed chain forward; previous flag restored before backward'}
    with (out/'environment.json').open('x') as f:
        json.dump(environment, f, indent=2); f.write('\n')
    assert environment['cudnn_allow_tf32'] and not environment['matmul_allow_tf32']
    draws = list(rows(M.C.ROOT/'outputs/M28-CORE-LABEL-SMOKE-R2/result/paired_training_draws.jsonl'))[:12]
    assert Counter(d['species'] for d in draws) == {'arabidopsis_thaliana': 6, 'oryza_sativa': 6}
    ids = {d['window_id'] for d in draws}
    index = {r['window_id']: r for r in rows(M.C.ROOT/'outputs/M28-PILOT-R4/features/index.jsonl')
             if r['window_id'] in ids}
    assert set(index) == ids and all(r['split'] == 'train' for r in index.values())
    base = M.C.ROOT/'outputs/M28-LABEL-POLICY-R4'
    windows, catalogs = {}, {}
    for species in ('arabidopsis_thaliana', 'oryza_sativa'):
        selected = [w for w in rows(base/species/'windows_train.jsonl') if w['id'] in ids]
        windows.update({w['id']: w for w in selected})
        needed = {ident for w in selected for ident in w['overlapping_transcript_ids']}
        catalogs[species] = {t['id']: t for t in rows(base/species/'transcripts.jsonl') if t['id'] in needed}
    assert set(windows) == ids
    state = torch.load(CHECKPOINT, map_location='cpu', weights_only=True)
    assert state['arm'] == 'B1' and state['step'] == 1536 and state['revision'] == F.REVISION
    torch.manual_seed(0)
    torch.set_num_threads(2)
    model = M28Core('B1', 2048, 128).to('cuda').train()
    optimizer = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=.01,
                                 betas=(.9, .999), eps=1e-8)
    variants = {'baseline': ChainHead.forward, 'packed': packed_chain_forward_no_tf32}
    equivalence, timing, diagnostics = [], [], []
    begin = time.perf_counter()
    with (out/'free_before_reference.jsonl').open('x') as free_log, (out/'steps.jsonl').open('x') as step_log:
        def step(variant, draw, phase, repeat=0, diagnostic=False):
            # Reset both weights and Adam state. No continuation across examples.
            optimizer.zero_grad(set_to_none=True)
            model.load_state_dict(state['model'])
            optimizer.load_state_dict(copy.deepcopy(state['optimizer']))
            component = {}
            forward = variants[variant]
            if diagnostic:
                def measured(head, *args, **kwargs):
                    torch.cuda.synchronize()
                    start = time.perf_counter()
                    result = forward(head, *args, **kwargs)
                    torch.cuda.synchronize()
                    component['instrumented_chain_forward_seconds'] = time.perf_counter()-start
                    return result
                model.chain.forward = MethodType(measured, model.chain)
            else:
                model.chain.forward = MethodType(forward, model.chain)
            torch.manual_seed(0)
            torch.cuda.synchronize()
            torch.cuda.reset_peak_memory_stats()
            started = time.perf_counter()
            ident = draw['window_id']
            w = windows[ident]
            cached = torch.load(M.C.ROOT/index[ident]['artifact'], map_location='cpu', weights_only=True)
            assert cached['revision'] == F.REVISION
            assert all(cached['window'][k] == w[k] for k in ('id', 'start', 'end', 'strand', 'valid_bases'))
            n, sequence = w['valid_bases'], cached['sequence']
            assert len(sequence) == n and tuple(cached['features'].shape) == (4096, 2048)
            features = cached['features'].float()[None].to('cuda')
            dna = one_hot(sequence, w['width'])[None].to('cuda')
            torch.cuda.synchronize()
            load_seconds = time.perf_counter()-started
            forward_start = time.perf_counter()
            outputs = model(features, dna, [n])
            if diagnostic:
                outputs['features'].retain_grad()
            torch.cuda.synchronize()
            shared_seconds = time.perf_counter()-forward_start
            proposal_start = time.perf_counter()
            with torch.no_grad():
                h = outputs['features'][0, :n].detach()
                p = outputs['segmentation_logits'][0, :n].softmax(-1)
                free = generate(sequence, outputs['endpoint_logits'][0, :n].detach().cpu().numpy(),
                                p[:, list(CDS_COLUMNS)].sum(-1).cpu().numpy(),
                                link_score_fn=lambda pairs: model.chain.link_logits(h, pairs).cpu().numpy(),
                                budget=Budget())
            torch.cuda.synchronize()
            proposal_seconds = time.perf_counter()-proposal_start
            free_log.write(json.dumps({'phase': phase, 'variant': variant, 'repeat': repeat,
                                       'window_id': ident, 'free': free})+'\n')
            free_log.flush()
            targets = build_targets(w, catalogs[w['species']], sequence)
            joint = joint_losses(model, outputs, free, w, catalogs[w['species']], targets,
                                 window_weight=draw['importance_weight_uniform_within_stratum'])
            assert torch.backends.cudnn.allow_tf32 and not torch.backends.cuda.matmul.allow_tf32
            loss = joint['total']
            if not bool(torch.isfinite(loss)):
                raise ValueError('Nonfinite diagnostic loss')
            loss.backward()
            snapshot = None
            if diagnostic:
                snapshot = {'free': free, 'pool': joint['training_chains'], 'labels': joint['chain_labels'],
                            'scores': {k: v.detach().cpu().clone() for k, v in joint['chain_scores'].items()},
                            'losses': {k: v.detach().cpu().clone() for k, v in joint['losses'].items()},
                            'total': loss.detach().cpu().clone(),
                            'shared_features': outputs['features'].detach().cpu().clone(),
                            'feature_gradient': outputs['features'].grad.detach().cpu().clone(),
                            'parameter_gradients': {k: None if p.grad is None else p.grad.detach().cpu().clone()
                                                    for k, p in model.named_parameters()}}
            norm = torch.nn.utils.clip_grad_norm_(model.parameters(), 1., error_if_nonfinite=True)
            optimizer.step()
            assert torch.backends.cudnn.allow_tf32 and not torch.backends.cuda.matmul.allow_tf32
            torch.cuda.synchronize()
            values = {k: float(v.detach()) for k, v in joint['losses'].items()}
            record = {'phase': phase, 'repeat': repeat, 'variant': variant, 'window_id': ident,
                      'losses': values, 'weighted_total': float(loss.detach()),
                      'gradient_norm_before_clip': float(norm), 'free_candidates': len(free['chains']),
                      'training_candidates': len(joint['training_chains']), 'origins': joint['origins']}
            step_log.write(json.dumps(record)+'\n')
            step_log.flush()
            # Diagnostic transfers/synchronizations are deliberately excluded
            # from performance conclusions by using only phase=timing records.
            elapsed = time.perf_counter()-started
            metrics = {'variant': variant, 'window_id': ident, 'repeat': repeat,
                       'full_step_seconds': elapsed, 'load_seconds': load_seconds,
                       'shared_forward_seconds': shared_seconds, 'proposal_seconds': proposal_seconds,
                       'free_candidates': len(free['chains']), 'training_candidates': len(joint['training_chains']),
                       'peak_GPU_allocated_bytes': torch.cuda.max_memory_allocated(),
                       'peak_GPU_reserved_bytes': torch.cuda.max_memory_reserved(), **component}
            if diagnostic:
                snapshot['norm'] = norm.detach().cpu().clone()
                snapshot['increments'] = {k: p.detach().cpu()-state['model'][k] for k, p in model.named_parameters()}
                snapshot['optimizer'] = {k: {s: v.detach().cpu().clone() for s, v in optimizer.state[p].items()}
                                         for k, p in model.named_parameters()}
                keys = {tuple(map(tuple, chain)): i for i, chain in enumerate(joint['training_chains'])}
                free_scores = snapshot['scores']['logits'][
                    [keys[tuple(map(tuple, chain))] for chain in free['chains']]]
                snapshot['free_positive'] = free_scores > 0
                snapshot['free_selected'] = select_nonoverlapping(free['chains'], free_scores.tolist())
            return metrics, snapshot

        # Both implementations are initialized before numerical comparisons.
        for variant in variants:
            step(variant, draws[0], 'warmup')
        with (out/'equivalence.jsonl').open('x') as eq_log:
            for i, draw in enumerate(draws):
                results = {}
                for variant in (('baseline', 'packed') if i % 2 == 0 else ('packed', 'baseline')):
                    metrics, snapshot = step(variant, draw, 'equivalence', diagnostic=True)
                    diagnostics.append(metrics)
                    results[variant] = snapshot
                old, new = results['baseline'], results['packed']
                assert old['free'] == new['free'] and old['pool'] == new['pool'] and old['labels'] == new['labels']
                errors = {}
                for field in ('scores', 'losses', 'parameter_gradients', 'increments'):
                    for key in old[field]:
                        a, b = new[field][key], old[field][key]
                        assert (a is None) == (b is None), (field, key)
                        if a is not None:
                            compare_tensors(field+'.'+key, a, b, errors)
                for field in ('total', 'norm', 'feature_gradient'):
                    compare_tensors(field, new[field], old[field], errors)
                for key in old['optimizer']:
                    assert old['optimizer'][key].keys() == new['optimizer'][key].keys()
                    for subkey in old['optimizer'][key]:
                        compare_tensors('optimizer.'+key+'.'+subkey, new['optimizer'][key][subkey],
                                        old['optimizer'][key][subkey], errors)
                shared_equal = torch.equal(new['shared_features'], old['shared_features'])
                result = {'window_id': draw['window_id'], 'shared_features_bitwise_equal': shared_equal, 'atol': ATOL, 'rtol': RTOL,
                          'free_candidates': len(old['free']['chains']), 'errors': errors,
                          'positive_sign_flips': int((old['free_positive'] != new['free_positive']).sum()),
                          'selected_indices_equal': old['free_selected'] == new['free_selected'],
                          'all_tensors_within_tolerance': all(v['within_tolerance'] for v in errors.values()),
                          'failed_tensor_keys': [k for k, v in errors.items() if not v['within_tolerance']]}
                result['accepted'] = (shared_equal and result['all_tensors_within_tolerance'] and
                                      result['positive_sign_flips'] == 0 and result['selected_indices_equal'])
                equivalence.append(result)
                eq_log.write(json.dumps(result)+'\n')
                eq_log.flush()
                if not result['accepted']:
                    raise AssertionError('R2 rejected on '+draw['window_id']+': '+json.dumps(result['failed_tensor_keys']))
                print(json.dumps({'equivalence_completed': i+1, 'window_id': draw['window_id']}), flush=True)
                del results, old, new, snapshot

        # Same GPU, alternating AB/BA per example and repeat. Inputs are warm.
        for repeat in range(2):
            for i, draw in enumerate(draws):
                for variant in (('baseline', 'packed') if (i+repeat) % 2 == 0 else ('packed', 'baseline')):
                    metrics, _ = step(variant, draw, 'timing', repeat)
                    timing.append(metrics)
            print(json.dumps({'timing_repeat_completed': repeat+1}), flush=True)
    by_variant = {v: [x['full_step_seconds'] for x in timing if x['variant'] == v] for v in variants}
    paired = {(x['repeat'], x['window_id'], x['variant']): x['full_step_seconds'] for x in timing}
    ratios = [paired[r, d['window_id'], 'baseline']/paired[r, d['window_id'], 'packed']
              for r in range(2) for d in draws]
    report = {'experiment': 'M28-PACKED-GRU-GPU-R2', 'environment': environment,
              'job_id': os.environ.get('SLURM_JOB_ID'), 'checkpoint': str(CHECKPOINT.relative_to(M.C.ROOT)),
              'checkpoint_step': 1536, 'feature_revision': F.REVISION, 'draws': draws,
              'GPU': torch.cuda.get_device_name(), 'torch': str(torch.__version__), 'CUDA': torch.version.cuda,
              'cudnn': torch.backends.cudnn.version(), 'matmul_allow_tf32': torch.backends.cuda.matmul.allow_tf32,
              'cudnn_allow_tf32': torch.backends.cudnn.allow_tf32,
              'equivalence': equivalence, 'instrumented_forward_diagnostics': diagnostics,
              'timed_steps': timing, 'full_step_seconds': {k: distribution(v) for k, v in by_variant.items()},
              'paired_baseline_over_packed': distribution(ratios),
              'all_free_selections_equal': all(x['selected_indices_equal'] and x['positive_sign_flips'] == 0
                                               for x in equivalence),
              'repeat_paired_median_ratio': {str(r): statistics.median(
                  paired[r, d['window_id'], 'baseline']/paired[r, d['window_id'], 'packed'] for d in draws)
                  for r in range(2)},
              'loop_seconds': time.perf_counter()-begin, 'discarded_one_step_updates': 74,
              'saved_model_or_optimizer': False, 'R4_modified_or_resumed': False,
              'DEV_test_setaria_used': False, 'timing_cache_condition': 'warm, paired alternating AB/BA',
              'timing_excludes': 'model/optimizer reset and initialization; full step includes load, proposal, free JSON flush, supervision, backward, optimizer and training trace write',
              'interpretation': 'Engineering equivalence and bounded timing only, not accuracy or complete-fit evidence'}
    with (out/'summary.json').open('x') as f:
        json.dump(report, f, indent=2)
        f.write('\n')
    print(json.dumps({k: report[k] for k in ('full_step_seconds', 'paired_baseline_over_packed',
                                            'all_free_selections_equal', 'loop_seconds')}, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, required=True)
    main(parser.parse_args().output_dir)
