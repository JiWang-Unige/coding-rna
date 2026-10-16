#!/usr/bin/env python3
"""One fixed TRAIN input, chain-only 2x2 precision diagnosis; no optimizer."""
import argparse
import copy
import json
import os
import time
from pathlib import Path
import torch
import feature_smoke as F
import manifest as M
from src.m28.core import M28Core, ChainHead, select_nonoverlapping
from src.m28.packed_chain_gru import packed_chain_forward
from src.m28.training import one_hot

ATOL, RTOL = 2e-6, 2e-4
CHECKPOINT = M.C.ROOT/'outputs/M28-PILOT-R4/B1/step_001536.pt'
R1 = M.C.ROOT/'outputs/M28-PACKED-GRU-GPU-R1/result/free_before_reference.jsonl'


def flags():
    return dict(cudnn_allow_tf32=torch.backends.cudnn.allow_tf32,
                matmul_allow_tf32=torch.backends.cuda.matmul.allow_tf32,
                cudnn_enabled=torch.backends.cudnn.enabled,
                cudnn_deterministic=torch.backends.cudnn.deterministic,
                cudnn_benchmark=torch.backends.cudnn.benchmark,
                deterministic_algorithms=torch.are_deterministic_algorithms_enabled(),
                float32_matmul_precision=torch.get_float32_matmul_precision(),
                autocast_cuda=torch.is_autocast_enabled())


def compare(new, old):
    result = {}
    for key in old:
        a, b = new[key].double(), old[key].double()
        delta = (a-b).abs()
        result[key] = dict(elements=b.numel(), max_abs=float(delta.max()),
                           max_relative_atol_floor=float((delta/b.abs().clamp_min(ATOL)).max()),
                           outside_original_tolerance=int((delta > ATOL+RTOL*b.abs()).sum()),
                           bitwise_equal=torch.equal(a, b))
    return result


def choices(scores, chains):
    gains = scores['logits']
    return dict(positive_indices=torch.where(gains > 0)[0].tolist(),
                selected_indices=select_nonoverlapping(chains, gains.tolist()),
                min_absolute_gain=float(gains.abs().min()))


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    environment = dict(job_id=os.environ.get('SLURM_JOB_ID'), GPU=torch.cuda.get_device_name(),
                       torch=str(torch.__version__), CUDA=torch.version.cuda,
                       cudnn=torch.backends.cudnn.version(), initial_flags=flags(),
                       NVIDIA_TF32_OVERRIDE=os.environ.get('NVIDIA_TF32_OVERRIDE'),
                       CUBLAS_WORKSPACE_CONFIG=os.environ.get('CUBLAS_WORKSPACE_CONFIG'))
    with (out/'environment.json').open('x') as f:
        json.dump(environment, f, indent=2); f.write('\n')
    # Changing these other flags would confound the intended single-flag contrast.
    assert torch.backends.cudnn.enabled and not torch.backends.cuda.matmul.allow_tf32
    torch.manual_seed(0); torch.set_num_threads(2)
    began = time.perf_counter()
    with R1.open() as f:
        original = json.loads(next(f))
    assert original['phase'] == 'warmup' and original['variant'] == 'baseline'
    ident = original['window_id']
    assert ident == 'arabidopsis_thaliana|NC_003070.9|29786112|-'
    chains = original['free']['chains']
    assert len(chains) == 192 and original['free']['reference_used'] is False
    with (M.C.ROOT/'outputs/M28-PILOT-R4/features/index.jsonl').open() as f:
        matches = [r for line in f if (r := json.loads(line))['window_id'] == ident]
    assert len(matches) == 1 and matches[0]['split'] == 'train'
    state = torch.load(CHECKPOINT, map_location='cpu', weights_only=True)
    assert state['arm'] == 'B1' and state['step'] == 1536 and state['revision'] == F.REVISION
    cached = torch.load(M.C.ROOT/matches[0]['artifact'], map_location='cpu', weights_only=True)
    w = cached['window']; n = w['valid_bases']
    assert w['id'] == ident and cached['revision'] == F.REVISION
    assert tuple(cached['features'].shape) == (4096, 2048) and len(cached['sequence']) == n
    model = M28Core('B1', 2048, 128).to('cuda').train()
    model.load_state_dict(state['model'], strict=True)
    # Produce h once under the explicitly fixed original-style setting. All four
    # chain-only cells receive exactly this same tensor, not recomputed contexts.
    torch.backends.cudnn.allow_tf32 = True
    outputs = model(cached['features'].float()[None].to('cuda'),
                    one_hot(cached['sequence'], w['width'])[None].to('cuda'), [n])
    h = outputs['features'][0, :n].detach().clone()
    del outputs, cached, state
    variants = dict(baseline=ChainHead.forward, packed=packed_chain_forward)
    snapshots = {}; rows = []
    with (out/'cells.jsonl').open('x') as log:
        for repeat in range(2):
            order = [(True, 'baseline'), (True, 'packed'), (False, 'baseline'), (False, 'packed')]
            if repeat:
                order.reverse()
            for tf32, variant in order:
                torch.backends.cudnn.allow_tf32 = tf32
                head = copy.deepcopy(model.chain)
                # Grad-enabled train-mode forward mirrors the failing chain call.
                result = variants[variant](head, h, chains)
                score = {k: v.detach().cpu().clone() for k, v in result.items()}
                if not all(torch.isfinite(v).all() for v in score.values()):
                    raise ValueError('Nonfinite chain scores')
                key = variant+'_'+str(tf32)+'_'+str(repeat)
                snapshots[key] = score
                row = dict(key=key, variant=variant, cudnn_allow_tf32=tf32, repeat=repeat,
                           flags=flags(), scores={k: v.tolist() for k, v in score.items()},
                           **choices(score, chains))
                rows.append(row); log.write(json.dumps(row)+'\n'); log.flush()
                del result, head
                print(json.dumps(dict(cell=key, complete=True)), flush=True)
    # Float64 is only a chain-head numerical reference on the same FP32-origin h;
    # not biological truth, full-model FP64, or a new trained checkpoint.
    reference_head = copy.deepcopy(model.chain).cpu().double()
    result = ChainHead.forward(reference_head, h.cpu().double(), chains)
    reference = {k: v.detach().clone() for k, v in result.items()}
    if not all(torch.isfinite(v).all() for v in reference.values()):
        raise ValueError('Nonfinite CPU float64 reference')
    with (out/'cpu_fp64_reference.json').open('x') as f:
        json.dump(dict(scores={k: v.tolist() for k, v in reference.items()},
                       **choices(reference, chains)), f); f.write('\n')
    pairs = {
        'within_original_setting': ('packed_True_0', 'baseline_True_0'),
        'within_tf32_disabled': ('packed_False_0', 'baseline_False_0'),
        'baseline_precision_shift': ('baseline_False_0', 'baseline_True_0'),
        'packed_precision_shift': ('packed_False_0', 'packed_True_0')}
    comparisons = {name: compare(snapshots[a], snapshots[b]) for name, (a, b) in pairs.items()}
    for variant in variants:
        for tf32 in (True, False):
            prefix = variant+'_'+str(tf32)
            comparisons['repeat_'+prefix] = compare(snapshots[prefix+'_1'], snapshots[prefix+'_0'])
            comparisons['vs_fp64_'+prefix] = compare(snapshots[prefix+'_0'], reference)
    selection = {}
    for name, (a, b) in pairs.items():
        sa, sb = snapshots[a]['logits'], snapshots[b]['logits']
        selection[name] = dict(positive_sign_flips=int(((sa > 0) != (sb > 0)).sum()),
                               selected_indices_equal=choices(snapshots[a], chains)['selected_indices'] ==
                                                      choices(snapshots[b], chains)['selected_indices'])
    report = dict(experiment='M28-CHAIN-PRECISION-D1', environment=environment,
                  checkpoint=str(CHECKPOINT.relative_to(M.C.ROOT)), checkpoint_step=1536,
                  window_id=ident, free_candidate_count=len(chains),
                  input_feature_source='One original-style shared forward; h fixed across all chain-only cells',
                  shared_input_flags='cudnn.allow_tf32=True, cuda.matmul.allow_tf32=False',
                  candidate_source=str(R1.relative_to(M.C.ROOT))+': first baseline warmup record',
                  new_reference_labels_read=False, optimizer_updates=0, new_weights_saved=False,
                  DEV_test_setaria_used=False, R4_modified_or_resumed=False,
                  atol=ATOL, rtol=RTOL, comparisons=comparisons, selection_comparisons=selection,
                  CPU_fp64_reference_scope='chain head only, same FP32-origin h and checkpoint values',
                  wall_loop_seconds=time.perf_counter()-began,
                  interpretation='Numerical diagnosis only; no equivalence admission, optimizer validation, throughput or accuracy claim')
    with (out/'summary.json').open('x') as f:
        json.dump(report, f, indent=2); f.write('\n')
    print(json.dumps(dict(comparisons=comparisons, selection=selection)), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, required=True)
    main(parser.parse_args().output_dir)
