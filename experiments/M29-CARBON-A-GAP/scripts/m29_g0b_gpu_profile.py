"""G0b (GPU part): per-step time of a real head-only training step of the fresh-joint recipe on the 4090, and a LoRA(E)-style
encoder forward+backward micro-benchmark.  <= K steps from the fixed fresh-draw plan, disposable AdamW state restored from
the original final4608 checkpoint; NOTHING is written except the timing json.  No Baobab/sealed data touched.
Usage: m29_g0b_gpu_profile.py steps K     |     m29_g0b_gpu_profile.py lora K
"""
import json
import os
import re
import sys
import time
from pathlib import Path
import torch
import m29_common as M
sys.path.insert(0, str(M.ROOT / 'pipeline/m28/scripts'))
import fit_full_grammar_objective_contrast as C    # noqa: E402
import fit_full_grammar_joint24 as J               # noqa: E402
import prepare_full_grammar_fresh_draws as PP      # noqa: E402

BAOBAB = '/home/users/j/jwang/coding-rna'
for _m in list(sys.modules.values()):
    for _k, _v in list(getattr(_m, '__dict__', {}).items()):
        if isinstance(_v, Path) and str(_v).startswith(BAOBAB):
            setattr(_m, _k, M.ROOT / Path(str(_v)).relative_to(BAOBAB))
G, A = C.G, J.A
OUT = M.M29 / 'G0b'
OUT.mkdir(parents=True, exist_ok=True)
T = {'enc': 0.0, 'enc_calls': 0}
_orig_encode = A.encode_window


def timed_encode(*a, **k):
    torch.cuda.synchronize(); t = time.perf_counter()
    r = _orig_encode(*a, **k)
    torch.cuda.synchronize(); T['enc'] += time.perf_counter() - t; T['enc_calls'] += 1
    return r


A.encode_window = timed_encode


def gpu_mem():
    return dict(alloc_gb=torch.cuda.max_memory_allocated() / 2**30, reserved_gb=torch.cuda.max_memory_reserved() / 2**30)


def steps(K):
    plan = [json.loads(x) for x in (M.OUTS / 'M28-FULL-GRAMMAR-FRESH-DRAWS-PREPARATION/training_plan.jsonl').read_text().splitlines()]
    windows, strata, catalogs, original = PP.load_population()
    sequences, _ = C.P.read_train_DNA()
    state = C.starting_state()
    model, optimizer = C.load_model_optimizer(state)
    backbone, tokenizer = G.load_native()
    torch.cuda.reset_peak_memory_stats()
    rows = []
    for ordinal in range(K):
        expected = plan[ordinal]
        w = windows[expected['window_id']]
        chromosome, geometry, context, offset, targets = C.example(w, expected, sequences, catalogs)
        optimizer.zero_grad(set_to_none=True)
        T['enc'] = 0.0; T['enc_calls'] = 0
        torch.cuda.synchronize(); t0 = time.perf_counter()
        joined = J.potentials(model, backbone, tokenizer, chromosome, geometry, w['strand'])
        torch.cuda.synchronize(); t_pot = time.perf_counter() - t0
        t1 = time.perf_counter()
        raw, _ = C.objective('STRUCTURED_CONTINUE', context, joined, offset, w['valid_bases'], targets)
        t_fwd = time.perf_counter() - t1
        t2 = time.perf_counter()
        total = raw * expected['objective_multiplier']
        total.backward()
        torch.cuda.synchronize(); t_bwd = time.perf_counter() - t2
        t3 = time.perf_counter()
        norm = torch.nn.utils.clip_grad_norm_(list(model.parameters()), 1.)
        optimizer.step(); torch.cuda.synchronize()
        t_opt = time.perf_counter() - t3
        row = dict(step=ordinal, species=w['species'], context_bp=len(context), windows_in_context=len(geometry['genomic_windows']),
                   encoder_calls=T['enc_calls'], encoder_forward_s=T['enc'], potentials_total_s=t_pot,
                   head_forward_and_transfer_s=t_pot - T['enc'], dp_forward_cpu_s=t_fwd, backward_s=t_bwd, optimizer_s=t_opt,
                   step_total_s=time.perf_counter() - t0, raw_NLL=float(raw.detach()), grad_norm=float(norm), loadavg=os.getloadavg()[0], **gpu_mem())
        rows.append(row)
        print(json.dumps(row), flush=True)
        del joined, raw, total
    (OUT / f'gpu_steps_{os.getpid()}.json').write_text(json.dumps(rows, indent=1))


def lora(K):
    """Encoder forward+backward with the LoRA(E) adapters (last-4-layer q/v rank 4) on the two joint24576 windows of a real context."""
    from src.m28 import joint_encoder_lora as LR
    plan = [json.loads(x) for x in (M.OUTS / 'M28-FULL-GRAMMAR-FRESH-DRAWS-PREPARATION/training_plan.jsonl').read_text().splitlines()]
    windows, strata, catalogs, original = PP.load_population()
    sequences, _ = C.P.read_train_DNA()
    backbone, tokenizer = G.load_native()
    adapters = LR.attach_adapters(backbone)
    params = LR.adapter_parameters(adapters)
    for p in params:
        p.requires_grad_(True)
    torch.cuda.reset_peak_memory_stats()
    rows = []
    for ordinal in range(K):
        expected = plan[ordinal]
        w = windows[expected['window_id']]
        chromosome, geometry, context, offset, targets = C.example(w, expected, sequences, catalogs)
        torch.cuda.synchronize(); t0 = time.perf_counter()
        loss = 0.
        n = 0
        for left, right in geometry['genomic_windows']:
            seq = re.sub('[^ACGT]', 'N', chromosome[left:right])
            if w['strand'] == '-':
                seq = C.Q.M.rc(seq)
            h = LR.encode_joint_window(backbone, tokenizer, seq, 'cuda')
            loss = loss + h.float().pow(2).mean()
            n += 1
        torch.cuda.synchronize(); t_fwd = time.perf_counter() - t0
        t1 = time.perf_counter()
        loss.backward()
        torch.cuda.synchronize(); t_bwd = time.perf_counter() - t1
        for p in params:
            p.grad = None
        rows.append(dict(step=ordinal, windows=n, encoder_forward_s=t_fwd, encoder_backward_s=t_bwd, **gpu_mem()))
        print(json.dumps(rows[-1]), flush=True)
    (OUT / f'gpu_lora_{os.getpid()}.json').write_text(json.dumps(rows, indent=1))


def enc(K):
    """Inference-path encoder benchmark: frozen BF16 GENERanno forward over 24,576-bp windows (real TRAIN DNA), batch 1/2/4/8."""
    import re as _re
    from torch.nn.attention import SDPBackend, sdpa_kernel
    backbone, tokenizer = G.load_native()
    seq = M.read_sequence('NC_003070.9')[5_000_000:5_000_000 + 24576 * 9]
    seq = _re.sub('[^ACGT]', 'N', seq)
    rows = []
    for bs in (1, 2, 4, 8):
        ids = tokenizer([seq[i * 24576:(i + 1) * 24576] for i in range(bs)], add_special_tokens=False, return_tensors='pt',
                        truncation=False)['input_ids'].cuda()
        torch.cuda.reset_peak_memory_stats()
        with torch.inference_mode(), sdpa_kernel(SDPBackend.FLASH_ATTENTION):
            backbone(input_ids=ids, attention_mask=torch.ones_like(ids)).last_hidden_state   # warm-up
            torch.cuda.synchronize(); t = time.perf_counter()
            for _ in range(K):
                h = backbone(input_ids=ids, attention_mask=torch.ones_like(ids)).last_hidden_state
            torch.cuda.synchronize(); dt = (time.perf_counter() - t) / K
        rows.append(dict(batch=bs, seconds_per_call=dt, seconds_per_window=dt / bs, kbp_per_s_encoder_only=bs * 24.576 / dt, **gpu_mem()))
        print(json.dumps(rows[-1]), flush=True)
    (OUT / f'gpu_enc_{os.getpid()}.json').write_text(json.dumps(rows, indent=1))


if __name__ == '__main__':
    mode, K = sys.argv[1], int(sys.argv[2])
    if mode == 'steps':
        steps(K)
    elif mode == 'lora':
        lora(K)
    elif mode == 'enc':
        enc(K)
