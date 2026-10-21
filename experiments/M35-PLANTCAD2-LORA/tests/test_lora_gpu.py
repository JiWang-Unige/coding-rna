"""GPU tests of the PlantCAD2 LoRA on the REAL PlantCAD2-Small (fused Triton Mamba2 kernels).  Short (a few minutes), peak memory a few GB.
Run only when the GPU is free (queue rule).  Writes ../results/lora_gpu_stats.json."""
import json
import time
from pathlib import Path

import numpy as np
import pytest
import torch

import pc2_adapter as A
import pc2_lora as L

pytestmark = pytest.mark.skipif(not torch.cuda.is_available(), reason='needs a GPU')
STATS = {}
RES = Path(__file__).resolve().parents[1] / 'results'


def rand_seq(n, seed):
    return ''.join(np.random.default_rng(seed).choice(list('ACGT'), n))


def rc(s):
    return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def cos(a, b):
    a, b = a.float().flatten(), b.float().flatten()
    return float(torch.dot(a, b) / (a.norm() * b.norm()))


@pytest.fixture(scope='module')
def env():
    enc = A.load_pc2()
    seqs = [rand_seq(2048, 1), rand_seq(2048, 2)]
    ref = A.encode_windows(enc, seqs, -1, batch=2)
    params = L.install_lora(enc, L.LORA_LAYERS, L.RANK, L.ALPHA)
    yield dict(enc=enc, seqs=seqs, ref=ref, params=params)
    RES.mkdir(exist_ok=True)
    (RES / 'lora_gpu_stats.json').write_text(json.dumps(STATS, indent=1))


def randomise_B(params, scale, seed=7):
    g = torch.Generator().manual_seed(seed)
    with torch.no_grad():
        for p in params[1::2]:                       # B matrices (A, B, A, B, ...)
            p.copy_((torch.randn(p.shape, generator=g) * scale).to(p.device))


def zero_B(params):
    with torch.no_grad():
        for p in params[1::2]:
            p.zero_()


def test_00_ties_and_param_count(env):
    n = L.tie_check_and_retie(env['enc'])
    assert n == 2 * L.N_LAYERS
    assert sum(p.numel() for p in env['params']) == 604_416 and len(env['params']) == 4 * len(L.LORA_LAYERS)
    STATS['lora_parameters'] = 604_416


def test_01_zero_lora_equals_frozen(env):
    out = A.encode_windows(env['enc'], env['seqs'], -1, batch=2)
    for a, b in zip(out, env['ref']):
        d = float((a.float() - b.float()).abs().max())
        STATS.setdefault('identity_at_init', []).append(dict(max_abs_diff=d, cosine=cos(a, b)))
        assert cos(a, b) > 0.99999 and d < 0.1


def test_02_rc_equivariance_with_lora(env):
    enc, params = env['enc'], env['params']
    randomise_B(params, .3)                          # large enough to move the features visibly (A ~ U(+-1/sqrt(768)), scale 2)
    s = env['seqs'][0]
    f = A.encode_windows(enc, [s], -1)[0]
    g = A.encode_windows(enc, [rc(s)], -1)[0]
    d = float((g.float() - f.flip(0).flip(-1).float()).abs().max())
    STATS['rc_equivariance_with_lora'] = dict(max_abs_diff=d, cosine=cos(g, f.flip(0).flip(-1)))
    base = env['ref'][0]
    STATS['lora_changes_features_cosine_to_frozen'] = cos(f, base)
    zero_B(params)
    assert STATS['rc_equivariance_with_lora']['cosine'] > 0.99999 and d < 0.1
    assert STATS['lora_changes_features_cosine_to_frozen'] < 0.9999                       # the randomised LoRA really acts


def test_03_encode_grad_matches_inference_features(env):
    enc, params = env['enc'], env['params']
    randomise_B(params, .02)
    inf = A.encode_windows(enc, env['seqs'], -1, batch=2)
    gr = L.encode_grad(enc, env['seqs'])
    for a, b in zip(inf, gr):
        assert cos(a, b.detach()) > 0.99999
    zero_B(params)


def grads(enc, params, seqs, scale_w=1.):
    for p in params:
        p.grad = None
    h = L.encode_grad(enc, seqs)
    w = torch.Generator(device='cpu').manual_seed(3)
    loss = sum((x.float() * torch.randn(x.shape, generator=w).to(x.device)).sum() for x in h) * 1e-3
    loss.backward()
    return [p.grad.clone() for p in params], float(loss)


def test_04_all_lora_params_get_grad_and_checkpoint_equals_plain(env):
    enc, params, seqs = env['enc'], env['params'], env['seqs']
    randomise_B(params, .02)
    g0, l0 = grads(enc, params, seqs)
    assert all(g is not None and bool(torch.isfinite(g).all()) for g in g0)
    n_zero = sum(1 for g in g0 if float(g.abs().sum()) == 0)
    L.enable_layer_checkpoint(enc, L.LORA_LAYERS)
    g1, l1 = grads(enc, params, seqs)
    cs = [cos(a, b) for a, b in zip(g0, g1)]
    STATS['checkpoint_vs_plain_grad'] = dict(min_cosine=min(cs), loss_plain=l0, loss_ckpt=l1, n_zero_grad_tensors=n_zero)
    assert min(cs) > 0.999 and n_zero == 0
    others = [n for n, p in enc.named_parameters() if p.grad is not None and id(p) not in {id(q) for q in params}]
    assert not others
    zero_B(params)


def test_05_directional_derivative_fp32(env):
    """Gradient correctness of the explicit adapter path through the fused kernels (fp32, no autocast, short window)."""
    enc, params = env['enc'], env['params']
    randomise_B(params, .01)
    s = [env['seqs'][0][:512]]
    ids = torch.stack([A.token_ids(x) for x in s])
    w = torch.randn(1, 512, A.HIDDEN, generator=torch.Generator().manual_seed(4)).cuda()

    def loss():
        return (enc(ids).last_hidden_state * w).sum() * 1e-2
    for p in params:
        p.grad = None
    loss().backward()
    gs = [p.grad.clone() for p in params]
    d = [torch.randn(p.shape, generator=torch.Generator().manual_seed(100 + i)).cuda() * (p.abs().mean()) for i, p in enumerate(params)]
    ana = sum(float((g * x).sum()) for g, x in zip(gs, d))
    eps = .05
    with torch.no_grad():
        for p, x in zip(params, d):
            p.add_(x * eps)
        lp = float(loss())
        for p, x in zip(params, d):
            p.sub_(x * 2 * eps)
        lm = float(loss())
        for p, x in zip(params, d):
            p.add_(x * eps)
    num = (lp - lm) / (2 * eps)
    STATS['directional_derivative_fp32'] = dict(analytic=ana, numeric=num, rel_err=abs(ana - num) / max(abs(num), 1e-12))
    zero_B(params)
    assert STATS['directional_derivative_fp32']['rel_err'] < .15


def test_06_state_roundtrip_and_inference_path(env):
    enc, params = env['enc'], env['params']
    randomise_B(params, .02, seed=11)
    st = L.lora_state(enc)
    a = A.encode_windows(enc, env['seqs'], -1, batch=2)
    enc2 = A.load_pc2()
    L.load_lora_state(enc2, st)
    b = A.encode_windows(enc2, env['seqs'], -1, batch=2)
    for x, y in zip(a, b):
        assert torch.equal(x, y)
    del enc2
    zero_B(params)


def test_07_memory_and_time_two_8k_windows(env):
    enc, params = env['enc'], env['params']
    randomise_B(params, .02)
    seqs = [rand_seq(A.PC_WIDTH, 5), rand_seq(A.PC_WIDTH, 6)]
    L.enable_layer_checkpoint(enc, L.LORA_LAYERS)
    torch.cuda.synchronize(); torch.cuda.reset_peak_memory_stats()
    t0 = time.perf_counter()
    h = L.encode_grad(enc, seqs)
    (h[0].float().pow(2).mean() + h[1].float().pow(2).mean()).backward()
    torch.cuda.synchronize()
    STATS['two_8192bp_windows_fwd_bwd'] = dict(seconds=time.perf_counter() - t0, peak_alloc_GB=torch.cuda.max_memory_allocated() / 2**30)
    t0 = time.perf_counter()
    with torch.no_grad():
        A.encode_windows(enc, seqs + seqs, -1, batch=4)
    torch.cuda.synchronize()
    STATS['four_8192bp_windows_forward_only_s'] = time.perf_counter() - t0
    zero_B(params)
    assert STATS['two_8192bp_windows_fwd_bwd']['peak_alloc_GB'] < 12
