"""CPU tests of the PlantCAD2 LoRA (pc2_lora) on a TOY Caduceus (real HF class, real RCPS/BiMamba/Mamba2 modules, tiny sizes, pure-torch scan injected
because the fused Triton kernel needs a GPU).  They check the STRUCTURE the real model relies on: re-tying, shared LoRA pair for mamba_fwd/mamba_rev,
identity at initialisation, exact reverse-complement equivariance with non-zero LoRA, gradient flow only to the LoRA, checkpointed == plain gradients, state
round-trip.  Equivalence with the real fused kernels is checked on the GPU by test_lora_gpu.py."""
import copy

import pytest
import torch
import torch.nn.functional as F

import pc2_adapter as A
import pc2_lora as L

torch.set_num_threads(2)
D_MODEL, NL = 64, 4


def toy_scan(zxbcdt, conv1d_weight, conv1d_bias, dt_bias, A_, D, chunk_size=None, activation='silu', rmsnorm_weight=None, rmsnorm_eps=1e-6,
             outproj_weight=None, outproj_bias=None, headdim=None, ngroups=1, norm_before_gate=True, **kw):
    """Pure-torch selective scan with the interface of mamba_split_conv1d_scan_combined (no out-projection used here)."""
    assert outproj_weight is None and ngroups == 1 and D.dim() == 1
    b, l, _ = zxbcdt.shape
    nh = A_.shape[0]
    dim = nh * headdim
    ds = (zxbcdt.shape[-1] - 2 * dim - nh) // 2
    z, xbc, dt = torch.split(zxbcdt, [dim, dim + 2 * ds, nh], -1)
    k = conv1d_weight.shape[1]
    xbc = F.silu(F.conv1d(xbc.transpose(1, 2), conv1d_weight[:, None, :], conv1d_bias, padding=k - 1, groups=conv1d_weight.shape[0])[..., :l]).transpose(1, 2)
    x, Bm, Cm = torch.split(xbc, [dim, ds, ds], -1)
    x = x.reshape(b, l, nh, headdim)
    dtp = F.softplus(dt + dt_bias)
    h = zxbcdt.new_zeros(b, nh, headdim, ds)
    ys = []
    for t in range(l):
        dA = torch.exp(dtp[:, t] * A_)
        h = h * dA[..., None, None] + dtp[:, t, :, None, None] * x[:, t, :, :, None] * Bm[:, t, None, None, :]
        ys.append((h * Cm[:, t, None, None, :]).sum(-1) + D[None, :, None] * x[:, t])
    y = torch.stack(ys, 1).reshape(b, l, dim)
    if norm_before_gate:                                  # Mamba2 default is norm_before_gate=False: RMSNorm(y * silu(z))
        return y * torch.rsqrt(y.pow(2).mean(-1, keepdim=True) + rmsnorm_eps) * rmsnorm_weight * F.silu(z)
    y = y * F.silu(z)
    return y * torch.rsqrt(y.pow(2).mean(-1, keepdim=True) + rmsnorm_eps) * rmsnorm_weight


@pytest.fixture(scope='module')
def base_cfg():
    from transformers import AutoConfig
    cfg = AutoConfig.from_pretrained(A.MODEL_DIR, trust_remote_code=True)
    cfg.d_model, cfg.n_layer = D_MODEL, NL
    cfg.fused_add_norm, cfg.rms_norm = False, False
    cfg.ssm_cfg = dict(cfg.ssm_cfg, d_state=8)
    return cfg


def make_model(cfg, seed=0, untie=True):
    from transformers import AutoModel
    torch.manual_seed(seed)
    m = AutoModel.from_config(cfg, trust_remote_code=True).float().eval()
    if untie:                                         # emulate the HF loader, which breaks the tie but keeps the values
        for _, bi in L.mixers(m):
            for name in L.NAMES:
                getattr(bi.mamba_rev, name).weight = torch.nn.Parameter(getattr(bi.mamba_fwd, name).weight.detach().clone())
    m.requires_grad_(False)
    plain_patched(m)                                  # CPU: every block uses the patched (injectable-scan) Mamba2 forward, not only the LoRA blocks
    return m


@pytest.fixture(autouse=True)
def inject_scan():
    L.scan_fn = toy_scan
    yield
    L.scan_fn = None


def rc_ids(ids, comp=(0, 1, 2, 6, 5, 4, 3, 7)):
    return torch.tensor(comp)[ids.flip(-1)]


def ids_batch(b=2, n=37, seed=1):
    return torch.randint(3, 7, (b, n), generator=torch.Generator().manual_seed(seed))


def plain_patched(model):
    """Baseline for 'identity at init': the same patched forward on plain Linear layers (no LoRA wrapper)."""
    import types
    for _, bi in L.mixers(model):
        for m in (bi.mamba_fwd, bi.mamba_rev):
            m.forward = types.MethodType(L._mamba2_forward, m)
    return model


def randomise(params, scale=.05, seed=5):
    g = torch.Generator().manual_seed(seed)
    with torch.no_grad():
        for p in params:
            p.copy_(torch.randn(p.shape, generator=g) * scale)


def test_retie_requires_identical_copies_and_ties(base_cfg):
    m = make_model(base_cfg)
    assert any(bi.mamba_rev.in_proj.weight is not bi.mamba_fwd.in_proj.weight for _, bi in L.mixers(m))
    n = L.tie_check_and_retie(m)
    assert n == NL * 2
    assert all(bi.mamba_rev.in_proj.weight is bi.mamba_fwd.in_proj.weight and bi.mamba_rev.out_proj.weight is bi.mamba_fwd.out_proj.weight for _, bi in L.mixers(m))
    m2 = make_model(base_cfg)
    with torch.no_grad():
        L.mixers(m2)[1][1].mamba_rev.out_proj.weight.add_(1e-3)
    with pytest.raises(AssertionError):
        L.tie_check_and_retie(m2)


def test_param_count_shared_pairs_and_real_size_arithmetic(base_cfg):
    m = make_model(base_cfg)
    layers = (2, 3)
    params = L.install_lora(m, layers, r=8, alpha=16)
    d_in_proj = m.backbone.layers[0].mixer.submodule.mamba_fwd.in_proj.weight.shape[0]
    d_inner = m.backbone.layers[0].mixer.submodule.mamba_fwd.out_proj.weight.shape[1]
    want = len(layers) * (8 * (D_MODEL + d_in_proj) + 8 * (d_inner + D_MODEL))
    assert len({id(p) for p in params}) == len(params) == 4 * len(layers)
    assert sum(p.numel() for p in params) == want
    assert [id(p) for p in params] == [id(p) for p in L.lora_params(m)]
    for i in layers:                                                    # one pair shared by mamba_fwd and mamba_rev
        bi = m.backbone.layers[i].mixer.submodule
        for name in L.NAMES:
            f, r = getattr(bi.mamba_fwd, name), getattr(bi.mamba_rev, name)
            assert f.A is r.A and f.B is r.B and f.base.weight is r.base.weight
    assert L.lora_layers(m) == [2, 3]
    # real PlantCAD2-Small sizes: 12 blocks x (in_proj 768->3224, out_proj 1536->768) at rank 8
    assert 12 * (8 * (768 + 3224) + 8 * (1536 + 768)) == 604_416
    assert 15 * 4 * 8 * 1024 == 491_520                                  # M31 B for comparison


def test_identity_at_init(base_cfg):
    ids = ids_batch()
    ref = plain_patched(make_model(base_cfg))
    m = make_model(base_cfg)
    L.install_lora(m, (1, 2, 3))
    with torch.no_grad():
        a, b = ref(ids).last_hidden_state, m(ids).last_hidden_state
    assert torch.equal(a, b) or float((a - b).abs().max()) < 1e-6


def test_lora_changes_output_when_nonzero(base_cfg):
    ids = ids_batch()
    m = make_model(base_cfg)
    params = L.install_lora(m, (1, 2, 3))
    with torch.no_grad():
        a = m(ids).last_hidden_state.clone()
    randomise(params)
    with torch.no_grad():
        b = m(ids).last_hidden_state
    assert float((a - b).abs().max()) > 1e-3


def test_rc_equivariance_with_nonzero_lora(base_cfg):
    ids = ids_batch(n=41)
    m = make_model(base_cfg)
    params = L.install_lora(m, (1, 2, 3))
    randomise(params, .08)
    with torch.no_grad():
        f = m(ids).last_hidden_state
        g = m(rc_ids(ids)).last_hidden_state
    assert f.shape[-1] == 2 * D_MODEL
    err = float((g - f.flip(-2).flip(-1)).abs().max())
    assert err < 2e-5, err
    # the control: a swap of the two channel halves instead of the channel flip is NOT equivariant
    assert float((g - f.flip(-2).roll(D_MODEL, -1)).abs().max()) > 1e-2


def test_untied_unshared_lora_would_break_nothing_but_tied_is_required(base_cfg):
    """Sanity of the claim 'the shared pair keeps the effective in/out projections tied': effective weights of mamba_fwd and mamba_rev are identical."""
    m = make_model(base_cfg)
    params = L.install_lora(m, (2, 3))
    randomise(params, .1)
    for i in (2, 3):
        bi = m.backbone.layers[i].mixer.submodule
        for name in L.NAMES:
            eff = [getattr(mm, name).base.weight + getattr(mm, name).scale * getattr(mm, name).B @ getattr(mm, name).A for mm in (bi.mamba_fwd, bi.mamba_rev)]
            assert torch.equal(eff[0], eff[1])


def loss_fn(m, ids):
    h = m(ids).last_hidden_state
    return (h.pow(2).mean() + h[:, :, :5].sum() * 1e-2)


def test_grads_only_to_lora_and_checkpoint_equals_plain(base_cfg):
    ids = ids_batch(n=33)
    m = make_model(base_cfg)
    params = L.install_lora(m, (2, 3))
    randomise(params, .05)
    loss_fn(m, ids).backward()
    plain = [p.grad.clone() for p in params]
    assert all(p.grad is not None and bool(torch.isfinite(p.grad).all()) and float(p.grad.abs().sum()) > 0 for p in params)
    others = [n for n, p in m.named_parameters() if p.grad is not None and id(p) not in {id(q) for q in params}]
    assert not others, others
    for p in params:
        p.grad = None
    L.enable_layer_checkpoint(m, (2, 3))
    loss_fn(m, ids).backward()
    for a, p in zip(plain, params):
        assert torch.allclose(a, p.grad, atol=1e-6, rtol=1e-4), float((a - p.grad).abs().max())
    # checkpointing is a no-op without autograd
    with torch.no_grad():
        m(ids)


def test_shared_pair_gets_gradient_from_both_directions(base_cfg):
    ids = ids_batch(n=33)
    m = make_model(base_cfg)
    params = L.install_lora(m, (3,))
    randomise(params, .05)
    loss_fn(m, ids).backward()
    ref = [p.grad.clone() for p in params]
    # break the sharing artificially: only mamba_fwd's wrapper keeps the pair, mamba_rev gets frozen copies -> gradient must change
    bi = m.backbone.layers[3].mixer.submodule
    for name in L.NAMES:
        r = getattr(bi.mamba_rev, name)
        r.A, r.B = torch.nn.Parameter(r.A.detach().clone(), requires_grad=False), torch.nn.Parameter(r.B.detach().clone(), requires_grad=False)
    for p in params:
        p.grad = None
    loss_fn(m, ids).backward()
    assert any(not torch.allclose(a, p.grad, atol=1e-7) for a, p in zip(ref, params))


def test_state_roundtrip_and_subset_loading(base_cfg):
    ids = ids_batch()
    m = make_model(base_cfg, seed=3)
    params = L.install_lora(m, (1, 3))
    randomise(params, .07, seed=9)
    st = L.lora_state(m)
    assert sorted({int(k.split('.')[0]) for k in st}) == [1, 3] and len(st) == 8
    with torch.no_grad():
        a = m(ids).last_hidden_state
    m2 = make_model(base_cfg, seed=3)
    p2 = L.load_lora_state(m2, st)
    assert L.lora_layers(m2) == [1, 3] and len(p2) == 8
    with torch.no_grad():
        b = m2(ids).last_hidden_state
    assert torch.equal(a, b)
    st2 = L.lora_state(m2)
    assert all(torch.equal(st[k], st2[k]) for k in st)
