"""M40 CPU unit tests (ont, env ~/envs/plantcad2 + m32_extra; no GPU, no M34 training data):
  * MLM pieces: LUT, lambda schedule, scope mapping on both strands, label-blind masking, 80/10/10 corruption, consistency across overlapping windows,
    each in-scope selected base counted exactly once (ownership), repeat weights;
  * the reproduced RCPS LM head == CaduceusForMaskedLM's own logits on a toy Caduceus (real HF classes), exact RC symmetry of the MLM logits,
    MLM gradient reaches only the M35 LoRA (embedding / tied LM head frozen);
  * m40_train.accumulate (chunked two-pass accumulation, fake encoder + fake DP): lambda = 0 -> same gradients as the M34 per-draw formula; lambda > 0 ->
    same gradients as ONE autograd call on the total loss; treatment-X draws contribute only the MLM term; diagnostic split norms add up."""
import importlib.util
import types
from pathlib import Path

import numpy as np
import pytest
import torch

import m40_mlm as MLM

torch.set_num_threads(2)
M35_TESTS = Path(__file__).resolve().parents[2] / 'M35-PLANTCAD2-LORA' / 'tests'


# ---- 1. pure pieces --------------------------------------------------------------------------------------------------------------------------
def test_lut_matches_pc2_adapter():
    import m40_common as W
    assert np.array_equal(MLM.LUT, W.T.A.LUT)
    assert list(MLM.ids_of('ACGTNacgtR')) == [3, 4, 5, 6, 2, 3, 4, 5, 6, 2]


def test_lambda_schedule():
    f = lambda s: MLM.lam_at(s, 3072, .1, .1, .1)
    assert f(0) == 0 and f(306) == 0
    assert abs(f(307) - .1 / 307) < 1e-15 and f(307) < f(400) < .1
    assert f(612) < .1 and f(613) == .1 and f(614) == .1 and f(3071) == .1
    assert max(f(s) for s in range(3072)) == .1
    assert MLM.lam_at(0, 3072, .1, 0, 0) == .1                                   # smoke setting: on from the first step


@pytest.mark.parametrize('strand', '+-')
def test_scope_mapping(strand):
    n_ctx, lo, hi = 100, 30, 70
    unk = np.zeros(40, bool); unk[:5] = True                                     # first five ORIENTED core positions unknown
    sc = MLM.scope_genomic(n_ctx, lo, hi, unk, strand)
    assert sc[:30].all() and sc[70:].all()                                       # flanks: no CRF label
    want = np.arange(30, 35) if strand == '+' else np.arange(65, 70)
    assert np.array_equal(np.flatnonzero(sc[30:70]) + 30, want)


def _geo(n):
    import m40_common as W
    import m34_backend as BK
    return BK.geometry(0, n)


def test_draw_masks_properties():
    rng0 = np.random.default_rng(3)
    n = 49152
    s = ''.join(rng0.choice(list('ACGTN'), n, p=[.2475] * 4 + [.01]))
    ids = MLM.ids_of(s)
    geo = _geo(n)
    win, own = geo['local_windows'], MLM.owned_ranges(geo)
    assert own[0][0] == 0 and own[-1][1] == n and all(own[i][1] == own[i + 1][0] for i in range(len(own) - 1))
    scope = np.ones(n, bool); scope[10000:30000] = False
    rep = np.zeros(n, bool); rep[:5000] = True
    m = MLM.draw_masks(ids, scope, rep, win, own, np.random.default_rng([40, 7]), .05)
    m_all = MLM.draw_masks(ids, np.ones(n, bool), rep, win, own, np.random.default_rng([40, 7]), .05)
    # label-blind: the corrupted input does not depend on the scope (= the labels)
    assert all(np.array_equal(x, y) for x, y in zip(m['ids'], m_all['ids'])) and m['n_selected'] == m_all['n_selected']
    # rate ~5 % of ACGT, never N; corruption 80/10/10
    cor = np.full(n, -1)
    for (a, b), x in zip(win, m['ids']):
        seg = cor[a:b]
        assert np.all((seg == -1) | (seg == x)), 'overlapping windows disagree on the corrupted input'
        cor[a:b] = x
    changed = cor != ids
    assert not changed[ids == MLM.UNK].any()
    frac = m['n_selected'] / (ids >= 3).sum()
    assert .045 < frac < .055
    n_mask = (cor == MLM.MASK).sum()
    assert .76 < n_mask / m['n_selected'] < .84
    # loss positions: in scope, counted once, true targets, repeat weights
    glob = np.concatenate([p + a for (a, _), p in zip(win, m['pos'])])
    assert len(glob) == len(set(glob.tolist())) == m['n_scope_selected']
    assert scope[glob].all()
    assert np.array_equal(np.concatenate(m['tgt']), ids[glob])
    w = np.concatenate(m['w'])
    assert np.allclose(w[rep[glob]], 1e-4) and np.all(w[~rep[glob]] == 1)
    assert abs(m['wsum'] - w.sum()) < 1e-3
    for (a, b), (lo, hi), p in zip(win, own, m['pos']):
        assert np.all((p + a >= lo) & (p + a < hi))


# ---- 2. LM head on a toy Caduceus ----------------------------------------------------------------------------------------------------------------
def _m35():
    spec = importlib.util.spec_from_file_location('m35_test_lora_cpu', M35_TESTS / 'test_lora_cpu.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope='module')
def toy():
    t = _m35()
    import pc2_lora as L
    from transformers import AutoConfig
    import pc2_adapter as A
    cfg = AutoConfig.from_pretrained(A.MODEL_DIR, trust_remote_code=True)
    cfg.d_model, cfg.n_layer = 32, 3
    cfg.fused_add_norm, cfg.rms_norm = False, False
    cfg.ssm_cfg = dict(cfg.ssm_cfg, d_state=8)
    return t, L, cfg


def _mlm_model(toy, seed=0):
    t, L, cfg = toy
    from transformers import AutoModelForMaskedLM
    torch.manual_seed(seed)
    m = AutoModelForMaskedLM.from_config(cfg, trust_remote_code=True).float().eval()
    m.requires_grad_(False)
    for _, bi in L.mixers(m.caduceus):
        for mm in (bi.mamba_fwd, bi.mamba_rev):
            mm.forward = types.MethodType(L._mamba2_forward, mm)
    return m


def test_lm_head_reproduces_caduceus_mlm_and_rc_symmetry(toy):
    t, L, cfg = toy
    L.scan_fn = t.toy_scan
    try:
        m = _mlm_model(toy)
        ids = torch.randint(3, 7, (2, 29), generator=torch.Generator().manual_seed(1))
        ids[0, 5] = MLM.MASK
        with torch.no_grad():
            want = m(input_ids=ids).logits
            h = m.caduceus(input_ids=ids).last_hidden_state                  # = what AutoModel (Caduceus) returns
            E, comp = MLM.lm_params(m.caduceus)
            got = MLM.lm_logits(h, E, comp)
            assert torch.equal(m.lm_head.weight, E)                          # tied
            assert list(comp.tolist()) == list(MLM.COMPLEMENT)
            assert float((got - want).abs().max()) < 1e-5
            # RC symmetry: logits(RC x) at position p = complement-permuted logits(x) at L-1-p
            rc = torch.tensor(MLM.COMPLEMENT)[ids.flip(-1)]
            g2 = MLM.lm_logits(m.caduceus(input_ids=rc).last_hidden_state, E, comp)
            assert float((g2 - got.flip(1)[..., list(MLM.COMPLEMENT)]).abs().max()) < 1e-4
    finally:
        L.scan_fn = None


def test_mlm_gradient_reaches_only_lora(toy):
    t, L, cfg = toy
    L.scan_fn = t.toy_scan
    try:
        m = _mlm_model(toy)
        enc = m.caduceus
        params = L.install_lora(enc, (1, 2))
        t.randomise(params[1::2], .05)                                       # B != 0 so that every LoRA matrix gets a gradient
        ids = torch.randint(3, 7, (1, 40), generator=torch.Generator().manual_seed(2))
        tgt = ids[0, [3, 9, 17]].clone()
        ids[0, [3, 9, 17]] = MLM.MASK
        h = enc(input_ids=ids).last_hidden_state[0]
        E, comp = MLM.lm_params(enc)
        s, ws, c = MLM.window_mlm(h, torch.tensor([3, 9, 17]), tgt, torch.ones(3), E, comp)
        s.backward()
        assert all(p.grad is not None and float(p.grad.abs().sum()) > 0 for p in params)
        lora_ids = {id(p) for p in params}
        assert E.grad is None and all(p.grad is None for p in enc.parameters() if id(p) not in lora_ids)
    finally:
        L.scan_fn = None


# ---- 3. two-pass accumulation of the trainer ---------------------------------------------------------------------------------------------------
class FakeEnc(torch.nn.Module):
    """Stand-in encoder: token embedding (frozen, also the 'LM head') -> trainable conv over positions -> [b, n, 1536]."""

    def __init__(self):
        super().__init__()
        g = torch.Generator().manual_seed(0)
        emb = torch.nn.Embedding(8, 768)
        with torch.no_grad():
            emb.weight.copy_(torch.randn(8, 768, generator=g) * .1)
        emb.weight.requires_grad_(False)
        self.backbone = types.SimpleNamespace(embeddings=types.SimpleNamespace(word_embeddings=types.SimpleNamespace(
            embedding=emb, complement_map=torch.tensor(MLM.COMPLEMENT))))
        self.emb = emb
        self.conv = torch.nn.Conv1d(768, 1536, 5, padding=2)

    def forward(self, ids):
        x = self.emb(ids).transpose(1, 2)
        return types.SimpleNamespace(last_hidden_state=torch.tanh(self.conv(x)).transpose(1, 2))


class FakeFut:
    def __init__(self, r):
        self.r = r

    def get(self):
        return self.r


def fake_submit(context, logits, cs, ce, targets, eps):
    g = np.random.default_rng(len(context) + cs).standard_normal(logits.shape).astype(np.float32) * .01
    return FakeFut(dict(grad=g, loss=float((logits * g).sum()), seconds=0.))


def _items(lam_on=True, x_second=False):
    import m40_common as W
    import m34_backend as BK
    rng = np.random.default_rng(5)
    chrom = ''.join(rng.choice(list('ACGT'), 60000))
    items = []
    for d, (strand, ca, cb, cs) in enumerate((('+', 1000, 26000, 4000), ('-', 3000, 21000, 2500))):
        geo = BK.geometry(ca, cb)
        n_core = 9000
        seg = np.zeros(n_core, np.int64); region = np.zeros(n_core, np.int64)
        seg[2000:5000] = -100; region[2000:5000] = -100
        counts = (n_core - 3000, 0, 3000)
        if x_second and d == 1:
            seg[:] = -100; region[:] = -100; counts = (0, 0, n_core)
        a = ca + (cs if strand == '+' else (cb - ca) - cs - n_core)
        it = dict(draw=dict(species='populus_trichocarpa', draw=d), w=dict(strand=strand), strand=strand, geo=geo, context=None, core=(cs, cs + n_core),
                  counts=counts, targets=dict(seg=seg, region=region), **BK.item_windows(chrom, geo, strand))
        ctx_seq = chrom[ca:cb]
        it['context'] = ctx_seq if strand == '+' else BK.T.Q.M.rc(ctx_seq)
        unk = (seg < 0) & (region < 0)
        scope = MLM.scope_genomic(cb - ca, a - ca, a - ca + n_core, unk, strand)
        it['mlm'] = MLM.draw_masks(MLM.ids_of(ctx_seq), scope, np.zeros(cb - ca, bool), geo['local_windows'], MLM.owned_ranges(geo),
                                   np.random.default_rng([40, d]), .05)
        items.append(it)
    return items


@pytest.fixture
def setup(monkeypatch):
    import m40_train as TR
    import m34_backend as BK
    enc = FakeEnc()
    monkeypatch.setattr(BK, '_pc2_encode', lambda e, strings, layer=-1, batch=8, grad=False, device='cpu':
                        TR.encode_ids(e, [MLM.ids_of(s) for s in strings], grad, 'cpu'))
    head = BK.MD.new_head(0).train()
    ctx = types.SimpleNamespace(enc=enc, layer=-1)
    params = list(head.parameters()) + list(enc.conv.parameters())
    return TR, BK, enc, head, ctx, params


def _grads(params):
    out = [p.grad.detach().clone() if p.grad is not None else torch.zeros_like(p) for p in params]
    for p in params:
        p.grad = None
    return out


def _close(a, b, tol=2e-4):
    num = sum(float((x - y).abs().max()) for x, y in zip(a, b))
    den = sum(float(y.abs().max()) for y in b) + 1e-12
    return num / den < tol


def test_lambda_zero_equals_m34_formula(setup):
    TR, BK, enc, head, ctx, params = setup
    items = _items()
    fwd = TR.make_fwd(head, ctx, 'cpu')
    st = TR.accumulate(items, fwd, fake_submit, 0., 8, 2, .1, 'cpu')
    got = _grads(params)
    assert st['mlm_draws'] == 0 and st['crf_draws'] == 2
    for it in items:                                                         # the M34 trainer's per-draw formula, chunked as in m34_train_lora
        n = len(it['seqs'])
        with torch.no_grad():
            p = [x.detach() for x in BK.window_pots(head, ctx, it, list(range(n)), False, device='cpu')]
        cl = BK.T.G.owned_context_logits(p, it['geo'], it['strand']).cpu().numpy()
        r = fake_submit(it['context'], cl.astype(np.float32), *it['core'], it['targets'], .1).get()
        g = torch.from_numpy(r['grad']) * (1. / (it['counts'][0] + it['counts'][1]) / 8)
        for a in range(0, n, 2):
            idx = list(range(a, min(a + 2, n)))
            new = BK.window_pots(head, ctx, it, idx, True, device='cpu')
            pieces = list(p)
            for i, x in zip(idx, new):
                pieces[i] = x
            (BK.T.G.owned_context_logits(pieces, it['geo'], it['strand']) * g).sum().backward()
    ref = _grads(params)
    assert _close(got, ref, 1e-6)


@pytest.mark.parametrize('x_second', [False, True])
def test_lambda_on_equals_single_autograd(setup, x_second):
    TR, BK, enc, head, ctx, params = setup
    items = _items(x_second=x_second)
    fwd = TR.make_fwd(head, ctx, 'cpu')
    lam, batch = .3, 8
    st = TR.accumulate(items, fwd, fake_submit, lam, batch, 2, .1, 'cpu', diag_params=list(enc.conv.parameters()))
    got = _grads(params)
    assert st['mlm_draws'] == 2 and st['crf_draws'] == (1 if x_second else 2) and st['mlm_n'] > 0
    # reference: one autograd call on the whole loss, all windows at once, same corrupted input in both passes
    E, comp = MLM.lm_params(enc)
    total = 0.
    for it in items:
        m = it['mlm']
        n = len(it['seqs'])
        h = TR.encode_ids(enc, m['ids'], True, 'cpu')
        if it['counts'][0] + it['counts'][1]:
            with torch.no_grad():
                p0 = BK.MD.window_potentials(head, BK.A.oriented_tokens([x.detach() for x in h], it['strand']), it['seqs'], device='cpu')
            cl = BK.T.G.owned_context_logits(p0, it['geo'], it['strand']).cpu().numpy()
            r = fake_submit(it['context'], cl.astype(np.float32), *it['core'], it['targets'], .1).get()
            g = torch.from_numpy(r['grad']) * (1. / (it['counts'][0] + it['counts'][1]) / batch)
            pots = BK.MD.window_potentials(head, BK.A.oriented_tokens(h, it['strand']), it['seqs'], device='cpu')
            total = total + (BK.T.G.owned_context_logits(pots, it['geo'], it['strand']) * g).sum()
        num = sum(MLM.window_mlm(h[i], torch.from_numpy(m['pos'][i]), torch.from_numpy(m['tgt'][i]), torch.from_numpy(m['w'][i]), E, comp)[0] for i in range(n))
        total = total + lam * num / m['wsum'] / batch
    total.backward()
    ref = _grads(params)
    assert _close(got, ref)
    if x_second:                                                             # the X draw has no CRF term: head gradient only from draw 0
        assert st['skipped'] == 1
    assert st['diag_grad_norm_crf'] > 0 and st['diag_grad_norm_mlm'] > 0
    conv_g = got[len(list(head.parameters())):]
    full = float(torch.sqrt(sum((g ** 2).sum() for g in conv_g)))
    assert full <= st['diag_grad_norm_crf'] + st['diag_grad_norm_mlm'] + 1e-6


def test_head_gets_no_mlm_gradient(setup):
    TR, BK, enc, head, ctx, params = setup
    items = _items(x_second=True)[1:]                                        # one X draw: MLM only
    st = TR.accumulate(items, TR.make_fwd(head, ctx, 'cpu'), fake_submit, .5, 8, 2, .1, 'cpu')
    assert st['crf_draws'] == 0 and st['mlm_draws'] == 1
    assert all(p.grad is None for p in head.parameters())
    assert all(p.grad is not None for p in enc.conv.parameters())
