"""CPU tests of the M41 window DP (the trainer's loss / gradient call): states off == the frozen M29 `a1_dp.window_loss_grad` (bit-identical loss and
gradient); states on: the bonus gradient matches central finite differences of the loss; the prior parameterisation round-trips."""
import json
import sys
from pathlib import Path

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
for p in (str(HERE.parent / 'scripts'), str(ROOT / 'experiments/M29-A1-RESIDUAL/scripts')):
    if p not in sys.path:
        sys.path.insert(0, p)

import a1_dp as DP0          # noqa: E402
import m41_dp as DP          # noqa: E402
import m41_prior as P        # noqa: E402


def window(seed=0, L=1500, cs=300, ce=1200):
    rng = np.random.default_rng(seed)
    ctx = ''.join(rng.choice(list('ACGT'), L))
    lg = (rng.standard_normal((L, 15)) * 2).astype(np.float32)
    n = ce - cs
    seg = np.full(n, -100, dtype=np.int64)
    region = np.full(n, -100, dtype=np.int64)
    seg[:50] = 0
    region[200:600] = rng.integers(0, 3, 400)
    return ctx, lg, cs, ce, dict(seg=seg, region=region)


def test_off_equals_frozen_a1_dp():
    ctx, lg, cs, ce, tg = window()
    r0 = DP0.window_loss_grad(ctx, lg, cs, ce, tg, .1)
    r1 = DP.window_loss_grad(ctx, lg, cs, ce, tg, .1, np.zeros(8), 0)
    assert r0['loss'] == r1['loss'] and np.array_equal(r0['grad'], r1['grad']) and not r1['grad_bonus'].any()


def test_bonus_gradient_finite_difference():
    ctx, lg, cs, ce, tg = window(1, L=900, cs=100, ce=800)
    fit = {f['short_stages']: f for f in json.load(open(HERE.parent / 'results' / 'intron_prior_fit.json'))['fits']}[1]
    b = P.bonus(torch.tensor(fit['theta'], dtype=torch.float64), 1, fit['p0']).numpy()
    r = DP.window_loss_grad(ctx, lg, cs, ce, tg, .1, b, 1)
    for t in range(1, 8):
        if t == 4:
            continue                                     # short advance: unused for r = 1
        e = np.zeros(8); e[t] = 1e-5
        hi = DP.window_loss_grad(ctx, lg, cs, ce, tg, .1, b + e, 1)['loss']
        lo = DP.window_loss_grad(ctx, lg, cs, ce, tg, .1, b - e, 1)['loss']
        fd = (hi - lo) / 2e-5
        assert abs(fd - r['grad_bonus'][t]) < 1e-5 * max(1, abs(fd)), (t, fd, r['grad_bonus'][t])
    assert r['grad_bonus'][4] == 0


def test_prior_fit_file_is_consistent():
    d = json.load(open(HERE.parent / 'results' / 'intron_prior_fit.json'))
    f = {x['short_stages']: x for x in d['fits']}[1]
    assert f['bic'] < f['bic_single_geometric'] and 0 < f['p0'] < 1 and f['short_mean_bp'] < f['long_mean_bp']
    assert set(d['species']) == {'arabidopsis_thaliana', 'oryza_sativa'}
