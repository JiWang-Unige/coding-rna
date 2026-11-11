"""M41 intron-length prior: parameterisation of the edge-bonus vector, closed-form length scores (for tests), and the maximum-likelihood fit of the
initial values on TRAINING-chromosome intron lengths (Ath NC_003070.9 / NC_003076.8, Osa NC_089043.1 / NC_089046.1 eligible transcripts; never DEV,
never validation chromosomes).

Length model of one intron of length n = 5 + k (k >= 0 bases in the I3* stages):
    P(k) = pi_S * NB(k - (r-1); r, p_S) + (1 - pi_S) * Geom(k; p_L)        [r = short_stages; r = 1: mixture of two geometrics]
The bonus is SHAPE-ONLY: log P(k) - log Geom(k; p0), p0 = the single-geometric MLE (fixed constant, stored with the arm).  The frozen grammar has no
length term, i.e. it is (up to what the head's per-base intron logit absorbs) a single geometric; subtracting Geom(k; p0) keeps the model's starting
point unshifted on average and leaves only the non-geometric part to the new states.
Learnable scalars theta = (a, b_S, b_L):  pi_S = sigmoid(a), p_S = sigmoid(b_S), p_L = sigmoid(b_L).
"""
import json
import math
import sys

import torch
import torch.nn.functional as F

import m41_grammar as GR

LS = F.logsigmoid


def bonus(theta, r, p0, gauge='shape'):
    """theta: float64 tensor (3,) (may require grad) -> bonus vector (N_TYPES,) for `m41_grammar` (differentiable).
    gauge 'shape' (pre-registered): log P(k) - log Geom(k; p0);  'absolute' (dry-run diagnostic only): log P(k)."""
    a, bs, bl = theta[0], theta[1], theta[2]
    lp0, lq0 = (math.log(p0), math.log1p(-p0)) if gauge == 'shape' else (0., 0.)
    z = theta.new_zeros(())
    lps, lqs, lpl, lql = LS(bs), LS(-bs), LS(bl), LS(-bl)
    vec = [z, LS(a), LS(-a), lps - lp0, (lqs - lp0) if r > 1 else z, lqs - lq0, lpl - lp0, lql - lq0]
    return torch.stack(vec)


def log_len_sum(k, theta, r, p0):
    """Closed form: sum over latent branch paths of the bonus for an intron with k interior bases = log P(k) - log Geom(k; p0)."""
    a, bs, bl = (float(x) for x in theta)
    lps, lqs, lpl, lql = (float(LS(torch.tensor(x, dtype=torch.float64))) for x in (bs, -bs, bl, -bl))
    s = k - (r - 1)
    terms = [float(LS(torch.tensor(-a, dtype=torch.float64))) + lql + k * lpl]
    if s >= 0:
        terms.append(float(LS(torch.tensor(a, dtype=torch.float64))) + math.lgamma(s + r) - math.lgamma(s + 1) - math.lgamma(r) + r * lqs + s * lps)
    m = max(terms)
    return m + math.log(sum(math.exp(t - m) for t in terms)) - (k * math.log(p0) + math.log1p(-p0))


def log_len_max(k, theta, r, p0):
    """Closed form: the best single latent path (what the max decoder uses) for k interior bases, minus log Geom(k; p0)."""
    a, bs, bl = (float(x) for x in theta)
    lps, lqs, lpl, lql = (float(LS(torch.tensor(x, dtype=torch.float64))) for x in (bs, -bs, bl, -bl))
    s = k - (r - 1)
    terms = [float(LS(torch.tensor(-a, dtype=torch.float64))) + lql + k * lpl]
    if s >= 0:
        terms.append(float(LS(torch.tensor(a, dtype=torch.float64))) + r * lqs + s * lps)
    return max(terms) - (k * math.log(p0) + math.log1p(-p0))


# ---- fit -------------------------------------------------------------------------------------------------------------------------
def _nll(theta, k, r):
    a, bs, bl = theta
    s = k - (r - 1)
    lgs = torch.lgamma(s.clamp(min=0) + r) - torch.lgamma(s.clamp(min=0) + 1) - math.lgamma(r) + r * LS(-bs) + s.clamp(min=0) * LS(bs)
    lgs = torch.where(s >= 0, lgs, torch.full_like(lgs, -1e30))
    lgl = LS(-bl) + k * LS(bl)
    return -torch.logaddexp(LS(a) + lgs, LS(-a) + lgl).sum()


def fit(lengths, r):
    """MLE of theta for intron lengths (bp, each >= 5) under short_stages r; several starts, best kept.  Returns a JSON-able dict."""
    k = torch.tensor([x - 5 for x in lengths if x >= 5], dtype=torch.float64)
    kbar = float(k.mean())
    p0 = kbar / (kbar + 1.)
    geo = -float((math.log1p(-p0) + k * math.log(p0)).sum())
    best = None
    for pi, ms, ml in ((.8, 80., 800.), (.9, 100., 1500.), (.6, 60., 400.), (.95, 120., 3000.)):
        ps = 1 - r / (ms + 1)
        th = torch.tensor([math.log(pi / (1 - pi)), math.log(ps / (1 - ps)), math.log(ml / 1.)], dtype=torch.float64, requires_grad=True)
        opt = torch.optim.LBFGS([th], max_iter=200, line_search_fn='strong_wolfe', tolerance_grad=1e-10, tolerance_change=1e-12)

        def closure():
            opt.zero_grad()
            v = _nll(th, k, r)
            v.backward()
            return v
        for _ in range(10):
            opt.step(closure)
        v = float(_nll(th.detach(), k, r))
        if best is None or v < best[0]:
            best = (v, th.detach().clone())
    v, th = best
    pi, ps, pl = (float(torch.sigmoid(x)) for x in th)
    n = len(k)
    return dict(short_stages=r, theta=[float(x) for x in th], p0=p0, n_introns=n, nll=v, nll_single_geometric=geo,
                bic=2 * v + 3 * math.log(n), bic_single_geometric=2 * geo + math.log(n), pi_short=pi,
                short_mean_bp=5 + (r - 1) + r * ps / (1 - ps), long_mean_bp=5 + pl / (1 - pl), single_geometric_mean_bp=5 + kbar)


def main():
    """m41_prior.py fit <introns.json> <out.json> [r ...]   introns.json: {species: {"introns": [bp, ...]}} (training chromosomes only)."""
    src, out = sys.argv[2], sys.argv[3]
    rs = [int(x) for x in sys.argv[4:]] or [1]
    d = json.load(open(src))
    lengths = [x for v in d.values() for x in v['introns']]
    res = dict(source=src, species={k: len(v['introns']) for k, v in d.items()}, fits=[fit(lengths, r) for r in rs])
    json.dump(res, open(out, 'w'), indent=1)
    for f in res['fits']:
        print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in f.items()}))


if __name__ == '__main__':
    main()
