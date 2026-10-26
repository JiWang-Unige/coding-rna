"""P2b scorer (m34_p2b_score.py), CPU only, on ont: (b) reproduces the M33 NT-v2 B unfiltered exact-chain F1 on evaluation set B (Dryas 49.92, Eschscholzia 44.09; report 18 section 3)
with the M30 ruler, (c) the paired 1-Mb block bootstrap: M33 vs itself is exactly 0 with CI [0, 0]; M33 vs the M30 Gate arm is only reported.  Plus synthetic checks of the bootstrap / verdict logic.
Run on ont:  cd .../M34-NTV2-SCALEUP/tests && taskset -c 25-31 nice -n 10 ~/envs/generanno/bin/pytest -q -s test_m34_p2b_score.py
"""
import json
import os
from pathlib import Path

import numpy as np
import pytest

import m34_p2b_score as S

OUT = Path(os.path.expanduser('~/coding-rna/outputs'))
M33 = OUT / 'M33-NTV2-FOLLOWUP'
M33_DIRS = [str(M33 / 'evalB_dryas'), str(M33 / 'evalB_esch')]
GATE = str(OUT / 'M30-HOLDOUT-ADAPT' / 'evalB' / 'gate')
needs_ont = pytest.mark.skipif(not (M33 / 'evalB_dryas').exists(), reason='needs the M33 / M30 outputs on ont')


@pytest.fixture(scope='module')
def res():
    r = S.run_scoring(('NTv2B', M33_DIRS), ('NTv2B_again', M33_DIRS), [('Gate', [GATE])], label='sanity')
    d = OUT / 'M34-NTV2-SCALEUP' / 'p2b_sanity'
    d.mkdir(parents=True, exist_ok=True)
    (d / 'sanity_M33_NTv2B_vs_self_and_Gate.json').write_text(json.dumps(r, indent=1))
    return r


@needs_ont
def test_scorer_reproduces_M33_ntv2b_unfiltered(res):
    t = res['runs']['NTv2B']['table']
    assert round(t['dryas_octopetala']['F1'], 2) == 49.92
    assert round(t['eschscholzia_californica']['F1'], 2) == 44.09
    for sp in S.SPECIES:                                         # counts equal those of M33's own result files (same Evaluator, same GFF3 files)
        m33 = json.load(open(M33 / f'results_{sp}.json'))['tables']
        for arm, name in (('NTv2B', 'NTv2B'), ('Gate', 'Gate')):
            p = m33[arm]['pooled']
            assert (t if arm == 'NTv2B' else res['runs']['Gate']['table'])[sp]['tp'] == p['tp']
            assert (t if arm == 'NTv2B' else res['runs']['Gate']['table'])[sp]['pred'] == p['pred']
            assert (t if arm == 'NTv2B' else res['runs']['Gate']['table'])[sp]['ref'] == p['ref']
            assert abs((t if arm == 'NTv2B' else res['runs']['Gate']['table'])[sp]['F1'] - p['F1']) < 1e-9
    print('pooled NTv2B', res['runs']['NTv2B']['table']['pooled'], 'blocks', res['bootstrap']['blocks'])


@needs_ont
def test_bootstrap_self_difference_is_exactly_zero(res):
    d = res['paired_new_minus']['NTv2B_again']
    for sc in S.SCOPES:
        assert d[sc]['point'] == 0.0 and d[sc]['CI95'] == [0.0, 0.0] and d[sc]['frac_positive'] == 0.0
    assert res['verdict']['met'] is False


@needs_ont
def test_bootstrap_vs_gate_is_consistent(res):
    d = res['paired_new_minus']['Gate']
    for sc in S.SCOPES:
        a, b = res['runs']['NTv2B']['bootstrap_F1'][sc], res['runs']['Gate']['bootstrap_F1'][sc]
        assert abs(d[sc]['point'] - (a['F1'] - b['F1'])) < 1e-9
        assert d[sc]['CI95'][0] <= d[sc]['point'] <= d[sc]['CI95'][1]
        assert a['CI95'][0] <= a['F1'] <= a['CI95'][1]
    print('NTv2B - Gate', json.dumps({sc: dict(point=round(d[sc]['point'], 2), CI95=[round(x, 2) for x in d[sc]['CI95']]) for sc in S.SCOPES}))


def test_bootstrap_logic_synthetic():
    rng = np.random.default_rng(3)
    nb = 60
    bl = [(c, i) for c, n in (('NC_141367.1', 25), ('NW_028238309.1', 35)) for i in range(n)]
    ref = rng.integers(5, 15, nb).astype(float)
    pred_old = ref + rng.integers(-2, 3, nb)
    tp_old = np.minimum(pred_old, ref) - rng.integers(0, 3, nb).clip(max=1)
    tp_new = np.minimum(tp_old + 2, np.minimum(pred_old, ref))                     # new = old with more hits in every block
    C = {'new': (tp_new, pred_old, ref), 'old': (tp_old, pred_old, ref), 'same': (tp_old, pred_old, ref)}
    runs, diffs = S.bootstrap(C, bl, 'new', draws=300, seed=0)
    assert diffs['same']['pooled']['point'] == diffs['old']['pooled']['point'] > 0
    assert diffs['old']['pooled']['CI95'] == diffs['same']['pooled']['CI95'] and diffs['old']['pooled']['CI95'][0] > 0
    rt = S.bootstrap({'new': C['new'], 'old': C['old']}, bl, 'new', draws=300, seed=0)                      # deterministic for a seed
    assert rt[1]['old']['pooled'] == diffs['old']['pooled']
    # pooled F1 of the full data
    assert abs(runs['new']['pooled']['F1'] - 200 * tp_new.sum() / (pred_old.sum() + ref.sum())) < 1e-9


def test_verdict_rule():
    mk = lambda p, lo, a, b: {'pooled': dict(point=p, CI95=[lo, p + 1]), 'dryas_octopetala': dict(point=a, CI95=[0, 1]), 'eschscholzia_californica': dict(point=b, CI95=[0, 1])}
    assert S.verdict(mk(2.0, 0.1, 1.0, 3.0))['met'] is True
    assert S.verdict(mk(1.99, 0.1, 1.0, 3.0))['met'] is False                      # pooled diff below +2.0
    assert S.verdict(mk(3.0, 0.0, 1.0, 3.0))['met'] is False                       # CI lower bound not > 0
    assert S.verdict(mk(3.0, 0.5, -0.1, 6.0))['met'] is False                      # one species negative
    assert S.verdict(mk(3.0, 0.5, 0.0, 6.0))['met'] is False                       # one species exactly zero is not positive
