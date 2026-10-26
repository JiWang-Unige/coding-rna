"""M36 scorer (m36_score.py).  (1) reproduces the M30 / M33 Carbon-A filtered (and unfiltered) B numbers with the project ruler: Dryas 73.98, Eschscholzia 73.32 (reads the
Carbon-A predictions and the B reference of an already-scored method only), (2) the per-species block bootstrap equals the M34 P2b implementation's, (3) closure ratio, (4) every branch of
the protocol 22 section 5 verdict, (5) one-shot sentinels, (6) an end-to-end run of `m36_score.py main` on a FAKE root whose "arms" are existing M33 / M30 outputs (no M36 output exists):
numbers, closure arithmetic, verdict, second call refused (exit 10), changed inputs refused (exit 11).
Run (ont): cd experiments/M36-BOOTSTRAP-ADAPT/tests && nice -n 10 taskset -c 25-31 ~/envs/generanno/bin/pytest -q -s test_m36_score.py
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import m36_common as W
import m36_score as SCR
import m36_sentinel as SN
import m34_p2b_score as PS

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
OUTS = Path(os.path.expanduser('~/coding-rna/outputs'))
M33 = OUTS / 'M33-NTV2-FOLLOWUP'
GATE = OUTS / 'M30-HOLDOUT-ADAPT' / 'evalB' / 'gate'
needs_ont = pytest.mark.skipif(not (M33 / 'evalB_dryas').exists() or not SCR.CARBON_A_DIR.exists(), reason='needs the M30 / M33 outputs on ont')
M30_COUNTS = {'dryas_octopetala': dict(CarbonA_filt=(6807, 8895, 9508), CarbonA_unf=(6885, 9532, 9508)),
              'eschscholzia_californica': dict(CarbonA_filt=(8584, 11048, 12367), CarbonA_unf=(8683, 11838, 12367))}


@needs_ont
def test_carbon_a_numbers_equal_the_m30_results():
    out = SCR.carbon_a_numbers()
    for sp, exp in M30_COUNTS.items():
        for name, (tp, pred, ref) in exp.items():
            r = out[sp][name]
            assert (r['tp'], r['pred'], r['ref']) == (tp, pred, ref), (sp, name, r)
        assert round(out[sp]['CarbonA_filt']['F1'], 2) == SCR.CARBON_A_FILT_EXPECTED[sp]
        m33 = json.load(open(W.EXPS / 'M33-NTV2-FOLLOWUP' / 'results' / f'results_{sp}.json'))['tables']['CarbonA_filt']['pooled'] if (W.EXPS / 'M33-NTV2-FOLLOWUP' / 'results' / f'results_{sp}.json').exists() else None
        if m33:
            assert abs(out[sp]['CarbonA_filt']['F1'] - m33['F1']) < 1e-9
    print({sp: {n: round(v['F1'], 3) for n, v in d.items()} for sp, d in out.items()})
    assert SCR.CARBON_A_FILT_EXPECTED == {'dryas_octopetala': 73.98, 'eschscholzia_californica': 73.32}


def test_block_bootstrap_equals_the_m34_p2b_implementation():
    rng = np.random.default_rng(3)
    n = 80
    bl = [('NC_141367.1', i) for i in range(n)]                                            # one species only: the species scope of the P2b bootstrap == all blocks
    ref = rng.integers(5, 15, n).astype(float)
    pred_a = ref + rng.integers(-2, 3, n)
    tp_a = np.minimum(pred_a, ref) - rng.integers(0, 3, n).clip(max=1)
    tp_b = np.minimum(tp_a + 2, np.minimum(pred_a, ref))
    counts = {'a': (tp_a, pred_a, ref), 'b': (tp_b, pred_a, ref), 'c': (tp_a, pred_a + 1, ref)}
    point, dr = SCR.bootstrap_draws(counts, n, draws=300, seed=0)
    runs, diffs = PS.bootstrap(counts, bl, 'b', draws=300, seed=0)
    for x in ('a', 'c'):
        mine = SCR.pair(point, dr, 'b', x)
        theirs = diffs[x]['dryas_octopetala']
        assert abs(mine['point'] - theirs['point']) < 1e-12 and np.allclose(mine['CI95'], theirs['CI95'], atol=1e-12) and abs(mine['frac_positive'] - theirs['frac_positive']) < 1e-12
    for x in counts:
        assert abs(point[x] - runs[x]['dryas_octopetala']['F1']) < 1e-12 and np.allclose(SCR.ci(dr[x]), runs[x]['dryas_octopetala']['CI95'], atol=1e-12)
    assert abs(point['a'] - 200 * tp_a.sum() / (pred_a.sum() + ref.sum())) < 1e-12
    p1, d1 = SCR.bootstrap_draws(counts, n, draws=50, seed=0)
    p2, d2 = SCR.bootstrap_draws(counts, n, draws=50, seed=0)
    assert all((d1[k] == d2[k]).all() for k in counts)                                       # deterministic
    z = SCR.pair(p1, d1, 'a', 'a')
    assert z['point'] == 0 and z['CI95'] == [0.0, 0.0] and z['frac_positive'] == 0.0           # a run against itself: exactly zero


def test_closure_formula():
    assert abs(SCR.closure(60., 50., 70.) - .5) < 1e-12 and SCR.closure(50., 50., 70.) == 0 and abs(SCR.closure(70., 50., 70.) - 1) < 1e-12 and SCR.closure(80., 50., 70.) > 1
    assert abs(SCR.closure(40., 50., 70.) + .5) < 1e-12
    assert np.isnan(SCR.closure(60., 70., 70.)) and np.isnan(SCR.closure(60., 75., 70.))      # Z at or above Carbon-A: ratio undefined


def mk(closure, d_z=(1., .5, 1.5), d_ca=(-10., -11., -9.)):
    """per-arm summary: closure, (point, lo, hi) of arm - Z and of arm - CarbonA_filt"""
    return dict(closure=closure, diff_vs_Z=dict(point=d_z[0], CI95=[d_z[1], d_z[2]]), diff_vs_CA=dict(point=d_ca[0], CI95=[d_ca[1], d_ca[2]]))


def both(**arms):
    """same summary in both species: arms given as {arm: summary}"""
    return {sp: dict(arms) for sp in W.SPECIES}


def test_verdict_go_strong():
    ok = mk(.95, (23., 22., 24.), (-1.5, -2.5, .4))                                              # within 2.0 points, CI upper bound >= 0
    v = SCR.verdict(both(S1=mk(.2), S2=ok, S3=mk(.4)))
    assert v['go_strong'] and v['label'] == 'GO-strong'
    for bad in (mk(.95, (23., 22., 24.), (-2.01, -3., -1.)),                                    # more than 2.0 below
                mk(.95, (23., 22., 24.), (-1.5, -2.5, -0.01))):                                  # CI excludes parity (upper bound < 0)
        assert not SCR.verdict(both(S1=mk(.2), S2=bad, S3=mk(.2)))['go_strong']
    edge = mk(.95, (23., 22., 24.), (-2.0, -3., 0.0))                                            # exactly -2.0 and an upper bound exactly 0: met
    assert SCR.verdict(both(S2=edge, S1=mk(.1), S3=mk(.1)))['go_strong']
    v = SCR.verdict({W.SPECIES[0]: dict(S2=ok, S1=mk(.1), S3=mk(.1)), W.SPECIES[1]: dict(S2=mk(.5, (12., 11., 13.), (-12., -13., -11.)), S1=mk(.1), S3=mk(.1))})
    assert not v['go_strong']                                                                    # must hold in BOTH species
    v = SCR.verdict({W.SPECIES[0]: dict(S2=ok, S3=ok, S1=mk(.1)), W.SPECIES[1]: dict(S2=mk(.5, (12., 11., 13.), (-12., -13., -11.)), S3=ok, S1=mk(.1))})
    assert v['go_strong'] and v['go_strong_by_arm']['S3'] == {sp: True for sp in W.SPECIES}      # S3 carries both species although S2 does not


def test_verdict_go_and_no_go():
    v = SCR.verdict(both(S1=mk(.31, (7., .1, 14.)), S2=mk(.1), S3=mk(.1)))
    assert v['go'] and v['label'] == 'GO' and not v['go_strong'] and not v['no_go']
    assert SCR.verdict(both(S1=mk(.30, (7., .01, 14.)), S2=mk(.1), S3=mk(.1)))['go']             # closure exactly 30% counts
    assert not SCR.verdict(both(S1=mk(.2999, (7., .1, 14.)), S2=mk(.1), S3=mk(.1)))['go']
    assert not SCR.verdict(both(S1=mk(.5, (12., 0.0, 14.)), S2=mk(.1), S3=mk(.1)))['go']         # CI lower bound exactly 0 is not > 0
    v = SCR.verdict({W.SPECIES[0]: dict(S1=mk(.5, (12., .5, 14.)), S2=mk(.1), S3=mk(.5, (12., .5, 14.))), W.SPECIES[1]: dict(S1=mk(.5, (12., -.5, 14.)), S2=mk(.1), S3=mk(.5, (12., .5, 14.)))})
    assert v['go'] and v['go_by_arm']['S1'][W.SPECIES[1]] is False and v['go_by_arm']['S3'] == {sp: True for sp in W.SPECIES}   # S3 qualifies in both species
    # NO-GO: some species where no arm reaches 30%
    v = SCR.verdict({W.SPECIES[0]: dict(S1=mk(.5, (12., .5, 14.)), S2=mk(.5), S3=mk(.5, (12., .5, 14.))), W.SPECIES[1]: dict(S1=mk(.29), S2=mk(.1), S3=mk(.05))})
    assert v['no_go'] and v['label'] == 'NO-GO' and not v['go']
    # S2 alone reaches 30% in the second species: not NO-GO under the main reading, NO-GO under the stricter S1/S3-only reading; GO fails (S1 / S3 miss it)
    v = SCR.verdict({W.SPECIES[0]: dict(S1=mk(.5, (12., .5, 14.)), S2=mk(.5), S3=mk(.5, (12., .5, 14.))), W.SPECIES[1]: dict(S1=mk(.2), S2=mk(.4), S3=mk(.2))})
    assert not v['no_go'] and v['no_go_reading_S1_S3_only'] and not v['go'] and v['label'].startswith('INCONCLUSIVE')
    # closure >= 30% in both species but a CI that includes 0 -> neither GO nor NO-GO
    v = SCR.verdict(both(S1=mk(.4, (8., -.2, 16.)), S2=mk(.4, (8., -.2, 16.)), S3=mk(.4, (8., -.2, 16.))))
    assert not v['go'] and not v['no_go'] and v['label'].startswith('INCONCLUSIVE')
    # a missing arm / NaN closure never counts as meeting a criterion
    v = SCR.verdict({W.SPECIES[0]: dict(S2=mk(.9)), W.SPECIES[1]: dict(S1=mk(float('nan')))})
    assert not v['go'] and not v['go_strong'] and v['no_go']


def test_verdict_c2_rule():
    base = both(S1=mk(.4, (8., .5, 16.)), S2=mk(.4), S3=mk(.4, (8., .5, 16.)))
    assert SCR.verdict(base, {W.SPECIES[0]: dict(point=-1.99, CI95=[-3., 0.])})['c2_retracted']        # F1 - S1 >= -2: filling background is not worse -> C2 retracted
    assert SCR.verdict(base, {W.SPECIES[0]: dict(point=-2.0, CI95=[-3., 0.])})['c2_retracted']
    assert not SCR.verdict(base, {W.SPECIES[0]: dict(point=-2.01, CI95=[-3., 0.])})['c2_retracted']
    assert not SCR.verdict(base, None)['c2_retracted'] and SCR.verdict(base, {})['c2_by_species'] == {}
    assert SCR.verdict(base, {W.SPECIES[0]: dict(point=-9., CI95=[-10., -8.]), W.SPECIES[1]: dict(point=0.5, CI95=[0., 1.])})['c2_by_species'] == {W.SPECIES[0]: False, W.SPECIES[1]: True}


def test_one_shot_sentinels(tmp_path):
    p = tmp_path / 'sub' / 'S.json'
    ident = dict(models={'a': 'x'}, runs=[1, 2])
    assert SN.start(p, ident) == 0 and p.exists()
    assert SN.start(p, json.loads(json.dumps(ident))) == 0                                          # identical identity: the (resumable) step may restart
    assert SN.start(p, dict(ident, models={'a': 'y'})) == 11                                          # another model / input: refused
    SN.done(p)
    assert SN.start(p, ident) == 10                                                                   # finished: never again
    r = subprocess.run([sys.executable, str(SCRIPTS / 'm36_sentinel.py'), 'start', str(tmp_path / 'T.json'), '{"m": 1}'])
    assert r.returncode == 0
    r = subprocess.run([sys.executable, str(SCRIPTS / 'm36_sentinel.py'), 'start', str(tmp_path / 'T.json'), '{"m": 2}'])
    assert r.returncode == 11


# ---- end to end on a fake root -----------------------------------------------------------------------------------------------------------------------------
def build_fake_root(tmp):
    """fake M36 + M34 output roots: Z = M33 NT-v2 B output (both species), arms = M30 Gate output (S2, S3), M33 NT-v2 output (S1, F1) -- real files, so the ruler runs for real"""
    o34, o36 = tmp / 'm34', tmp / 'm36'
    contigs = {sp: W.TARGETS[sp]['B_never_touched'] for sp in W.SPECIES}
    src = {'ntv2': {'dryas_octopetala': M33 / 'evalB_dryas', 'eschscholzia_californica': M33 / 'evalB_esch'}}

    def link_run(run, sp, kind):
        run.mkdir(parents=True, exist_ok=True)
        for c in contigs[sp]:
            d = run / c / 'prediction'
            d.mkdir(parents=True, exist_ok=True)
            s = (src['ntv2'][sp] / c / 'prediction' / 'full_grammar.gff3') if kind == 'ntv2' else (GATE / c / 'prediction' / 'full_grammar.gff3')
            (d / 'full_grammar.gff3').symlink_to(s)

    z = o34 / 'p2b' / 'T10'
    for sp in W.SPECIES:
        link_run(z, sp, 'ntv2')
    (z / 'DONE_aaaaaaaa.json').write_text(json.dumps(dict(completed=True, model=str(W.BASE_MODEL), regions='all-B', step=3072)))
    kinds = {'S1': 'ntv2', 'S2': 'gate', 'S3': 'gate', 'F1': 'ntv2'}
    for arm, sp in W.order():
        key = W.arm_key(arm, sp)
        link_run(o36 / 'evalB' / key, sp, kinds[arm])
        (o36 / 'evalB' / key / 'DONE_bbbbbbbb.json').write_text(json.dumps(dict(completed=True, model=f'train/{key}/model_000512.pt', regions=','.join(contigs[sp]), step=512)))
        (o36 / 'train' / key).mkdir(parents=True, exist_ok=True)
        (o36 / 'train' / key / 'model_000512.pt').write_bytes(key.encode() * 10)
    return o36, o34


def run_main(o36, o34, cmd='main'):
    env = dict(os.environ, M36_OUT=str(o36), M36_M34_OUT=str(o34))
    return subprocess.run(['nice', '-n', '10', sys.executable, str(SCRIPTS / 'm36_score.py'), cmd, '--root', str(o36)], cwd=SCRIPTS, env=env, capture_output=True, text=True, timeout=7200)


@pytest.fixture(scope='module')
def scored(tmp_path_factory):
    tmp = tmp_path_factory.mktemp('fake')
    o36, o34 = build_fake_root(tmp)
    r = run_main(o36, o34)
    assert r.returncode == 0, r.stderr[-3000:]
    print(r.stdout)
    return o36, o34, tmp, json.load(open(o36 / 'results' / 'm36_result.json'))


@needs_ont
def test_main_end_to_end_on_a_fake_root(scored):
    o36, o34, tmp, res = scored
    assert (o36 / 'results' / 'SCORE_MAIN_STARTED.json.done').exists()
    for sp in W.SPECIES:
        s = res['species'][sp]
        assert s['CarbonA_filt_matches_M30'] and round(s['F1_CarbonA_filt'], 2) == SCR.CARBON_A_FILT_EXPECTED[sp]
        fz, fca = s['F1_Z'], s['F1_CarbonA_filt']
        assert round(fz, 2) == {'dryas_octopetala': 49.92, 'eschscholzia_californica': 44.09}[sp]            # Z = the M33 NT-v2 B output here
        for arm, row in s['arms'].items():
            assert abs(row['closure'] - (row['F1'] - fz) / (fca - fz)) < 1e-12
            assert abs(row['diff_vs_Z']['point'] - (row['F1'] - fz)) < 1e-9 and abs(row['diff_vs_CA']['point'] - (row['F1'] - fca)) < 1e-9
            assert row['diff_vs_Z']['CI95'][0] <= row['diff_vs_Z']['point'] <= row['diff_vs_Z']['CI95'][1]
        assert abs(s['arms']['S1']['closure']) < 1e-12 and s['arms']['S1']['diff_vs_Z']['CI95'] == [0.0, 0.0]        # S1 output == Z output here
        assert 0.89 < s['arms']['S2']['closure'] < 0.95                                                           # Gate output (71.97 / 71.38): ~92% of the Z -> Carbon-A gap
        assert len(s['arms']) == len([1 for a, sp2 in W.order() if sp2 == sp])
        assert isinstance(s['m30_reference_rows'], dict) and round(s['m30_reference_rows']['k25_U'], 2) == {'dryas_octopetala': 71.65, 'eschscholzia_californica': 70.59}[sp]
        assert s['gene_unit_bootstrap']['diffs']['S2 - Z']['lo'] > 0                                              # the M30 gene-unit CI is reported next to the block CI
    v = res['verdict']
    assert v['go_strong'] is False                                                                              # S2 / S3 (Gate output) are 2.0035 points below Carbon-A filtered in Dryas: just outside the 2.0 margin
    assert v['go'] is True and v['label'] == 'GO' and v['no_go'] is False                                       # S3 closes ~92% of the gap with a CI above 0 in both species
    assert res['species']['dryas_octopetala']['pairs']['F1 - S1']['point'] == 0.0 and v['c2_by_species'] == {'dryas_octopetala': True} and v['c2_retracted'] is True
    assert res['verdict_block_bootstrap']['label'] == 'GO' and res['verdict_statistic_sensitive'] is False
    assert res['bootstrap']['kind'].startswith('M30 / M33 gene-unit')                                          # primary statistic per decision m36-ambiguities
    # one-shot: the second call is refused with exit code 10
    r2 = run_main(o36, o34)
    assert r2.returncode == 10 and 'already finished' in r2.stderr


@needs_ont
def test_changed_inputs_and_missing_arms_are_refused(tmp_path):
    o36, o34 = build_fake_root(tmp_path / 'a')
    key = W.arm_key(*W.order()[0])
    env = dict(os.environ, M36_OUT=str(o36), M36_M34_OUT=str(o34))
    code = f"import m36_score as S, m36_sentinel as SN, pathlib; r = pathlib.Path({str(o36)!r}); SN.start(r / 'results' / 'SCORE_MAIN_STARTED.json', S.sentinel_identity(r))"
    assert subprocess.run([sys.executable, '-c', code], cwd=SCRIPTS, env=env).returncode == 0                  # an earlier (interrupted) start with the same inputs
    (o36 / 'train' / key / 'model_000512.pt').write_bytes(b'changed')                                          # ... and now a different model
    r = run_main(o36, o34)
    assert r.returncode == 11 and 'refused' in r.stderr
    # missing arms: scoring does not start (exit 2) and writes no sentinel
    o36b, o34b = build_fake_root(tmp_path / 'b')
    (o36b / 'evalB' / W.arm_key(*W.order()[0]) / 'DONE_bbbbbbbb.json').unlink()
    r4 = run_main(o36b, o34b)
    assert r4.returncode == 2 and not (o36b / 'results' / 'SCORE_MAIN_STARTED.json').exists()


@needs_ont
def test_posterior_filtered_descriptive_scoring(scored):
    """fake posteriors (true chains log P = -0.5, false chains -3.0, computed from the same fake runs): the filter keeps exactly the true chains -> precision 100%, AUROC 1"""
    o36, o34, tmp, res = scored
    post = o36 / 'postB'
    for sp in W.SPECIES:
        S = SCR.SpeciesScorer(sp)
        runs = {f'Z_{sp}': o34 / 'p2b' / 'T10'} | {W.arm_key(a, s): o36 / 'evalB' / W.arm_key(a, s) for a, s in W.order() if s == sp}
        for name, run in runs.items():
            ch = SCR.gff_chains(run, sp)
            for c in SCR.B[sp]:
                for st in '+-':
                    genes = [dict(exons=[list(e) for e in k[2]], gain=0., log_p=-0.5 if k in S.ev.refkeys[c] else -3.0) for k in ch[c] if k[1] == st]
                    d = post / name / 'tau1'
                    d.mkdir(parents=True, exist_ok=True)
                    json.dump(dict(seq=c, strand=st, tau=1.0, genes=genes, log_partition=0., seconds=0.), open(d / f'{c}_{st}.json', 'w'))
    r = run_main(o36, o34, 'post')
    assert r.returncode == 0, r.stderr[-3000:]
    out = json.load(open(o36 / 'results' / 'm36_result_posterior_filtered.json'))
    for sp in W.SPECIES:
        o = out['species'][sp]
        for n, cal in o['calibration'].items():
            assert cal['AUROC'] == 1.0 and cal['chains_without_posterior'] == 0
        z, zp = o['F1']['Z']['F1'], o['F1']['Z_post']['F1']
        assert zp > z + 5                                                                                   # keeping only true chains: F1 = 200 tp / (tp + ref) > unfiltered
        for n, c in o['closure_post'].items():
            assert abs(c - (o['F1'][n + '_post']['F1'] - zp) / (o['F1']['CarbonA_filt']['F1'] - zp)) < 1e-12
        assert set(o['closure_post']) == {a for a, s in W.order() if s == sp}
    assert out['thr'] == W.S1_LOGP_THR and (o36 / 'results' / 'SCORE_POST_STARTED.json.done').exists()
    assert run_main(o36, o34, 'post').returncode == 10                                                       # one-shot as well
