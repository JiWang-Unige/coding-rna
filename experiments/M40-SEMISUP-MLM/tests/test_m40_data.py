"""M40 CPU dry run on the real M34 data (ont; skipped where the M34 catalogs / schedules / M40 CPU prep are absent):
  * k0 schedule: same species per draw as the M34 T10_k25U schedule (same sampler, same seed), no gene-first draw for any new species, no forbidden contig;
  * repeat index == lowercase of the soft-masked legacy genome on sampled slices; repeat fractions of all 10 species plausible;
  * real draws through m40_train.build_item: CRF targets identical to the M34 example builder for U arms; treatment X wipes new-species targets only;
    clean window ids == the ids the M34 encoder path tokenises; masking rate; loss positions only at unlabeled bases; legacy draws without MLM."""
import json

import numpy as np
import pytest

import m40_common as W
import m40_mlm as MLM

C, T, S = W.C, W.T, W.S
HAVE = (T.SCHEDULES / 'schedule_T10_k25U.jsonl').exists() and (W.SCHEDULES40 / 'schedule_M40_k0.jsonl').exists() and (W.REPEATS / 'summary.json').exists()
pytestmark = pytest.mark.skipif(not HAVE, reason='M34 data / M40 CPU prep not present')


def test_k0_schedule():
    k0 = T.read_schedule(W.schedule_path('M40_k0'))
    k25 = T.read_schedule(W.schedule_path('M40_k25U_MLM'))
    assert len(k0) == len(k25) == 3072 * 8
    assert [r['species'] for r in k0] == [r['species'] for r in k25]
    assert not any(r['mode'] == 'G' for r in k0 if r['species'] in W.NEW)
    assert any(r['mode'] == 'G' for r in k0 if r['species'] in W.LEGACY)
    forb = T.forbidden_contigs()
    assert not any(r['seqid'] in forb for r in k0)
    assert W.schedule_path('M40_k50U_MLM').name == 'schedule_T10_k50U.jsonl'


def test_repeat_index():
    import m40_repeats as RP
    summ = json.load(open(W.REPEATS / 'summary.json'))
    assert set(summ) == set(C.species_table())
    assert all(.05 < v['repeat_fraction'] < .9 for v in summ.values()), {k: v['repeat_fraction'] for k, v in summ.items()}
    key = 'arabidopsis_thaliana'
    sp = C.species_table()[key]
    c = sp['train'][0]
    soft = RP.read_fasta(sp['dir'] / 'genome.fa', {c})[c]
    idx = RP.RepeatIndex([key])
    rng = np.random.default_rng(0)
    for _ in range(20):
        a = int(rng.integers(0, len(soft) - 60000)); b = a + int(rng.integers(1, 60000))
        want = np.array([ch.islower() for ch in soft[a:b]])
        assert np.array_equal(idx.mask(key, c, a, b), want)


@pytest.fixture(scope='module')
def built():
    import m40_repeats as RP
    import m40_train as TR
    out = {}
    for arm in ('M40_k25U_MLM', 'M40_k0'):
        rows = T.read_schedule(W.schedule_path(arm))[:48]
        labs, _ = T.arm_labelers(arm)
        rep = RP.RepeatIndex([k for k in S.ARMS[arm]['species'] if W.mlm_on_species(k)])
        out[arm] = (rows, labs, [TR.build_item(r, labs, arm, rep, TR.CFG) for r in rows])
    return out


def test_targets_and_treatments(built):
    rows, labs, items = built['M40_k25U_MLM']
    for r, it in zip(rows, items):
        ref = T.example(labs[r['species']], 'U', r)
        assert np.array_equal(it['targets']['seg'], ref[4]['seg']) and np.array_equal(it['targets']['region'], ref[4]['region']) and it['counts'] == ref[5]
    rows, labs, items = built['M40_k0']
    assert {r['species'] for r in rows} & set(W.NEW) and {r['species'] for r in rows} & set(W.LEGACY)
    for r, it in zip(rows, items):
        if r['species'] in W.NEW:
            assert it['counts'][:2] == (0, 0) and (it['targets']['seg'] < 0).all() and it['mlm'] is not None
        else:
            ref = T.example(labs[r['species']], 'U', r)
            assert np.array_equal(it['targets']['seg'], ref[4]['seg']) and it['mlm'] is None


def test_masks_on_real_draws(built):
    sel = acgt = 0
    for arm in built:
        rows, labs, items = built[arm]
        for r, it in zip(rows, items):
            if r['species'] in W.LEGACY:
                assert it['mlm'] is None
                continue
            m = it['mlm']
            clean = [MLM.ids_of(s) for s in it['gen']]
            for (a, b), x, y in zip(it['geo']['local_windows'], m['ids'], clean):
                assert len(x) == len(y) == b - a
                assert ((x == y) | (x == MLM.MASK) | (y >= 3)).all()
            ca, cb = it['geo']['genomic_context']
            ctx = MLM.ids_of(labs[r['species']].seqs[r['seqid']][ca:cb])
            sel += m['n_selected']; acgt += int((ctx >= 3).sum())
            # loss positions never at a labeled core base
            known = (it['targets']['seg'] >= 0) | (it['targets']['region'] >= 0)
            lo, hi = r['start'] - ca, r['end'] - ca
            kg = np.zeros(cb - ca, bool)
            kg[lo:hi] = known if r['strand'] == '+' else known[::-1]
            glob = np.concatenate([p + a for (a, _), p in zip(it['geo']['local_windows'], m['pos'])])
            assert not kg[glob].any()
            if arm == 'M40_k0':
                assert m['n_scope_selected'] == m['n_selected']
    assert .045 < sel / acgt < .055
