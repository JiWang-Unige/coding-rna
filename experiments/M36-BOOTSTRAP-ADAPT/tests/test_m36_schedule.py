"""M36 draw schedule (m36_schedule.py): the 'seed-first + uniform' half/half law of protocol 22 §3, isomorphic to the M34 gene-prior sampler.
Synthetic part (no data): determinism, 4,096 draws = 512 steps x 8, half / half modes, seed-first windows cover their seed, only training contigs, S1 and F1 share one draw
sequence, S3 differs, a B contig is refused, the law reproduces the M34 sampler's draws when fed the same pools.  Data part (ont): the real S2 seeds of both species.
Run (ont): cd experiments/M36-BOOTSTRAP-ADAPT/tests && nice -n 10 taskset -c 25-31 ~/envs/generanno/bin/pytest -q test_m36_schedule.py
"""
import json

import numpy as np
import pytest

import m36_common as W
import m36_schedule as SCH
import m36_seed_labels as SL
import m34_common as C
import m34_schedule as S34

LENS = {'cA': 90_000, 'cB': 61_000}


def seeds(n=60, seed=0, lens=LENS):
    rng = np.random.default_rng(seed)
    rows = []
    for k in range(n):
        c = 'cA' if k % 3 else 'cB'
        a = int(rng.integers(1000, lens[c] - 30000))
        L = int(rng.integers(300, 20000))                                   # includes spans longer than the 12,288 stride
        rows.append(SL.seed_row(c, '+-'[k % 2], [(a, a + L, '0')], f's{k:03d}'))
    return rows


def test_shape_modes_and_determinism():
    rows = SCH.draw_schedule('syn', seeds(), LENS)
    assert len(rows) == 512 * 8 and rows == SCH.draw_schedule('syn', seeds(), LENS)
    assert [r['draw'] for r in rows] == list(range(4096)) and [r['step'] for r in rows] == [i // 8 + 1 for i in range(4096)]
    g = sum(r['mode'] == 'G' for r in rows) / len(rows)
    assert abs(g - .5) < 5 * (.25 / len(rows)) ** .5
    for r in rows[:3000]:
        assert r['seqid'] in LENS and r['strand'] in '+-' and 0 < r['end'] - r['start'] <= 24576 and r['id'] == f"syn|{r['seqid']}|{r['start']}|{r['strand']}"
        assert (r['mode'] == 'G') == ('gene' in r)


def test_seed_first_windows_cover_their_seed_on_its_strand():
    S = seeds()
    by = {r['id']: r for r in S}
    rows = SCH.draw_schedule('syn', S, LENS)
    n = 0
    for r in rows:
        if r['mode'] != 'G':
            continue
        t = by[r['gene']]
        assert r['seqid'] == t['seqid'] and r['strand'] == t['strand'] and r['start'] < t['span'][1] and r['end'] > t['span'][0]
        assert (r['start'], r['end']) in C.gene_cells(t['span'], LENS[r['seqid']])
        n += 1
    assert n > 1800
    assert len({r['gene'] for r in rows if r['mode'] == 'G'}) > 50                        # nearly every seed is drawn at least once


def test_uniform_half_matches_the_m34_law():
    """with no seeds the sampler must equal the M34 uniform law on the same (contig, strand, window) pool: same rng stream, same cell indexing"""
    rows = SCH.draw_schedule('syn', [], LENS)
    assert all(r['mode'] == 'U' for r in rows)
    rng = np.random.default_rng(1)
    names, nwin = list(LENS), [len(C.grid(L)) for L in LENS.values()]
    for r in rows[:500]:
        rng.random()                                                                       # the `rng.random() < uniform_fraction` draw of M34
        tot = 2 * sum(nwin)
        u = int(rng.integers(tot))
        strand = '+' if u < tot // 2 else '-'
        u %= tot // 2
        for c, n in zip(names, nwin):
            if u < n:
                a, b = C.grid(LENS[c])[u]
                break
            u -= n
        assert (r['seqid'], r['strand'], r['start'], r['end']) == (c, strand, a, b)


def test_arms_share_or_differ_as_specified():
    s1 = seeds(60, 0)
    s2 = seeds(40, 1)
    f1 = SCH.draw_schedule('syn', SL.arm_seeds('F1', s1, s2)[0], LENS)
    a1 = SCH.draw_schedule('syn', SL.arm_seeds('S1', s1, s2)[0], LENS)
    a2 = SCH.draw_schedule('syn', SL.arm_seeds('S2', s1, s2)[0], LENS)
    a3 = SCH.draw_schedule('syn', SL.arm_seeds('S3', s1, s2)[0], LENS)
    assert f1 == a1                                                                          # F1 uses S1's seed set -> the same draw sequence (as M34's U / F arms)
    assert a2 != a1 and a3 != a1 and a3 != a2
    assert [r['mode'] for r in a2] != [r['mode'] for r in a1] or [r.get('gene') for r in a2] != [r.get('gene') for r in a1]


def test_pool_is_the_given_training_contigs_only():
    rows = SCH.draw_schedule('syn', seeds(), LENS, contigs=['cA'])
    assert {r['seqid'] for r in rows} == {'cA'}
    rows = SCH.draw_schedule('syn', seeds(), {'cA': LENS['cA']})
    assert {r['seqid'] for r in rows} == {'cA'}                                          # seeds on contigs outside the pool are ignored, never sampled


def test_B_contigs_are_refused():
    sp = 'dryas_octopetala'
    with pytest.raises(RuntimeError, match='evaluation-set B'):
        SCH.draw_schedule(sp, [], {'NC_141367.1': 33_033_793})
    L = W.a_lengths(sp)
    assert not set(L) & set(W.TARGETS[sp]['B_never_touched'])                            # label-set A lengths contain no B contig


def test_write_and_reread(tmp_path):
    rows = SCH.draw_schedule('syn', seeds(), LENS)
    SCH.write_schedule(rows, tmp_path / 's.jsonl')
    assert [json.loads(x) for x in open(tmp_path / 's.jsonl')] == rows


def have_s2():
    import m36_train_inputs as TI
    return all(TI.s2_path(sp).exists() for sp in W.SPECIES)


@pytest.mark.skipif(not have_s2(), reason='needs the S2 seeds on ont')
@pytest.mark.parametrize('species', list(W.SPECIES))
def test_real_s2_schedule(species):
    import m36_train_inputs as TI
    inp = TI.load_arm_inputs('S2', species)
    train = W.TARGETS[species]['train']
    lens = {c: W.a_lengths(species)[c] for c in train}
    rows = SCH.draw_schedule(species, inp['rows'], lens)
    assert {r['seqid'] for r in rows} <= set(train) and not {r['seqid'] for r in rows} & set(W.TARGETS[species]['B_never_touched'])
    assert not any(r['seqid'] == W.TARGETS[species]['val'][0] for r in rows)
    n_train = sum(r['seqid'] in train for r in inp['rows'])
    print(species, 'train seeds', n_train, 'unique windows', len({r['id'] for r in rows}), 'seed-first', sum(r['mode'] == 'G' for r in rows))
    assert n_train > 1900 and sum(r['mode'] == 'G' for r in rows) > 1900
    assert rows == SCH.draw_schedule(species, list(reversed(inp['rows'])), lens)           # independent of the order of the seed list
