"""M39 CPU tests (ont data required; skipped otherwise), adapted from the M34 tests: data rules of protocol 29 §3.2 (main chromosomes incl. linkage groups,
S / D / validation-crop choice, trans-spliced transcripts dropped), label masks (U / F / FULL, catalog labeler == parse-based labeler, old species always full),
deterministic schedules (U and F share one draw sequence, no S / D / validation contig), registration of every evaluation contig in the frozen M29 helpers
(references == catalog counts, no plant contig shadowed), resumable-checkpoint round trip, smoke schedule.  Run: ~/envs/generanno/bin/pytest -q tests"""
import hashlib
import json
import random

import numpy as np
import pytest
import torch

import m39_train_common as T
from src.m28.labels import IGNORE
import m39_common as C
import m39_schedule as S

HAVE = (C.OUT / 'schedules' / 'schedule_I10_k50U.jsonl').exists() and all((C.CATALOG / k / 'transcripts.jsonl').exists() for k in C.roster())
pytestmark = pytest.mark.skipif(not HAVE, reason='M39 prepared data / schedules not present')

OLD, NEW = S.LEGACY, S.NEW
SMALL = 'ips_typographus'                       # new species with the smallest validation crop (8 Mb)


def sched(arm):
    return T.read_schedule(T.SCHEDULES / f'schedule_{arm}.jsonl')


# ---- 1. data rules ------------------------------------------------------------------------------------------------------------------
def test_roster_and_main_chromosome_rule():
    assert len(OLD) == 2 and len(NEW) == 8
    fam = [C.roster()[k]['family'] for k in OLD + NEW]
    assert len(set(fam)) == len(fam)                                                     # E7: all families distinct
    assert len({C.roster()[k]['insect_order'] for k in NEW}) >= 4
    rows = [dict(role='assembled-molecule', location=loc, unit='Primary Assembly', length=n, refseq=r) for loc, n, r in
            (('Chromosome', 6e6, 'A'), ('Linkage Group', 7e6, 'B'), ('Linkage Group', 7e6, 'na'), ('Mitochondrion', 9e6, 'C'), ('Chromosome', 4e6, 'D'))]
    assert [r['refseq'] for r in C.main_chromosomes(rows)] == ['A', 'B']
    ch = [dict(refseq=f'c{i}', length=n) for i, n in enumerate((10, 30, 20, 50, 40))]
    assert C.choose_heldout(ch)['refseq'] == 'c1' and C.choose_d(ch)['refseq'] == 'c2'     # lower median of 10,20,30,40,50 = 30; D = next shorter = 20


def test_prepared_splits():
    tab = C.species_table()
    for key in OLD + NEW:
        s = json.load(open(C.PREP / key / 'splits.json'))
        lens = s['lengths']
        assert set(lens) == set(s['train']) | {s['validation_contig']} | set(s['heldout_d'])
        order = sorted(s['main_chromosomes'].items(), key=lambda x: (x[1], x[0]))
        assert s['heldout_chromosome'] == order[(len(order) - 1) // 2][0]
        if key in OLD:
            assert s['heldout_d'] == [order[(len(order) - 1) // 2 - 1][0]] and lens[s['validation_contig']] == s['main_chromosomes'][s['validation_contig']]
        else:
            crop = s['validation_crop']
            assert s['heldout_chromosome'] not in lens and lens[s['validation_contig']] == crop['end'] >= min(C.VAL_MIN_BP, s['heldout_length'])
            assert crop['whole_chromosome'] or crop['end'] % C.VAL_STEP_BP == 0
            assert crop['whole_chromosome'] or crop['references_in_crop_estimate'] >= C.VAL_TARGET_REFS
        assert tab[key]['val'] == [s['validation_contig']]


def test_trans_spliced_dropped_only_in_drosophila():
    for key in OLD + NEW:
        s = json.load(open(C.PREP / key / 'splits.json'))
        n = len(s['dropped_trans_spliced_transcripts']) if 'dropped_trans_spliced_transcripts' in s else 0
        assert (n == 23) if key == 'drosophila_melanogaster' else (n == 0), key
        gff = open(C.PREP / key / 'reference.gff3').read()
        assert 'exception=trans-splicing' not in '\n'.join(l for l in gff.split('\n') if '\tgene\t' not in l)


# ---- 2. schedules -------------------------------------------------------------------------------------------------------------------
def test_schedule_regeneration_is_bitwise_deterministic():
    spec = S.ARMS['I10_k50U']
    pools = {k: S.load_species(k, 1.0 if k in OLD else .5) for k in spec['species']}
    a = S.draw_schedule('I10_k50U', spec, pools)
    assert [json.dumps(r) for r in a] == [json.dumps(r) for r in sched('I10_k50U')]


@pytest.mark.parametrize('arm', list(S.ARMS))
def test_schedule_invariants(arm):
    rows, spec = sched(arm), S.ARMS[arm]
    assert len(rows) == S.STEPS * S.BATCH
    assert [r['draw'] for r in rows] == list(range(len(rows))) and [r['step'] for r in rows] == [i // 8 + 1 for i in range(len(rows))]
    forb = T.forbidden_contigs()
    assert not any(r['seqid'] in forb or '__val' in r['seqid'] for r in rows)
    tab = C.species_table()
    assert all(r['seqid'] in tab[r['species']]['train'] for r in rows)
    assert {r['species'] for r in rows} == set(spec['species'])
    k = len(spec['species'])
    sd = (len(rows) * (1 / k) * (1 - 1 / k)) ** .5
    assert all(abs(sum(r['species'] == s for r in rows) - len(rows) / k) < 5 * sd for s in spec['species'])
    for r in rows[:2000]:
        assert 0 < r['end'] - r['start'] <= 24576 and r['id'] == f"{r['species']}|{r['seqid']}|{r['start']}|{r['strand']}" and r['strand'] in '+-'


def test_U_and_F_share_one_draw_sequence_and_gene_first_draws_are_labeled():
    a, b = sched('I10_k50U'), sched('I10_k50F')
    assert hashlib.md5(json.dumps(a).encode()).hexdigest() == hashlib.md5(json.dumps(b).encode()).hexdigest()
    for key in NEW[:2]:
        labeled = set(json.load(open(C.CATALOG / key / 'labeled_k50.json')))
        genes = {r['gene'] for r in a if r['species'] == key and r['mode'] == 'G'}
        assert genes and genes <= labeled


# ---- 3. label masks -----------------------------------------------------------------------------------------------------------------
def test_catalog_labeler_equals_parsed_labeler():
    import m39_labels as LB
    contig = C.species_table()[SMALL]['val'][0]
    parsed = LB.Labeler(SMALL, [contig], 1.0)
    cat = T.CatalogLabeler(SMALL, [contig], 1.0)
    assert set(parsed.catalog) == set(cat.catalog) and len(cat.catalog) > 300
    rng = random.Random(0)
    for treat in ('FULL', 'U', 'F'):
        for _ in range(12):
            s = rng.choice('+-')
            a, b = rng.choice(parsed.grid(contig))
            x, _, _ = parsed.targets(contig, s, a, b, treat)
            y, _, _ = cat.targets(contig, s, a, b, treat)
            assert np.array_equal(x['seg'], y['seg']) and np.array_equal(x['region'], y['region'])


@pytest.mark.parametrize('key', ['anopheles_gambiae', 'apis_cerana'])
def test_partial_label_masks_U_F_FULL(key):
    lab = T.CatalogLabeler(key, C.species_table()[key]['train'][:2], .5)
    assert abs(lab.n_labeled / lab.n_eligible - .5) < .002
    un = sorted((t for t in lab.catalog.values() if t['eligible'] and t['id'] not in lab.labeled and 1500 < t['span'][1] - t['span'][0] < 15000),
                key=lambda t: (t['seqid'], t['span'][0]))
    checked = 0
    for t in un[:300]:
        c, s = t['seqid'], t['strand']
        a = max(0, t['span'][0] - 3000)
        b = min(lab.lengths[c], a + 24576)
        if any(u is not t and u['span'][0] < t['span'][1] + 200 and u['span'][1] > t['span'][0] - 200 for u in lab.local(c, s, a, b)):
            continue
        U, _, _ = lab.targets(c, s, a, b, 'U')
        F, _, _ = lab.targets(c, s, a, b, 'F')
        FU, _, _ = lab.targets(c, s, a, b, 'FULL')
        x, y = (t['span'][0] - a, t['span'][1] - a) if s == '+' else (b - t['span'][1], b - t['span'][0])
        x, y = max(0, x), min(b - a, y)
        assert (U['seg'][x:y] == IGNORE).all() and (U['region'][x:y] == IGNORE).all()      # U: unlabeled locus marginalised
        assert (F['seg'][x:y] == 0).all() and (F['region'][x:y] == 0).all()                # F: background
        assert (FU['seg'][x:y] > 0).any()                                                  # FULL: labeled
        roles = F['seg'] > 0
        assert (U['seg'][roles] == F['seg'][roles]).all() and (FU['seg'][roles] == F['seg'][roles]).all()
        checked += 1
        if checked == 6:
            break
    assert checked >= 3


def test_old_species_ignore_label_fraction():
    labs, treat = T.arm_labelers('I10_k50F', species=['drosophila_melanogaster', SMALL])
    assert treat == 'F'
    assert labs['drosophila_melanogaster'].n_labeled == labs['drosophila_melanogaster'].n_eligible
    assert abs(labs[SMALL].n_labeled / labs[SMALL].n_eligible - .5) < .002


# ---- 4. registration of the evaluation contigs ---------------------------------------------------------------------------------------
def test_registered_references_match_catalog_counts():
    import m29_common as MC0
    plant = set(MC0.SEQ2SPECIES)
    MC = T.register()
    tab = C.species_table()
    for key in OLD + NEW:
        contig = tab[key]['val'][0]
        assert contig not in plant
        ref = MC.build_reference(contig)
        summ = json.load(open(C.CATALOG / key / 'summary.json'))
        assert len(ref['refs']) == summ['stats']['val']['cds_assessable_references'], key
        assert len(MC.read_sequence(contig)) == MC.LENGTH[contig] == tab[key]['lengths'][contig]
    for c in T.d_regions():
        assert c in MC.LENGTH and len(MC.build_reference(c)['refs']) > 300
    assert T.select_regions() == ['NT_033778.4', 'NC_087399.1'] and T.d_regions() == ['NC_004354.4', 'NC_087400.1']


# ---- 5. checkpoint round trip, smoke schedule -----------------------------------------------------------------------------------------
def _toy(seed):
    import types
    torch.manual_seed(seed)
    import ntv2_adapter as A
    layers = [types.SimpleNamespace(attention=types.SimpleNamespace(self=types.SimpleNamespace(query=torch.nn.Linear(8, 8), value=torch.nn.Linear(8, 8)))) for _ in range(A.N_LAYERS)]
    esm = types.SimpleNamespace(encoder=types.SimpleNamespace(layer=layers))
    params = A.add_lora(esm, [12, 13])
    for p in params[1::2]:
        torch.nn.init.normal_(p, std=.01)
    head = torch.nn.Linear(8, 3)
    opt = torch.optim.AdamW([dict(params=list(head.parameters())), dict(params=params)], lr=1e-3, weight_decay=.01)
    return head, esm, params, opt


def _step(head, esm, params, opt, i):
    x = torch.from_numpy(np.random.default_rng(i).standard_normal((4, 8)).astype(np.float32))
    q = esm.encoder.layer[12].attention.self
    y = head(q.query(x) + q.value(x))
    opt.zero_grad()
    (y.pow(2).mean() + 1e-3 * torch.rand(1)[0]).backward()
    opt.step()


def test_checkpoint_resume_equals_uninterrupted(tmp_path):
    import types
    import m39_ckpt as CK
    import m39_backend as BK
    import ntv2_adapter as A
    cfg = dict(a=1, steps=10)
    h1, e1, p1, o1 = _toy(0)
    for i in range(10):
        _step(h1, e1, p1, o1, i)
    h2, e2, p2, o2 = _toy(0)
    for i in range(6):
        _step(h2, e2, p2, o2, i)
        if i + 1 == 4:
            CK.save_resumable(tmp_path / 'ck.pt', h2, A.lora_state(e2), o2, 4, cfg, 32)
    h3, e3, p3, o3 = _toy(0)
    with torch.no_grad():
        for q in list(h3.parameters()) + p3:
            q.add_(torch.randn_like(q))
    ctx3 = types.SimpleNamespace(esm=e3)
    assert CK.restore(tmp_path / 'ck.pt', h3, lambda st: BK.copy_lora(ctx3, st), o3, cfg) == 4
    for i in range(4, 10):
        _step(h3, e3, p3, o3, i)
    for a, b in zip(list(h1.parameters()) + p1, list(h3.parameters()) + p3):
        assert torch.equal(a, b)
    with pytest.raises(AssertionError):
        CK.restore(tmp_path / 'ck.pt', h3, lambda st: BK.copy_lora(ctx3, st), o3, dict(cfg, steps=11))


def test_smoke_schedule_builder(tmp_path):
    import subprocess, sys
    out = tmp_path / 's.jsonl'
    subprocess.run([sys.executable, str(T.HERE / 'm39_smoke_schedule.py'), 'I10_T', str(out)], check=True, cwd=T.HERE)
    rows = T.read_schedule(out)
    real = sched('I10_T')
    assert len(rows) == len(real) and rows[160:] == real[160:]
    assert len({r['species'] for r in rows[:20]}) == 10
    assert any(int(r['note'].split('span=')[-1]) > 24576 for r in rows[:20])
    assert not any(r['seqid'] in T.forbidden_contigs() for r in rows)


def test_init_head_is_m31_runA():
    import m39_backend as BK
    init = torch.load(T.INIT, map_location='cpu', weights_only=False)
    BK.check_init(init)                                   # step 2048, layers (21, 27)
    assert T.INIT.name == 'model_002048.pt' and T.INIT.parent.name == 'runA' and 'M31-C1-NTV2' in str(T.INIT)
