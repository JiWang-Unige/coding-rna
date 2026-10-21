"""M34 trainer CPU tests (ont data required; skipped otherwise): deterministic multi-species sampling, label masks (U / F / FULL, catalog labeler == parse-based
labeler), T2 example/loss/gradient identity with the M31 single-species path, resumable-checkpoint round trip, smoke schedule, GPU gate, pseudo-contig references."""
import hashlib
import json
import random

import numpy as np
import pytest
import torch

import m34_train_common as T
from src.m28.labels import IGNORE
import m34_common as C
import m34_schedule as S

HAVE = T.BACKEND == 'ntv2' and (C.OUT / 'schedules' / 'schedule_T10.jsonl').exists() and (C.CATALOG / 'populus_trichocarpa' / 'transcripts.jsonl').exists()
pytestmark = pytest.mark.skipif(not HAVE, reason='M34 prepared data / schedules not present')

NEW = S.NEW


def sched(arm):
    return T.read_schedule(T.SCHEDULES / f'schedule_{arm}.jsonl')


# ---- 1. multi-species sampling is deterministic and well formed -----------------------------------------------------------------------
def test_schedule_regeneration_is_bitwise_deterministic():
    for arm in ('T2', 'T4'):
        spec = S.ARMS[arm]
        pools = {k: S.load_species(k, 1.0) for k in spec['species']}
        a = S.draw_schedule(arm, spec, pools)
        b = S.draw_schedule(arm, spec, pools)
        assert a == b
        on_disk = sched(arm)
        assert [json.dumps(r) for r in a] == [json.dumps(r) for r in on_disk], arm


@pytest.mark.parametrize('arm', list(S.ARMS))
def test_schedule_invariants(arm):
    rows, spec = sched(arm), S.ARMS[arm]
    assert len(rows) == S.STEPS * S.BATCH
    assert [r['draw'] for r in rows] == list(range(len(rows))) and [r['step'] for r in rows] == [i // 8 + 1 for i in range(len(rows))]
    forb = T.forbidden_contigs()
    assert not any(r['seqid'] in forb for r in rows)
    tab = C.species_table()
    assert all(r['seqid'] in tab[r['species']]['train'] for r in rows)
    assert {r['species'] for r in rows} == set(spec['species'])
    cnt = {k: sum(r['species'] == k for r in rows) for k in spec['species']}
    exp = len(rows) / len(spec['species'])                                           # equal species weights: counts within 5 sigma of binomial
    sd = (len(rows) * (1 / len(spec['species'])) * (1 - 1 / len(spec['species']))) ** .5
    assert all(abs(v - exp) < 5 * sd for v in cnt.values()), cnt
    for r in rows[:2000]:
        assert 0 < r['end'] - r['start'] <= 24576 and r['id'] == f"{r['species']}|{r['seqid']}|{r['start']}|{r['strand']}" and r['strand'] in '+-'


def test_arm_species_nesting_and_same_draw_law():
    assert S.ARMS['T2']['species'] < S.ARMS['T4']['species'] < S.ARMS['T6']['species'] < S.ARMS['T10']['species']
    for arm in ('T10_k50U', 'T10_k50F'):
        assert S.ARMS[arm]['new_k'] == .5 and S.ARMS[arm]['species'] == S.ARMS['T10']['species']
    a, b = sched('T10_k50U'), sched('T10_k50F')
    assert hashlib.md5(json.dumps(a).encode()).hexdigest() == hashlib.md5(json.dumps(b).encode()).hexdigest()   # U and F share one draw sequence


def test_gene_first_draws_are_labeled_genes():
    elig = {}
    for key in NEW[:2]:
        tr = set(C.species_table()[key]['train'])
        elig[key] = {t['id'] for t in map(json.loads, open(C.CATALOG / key / 'transcripts.jsonl')) if t['eligible'] and t['seqid'] in tr}
    for arm, k in (('T10_k25U', .25), ('T10_k50U', .5), ('T10', 1.0)):
        rows = sched(arm)
        for key in NEW[:2]:
            labeled = set(json.load(open(C.CATALOG / key / f'labeled_k{int(k * 100)}.json'))) if k < 1 else elig[key]
            genes = {r['gene'] for r in rows if r['species'] == key and r['mode'] == 'G'}
            assert genes and genes <= labeled


# ---- 2. label masks ---------------------------------------------------------------------------------------------------------------
def test_catalog_labeler_equals_parsed_labeler():
    import m34_labels as LB
    key = 'brachypodium_distachyon'
    sp = C.species_table()[key]
    contig = sp['val'][0]                                                   # small (8 Mb); the catalog was built over train + val contigs together
    parsed = LB.Labeler(key, [contig], 1.0)
    cat = T.CatalogLabeler(key, [contig], 1.0)
    assert set(parsed.catalog) == set(cat.catalog) and len(cat.catalog) > 300
    for i, t in parsed.catalog.items():
        u = cat.catalog[i]
        for f in t:
            if f == 'split':
                continue
            assert json.dumps(t[f], default=list) == json.dumps(u[f], default=list) or t[f] == u[f], (i, f)
    c, rng = contig, random.Random(0)
    for lab, k in ((parsed, 1.0),):
        for treat in ('FULL', 'U', 'F'):
            for _ in range(12):
                s = rng.choice('+-')
                a, b = rng.choice(parsed.grid(c))
                x, _, _ = parsed.targets(c, s, a, b, treat)
                y, _, _ = cat.targets(c, s, a, b, treat)
                assert np.array_equal(x['seg'], y['seg']) and np.array_equal(x['region'], y['region'])


@pytest.mark.parametrize('key', ['brachypodium_distachyon', 'dioscorea_rotundata'])
def test_partial_label_masks_U_F_FULL(key):
    lab = T.CatalogLabeler(key, C.species_table()[key]['train'][:2], .5)
    assert abs(lab.n_labeled / lab.n_eligible - .5) < .002
    un = [t for t in lab.catalog.values() if t['eligible'] and t['id'] not in lab.labeled and 1500 < t['span'][1] - t['span'][0] < 15000]
    un.sort(key=lambda t: (t['seqid'], t['span'][0]))
    checked = 0
    for t in un[:200]:
        c, s = t['seqid'], t['strand']
        a = max(0, t['span'][0] - 3000)
        b = min(lab.lengths[c], a + 24576)
        # skip windows where another (labeled or unlabeled) gene overlaps the locus, so the locus mask is attributable to this gene alone
        if any(u is not t and u['span'][0] < t['span'][1] + 200 and u['span'][1] > t['span'][0] - 200 for u in lab.local(c, s, a, b)):
            continue
        U, _, _ = lab.targets(c, s, a, b, 'U')
        F, _, _ = lab.targets(c, s, a, b, 'F')
        FU, _, _ = lab.targets(c, s, a, b, 'FULL')
        n = b - a
        x, y = (t['span'][0] - a, t['span'][1] - a) if s == '+' else (b - t['span'][1], b - t['span'][0])
        x, y = max(0, x), min(n, y)
        assert (U['seg'][x:y] == IGNORE).all() and (U['region'][x:y] == IGNORE).all()   # U: the whole unlabeled locus is unknown (marginalised)
        assert (F['seg'][x:y] == 0).all() and (F['region'][x:y] == 0).all()             # F: filled with background
        assert (FU['seg'][x:y] > 0).any()                                               # FULL: the gene is labeled
        lab_roles = (F['seg'] > 0)
        assert (U['seg'][lab_roles] == F['seg'][lab_roles]).all() and (FU['seg'][lab_roles] == F['seg'][lab_roles]).all()   # labeled genes identical in all treatments
        checked += 1
        if checked == 6:
            break
    assert checked >= 3


def test_legacy_species_ignore_label_fraction():
    labs, treat = T.arm_labelers('T10_k25U', species=['arabidopsis_thaliana', 'brachypodium_distachyon'])
    assert treat == 'U'
    assert labs['arabidopsis_thaliana'].n_labeled == labs['arabidopsis_thaliana'].n_eligible          # Ath / Osa stay fully labeled
    assert abs(labs['brachypodium_distachyon'].n_labeled / labs['brachypodium_distachyon'].n_eligible - .25) < .002


# ---- 3. T2 path == M31 single-species path (same inputs -> same loss / gradient bitwise) ----------------------------------------------
def test_t2_examples_and_dp_loss_grad_equal_m31_path():
    import m31_common as M
    import a1_dp as DP
    windows, strata, catalogs, original = M.PP.load_population()
    sequences, _ = M.C.P.read_train_DNA()
    labs, treat = T.arm_labelers('T2')
    assert treat == 'FULL'
    rows = [r for r in sched('T2')[:6000]]
    rng = np.random.default_rng(3)
    picks, stats = [], dict(compared=0, windows_missing=0, target_diff_windows=0, target_diff_positions=0, total_positions=0)
    quota = {}
    for r in rows:
        key = (r['species'], r['strand'], r['mode'])
        if quota.get(key, 0) >= 3:
            continue
        w = windows.get(r['id'])
        if w is None:
            stats['windows_missing'] += 1
            continue
        quota[key] = quota.get(key, 0) + 1
        ch1, geo1, ctx1, off1, tg1, cnt1 = M.example(w, sequences, catalogs)
        ch2, geo2, ctx2, off2, tg2, cnt2, win2 = T.example(labs[r['species']], treat, r)
        assert ch1 is ch2 or ch1 == ch2
        assert geo1['genomic_context'] == geo2['genomic_context'] and ctx1 == ctx2 and off1 == off2
        n = w['valid_bases']
        assert n == win2['valid_bases']
        s1, r1 = np.asarray(tg1['seg'])[:n], np.asarray(tg1['region'])[:n]
        diff = (s1 != tg2['seg']) | (r1 != tg2['region'])
        # the M34 catalog labeler never contradicts a frozen label: every difference is a position the frozen R4 labels mask (IGNORE, alternative-isoform conflict) and M34 labels
        dseg, dreg = s1 != tg2['seg'], r1 != tg2['region']
        assert (s1[dseg] == IGNORE).all() and (r1[dreg] == IGNORE).all()
        d = int(diff.sum())
        stats['compared'] += 1; stats['total_positions'] += n
        stats['target_diff_windows'] += d > 0; stats['target_diff_positions'] += d
        picks.append((ctx1, off1, n, {k: np.asarray(tg1[k])[:n] for k in ('seg', 'region')}, {k: np.asarray(tg2[k]) for k in ('seg', 'region')}, d, cnt1, cnt2))
        if stats['compared'] >= 12:
            break
    print('T2 vs M31 examples:', stats)
    assert stats['compared'] >= 8 and stats['windows_missing'] == 0
    # identical targets -> identical loss and gradient with the shared DP code (same random potentials); differing targets are reported, not hidden
    pool = DP.DPPool(4)
    identical = 0
    for ctx, off, n, t1, t2, d, c1, c2 in picks:
        head = (rng.standard_normal((len(ctx), 15)) * 1.5).astype(np.float32)
        if d == 0:
            assert c1 == c2
            f1 = pool.submit(ctx, head, off, off + n, t1)
            f2 = pool.submit(ctx, head, off, off + n, t2)
            a, b = f1.get(), f2.get()
            assert a['loss'] == b['loss'] and np.array_equal(a['grad'], b['grad'])
            identical += 1
    pool.close()
    print('bitwise-identical loss/grad windows:', identical, 'of', len(picks))
    assert identical >= 1
    # documented label-policy difference (plan 19 section 2.3): the M34 catalog holds primary transcripts only, so the frozen alternative-isoform conflict mask is absent -> a few
    # percent of the gene positions that M31 marks unknown carry the primary-isoform label in M34 (never a contradiction, asserted above).  Only the size is bounded here.
    assert stats['target_diff_positions'] / stats['total_positions'] < 0.03


# ---- 4. checkpoint round trip -------------------------------------------------------------------------------------------------------
def _toy(seed):
    import types
    torch.manual_seed(seed)
    import ntv2_adapter as A
    layers = [types.SimpleNamespace(attention=types.SimpleNamespace(self=types.SimpleNamespace(query=torch.nn.Linear(8, 8), value=torch.nn.Linear(8, 8)))) for _ in range(A.N_LAYERS)]
    esm = types.SimpleNamespace(encoder=types.SimpleNamespace(layer=layers))
    params = A.add_lora(esm, [12, 13])
    for p in params[1::2]:
        torch.nn.init.normal_(p, std=.01)                  # non-zero B so the LoRA matters
    head = torch.nn.Linear(8, 3)
    opt = torch.optim.AdamW([dict(params=list(head.parameters())), dict(params=params)], lr=1e-3, weight_decay=.01)
    return head, esm, params, opt


def _step(head, esm, params, opt, i):
    x = torch.from_numpy(np.random.default_rng(i).standard_normal((4, 8)).astype(np.float32))
    q = esm.encoder.layer[12].attention.self
    y = head(q.query(x) + q.value(x))
    opt.zero_grad()
    (y.pow(2).mean() + 1e-3 * torch.rand(1)[0]).backward()          # torch RNG use -> RNG state must be restored too
    opt.step()


def test_checkpoint_resume_equals_uninterrupted(tmp_path):
    import types
    import m34_ckpt as CK
    import m34_backend as BK
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
    # "crash" after step 6; restart from the step-4 checkpoint in fresh objects and finish
    h3, e3, p3, o3 = _toy(0)                                     # same frozen base weights; trainable state is perturbed so that the restore is what makes it match
    with torch.no_grad():
        for q in list(h3.parameters()) + p3:
            q.add_(torch.randn_like(q))
    ctx3 = types.SimpleNamespace(esm=e3)
    step = CK.restore(tmp_path / 'ck.pt', h3, lambda st: BK.copy_lora(ctx3, st), o3, cfg)
    assert step == 4
    for i in range(4, 10):
        _step(h3, e3, p3, o3, i)
    for a, b in zip(list(h1.parameters()) + p1, list(h3.parameters()) + p3):
        assert torch.equal(a, b)
    assert A.lora_state(e1).keys() == A.lora_state(e3).keys()
    with pytest.raises(AssertionError):
        CK.restore(tmp_path / 'ck.pt', h3, lambda st: BK.copy_lora(ctx3, st), o3, dict(cfg, steps=11))            # config change is refused


# ---- 5. smoke schedule, gate, pseudo-contig references -----------------------------------------------------------------------------
def test_smoke_schedule_builder(tmp_path):
    import subprocess, sys
    out = tmp_path / 's.jsonl'
    subprocess.run([sys.executable, str(T.HERE / 'm34_smoke_schedule.py'), 'T10', str(out)], check=True, cwd=T.HERE)
    rows = T.read_schedule(out)
    real = sched('T10')
    assert len(rows) == len(real) and rows[160:] == real[160:]
    sp = {r['species'] for r in rows[:20]}
    assert len(sp) >= 8 and any(r.get('note', '').split('span=')[-1].isdigit() and int(r['note'].split('span=')[-1]) > 24576 for r in rows[:160])
    forb = T.forbidden_contigs()
    assert not any(r['seqid'] in forb for r in rows)


def test_gate_does_not_flag_itself_and_blocks_without_done_files():
    import m34_gate
    assert m34_gate.os.getpid() not in [p for p, _ in m34_gate.scan()]            # never the pytest process itself
    if not (m34_gate.os.path.exists(m34_gate.FILES[0]) and m34_gate.os.path.exists(m34_gate.FILES[1])):
        assert m34_gate.main() == 1


def test_pseudo_contig_references_match_catalog_counts():
    MC = T.register()
    for key in T.new_keys():
        contig = T.val_contig(key)
        ref = MC.build_reference(contig)
        summ = json.load(open(C.CATALOG / key / 'summary.json'))
        assert len(ref['refs']) == summ['stats']['val']['cds_assessable_references'], key
        assert len(MC.read_sequence(contig)) == MC.LENGTH[contig] == 8_000_000
