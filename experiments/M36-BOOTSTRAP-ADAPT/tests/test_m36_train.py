"""M36 trainer wrapper (m36_train.py) on CPU: hooks into the M34 trainer, protocol 22 §3 hyper-parameters, the draw -> example -> DP path of every arm, and a full CPU run of the REAL
`m34_train_lora.main` loop with fakes (scripts/m36_cpu_smoke.py, pc2 backend; no GPU, no real encoder).
  Part A (generanno env, no torch GPU): Patched hooks + restore, LoRA continuation call, base-model checks, lr schedule, fingerprint refusal, per-arm examples and window DP loss / gradient
          (synthetic contigs), per-arm real-data windows (ont, Dryas S2 seeds; F1 built from the same seeds).
  Part B (pc2 env as a subprocess): training loop with checkpoint / resume; LoRA of the base model reaches the first forward pass, moves under training, the resumed run starts from the saved state.
Run (ont): cd experiments/M36-BOOTSTRAP-ADAPT/tests && nice -n 10 taskset -c 25-31 ~/envs/generanno/bin/pytest -q -s test_m36_train.py
"""
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import torch

import m34_train_common as T
import m36_common as W
import m36_seed_labels as SL
import m36_schedule as SCH
import m36_train as TR
import m36_train_inputs as TI
import m34_backend as BK
import m34_train_lora as TL
from src.m28.labels import IGNORE
import test_m36_s1_seeds as TS1

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
PC2_PY = os.path.expanduser('~/envs/plantcad2/bin/python')
SP = 'dryas_octopetala'


class FakeCtx:
    def __init__(self):
        self.lora_loaded = None


def make_base(backend=None, lora=None, cfg_lora=None):
    backend = backend or T.BACKEND
    cfg = dict(backend=backend, layer=-1, head_width=1)
    if backend == 'pc2':
        cfg['lora'] = cfg_lora if cfg_lora is not None else dict(rank=8)
    else:
        cfg.update(lora_r=8, lora_alpha=16, lora_encoder_layers=list(range(12, 27)))
    return dict(model={}, lora={'k': torch.ones(2)} if lora is None else lora, step=3072, config=cfg)


# ---- Part A1: hooks ---------------------------------------------------------------------------------------------------------------------------
def test_patched_hooks_install_and_restore(monkeypatch):
    saved = dict(arm_spec=T.arm_spec, arm_labelers=T.arm_labelers, forbidden=T.forbidden_contigs, load=BK.load, check=BK.check_init, copy=BK.copy_lora, cfg=dict(TL.CFG))
    assert TL.CFG['lr_final'] == 1e-5                                                         # M34's value, untouched by importing m36_train
    copied = []
    ctx = FakeCtx()
    cfg_seen = {}

    def fake_load(cfg, lora=True, device='cuda'):
        cfg_seen.update(cfg)
        if T.BACKEND == 'pc2':
            cfg['lora'] = dict(rank=8)
        return ctx, ['p']

    monkeypatch.setattr(BK, 'load', fake_load)
    monkeypatch.setattr(BK, 'copy_lora', lambda c, st: copied.append((c, st)))
    base = make_base()
    lab = object()
    with TR.Patched('S2_dryas_octopetala', 'S2', SP, lab, base):
        assert T.arm_spec('S2_dryas_octopetala') == dict(species=[SP], new_k=None, treatment='S2')
        labs, treat = T.arm_labelers('S2_dryas_octopetala')
        assert labs == {SP: lab} and treat == 'S2'
        assert T.arm_spec('T10') == saved['arm_spec']('T10') and T.arm_spec('T10')['treatment'] == 'FULL'      # other arms still resolve to the M34 definitions
        forb = T.forbidden_contigs()
        assert forb >= saved['forbidden']() and forb >= set(W.TARGETS[SP]['B_never_touched']) and W.TARGETS[SP]['val'][0] in forb
        assert TL.CFG['lr_final'] == 1e-6                                                      # protocol 22 §3: cosine to 1e-6
        cfg = dict(lora_r=8, lora_alpha=16, lora_encoder_layers=list(range(12, 27)))
        c, params = BK.load(cfg, True)
        assert c is ctx and copied == [(ctx, base['lora'])]                                      # the base model's LoRA is loaded into the freshly adapted encoder
        with pytest.raises(RuntimeError, match='head \\+ LoRA'):
            BK.load(cfg, False)
    assert (T.arm_spec, T.arm_labelers, T.forbidden_contigs, BK.load, BK.check_init) == (saved['arm_spec'], saved['arm_labelers'], saved['forbidden'], fake_load, saved['check'])
    assert dict(TL.CFG) == saved['cfg']                                                        # everything restored: M34 behaviour is unchanged outside the context


def test_lora_config_mismatch_is_refused(monkeypatch):
    monkeypatch.setattr(BK, 'load', lambda cfg, lora=True, device='cuda': (FakeCtx(), []))
    monkeypatch.setattr(BK, 'copy_lora', lambda c, st: None)
    base = make_base()
    if T.BACKEND == 'pc2':
        cfg = dict(lora=dict(rank=16))
    else:
        cfg = dict(lora_r=16, lora_alpha=16, lora_encoder_layers=list(range(12, 27)))
    with TR.Patched('S2_dryas_octopetala', 'S2', SP, object(), base):
        with pytest.raises(RuntimeError, match='LoRA configuration'):
            BK.load(cfg, True)


def test_base_model_checks():
    TR.check_base(make_base())
    with pytest.raises(RuntimeError, match='backend'):
        TR.check_base(make_base(backend='pc2' if T.BACKEND == 'ntv2' else 'ntv2'))
    with pytest.raises(RuntimeError, match='no LoRA'):
        TR.check_base(make_base(lora={}))


def test_lr_schedule_is_protocol_22():
    saved = dict(TL.CFG)
    try:
        with TR.Patched('S2_dryas_octopetala', 'S2', SP, object(), make_base()):
            TL.CFG['lr_head_peak'], TL.CFG['lr_lora_peak'] = 5e-5, 1e-4                     # what m36_train passes as --lr-head / --lr-lora
            lr = [TL.lr_pair(i, 512) for i in range(512)]
            assert math.ceil(.05 * 512) == 26 and abs(lr[25][0] - 5e-5) < 1e-18 and abs(lr[25][1] - 1e-4) < 1e-18          # 5% warm-up reaches the peak at step 26
            assert abs(lr[0][0] - 5e-5 / 26) < 1e-18 and all(lr[i][0] < lr[i + 1][0] for i in range(25))
            assert all(lr[i][0] >= lr[i + 1][0] for i in range(26, 511))                       # cosine decay
            assert 1e-6 <= lr[511][0] < 1.05e-6 and 1e-6 <= lr[511][1] < 1.1e-6
            assert all(abs(h - (1e-6 + (5e-5 - 1e-6) * f)) < 1e-15 and abs(l - (1e-6 + (1e-4 - 1e-6) * f)) < 1e-15
                       for i, (h, l) in enumerate(lr[26:], 26) for f in [.5 * (1 + math.cos(math.pi * (i - 26) / (512 - 26)))])     # cosine to 1e-6 from the two peaks
            assert TL.CFG['clip'] == 1. and TL.CFG['weight_decay'] == .01 and TL.CFG['betas'] == [.9, .999] and TL.CFG['epsilon'] == .1 and TL.CFG['max_exons'] == 128
    finally:
        TL.CFG.clear(); TL.CFG.update(saved)
    assert W.TRAIN['steps'] == 512 and W.TRAIN['batch'] == 8 and W.TRAIN['keep_every'] == 128 and W.TRAIN['main_step'] == 512


def test_parse_defaults_match_protocol():
    a = TR.parse(['run', '--arm', 'S1', '--species', SP])
    assert (a.steps, a.batch, a.ckpt_every, a.keep_every, a.workers) == (512, 8, 64, 128, 12) and a.base == str(W.BASE_MODEL)
    with pytest.raises(SystemExit):
        TR.parse(['run', '--arm', 'S1', '--species', 'zea_mays'])
    with pytest.raises(SystemExit):
        TR.parse(['run', '--arm', 'T10', '--species', SP])


def test_fingerprint_refuses_changed_inputs(tmp_path):
    base = tmp_path / 'base.pt'
    base.write_bytes(b'abc')
    sched = tmp_path / 's.jsonl'
    sched.write_text('{}\n')
    args = TR.parse(['run', '--arm', 'S2', '--species', SP])
    inputs = dict(files={'S2': ('x.json', 'aaa')})
    fp = TR.fingerprint('S2', SP, inputs, base, sched, args)
    assert TR.check_fingerprint(tmp_path, fp) is None
    json.dump(fp, open(tmp_path / 'm36_arm.json', 'w'))
    assert TR.check_fingerprint(tmp_path, TR.fingerprint('S2', SP, inputs, base, sched, args))['arm'] == 'S2'             # same inputs: a restart is fine
    for what in ('seed', 'base', 'schedule'):
        base2, sched2, inputs2 = base, sched, inputs
        if what == 'seed':
            inputs2 = dict(files={'S2': ('x.json', 'bbb')})
        elif what == 'base':
            base2 = tmp_path / 'base2.pt'; base2.write_bytes(b'abd')
        else:
            sched2 = tmp_path / 's2.jsonl'; sched2.write_text('{"a":1}\n')
        with pytest.raises(SystemExit) as e:
            TR.check_fingerprint(tmp_path, TR.fingerprint('S2', SP, inputs2, base2, sched2, args))
        assert e.value.code == 11


# ---- Part A2: per-arm examples and DP on synthetic windows -------------------------------------------------------------------------------------
@pytest.fixture(scope='module')
def toy(tmp_path_factory):
    rng = np.random.default_rng(21)
    dna, rows = TS1.synth_contig(rng, length=7000, n_genes=4, seq='SYN', start=300, gap=(100, 400))
    d = tmp_path_factory.mktemp('toy')
    fa = d / 'g.fa'
    fa.write_text('>SYN\n' + '\n'.join(dna[i:i + 80] for i in range(0, len(dna), 80)) + '\n')
    s1, s2 = rows[:3], rows[2:]                                                                # overlapping on rows[2]
    s1 = [dict(r, id='a' + r['id']) for r in s1]
    s2 = [dict(r, id='b' + r['id']) for r in s2]
    bg = {('SYN', '+'): [(0, 150)], ('SYN', '-'): [(6800, 7000)]}
    return dict(fa=fa, dna=dna, lens={'SYN': len(dna)}, s1=s1, s2=s2, bg=bg)


def arm_lab(toy, arm):
    rows, dropped = SL.arm_seeds(arm, toy['s1'], toy['s2'])
    bg = toy['bg'] if arm in ('S1', 'S3') else None
    return SL.SeedLabeler('syn', toy['fa'], toy['lens'], rows, bg), rows, dropped


def draw_rows(rows, n=12):
    sched = SCH.draw_schedule('syn', rows, {'SYN': 7000}, steps=n, batch=1)
    for r in sched:
        r['species'] = 'syn'
    return sched


@pytest.mark.parametrize('arm', list(W.ARMS))
def test_examples_and_dp_for_each_arm(toy, arm):
    import a1_dp as DP
    lab, rows, dropped = arm_lab(toy, arm)
    assert (dropped == 1) == (arm == 'S3')                                                      # S3: the overlapping S1 seed is dropped, S2 wins
    DP._install()
    rng = np.random.default_rng(1)
    skipped = checked = 0
    for r in draw_rows(rows):
        chromosome, geometry, context, offset, targets, counts, win = T.example(lab, arm, r)
        n = r['end'] - r['start']
        seg, region = np.asarray(targets['seg']), np.asarray(targets['region'])
        assert seg.shape == region.shape == (n,) and sum(counts) == n and counts[0] == int((seg >= 0).sum())
        known = np.zeros(n, bool)
        for t in lab.local(r['seqid'], r['strand'], r['start'], r['end']):
            x, y = max(r['start'], t['span'][0]) - r['start'], min(r['end'], t['span'][1]) - r['start']
            known[x:y] = True                                                                  # (oriented bookkeeping below only needs the COUNT of seed bases)
        n_seed_bases = int(known.sum())
        if arm == 'F1':
            assert counts[2] == 0                                                              # Carbon-A-style control: nothing unknown (no N in the toy genome)
        elif arm == 'S2':
            assert counts[1] + counts[0] <= n and (counts[2] > 0 or n_seed_bases == n)         # unknown = everything outside the seeds (no pseudo-background in S2)
            assert counts[0] + counts[1] == n_seed_bases
        else:
            assert counts[0] + counts[1] >= n_seed_bases                                       # S1 / S3: seeds + pseudo-background are known
        if counts[0] + counts[1] == 0:                                                         # all-unknown window: the trainer skips it (zero gradient)
            skipped += 1
            continue
        logits = (rng.standard_normal((len(context), 15)) * 1.5).astype(np.float32)
        out = DP.window_loss_grad(context, logits, offset, offset + n, dict(seg=seg, region=region))
        assert np.isfinite(out['loss']) and out['grad'].shape == logits.shape and np.isfinite(out['grad']).all() and abs(out['grad']).sum() > 0
        checked += 1
    assert checked >= 6
    print(arm, 'windows with a loss', checked, 'skipped (all unknown)', skipped)


def test_f1_equals_s_on_seeds_and_background_elsewhere(toy):
    lab = SL.SeedLabeler('syn', toy['fa'], toy['lens'], toy['s1'], None)                    # one seed set, two treatments
    for strand in '+-':
        tf, _, _ = lab.targets('SYN', strand, 0, 7000, 'F1')
        ts, _, _ = lab.targets('SYN', strand, 0, 7000, 'S1')
        same = ts['seg'][:7000] != IGNORE
        assert (tf['seg'][:7000][same] == ts['seg'][:7000][same]).all()
        assert (tf['seg'][:7000][~same] == 0).all() and (tf['region'][:7000][~same] == 0).all()


# ---- Part A3: real windows (ont) ------------------------------------------------------------------------------------------------------------------
def have_real():
    return TI.s2_path(SP).exists() and (W.M1 / SP / 'genome.fa').exists()


@pytest.mark.skipif(not have_real(), reason='needs the Dryas S2 seeds and genome on ont')
def test_real_dryas_windows_for_each_arm():
    """F1 is built from the same (S2) seeds as a stand-in for the S1 seeds that do not exist yet; S1 / S3 differ from S2 only by the background / seed set (tests above)"""
    inp = TI.load_arm_inputs('S2', SP)
    lens = {c: W.a_lengths(SP)[c] for c in W.TARGETS[SP]['train']}
    lab = {'S2': TI.make_labeler(inp)}
    lab['F1'] = SL.SeedLabeler(SP, None, lens, inp['rows'], None)
    lab['F1'].seqs = lab['S2'].seqs
    sched = SCH.draw_schedule(SP, inp['rows'], lens, steps=8, batch=1)
    seen = []
    for r in sched:
        n = r['end'] - r['start']
        ex = {arm: T.example(lab[arm], arm, r) for arm in ('S2', 'F1')}
        c2, c1 = ex['S2'][5], ex['F1'][5]
        nN = lab['S2'].seqs[r['seqid']][r['start']:r['end']].count('N')
        assert sum(c2) == sum(c1) == n and ex['S2'][0] is lab['S2'].seqs[r['seqid']]
        assert c1[2] <= nN and c1[0] + c1[1] >= c2[0] + c2[1]                                      # F1 knows everything S2 knows, plus the background; only N bases stay unknown
        if r['mode'] == 'G':
            assert c2[0] > 0 and c2[2] > 0                                                         # seed-first windows contain labelled loci; most of the window is unknown in S2
        seen.append((r['mode'], c2, c1))
    print('real Dryas windows (mode, S2 fine/coarse/unknown, F1 fine/coarse/unknown):', seen)


# ---- Part B: real training loop on CPU with fakes (pc2 backend, separate process) ------------------------------------------------------------------
needs_pc2 = pytest.mark.skipif(not Path(PC2_PY).exists(), reason='needs ~/envs/plantcad2')


def run_smoke(work, arm, stop, resume=False):
    env = dict(os.environ, M34_BACKEND='pc2', M36_OUT=str(work / 'm36out'))
    cmd = ['nice', '-n', '10', PC2_PY, str(SCRIPTS / 'm36_cpu_smoke.py'), str(work), arm, str(stop)] + (['resume'] if resume else [])
    r = subprocess.run(cmd, cwd=SCRIPTS, env=env, capture_output=True, text=True, timeout=3000)
    assert r.returncode == 0, (r.stdout[-3000:], r.stderr[-3000:])
    return r


def rows_of(p):
    return [json.loads(x) for x in open(p)]


@needs_pc2
@pytest.mark.parametrize('arm', list(W.ARMS))
def test_cpu_training_loop_runs_for_each_arm(tmp_path_factory, arm):
    work = tmp_path_factory.mktemp(f'smoke_{arm}')
    run_smoke(work, arm, 2)
    out = work / f'train_{arm}'
    log = rows_of(out / 'training.jsonl')
    assert [r['step'] for r in log] == [1, 2] and all(math.isfinite(r['batch_loss']) and math.isfinite(r['grad_norm']) for r in log)
    assert abs(log[0]['lr_head'] - 5e-5) < 1e-12 and abs(log[0]['lr_lora'] - 1e-4) < 1e-12                          # steps = 3 -> warm-up of one step: peak lrs at step 1
    assert (out / 'DONE_2.json').exists() and (out / 'model_000002.pt').exists() and (out / 'ckpt_latest.pt').exists()
    fp = json.load(open(out / 'm36_arm.json'))
    assert fp['arm'] == arm and fp['species'] == SP and fp['backend'] == 'pc2' and fp['hyper']['lr_head'] == 5e-5
    base = torch.load(work / 'data' / 'base_model.pt', map_location='cpu', weights_only=False)
    first = json.load(open(work / f'first_lora_{arm}_start.json'))
    assert first == {k: v.tolist() for k, v in base['lora'].items()}                                                   # the base LoRA reached the first forward pass
    m = torch.load(out / 'model_000002.pt', map_location='cpu', weights_only=False)
    assert m['step'] == 2 and set(m) == {'model', 'lora', 'step', 'config'} and m['config']['lora'] == base['config']['lora'] and m['config']['treatment'] == arm
    assert m['config']['lr_final'] == 1e-6 and m['config']['init_step'] == 3072 and m['config']['species'] == [SP]
    assert any(not torch.equal(m['lora'][k], base['lora'][k]) for k in m['lora'])                                      # LoRA trained (gradient reached it)
    assert any(not torch.equal(m['model'][k], base['model'][k]) for k in m['model'])                                   # head trained


@needs_pc2
def test_cpu_training_resume_equals_uninterrupted(tmp_path_factory):
    a = tmp_path_factory.mktemp('resume_a')
    run_smoke(a, 'S1', 2)
    run_smoke(a, 'S1', 3, resume=True)                                                                                   # restart: ckpt_latest.pt (step 2) -> step 3
    b = tmp_path_factory.mktemp('resume_b')
    run_smoke(b, 'S1', 3)
    la, lb = rows_of(a / 'train_S1' / 'training.jsonl'), rows_of(b / 'train_S1' / 'training.jsonl')
    assert [r['step'] for r in la] == [1, 2, 3] == [r['step'] for r in lb]
    first = json.load(open(a / 'first_lora_S1_resume.json'))
    saved = torch.load(a / 'train_S1' / 'model_000002.pt', map_location='cpu', weights_only=False)['lora']
    assert first == {k: v.tolist() for k, v in saved.items()}                                                          # the resumed run starts from the saved LoRA, not from the base
    ma = torch.load(a / 'train_S1' / 'model_000003.pt', map_location='cpu', weights_only=False)
    mb = torch.load(b / 'train_S1' / 'model_000003.pt', map_location='cpu', weights_only=False)
    for k in ma['lora']:
        assert torch.allclose(ma['lora'][k], mb['lora'][k], atol=1e-6)
    for k in ma['model']:
        assert torch.allclose(ma['model'][k], mb['model'][k], atol=1e-6)
    assert abs(la[2]['batch_loss'] - lb[2]['batch_loss']) < 1e-4


@needs_pc2
def test_cpu_training_refuses_a_changed_seed_set(tmp_path_factory):
    w = tmp_path_factory.mktemp('refuse')
    run_smoke(w, 'S2', 1)
    out = w / 'train_S2'
    fp = json.load(open(out / 'm36_arm.json'))
    fp['schedule_md5'] = '0' * 32
    json.dump(fp, open(out / 'm36_arm.json', 'w'))
    env = dict(os.environ, M34_BACKEND='pc2', M36_OUT=str(w / 'm36out'))
    r = subprocess.run(['nice', '-n', '10', PC2_PY, str(SCRIPTS / 'm36_cpu_smoke.py'), str(w), 'S2', '2', 'resume'], cwd=SCRIPTS, env=env, capture_output=True, text=True, timeout=3000)
    assert r.returncode == 11 and 'fingerprint changed' in r.stdout
