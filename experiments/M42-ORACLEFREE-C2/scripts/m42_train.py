"""M42 gate-1 adaptation trainer (protocol 35 §4): the M34 trainer `m34_train_lora.main` (imported, NOT edited) on the training chromosomes A of ONE target species,
continuing head + LoRA of T2 (step 2304), with the M36 recipe (512 steps x batch 8, head lr 5e-5, LoRA lr 1e-4, cosine to 1e-6, 5% warm-up, clip 1; m36_common.TRAIN).
Targets come from the M36 seed labeler (`m36_seed_labels.SeedLabeler`, unchanged) with
  * seeds = the frozen OLD-release positives (m42_prep.py: positives.json),
  * known background = A minus the frozen unknown mask (background_<V>.json),
  * treatment 'S' (U, L, U3: everything that is neither a positive locus nor known background is IGNORE = marginalised) or 'F' (F: everything else background).
Arm L replaces the partition-function loss by a per-position cross-entropy over the same 15 role potentials at the known fine-label positions (the DP pool is swapped
inside this process only; the CRF decoder at inference is unchanged).  One draw schedule ('positive-first + uniform', m36_schedule.draw_schedule, seed 1) serves
every arm.  The M34 / M36 modules are only monkey-patched inside this process.
usage: m42_train.py run --arm U|L|F|U3 [--out DIR] [--stop-at N] [--resume] [--workers 8] [--cpus ...]
       m42_train.py prepare --arm A        CPU only: labeler + schedule + fingerprint (m42_arm.json) and label statistics on a sample of draws
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'M34-NTV2-SCALEUP' / 'scripts'))
import m34_train_common as T                 # FIRST: the backend adapter (pc2: the M32 `a1_common` shim) must load before the M29-A1 original can shadow it
import torch
import m34_backend as BK
import m34_train_lora as TL
sys.path.insert(1, str(HERE.parents[1] / 'M36-BOOTSTRAP-ADAPT' / 'scripts'))
import m36_common as W
import m36_seed_labels as SL
import m36_schedule as SCH
import m36_train as MT
for _p in (str(T.ENC_SCRIPTS), str(T.HERE)):             # keep the backend's script dirs first
    if _p in sys.path:
        sys.path.remove(_p)
    sys.path.insert(0, _p)
if str(HERE) not in sys.path:
    sys.path.append(str(HERE))
import m42_common as Z


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()


# ---- arm L: per-position cross-entropy on the known fine labels (same potentials, same normalisation by n_known as the CRF loss) ----------------
def ce_loss_grad(logits, core_start, core_end, targets):
    t0 = time.perf_counter()
    x = np.asarray(logits, dtype=np.float64)
    seg = np.asarray(targets['seg'])
    n = core_end - core_start
    assert seg.shape[0] == n, (seg.shape, n)
    idx = np.flatnonzero(seg >= 0)
    grad = np.zeros(x.shape, dtype=np.float32)
    if idx.size == 0:
        return dict(loss=0., grad=grad, seconds=time.perf_counter() - t0)
    z = x[core_start + idx]
    m = z.max(1, keepdims=True)
    lse = m[:, 0] + np.log(np.exp(z - m).sum(1))
    lab = seg[idx]
    loss = float((lse - z[np.arange(idx.size), lab]).sum())
    p = np.exp(z - lse[:, None])
    p[np.arange(idx.size), lab] -= 1.
    grad[core_start + idx] = p.astype(np.float32)
    return dict(loss=loss, grad=grad, seconds=time.perf_counter() - t0)


class _Ready:
    def __init__(self, r):
        self.r = r

    def get(self):
        return self.r


class CEPool:
    """drop-in for a1_dp.DPPool (same submit signature); computed in the calling thread (a few ms per window)"""

    def __init__(self, workers=0, cpus=None):
        pass

    def submit(self, context, logits, core_start, core_end, targets, epsilon=.1):
        return _Ready(ce_loss_grad(logits, core_start, core_end, targets))

    def close(self):
        pass


# ---- inputs -------------------------------------------------------------------------------------------------------------------------------------
def make_labeler(arm, genome_fa=None):
    spec = Z.ARMS[arm]
    lens = Z.a_lengths()
    rows = SL.load_seed_json(Z.DATA / 'positives.json', set(lens))
    bg = Z.load_background(spec['mask']) if spec['mask'] else None
    lab = SL.SeedLabeler(Z.SP, genome_fa or (Z.M1 / 'genome.fa'), lens, rows, bg)
    return lab, rows


def schedule(rows, lens, steps, batch):
    path = Z.OUT / 'schedules' / f'schedule_{Z.SP}.jsonl'          # one schedule for every arm (same positives)
    draws = SCH.draw_schedule(Z.SP, rows, lens, steps=steps, batch=batch)
    if path.exists():
        assert [json.loads(x) for x in open(path)] == draws, f'{path} differs from the schedule rebuilt from the frozen positives'
    else:
        SCH.write_schedule(draws, path)
    return path, draws


def label_stats(lab, treatment, draws, k=64):
    """fine-known / coarse-only / unknown base fractions over the first k draws (prepare mode)"""
    tot = np.zeros(3)
    for r in draws[:k]:
        _, _, _, _, _, counts, _ = T.example(lab, treatment, r)
        tot += counts
    return dict(draws=k, fine=tot[0] / tot.sum(), coarse_only=tot[1] / tot.sum(), unknown=tot[2] / tot.sum())


def fingerprint(arm, base, sched, args):
    spec = Z.ARMS[arm]
    files = {'positives': md5(Z.DATA / 'positives.json'), 'genome': md5(Z.M1 / 'genome.fa')}
    if spec['mask']:
        files[f"background_{spec['mask']}"] = md5(Z.DATA / f"background_{spec['mask']}.json")
    code = {p.name: md5(p) for p in sorted(HERE.glob('m42_*.py'))}
    return dict(arm=arm, species=Z.SP, spec=spec, base_model=str(base), base_md5=md5(base), inputs=files, schedule=str(sched), schedule_md5=md5(sched),
                steps=args.steps, batch=args.batch, hyper=W.TRAIN, backend=T.BACKEND, code=code, torch=torch.__version__)


def check_fingerprint(out, fp):
    p = Path(out) / 'm42_arm.json'
    keys = ('arm', 'species', 'spec', 'base_md5', 'inputs', 'schedule_md5', 'steps', 'batch', 'backend')
    if p.exists():
        old = json.load(open(p))
        diff = [k for k in keys if json.loads(json.dumps(old.get(k))) != json.loads(json.dumps(fp.get(k)))]
        if diff:
            print(json.dumps(dict(event='refused', reason='fingerprint changed since first run', keys=diff)), flush=True)
            sys.exit(11)
        return old
    return None


class Patched:
    """M34 hooks for one arm (restored on exit)."""

    def __init__(self, armkey, arm, lab, base, forbidden):
        self.armkey, self.arm, self.lab, self.base, self.forbidden = armkey, arm, lab, base, set(forbidden)
        self.saved = []

    def _set(self, obj, name, value):
        self.saved.append((obj, name, getattr(obj, name)))
        setattr(obj, name, value)

    def __enter__(self):
        spec = Z.ARMS[self.arm]
        key, sp, base = self.armkey, Z.SP, self.base
        orig_spec, orig_labelers, orig_load, orig_check = T.arm_spec, T.arm_labelers, BK.load, BK.check_init
        self._set(T, 'arm_spec', lambda a: dict(species=[sp], new_k=None, treatment=spec['treatment']) if a == key else orig_spec(a))
        self._set(T, 'arm_labelers', lambda a, species=None: ({sp: self.lab}, spec['treatment']) if a == key else orig_labelers(a, species))
        self._set(T, 'forbidden_contigs', lambda: self.forbidden)     # evaluation set E is never sampled (the M34 species table is not needed here)

        def load(cfg, lora=True, device='cuda'):
            ctx, params = orig_load(cfg, lora, device)
            if not lora:
                raise RuntimeError('M42 trains head + LoRA (protocol 35 §4.2)')
            if not MT.lora_config_matches(base['config'], cfg):
                raise RuntimeError(f"LoRA configuration differs from the base model's: {base['config'].get('lora')} vs {cfg.get('lora')}")
            BK.copy_lora(ctx, base['lora'])                    # continue from T2's LoRA
            return ctx, params

        def check_init(init):
            orig_check(init)
            MT.check_base(init)

        self._set(BK, 'load', load)
        self._set(BK, 'check_init', check_init)
        if spec['loss'] == 'ce':
            self._set(TL.DP, 'DPPool', CEPool)
        self.cfg_saved = dict(TL.CFG)
        TL.CFG['lr_final'] = W.TRAIN['lr_final']
        return self

    def __exit__(self, *exc):
        for obj, name, value in reversed(self.saved):
            setattr(obj, name, value)
        TL.CFG.clear()
        TL.CFG.update(self.cfg_saved)
        return False


def parse(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=('run', 'prepare'))
    ap.add_argument('--arm', required=True, choices=tuple(Z.ARMS))
    ap.add_argument('--out', default=None)
    ap.add_argument('--base', default=str(Z.BASE))
    ap.add_argument('--steps', type=int, default=W.TRAIN['steps'])
    ap.add_argument('--batch', type=int, default=W.TRAIN['batch'])
    ap.add_argument('--stop-at', type=int, default=None)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--cpus', default='')
    ap.add_argument('--ckpt-every', type=int, default=32)
    ap.add_argument('--keep-every', type=int, default=128)
    ap.add_argument('--resume', action='store_true')
    return ap.parse_args(argv)


def main(argv=None):
    args = parse(argv if argv is not None else sys.argv[1:])
    spec = Z.ARMS[args.arm]
    armkey = f'{args.arm}_{Z.SP}'
    out = Path(args.out) if args.out else Z.OUT / 'train' / armkey
    out.mkdir(parents=True, exist_ok=True)
    sp = Z.split()
    lab, rows = make_labeler(args.arm)
    sched, draws = schedule(rows, Z.a_lengths(), args.steps, args.batch)
    assert not any(r['seqid'] in set(sp['E']) for r in draws)
    fp = fingerprint(args.arm, args.base, sched, args)
    if check_fingerprint(out, fp) is None:
        fp = dict(fp, positives=len(rows), positives_in_labeler=len(lab.catalog), created=time.strftime('%Y-%m-%dT%H:%M:%S'))
        json.dump(fp, open(out / 'm42_arm.json', 'w'), indent=1)
    if args.mode == 'prepare':
        st = label_stats(lab, spec['treatment'], draws)
        print(json.dumps(dict(event='prepared', arm=armkey, schedule=str(sched), draws=len(draws), seed_first=sum(r['mode'] == 'G' for r in draws) / len(draws),
                              positives=len(lab.catalog), base_md5=fp['base_md5'], label_fractions=st)), flush=True)
        return 0
    base = torch.load(args.base, map_location='cpu', weights_only=False)
    MT.check_base(base)
    tl_argv = ['run', '--arm', armkey, '--out', str(out), '--schedule', str(sched), '--init', args.base, '--steps', str(args.steps), '--batch', str(args.batch),
               '--workers', str(args.workers), '--cpus', args.cpus, '--ckpt-every', str(args.ckpt_every), '--keep-every', str(args.keep_every),
               '--lr-head', repr(W.TRAIN['lr_head']), '--lr-lora', repr(W.TRAIN['lr_lora']), '--lora', 'on']
    if args.stop_at:
        tl_argv += ['--stop-at', str(args.stop_at)]
    if args.resume:
        tl_argv += ['--resume']
    sys.argv = ['m34_train_lora.py'] + tl_argv
    with Patched(armkey, args.arm, lab, base, sp['E']):
        TL.main()                                            # ends with os._exit(0) after DONE_<end>.json
    return 0


if __name__ == '__main__':
    sys.exit(main())
