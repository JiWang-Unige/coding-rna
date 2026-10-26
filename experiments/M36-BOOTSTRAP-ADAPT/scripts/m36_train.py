"""M36 adaptation trainer (protocol 22 §3): the M34 trainer `m34_train_lora.main` (imported, NOT edited) run on ONE target species with
  * the base model (default W.BASE_MODEL = M34 T10 step 3072): head AND LoRA are loaded from it and training CONTINUES (no new layers);
  * targets from the M36 seed labeler (`m36_seed_labels.SeedLabeler`: S1 / S2 / S3 seed loci labelled and everything else unknown = IGNORE (+ pseudo-background for S1/S3),
    F1: seed loci labelled and everything else background), through the unchanged `m34_train_common.example` path;
  * a precomputed 'seed-first + uniform' half/half draw schedule (m36_schedule.py), 512 steps x batch 8, head lr 5e-5, LoRA lr 1e-4, cosine to 1e-6, 5% warm-up, clip 1,
    AdamW (M34 defaults: weight decay .01, betas .9/.999), resumable snapshot every 64 steps (resume = M34 `ckpt_latest.pt` mechanism), model files every 128 steps.
The M34 modules are only monkey-patched inside this process (T.arm_spec / T.arm_labelers / T.forbidden_contigs / BK.load / BK.check_init / TL.CFG['lr_final']); no M34 file is edited,
so M34 behaviour is unchanged for every M34 invocation.  Backend: M34_BACKEND=pc2 (env ~/envs/plantcad2), asserted against the base model's config.
usage: m36_train.py run --arm S1|S2|S3|F1 --species dryas_octopetala|eschscholzia_californica --out DIR [--base MODEL] [--stop-at N] [--resume] [--workers 12] [--cpus ...]
       m36_train.py prepare --arm A --species S [--out DIR]     CPU only: seeds -> labeler -> schedule -> fingerprint (m36_arm.json); no model, no GPU
"""
import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'M34-NTV2-SCALEUP' / 'scripts'))
import m34_train_common as T                 # FIRST: the backend adapter (pc2: the M32 `a1_common` shim) must load before the M29-A1 original can shadow it
import torch
import m34_backend as BK
import m34_train_lora as TL
import m36_common as W
import m36_train_inputs as TI
import m36_schedule as SCH

for _p in (str(T.ENC_SCRIPTS), str(T.HERE)):             # keep the backend's script dirs first (m36_common put the M30 / M34 dirs at the head)
    if _p in sys.path:
        sys.path.remove(_p)
    sys.path.insert(0, _p)


def md5(path):
    return TI.md5(path)


def arm_spec_for(arm, species):
    return dict(species=[species], new_k=None, treatment=arm)


def check_base(base, backend=None):
    """the base model file (M34 model_N.pt: head + LoRA + config) must come from the backend this process runs"""
    backend = backend or T.BACKEND
    cfg = base['config']
    if cfg.get('backend') != backend:
        raise RuntimeError(f"base model backend {cfg.get('backend')!r} != process backend {backend!r} (set M34_BACKEND / use the matching env)")
    if not base.get('lora'):
        raise RuntimeError('base model has no LoRA state: M36 continues head + LoRA of the base (protocol 22 §3)')


def lora_config_matches(base_cfg, cfg):
    if 'lora' in cfg or 'lora' in base_cfg:                   # pc2: the M35 LoRA configuration dict
        return base_cfg.get('lora') == cfg.get('lora')
    return all(base_cfg.get(k) == cfg.get(k) for k in ('lora_r', 'lora_alpha', 'lora_encoder_layers'))


class Patched:
    """Context manager: M34 hooks for one (arm, species).  Restores everything on exit (tests)."""

    def __init__(self, armkey, arm, species, labeler, base):
        self.armkey, self.arm, self.species, self.labeler, self.base = armkey, arm, species, labeler, base
        self.saved = []

    def _set(self, obj, name, value):
        self.saved.append((obj, name, getattr(obj, name)))
        setattr(obj, name, value)

    def __enter__(self):
        orig_spec, orig_labelers, orig_forb, orig_load, orig_check = T.arm_spec, T.arm_labelers, T.forbidden_contigs, BK.load, BK.check_init
        key, arm, sp = self.armkey, self.arm, self.species
        self._set(T, 'arm_spec', lambda a: arm_spec_for(arm, sp) if a == key else orig_spec(a))
        self._set(T, 'arm_labelers', lambda a, species=None: ({sp: self.labeler}, arm) if a == key else orig_labelers(a, species))
        forb = set(orig_forb())
        if sp in W.TARGETS:                                  # evaluation-set B and the validation pseudo-contig are never sampled
            forb |= set(W.TARGETS[sp]['B_never_touched']) | set(W.TARGETS[sp]['val'])
        self._set(T, 'forbidden_contigs', lambda: forb)
        base = self.base

        def load(cfg, lora=True, device='cuda'):
            ctx, params = orig_load(cfg, lora, device)
            if not lora:
                raise RuntimeError('M36 trains head + LoRA (protocol 22 §3); --lora none is not allowed')
            if not lora_config_matches(base['config'], cfg):
                raise RuntimeError(f"LoRA configuration of the base model differs from the one the encoder was adapted with: {base['config'].get('lora')} vs {cfg.get('lora')}")
            BK.copy_lora(ctx, base['lora'])                  # continue from the base model's LoRA (the M34 trainer would start a fresh zero-initialised one)
            return ctx, params

        def check_init(init):
            orig_check(init)
            check_base(init)

        self._set(BK, 'load', load)
        self._set(BK, 'check_init', check_init)
        self.cfg_saved = dict(TL.CFG)
        TL.CFG['lr_final'] = W.TRAIN['lr_final']             # protocol 22 §3: cosine to 1e-6 (M34: 1e-5)
        return self

    def __exit__(self, *exc):
        for obj, name, value in reversed(self.saved):
            setattr(obj, name, value)
        TL.CFG.clear()
        TL.CFG.update(self.cfg_saved)
        return False


def fingerprint(arm, species, inputs, base_path, sched_path, args):
    return dict(arm=arm, species=species, base_model=str(base_path), base_md5=md5(base_path), seed_files=inputs['files'], schedule=str(sched_path),
                schedule_md5=md5(sched_path), steps=args.steps, batch=args.batch, hyper=W.TRAIN, backend=T.BACKEND)


def check_fingerprint(out, fp):
    """the first run records what it trains on; a restart with different seeds / base / schedule is refused (exit 11)"""
    p = Path(out) / 'm36_arm.json'
    keys = ('arm', 'species', 'base_md5', 'seed_files', 'schedule_md5', 'steps', 'batch', 'backend')
    if p.exists():
        old = json.load(open(p))
        diff = [k for k in keys if json.loads(json.dumps(old.get(k))) != json.loads(json.dumps(fp.get(k)))]
        if diff:
            print(json.dumps(dict(event='refused', reason='fingerprint changed since first run', keys=diff)), flush=True)
            sys.exit(11)
        return old
    return None


def parse(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=('run', 'prepare'))
    ap.add_argument('--arm', required=True, choices=W.ARMS)
    ap.add_argument('--species', required=True, choices=W.SPECIES)
    ap.add_argument('--out', default=None)
    ap.add_argument('--base', default=str(W.BASE_MODEL))
    ap.add_argument('--schedule', default=None)
    ap.add_argument('--steps', type=int, default=W.TRAIN['steps'])
    ap.add_argument('--batch', type=int, default=W.TRAIN['batch'])
    ap.add_argument('--stop-at', type=int, default=None)
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--cpus', default='')
    ap.add_argument('--ckpt-every', type=int, default=W.TRAIN['ckpt_every'])
    ap.add_argument('--keep-every', type=int, default=W.TRAIN['keep_every'])
    ap.add_argument('--resume', action='store_true')
    return ap.parse_args(argv)


def main(argv=None, hooks=None):
    args = parse(argv if argv is not None else sys.argv[1:])
    hooks = hooks or {}
    armkey = W.arm_key(args.arm, args.species)
    out = Path(args.out) if args.out else W.OUT / 'train' / armkey
    out.mkdir(parents=True, exist_ok=True)
    inputs = hooks['inputs'](args.arm, args.species) if 'inputs' in hooks else TI.load_arm_inputs(args.arm, args.species)
    lab = hooks['labeler'](inputs) if 'labeler' in hooks else TI.make_labeler(inputs)
    if args.schedule:
        sched = Path(args.schedule)
    else:
        sched, _ = SCH.make(args.arm, args.species)
    fp = fingerprint(args.arm, args.species, inputs, args.base, sched, args)
    old = check_fingerprint(out, fp)
    if old is None:
        fp = dict(fp, stats=TI.train_stats(inputs, lab), created=time.strftime('%Y-%m-%dT%H:%M:%S'))
        json.dump(fp, open(out / 'm36_arm.json', 'w'), indent=1)
    if args.mode == 'prepare':
        print(json.dumps(dict(event='prepared', arm=armkey, stats=TI.train_stats(inputs, lab), schedule=str(sched), base_md5=fp['base_md5'])), flush=True)
        return 0
    base = torch.load(args.base, map_location='cpu', weights_only=False)
    check_base(base)
    tl_argv = ['run', '--arm', armkey, '--out', str(out), '--schedule', str(sched), '--init', args.base, '--steps', str(args.steps), '--batch', str(args.batch),
               '--workers', str(args.workers), '--cpus', args.cpus, '--ckpt-every', str(args.ckpt_every), '--keep-every', str(args.keep_every),
               '--lr-head', repr(W.TRAIN['lr_head']), '--lr-lora', repr(W.TRAIN['lr_lora']), '--lora', 'on']
    if args.stop_at:
        tl_argv += ['--stop-at', str(args.stop_at)]
    if args.resume:
        tl_argv += ['--resume']
    sys.argv = ['m34_train_lora.py'] + tl_argv
    with Patched(armkey, args.arm, args.species, lab, base):
        TL.main()                                            # ends with os._exit(0) after DONE_<end>.json
    return 0


if __name__ == '__main__':
    sys.exit(main())
