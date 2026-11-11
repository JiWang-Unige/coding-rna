"""M41 trainer = the M34 trainer (`M34-NTV2-SCALEUP/scripts/m34_train_lora.py`, pc2 backend, T2 recipe) with ONE change: the window DP runs on the M41
grammar (short / long intron states, `m41_grammar`) and the three intron-length scalars theta are trained with their own AdamW (protocol
docs/review_20261016/33 section 3).  Everything else -- schedule file, draws, labeler, head, LoRA, optimizer / LR schedule / clipping of head + LoRA,
batch, epsilon, exon cap, checkpoint cadence, seeds -- is the M34 code path; the main loop below is the M34 loop with the marked [M41] lines added.

  [M41] theta: float64 CPU tensor (3,), initialised from results/intron_prior_fit.json (MLE on training-chromosome intron lengths), its own AdamW
        (betas / eps as the head, weight decay 0), cosine LR with peak --lr-intron (pre-registered 1e-3) and final 0.1 x peak (the head's ratio), warmup as the head;
        gradient = the same per-window weights as the head (1 / n_known / batch); clipped separately to norm 1 (the head + LoRA clipping is unchanged);
        --intron-learn off keeps theta fixed at its initial value.
  [M41] checkpoints carry `intron` = dict(theta, short_stages, p0, gauge, optimizer state); model files carry dict(theta, short_stages, p0, gauge, bonus).
Usage: M34_BACKEND=pc2 m41_train.py run --arm T2 --label T2_ILS --out DIR [--short-stages 1] [--steps 3072] [--stop-at N] [--resume]
"""
import argparse
import hashlib
import json
import os
import queue
import sys
import threading
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
M34S = HERE.parents[1] / 'M34-NTV2-SCALEUP' / 'scripts'
sys.path.insert(0, str(M34S))
import numpy as np                       # noqa: E402
import torch                             # noqa: E402

import m34_train_common as T             # noqa: E402  (backend paths)
import m34_backend as BK                 # noqa: E402
import m34_ckpt as CK                    # noqa: E402
import m34_train_lora as TL              # noqa: E402  (CFG, DataPrefetch, lr_pair: unchanged)
if str(HERE) not in sys.path:
    sys.path.insert(1, str(HERE))
import m41_dp as DP                      # noqa: E402
import m41_grammar as GR                 # noqa: E402
import m41_prior as P                    # noqa: E402
M = T.M
from src.m28.full_grammar_observed_partition import InfeasibleObservation   # noqa: E402

CFG = TL.CFG
OUT = T.C.K.OUTS / 'M41-INTRON-LENGTH'
PRIOR = HERE.parent / 'results' / 'intron_prior_fit.json'


def intron_init(r, path=PRIOR):
    fits = {f['short_stages']: f for f in json.load(open(path))['fits']}
    f = fits[r]
    return [float(x) for x in f['theta']], float(f['p0'])


def describe(theta, r):
    a, bs, bl = (float(x) for x in theta)
    sg = lambda x: 1 / (1 + np.exp(-x))
    ps, pl = sg(bs), sg(bl)
    return dict(pi_short=sg(a), short_mean_bp=5 + (r - 1) + r * ps / (1 - ps), long_mean_bp=5 + pl / (1 - pl))


def save_resumable(path, head, lora_state, opt, step, cfg, draw_position, theta, opt_il, il_meta):
    st = CK.snapshot(head, lora_state, opt, step, cfg, draw_position)
    st['intron'] = dict(il_meta, theta=theta.detach().clone(), optimizer=opt_il.state_dict())
    M.atomic_save(st, path)
    return st


def save_model(path, st, r, p0):
    th = st['intron']['theta']
    M.atomic_save(dict(model=st['model'], lora=st['lora'], step=st['step'], config=st['config'],
                       intron=dict(theta=[float(x) for x in th], short_stages=r, p0=p0, gauge='shape',
                                   bonus=[float(x) for x in P.bonus(th.double(), r, p0)])), path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=('run',))
    ap.add_argument('--arm', required=True, help='M34 arm whose data / schedule / seeds are used (T2)')
    ap.add_argument('--label', required=True, help='M41 arm name (T2_ILS)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--schedule', default=None)
    ap.add_argument('--init', default=str(T.INIT))
    ap.add_argument('--lora', choices=('on', 'none'), default='on')
    ap.add_argument('--steps', type=int, default=3072)
    ap.add_argument('--stop-at', type=int, default=None)
    ap.add_argument('--batch', type=int, default=8)
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--cpus', default='')
    ap.add_argument('--ckpt-every', type=int, default=64)
    ap.add_argument('--keep-every', type=int, default=768)
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--short-stages', type=int, default=1)                       # [M41]
    ap.add_argument('--lr-intron', type=float, default=1e-3)                    # [M41]
    ap.add_argument('--intron-learn', choices=('on', 'off'), default='on')      # [M41]
    ap.add_argument('--prior-fit', default=str(PRIOR))                          # [M41]
    args = ap.parse_args()
    out = M.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cpus = M.parse_cpus(args.cpus)
    if cpus:
        os.sched_setaffinity(0, set(cpus))
    torch.set_num_threads(2)
    init = torch.load(args.init, map_location='cpu', weights_only=False)
    BK.check_init(init)
    began = time.perf_counter()
    sched_path = M.Path(args.schedule) if args.schedule else T.SCHEDULES / f'schedule_{args.arm}.jsonl'
    rows = T.read_schedule(sched_path)
    assert len(rows) >= args.steps * args.batch, (len(rows), args.steps, args.batch)
    rows = rows[:args.steps * args.batch]
    forbidden = T.forbidden_contigs()
    assert not any(r['seqid'] in forbidden for r in rows), 'schedule contains DEV / held-out / validation contigs'
    spec = T.arm_spec(args.arm)
    assert {r['species'] for r in rows} <= set(spec['species'])
    for key in ('sampling_seed', 'lora_seed'):
        if key in spec:
            CFG[key] = spec[key]
    sched_md5 = hashlib.md5(open(sched_path, 'rb').read()).hexdigest()
    r_il = args.short_stages                                                     # [M41]
    theta0, p0 = intron_init(r_il, args.prior_fit)                              # [M41]
    il = GR.IntronLength(True, r_il)                                            # [M41]
    CFG.update(BK.encoder_cfg(init['config']), lora_mode=args.lora, init_file=str(args.init), steps=args.steps, batch=args.batch, init_step=init['step'], arm=args.arm,
               species=list(spec['species']), new_species_label_fraction=spec['new_k'], treatment=spec['treatment'], schedule_md5=sched_md5,
               schedule=sched_path.name, heldout=list(T.VAL_REGIONS_LEGACY),
               m41=dict(label=args.label, grammar='short/long intron states', short_stages=r_il, theta_init=theta0, p0=p0, gauge='shape',   # [M41]
                        lr_intron_peak=args.lr_intron, lr_intron_final=.1 * args.lr_intron, intron_learn=args.intron_learn, prior_fit=str(args.prior_fit)))
    labs, treatment = T.arm_labelers(args.arm)
    assert all(r['seqid'] in labs[r['species']].lengths for r in rows)
    label_load_s = time.perf_counter() - began
    head = BK.new_head(CFG, CFG['seed'])
    head.load_state_dict(init['model'])
    head = head.to('cuda').train()
    ctx, lora_params = BK.load(CFG, lora=args.lora == 'on')
    CFG['head_parameters'] = sum(p.numel() for p in head.parameters())
    CFG['lora_parameters'] = sum(p.numel() for p in lora_params)
    groups = [dict(params=list(head.parameters()), name='head')] + ([dict(params=lora_params, name='lora')] if lora_params else [])
    opt = torch.optim.AdamW(groups, lr=CFG['lr_head_peak'], weight_decay=CFG['weight_decay'], betas=tuple(CFG['betas']), eps=CFG['eps'])
    theta = torch.tensor(theta0, dtype=torch.float64, requires_grad=True)                                                     # [M41]
    opt_il = torch.optim.AdamW([theta], lr=args.lr_intron, weight_decay=0., betas=tuple(CFG['betas']), eps=CFG['eps'])        # [M41]
    il_meta = dict(short_stages=r_il, p0=p0, gauge='shape')
    start = 0
    ck = out / 'ckpt_latest.pt'
    if args.resume and ck.exists():
        start = CK.restore(ck, head, lambda st: BK.copy_lora(ctx, st), opt, CFG)
        st_il = torch.load(ck, map_location='cpu', weights_only=False)['intron']                                                # [M41]
        assert {k: st_il[k] for k in il_meta} == il_meta
        with torch.no_grad():
            theta.copy_(st_il['theta'])
        opt_il.load_state_dict(st_il['optimizer'])
        print(json.dumps(dict(event='resumed', from_step=start, time=time.time())), flush=True)
    pool = DP.DPPool(args.workers, cpus)
    pf = TL.DataPrefetch(rows[start * args.batch:], labs, treatment)
    pf.start()
    log = (out / 'training.jsonl').open('a', buffering=1)
    print(json.dumps(dict(event='start', arm=args.arm, label=args.label, start_step=start, steps=args.steps, head_parameters=CFG['head_parameters'],
                          lora_parameters=CFG['lora_parameters'], setup_s=time.perf_counter() - began, label_load_s=label_load_s,
                          state_count=len(GR.state_table(128, il.key())[0]), theta=[float(x) for x in theta], **describe(theta, r_il),
                          eligible={k: v.n_eligible for k, v in labs.items()}, labeled={k: v.n_labeled for k, v in labs.items()})), flush=True)
    params = list(head.parameters()) + lora_params
    end = args.stop_at or args.steps
    ch = CFG['chunk']
    G = M.G

    def fwd(it, idx, grad):
        return BK.window_pots(head, ctx, it, idx, grad)

    for step in range(start, end):
        t_step = time.perf_counter()
        lr_h, lr_l = TL.lr_pair(step, args.steps)
        opt.param_groups[0]['lr'] = lr_h
        if len(opt.param_groups) > 1:
            opt.param_groups[1]['lr'] = lr_l
        lr_i = M.lr_at(step, args.steps, args.lr_intron, .1 * args.lr_intron, CFG['warmup_frac'])                            # [M41]
        opt_il.param_groups[0]['lr'] = lr_i
        bonus_np = P.bonus(theta.detach(), r_il, p0).numpy().copy()                                                            # [M41]
        g_bonus = np.zeros_like(bonus_np)
        opt.zero_grad(set_to_none=True)
        opt_il.zero_grad(set_to_none=True)
        items, futs, pots = [], [], []
        t_wait = 0.
        for _ in range(args.batch):
            t0 = time.perf_counter()
            it = pf.q.get()
            t_wait += time.perf_counter() - t0
            if it is None:
                raise RuntimeError('prefetch failed') from pf.err
            n = len(it['seqs'])
            with torch.no_grad():
                p = []
                for a in range(0, n, 4):
                    p += [x.detach() for x in fwd(it, list(range(a, min(a + 4, n))), False)]
            ctx_logits = G.owned_context_logits(p, it['geo'], it['w']['strand']).cpu().numpy()
            cs, ce = it['core']
            futs.append(pool.submit(it['context'], ctx_logits.astype(np.float32), cs, ce, it['targets'], CFG['epsilon'], bonus_np, r_il)
                        if it['counts'][0] + it['counts'][1] > 0 else None)
            items.append(it); pots.append(p)
        total_loss, dp_s, t_bwd, skipped, infeasible = 0., 0., 0., 0, 0
        by_sp = {}
        for it, fut, p in zip(items, futs, pots):
            sp = it['draw']['species']
            if fut is None:
                skipped += 1; by_sp.setdefault(sp, []).append(0.); continue
            try:
                r = fut.get()
            except InfeasibleObservation:
                infeasible += 1; continue
            n_known = max(1, it['counts'][0] + it['counts'][1])
            grad = torch.from_numpy(r['grad']).to('cuda') * (1. / n_known / args.batch)
            g_bonus += r['grad_bonus'] * (1. / n_known / args.batch)                                                          # [M41]
            t0 = time.perf_counter()
            n = len(it['seqs'])
            for a in range(0, n, ch):
                idx = list(range(a, min(a + ch, n)))
                new = fwd(it, idx, True)
                pieces = list(p)
                for i, x in zip(idx, new):
                    pieces[i] = x
                context = G.owned_context_logits(pieces, it['geo'], it['w']['strand'])
                (context * grad).sum().backward()
                del new, pieces, context
            torch.cuda.synchronize(); t_bwd += time.perf_counter() - t0
            total_loss += r['loss'] / n_known / args.batch
            dp_s += r['seconds']
            by_sp.setdefault(sp, []).append(r['loss'] / n_known)
        norm_g = torch.nn.utils.clip_grad_norm_(params, CFG['clip'], error_if_nonfinite=True)
        opt.step()
        if any(not bool(torch.isfinite(q).all()) for q in params):
            raise RuntimeError('nonfinite parameter')
        (P.bonus(theta, r_il, p0) * torch.from_numpy(g_bonus)).sum().backward()                                             # [M41]
        norm_i = float(torch.nn.utils.clip_grad_norm_([theta], 1., error_if_nonfinite=True))
        if args.intron_learn == 'on':
            opt_il.step()
        if not bool(torch.isfinite(theta).all()):
            raise RuntimeError('nonfinite intron-length parameter')
        row = dict(step=step + 1, lr_head=lr_h, lr_lora=lr_l, batch_loss=total_loss, grad_norm=float(norm_g), step_s=time.perf_counter() - t_step,
                   wait_prefetch_s=t_wait, bwd_s=t_bwd, dp_sum_s=dp_s, skipped_all_unknown=skipped, infeasible=infeasible,
                   nll_by_species={k: float(np.mean(v)) for k, v in by_sp.items()}, time=time.time(),
                   lr_intron=lr_i, intron_grad_norm=norm_i, theta=[float(x) for x in theta], **describe(theta, r_il))            # [M41]
        log.write(json.dumps(row, allow_nan=False) + '\n')
        if (step + 1) % 16 == 0 or step < 3:
            print(json.dumps({k: row[k] for k in ('step', 'lr_head', 'lr_lora', 'batch_loss', 'grad_norm', 'step_s', 'bwd_s', 'dp_sum_s', 'infeasible',
                                                  'pi_short', 'short_mean_bp', 'long_mean_bp')}), flush=True)
        done = step + 1
        if done % args.ckpt_every == 0 or done == end:
            st = save_resumable(ck, head, BK.lora_state(ctx), opt, done, CFG, done * args.batch, theta, opt_il, il_meta)
            if done % args.keep_every == 0 or done == end:
                save_model(out / f'model_{done:06d}.pt', st, r_il, p0)
    json.dump(dict(completed=True, arm=args.arm, label=args.label, steps=args.steps, reached=end, seconds=time.perf_counter() - began,
                   theta=[float(x) for x in theta], **describe(theta, r_il),
                   peak_gpu_alloc_GB=torch.cuda.max_memory_allocated() / 2**30, peak_gpu_reserved_GB=torch.cuda.max_memory_reserved() / 2**30),
              open(out / f'DONE_{end}.json', 'w'))
    pool.close()
    os._exit(0)


if __name__ == '__main__':
    main()
