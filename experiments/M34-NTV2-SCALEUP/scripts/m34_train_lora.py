"""M34-NTV2-SCALEUP multi-species trainer = the M31 training-B recipe (NT-v2-500M + q/v LoRA on blocks 13..27, hidden_states 21+27, C0 head, partition-function
CRF loss, parallel window DP) fed by a precomputed M34 draw schedule over 2/4/6/10 species.  Differences from `M31-C1-NTV2/scripts/c1_train_lora.py`
(everything else is the same code path, imported from M31 unchanged):
  * draws come from `outputs/M34-NTV2-SCALEUP/schedules/schedule_<arm>.jsonl` (protocol 19 §3.4) starting at draw 0 (the M31 A-stage offset does not apply);
  * targets are built by the M34 catalog labeler for every species with the arm's treatment: FULL (all assessable genes labeled), U (unselected eligible-gene
    loci unknown = marginalised), F (unselected genes filled with background);
  * 3,072-step cosine schedule, checkpoints every `--ckpt-every` (default 64) steps, model files every `--keep-every` (default 768).
Usage: m34_train_lora.py run --arm T10 --out DIR [--steps 3072] [--stop-at N] [--resume]
       m34_train_lora.py run --arm T10 --schedule FILE --out DIR --steps 3072 --stop-at 20      (smoke: explicit schedule file)
"""
import argparse
import hashlib
import json
import os
import queue
import threading
import time

import numpy as np
import torch

import m34_train_common as T
import m34_backend as BK
import m34_ckpt as CK
import a1_dp as DP
M = T.M
from src.m28.full_grammar_observed_partition import InfeasibleObservation   # noqa: E402

G, Q = M.G, M.Q
CFG = dict(seed=0, sampling_seed=1, lr_head_peak=1e-4, lr_lora_peak=3e-4, lr_final=1e-5, warmup_frac=.05, weight_decay=.01, betas=[.9, .999],
           eps=1e-8, clip=1., epsilon=.1, max_exons=128, checkpoint_bp=256, chunk=2)


class DataPrefetch(threading.Thread):
    """CPU-side example construction (strings, targets, geometry); the encoder runs in the main thread."""

    def __init__(self, rows, labelers, treatment, depth=16):
        super().__init__(daemon=True)
        self.q = queue.Queue(depth)
        self.args = (rows, labelers, treatment)
        self.err = None

    def run(self):
        rows, labs, treatment = self.args
        try:
            for r in rows:
                lab = labs[r['species']]
                chromosome, geometry, context, offset, targets, counts, win = T.example(lab, treatment, r)
                a, b = geometry['genomic_context']
                geo = BK.geometry(a, b)
                self.q.put(dict(draw=r, w=win, strand=r['strand'], geo=geo, context=context, core=(offset, offset + win['valid_bases']), counts=counts,
                                targets={k: np.asarray(targets[k]) for k in ('seg', 'region')}, **BK.item_windows(chromosome, geo, r['strand'])))
        except BaseException as e:
            self.err = e
            self.q.put(None)


def lr_pair(step, total):
    return (M.lr_at(step, total, CFG['lr_head_peak'], CFG['lr_final'], CFG['warmup_frac']),
            M.lr_at(step, total, CFG['lr_lora_peak'], CFG['lr_final'], CFG['warmup_frac']))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=('run',))
    ap.add_argument('--arm', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--schedule', default=None, help='explicit draw schedule file (smoke); default schedules/schedule_<arm>.jsonl')
    ap.add_argument('--init', default=str(T.INIT))
    ap.add_argument('--lora', choices=('on', 'none'), default='on', help='none: frozen encoder, head only (pc2 backend until the M35 LoRA module exists)')
    ap.add_argument('--lr-head', type=float, default=None)
    ap.add_argument('--lr-lora', type=float, default=None)
    ap.add_argument('--steps', type=int, default=3072)
    ap.add_argument('--stop-at', type=int, default=None)
    ap.add_argument('--batch', type=int, default=8)
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--cpus', default='')
    ap.add_argument('--ckpt-every', type=int, default=64)
    ap.add_argument('--keep-every', type=int, default=768)
    ap.add_argument('--resume', action='store_true')
    args = ap.parse_args()
    out = M.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cpus = M.parse_cpus(args.cpus)
    if cpus:
        os.sched_setaffinity(0, set(cpus))
    torch.set_num_threads(2)
    init = torch.load(args.init, map_location='cpu', weights_only=False)
    BK.check_init(init)
    if args.lr_head is not None:
        CFG['lr_head_peak'] = args.lr_head
    if args.lr_lora is not None:
        CFG['lr_lora_peak'] = args.lr_lora
    began = time.perf_counter()
    sched_path = M.Path(args.schedule) if args.schedule else T.SCHEDULES / f'schedule_{args.arm}.jsonl'
    rows = T.read_schedule(sched_path)
    assert len(rows) >= args.steps * args.batch, (len(rows), args.steps, args.batch)
    rows = rows[:args.steps * args.batch]
    forbidden = T.forbidden_contigs()
    assert not any(r['seqid'] in forbidden for r in rows), 'schedule contains DEV / held-out / validation contigs'
    spec = T.arm_spec(args.arm)
    assert {r['species'] for r in rows} <= set(spec['species'])
    sched_md5 = hashlib.md5(open(sched_path, 'rb').read()).hexdigest()
    CFG.update(BK.encoder_cfg(init['config']), lora_mode=args.lora, init_file=str(args.init), steps=args.steps, batch=args.batch, init_step=init['step'], arm=args.arm,
               species=list(spec['species']), new_species_label_fraction=spec['new_k'], treatment=spec['treatment'], schedule_md5=sched_md5,
               schedule=sched_path.name, heldout=list(T.VAL_REGIONS_LEGACY))
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
    start = 0
    ck = out / 'ckpt_latest.pt'
    if args.resume and ck.exists():
        start = CK.restore(ck, head, lambda st: BK.copy_lora(ctx, st), opt, CFG)
        print(json.dumps(dict(event='resumed', from_step=start, time=time.time())), flush=True)
    pool = DP.DPPool(args.workers, cpus)
    pf = DataPrefetch(rows[start * args.batch:], labs, treatment)
    pf.start()
    log = (out / 'training.jsonl').open('a', buffering=1)
    print(json.dumps(dict(event='start', arm=args.arm, start_step=start, steps=args.steps, head_parameters=CFG['head_parameters'],
                          lora_parameters=CFG['lora_parameters'], setup_s=time.perf_counter() - began, label_load_s=label_load_s,
                          eligible={k: v.n_eligible for k, v in labs.items()}, labeled={k: v.n_labeled for k, v in labs.items()})), flush=True)
    params = list(head.parameters()) + lora_params
    end = args.stop_at or args.steps
    ch = CFG['chunk']

    def fwd(it, idx, grad):
        return BK.window_pots(head, ctx, it, idx, grad)

    for step in range(start, end):
        t_step = time.perf_counter()
        lr_h, lr_l = lr_pair(step, args.steps)
        opt.param_groups[0]['lr'] = lr_h
        if len(opt.param_groups) > 1:
            opt.param_groups[1]['lr'] = lr_l
        opt.zero_grad(set_to_none=True)
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
            ctx = G.owned_context_logits(p, it['geo'], it['w']['strand']).cpu().numpy()
            cs, ce = it['core']
            futs.append(pool.submit(it['context'], ctx.astype(np.float32), cs, ce, it['targets'], CFG['epsilon']) if it['counts'][0] + it['counts'][1] > 0 else None)
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
        row = dict(step=step + 1, lr_head=lr_h, lr_lora=lr_l, batch_loss=total_loss, grad_norm=float(norm_g), step_s=time.perf_counter() - t_step,
                   wait_prefetch_s=t_wait, bwd_s=t_bwd, dp_sum_s=dp_s, skipped_all_unknown=skipped, infeasible=infeasible,
                   nll_by_species={k: float(np.mean(v)) for k, v in by_sp.items()}, time=time.time())
        log.write(json.dumps(row, allow_nan=False) + '\n')
        if (step + 1) % 16 == 0 or step < 3:
            print(json.dumps({k: row[k] for k in ('step', 'lr_head', 'lr_lora', 'batch_loss', 'grad_norm', 'step_s', 'bwd_s', 'dp_sum_s', 'infeasible')}), flush=True)
        done = step + 1
        if done % args.ckpt_every == 0 or done == end:
            st = CK.save_resumable(ck, head, BK.lora_state(ctx), opt, done, CFG, done * args.batch)
            if done % args.keep_every == 0 or done == end:
                CK.save_model(out / f'model_{done:06d}.pt', st)
    json.dump(dict(completed=True, arm=args.arm, steps=args.steps, reached=end, seconds=time.perf_counter() - began,
                   peak_gpu_alloc_GB=torch.cuda.max_memory_allocated() / 2**30, peak_gpu_reserved_GB=torch.cuda.max_memory_reserved() / 2**30),
              open(out / f'DONE_{end}.json', 'w'))
    pool.close()
    os._exit(0)


if __name__ == '__main__':
    main()
