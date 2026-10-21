"""M31-C1-NTV2 training B: continue from the selected training-A head (step 2048) and add rank-r LoRA on query/value of NT-v2 blocks
(default blocks 13..27, i.e. 0-based encoder layers 12..26; the features are hidden_states 21 and 27, so blocks 28-29 are never run).
Same loss, sampler (draws continue after the 2048*8 draws training A consumed), DP workers and checkpointing as c1_train.py.
Because the encoder changes with the LoRA weights, features cannot be prefetched: per sample, phase 1 = no-grad encoder+head forward -> DP
(workers, in parallel); phase 2 = recompute with grad in chunks of 2 windows and back-propagate the DP gradient (gradient-checkpoint style).
Usage: c1_train_lora.py run --out DIR --init runA/model_002048.pt --steps N [--stop-at K] [--resume]
"""
import argparse, json, os, queue, threading, time
import numpy as np
import torch

import ntv2_adapter as A
import m31_common as M
import c1_model as MD
import a1_dp as DP
from src.m28.full_grammar_observed_partition import InfeasibleObservation   # noqa: E402

G, Q = M.G, M.Q
A_STEPS_USED = 2048                      # draws 0 .. 2048*8-1 were used by training A; B continues the same schedule
CFG = dict(seed=0, sampling_seed=1, lr_head_peak=1e-4, lr_lora_peak=3e-4, lr_final=1e-5, warmup_frac=.05, weight_decay=.01, betas=[.9, .999],
           eps=1e-8, clip=1., epsilon=.1, max_exons=128, checkpoint_bp=256, encoder=A.MODEL, nt_width=A.NT_WIDTH, nt_stride=A.NT_STRIDE,
           init_from='runA/model_002048.pt (selected by internal validation)', heldout=list(M.HELDOUT), lora_r=8, lora_alpha=16,
           lora_encoder_layers=list(range(12, 27)), chunk=2)


class DataPrefetch(threading.Thread):
    """CPU-side example construction (strings, targets, geometry); the encoder runs in the main thread."""
    def __init__(self, draws, windows, sequences, catalogs, depth=16):
        super().__init__(daemon=True)
        self.q = queue.Queue(depth)
        self.args = (draws, windows, sequences, catalogs)
        self.err = None

    def run(self):
        draws, windows, sequences, catalogs = self.args
        try:
            for d in draws:
                w = windows[d['window_id']]
                chromosome, geometry, context, offset, targets, counts = M.example(w, sequences, catalogs)
                a, b = geometry['genomic_context']
                geo = A.nt_geometry(a, b)
                self.q.put(dict(draw=d, w=w, geo=geo, context=context, core=(offset, offset + w['valid_bases']), counts=counts,
                                targets={k: np.asarray(targets[k]) for k in ('seg', 'region')},
                                seqs=MD.oriented_windows(chromosome, geo, w['strand'])))
        except BaseException as e:
            self.err = e
            self.q.put(None)


def lr_pair(step, total):
    return (M.lr_at(step, total, CFG['lr_head_peak'], CFG['lr_final'], CFG['warmup_frac']),
            M.lr_at(step, total, CFG['lr_lora_peak'], CFG['lr_final'], CFG['warmup_frac']))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=('run', 'smoke'))
    ap.add_argument('--out', required=True)
    ap.add_argument('--init', required=True)
    ap.add_argument('--steps', type=int, required=True)
    ap.add_argument('--stop-at', type=int, default=None)
    ap.add_argument('--batch', type=int, default=8)
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--cpus', default='')
    ap.add_argument('--ckpt-every', type=int, default=64)
    ap.add_argument('--keep-every', type=int, default=512)
    ap.add_argument('--resume', action='store_true')
    args = ap.parse_args()
    out = M.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cpus = M.parse_cpus(args.cpus)
    if cpus:
        os.sched_setaffinity(0, set(cpus))
    torch.set_num_threads(2)
    init = torch.load(args.init, map_location='cpu', weights_only=False)
    layers, norm = tuple(init['config']['layers']), init['config']['layer_norm_per_layer']
    CFG.update(layers=list(layers), layer_norm_per_layer=norm, steps=args.steps, batch=args.batch, init_step=init['step'])
    assert init['step'] == A_STEPS_USED and max(layers) <= 27
    began = time.perf_counter()
    windows, strata, catalogs, original = M.PP.load_population()
    sequences, _ = M.C.P.read_train_DNA()
    sched = M.schedule((A_STEPS_USED + args.steps) * args.batch, strata, original, CFG['sampling_seed'])[A_STEPS_USED * args.batch:]
    assert not any(windows[d['window_id']]['seqid'] in M.HELDOUT for d in sched)
    head = MD.new_head(layers, CFG['seed'])
    head.load_state_dict(init['model'])
    head = head.to('cuda').train()
    esm, vocab = A.load_ntv2()
    lora_params = A.add_lora(esm, CFG['lora_encoder_layers'], CFG['lora_r'], CFG['lora_alpha'])
    CFG['head_parameters'] = sum(p.numel() for p in head.parameters())
    CFG['lora_parameters'] = sum(p.numel() for p in lora_params)
    opt = torch.optim.AdamW([dict(params=list(head.parameters()), name='head'), dict(params=lora_params, name='lora')],
                            lr=CFG['lr_head_peak'], weight_decay=CFG['weight_decay'], betas=tuple(CFG['betas']), eps=CFG['eps'])
    start = 0
    ck = out / 'ckpt_latest.pt'
    if args.resume and ck.exists():
        st = torch.load(ck, map_location='cpu', weights_only=False)
        head.load_state_dict(st['model'])
        for k, v in st['lora'].items():
            i, n, w = k.split('.')
            getattr(getattr(esm.encoder.layer[int(i)].attention.self, n), w).data.copy_(v)
        opt.load_state_dict(st['optimizer'])
        torch.set_rng_state(st['rng_torch']); np.random.set_state(st['rng_numpy'])
        start = st['step']
        assert st['config'] == CFG, 'config changed since checkpoint'
        print(json.dumps(dict(event='resumed', from_step=start, time=time.time())), flush=True)
    pool = DP.DPPool(args.workers, cpus)
    pf = DataPrefetch(sched[start * args.batch:], windows, sequences, catalogs)
    pf.start()
    log = (out / 'training.jsonl').open('a', buffering=1)
    print(json.dumps(dict(event='start', start_step=start, steps=args.steps, head_parameters=CFG['head_parameters'], lora_parameters=CFG['lora_parameters'],
                          setup_s=time.perf_counter() - began)), flush=True)
    params = list(head.parameters()) + lora_params
    end = args.stop_at or args.steps
    ch = CFG['chunk']

    def fwd(seqs, idx, grad):
        feats = A.encode_windows(esm, vocab, [seqs[i] for i in idx], layers, batch=len(idx), norm=norm, custom=True, grad=grad)
        return MD.window_potentials(head, feats, [seqs[i] for i in idx])

    for step in range(start, end):
        t_step = time.perf_counter()
        lr_h, lr_l = lr_pair(step, args.steps)
        opt.param_groups[0]['lr'], opt.param_groups[1]['lr'] = lr_h, lr_l
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
                    p += [x.detach() for x in fwd(it['seqs'], list(range(a, min(a + 4, n))), False)]
            ctx = G.owned_context_logits(p, it['geo'], it['w']['strand']).cpu().numpy()
            cs, ce = it['core']
            futs.append(pool.submit(it['context'], ctx.astype(np.float32), cs, ce, it['targets'], CFG['epsilon']) if it['counts'][0] + it['counts'][1] > 0 else None)
            items.append(it); pots.append(p)
        total_loss, dp_s, t_bwd, skipped, infeasible = 0., 0., 0., 0, 0
        by_sp = {}
        for it, fut, p in zip(items, futs, pots):
            sp = it['w']['species']
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
                new = fwd(it['seqs'], idx, True)
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
            print(json.dumps({k: row[k] for k in ('step', 'lr_head', 'lr_lora', 'batch_loss', 'grad_norm', 'step_s', 'bwd_s', 'dp_sum_s')}), flush=True)
        done = step + 1
        if args.mode == 'run' and (done % args.ckpt_every == 0 or done == end):
            st = dict(model=head.state_dict(), lora=A.lora_state(esm), optimizer=opt.state_dict(), step=done, rng_torch=torch.get_rng_state(),
                      rng_numpy=np.random.get_state(), config=dict(CFG), draw_position=(A_STEPS_USED + done) * args.batch)
            M.atomic_save(st, ck)
            if done % args.keep_every == 0 or done == end:
                M.atomic_save(dict(model=head.state_dict(), lora=st['lora'], step=done, config=st['config']), out / f'model_{done:06d}.pt')
    json.dump(dict(completed=True, steps=args.steps, reached=end, seconds=time.perf_counter() - began,
                   peak_gpu_alloc_GB=torch.cuda.max_memory_allocated() / 2**30, peak_gpu_reserved_GB=torch.cuda.max_memory_reserved() / 2**30),
              open(out / f'DONE_{end}.json', 'w'))
    pool.close()
    os._exit(0)


if __name__ == '__main__':
    main()
