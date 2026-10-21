"""M32-C1-PLANTCAD2 training: frozen PlantCAD2-Small (self-supervised MLM only; windows 8,192 bp / stride 6,144) + C0 full-grammar role head, partition-function CRF loss
(log Zall - log Zobs, fine hard constraints, coarse observation channel eps=.1, unknown positions marginalised), no Carbon-A / GENERanno anything.
Same recipe, head, targets, sampler and internal validation as M31-C1-NTV2 (only the frozen encoder and its window geometry differ).

Same recipe as the M29 parallel trainer: parallel DP workers (bit-identical to the frozen serial DP, tested), prefetching frozen encoder,
batch 8, AdamW + cosine (+5 % linear warm-up), grad clip 1, resumable checkpoints (<= 256 steps), compact jsonl.
Internal-validation chromosomes (NC_003071.7, NC_089044.1) are never sampled.
Usage: c1_train.py {run|smoke} --out DIR --steps N [--stop-at K] [--layer -1] [--resume]
"""
import argparse, json, os, queue, threading, time
import numpy as np
import torch

import pc2_adapter as A
import m32_common as M
import c1_model as MD
import a1_dp as DP                                   # frozen-equivalent parallel DP (M29-A1), imported only
from src.m28.full_grammar_observed_partition import InfeasibleObservation   # noqa: E402

G, Q = M.G, M.Q
CFG = dict(seed=0, sampling_seed=1, lr_peak=3e-4, lr_final=1e-5, warmup_frac=.05, weight_decay=.01, betas=[.9, .999], eps=1e-8, clip=1.,
           epsilon=.1, max_exons=128, checkpoint_bp=256, encoder=A.MODEL, pc_width=A.PC_WIDTH, pc_stride=A.PC_STRIDE, head_width=A.HEAD_W, token_pool='mean over valid bases of aligned 6-base groups',
           init='C0 default PyTorch init, seed 0', loss='coarse_channel_nll/(max(1,fine+coarse)) per window / batch', heldout=list(M.HELDOUT))


class Prefetch(threading.Thread):
    """Frozen-encoder features do not depend on the head, so encode future draws (~8 PlantCAD2 windows of the 49 kb context) while the DP runs."""
    def __init__(self, draws, windows, sequences, catalogs, enc, layer, depth=8):
        super().__init__(daemon=True)
        self.q = queue.Queue(depth)
        self.args = (draws, windows, sequences, catalogs, enc, layer)
        self.err = None
        self.enc_s = 0.        # cumulative encoder wall time of the prefetch thread (GPU shared with the main thread's head work)
        self.n = 0

    def run(self):
        draws, windows, sequences, catalogs, enc, layer = self.args
        try:
            for d in draws:
                w = windows[d['window_id']]
                chromosome, geometry, context, offset, targets, counts = M.example(w, sequences, catalogs)
                a, b = geometry['genomic_context']
                geo = A.pc_geometry(a, b)
                gen = MD.genomic_window_strings(chromosome, geo)
                seqs = MD.orient_strings(gen, w['strand'])
                if counts[0] + counts[1] > 0:      # all-unknown draws have loss and gradient exactly 0 and are skipped by the trainer: no encoding
                    t0 = time.perf_counter()
                    feats = A.oriented_tokens(A.encode_windows(enc, gen, layer, batch=8), w['strand'])
                    torch.cuda.synchronize()
                    self.enc_s += time.perf_counter() - t0
                    self.n += 1
                else:
                    feats = None
                self.q.put(dict(draw=d, w=w, geo=geo, context=context, core=(offset, offset + w['valid_bases']), counts=counts,
                                targets={k: np.asarray(targets[k]) for k in ('seg', 'region')}, seqs=seqs, feats=feats))
        except BaseException as e:      # surface in the main thread
            self.err = e
            self.q.put(None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=('run', 'smoke'))
    ap.add_argument('--out', required=True)
    ap.add_argument('--steps', type=int, required=True)
    ap.add_argument('--stop-at', type=int, default=None)
    ap.add_argument('--batch', type=int, default=8)
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--cpus', default='')
    ap.add_argument('--ckpt-every', type=int, default=128)
    ap.add_argument('--keep-every', type=int, default=512)
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--lr-peak', type=float, default=None)
    ap.add_argument('--layer', type=int, default=-1, help='-1 = last_hidden_state (final norm), else hidden_states index')
    args = ap.parse_args()
    out = M.Path(args.out) if hasattr(M, 'Path') else __import__('pathlib').Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cpus = M.parse_cpus(args.cpus)
    if cpus:
        os.sched_setaffinity(0, set(cpus))
    torch.set_num_threads(2)
    layer = args.layer
    CFG.update(layer=layer, steps=args.steps, batch=args.batch)
    if args.lr_peak is not None:
        CFG['lr_peak'] = args.lr_peak
    total_draws = args.steps * args.batch
    began = time.perf_counter()
    windows, strata, catalogs, original = M.PP.load_population()
    sequences, _ = M.C.P.read_train_DNA()
    sched = M.schedule(total_draws, strata, original, CFG['sampling_seed'])
    assert not any(windows[d['window_id']]['seqid'] in M.HELDOUT for d in sched), 'held-out chromosomes must not be trained on'
    model = MD.new_head(CFG['seed']).to('cuda').train()
    CFG['head_parameters'] = sum(p.numel() for p in model.parameters())
    opt = torch.optim.AdamW(model.parameters(), lr=CFG['lr_peak'], weight_decay=CFG['weight_decay'], betas=tuple(CFG['betas']), eps=CFG['eps'])
    start = 0
    ck = out / 'ckpt_latest.pt'
    if args.resume and ck.exists():
        st = torch.load(ck, map_location='cpu', weights_only=False)
        model.load_state_dict(st['model']); opt.load_state_dict(st['optimizer'])
        torch.set_rng_state(st['rng_torch']); np.random.set_state(st['rng_numpy'])
        start = st['step']
        assert st['config'] == CFG, 'config changed since checkpoint'
        print(json.dumps(dict(event='resumed', from_step=start, time=time.time())), flush=True)
    enc = A.load_pc2()
    pool = DP.DPPool(args.workers, cpus)
    pf = Prefetch(sched[start * args.batch:], windows, sequences, catalogs, enc, layer)
    pf.start()
    log = (out / 'training.jsonl').open('a', buffering=1)
    print(json.dumps(dict(event='start', start_step=start, steps=args.steps, batch=args.batch, layer=layer, head_parameters=CFG['head_parameters'],
                          setup_s=time.perf_counter() - began)), flush=True)
    params = list(model.parameters())
    end = args.stop_at or args.steps
    bases_seen = 0
    t_run = time.perf_counter()
    for step in range(start, end):
        t_step = time.perf_counter()
        lr = M.lr_at(step, args.steps, CFG['lr_peak'], CFG['lr_final'], CFG['warmup_frac'])
        for g in opt.param_groups:
            g['lr'] = lr
        opt.zero_grad(set_to_none=True)
        items, futs = [], []
        t_wait = 0.
        for _ in range(args.batch):
            t0 = time.perf_counter()
            it = pf.q.get()
            t_wait += time.perf_counter() - t0
            if it is None:
                raise RuntimeError('prefetch failed') from pf.err
            if it['feats'] is not None:
                with torch.no_grad():
                    hl = MD.head_context_logits(model, it['feats'], it['seqs'], it['geo'], it['w']['strand']).cpu().numpy()
                futs.append(DP_submit(pool, it, hl))
            else:
                futs.append(None)
            items.append(it)
        total_loss, dp_s, t_head, skipped, infeasible = 0., 0., 0., 0, 0
        by_sp = {}
        for it, fut in zip(items, futs):
            sp = it['w']['species']
            if fut is None:      # no fine/coarse observation: loss and gradient are exactly 0 (frozen verification) -> skip the DP
                skipped += 1; by_sp.setdefault(sp, []).append(0.); continue
            try:
                r = fut.get()
            except InfeasibleObservation:
                infeasible += 1; continue
            n_known = max(1, it['counts'][0] + it['counts'][1])
            t0 = time.perf_counter()
            hl = MD.head_context_logits(model, it['feats'], it['seqs'], it['geo'], it['w']['strand'])
            hl.backward(torch.from_numpy(r['grad']).to('cuda') * (1. / n_known / args.batch))
            torch.cuda.synchronize(); t_head += time.perf_counter() - t0
            total_loss += r['loss'] / n_known / args.batch
            dp_s += r['seconds']
            by_sp.setdefault(sp, []).append(r['loss'] / n_known)
            bases_seen += it['w']['valid_bases']
        norm_g = torch.nn.utils.clip_grad_norm_(params, CFG['clip'], error_if_nonfinite=True)
        opt.step()
        if any(not bool(torch.isfinite(p).all()) for p in params):
            raise RuntimeError('nonfinite parameter')
        row = dict(step=step + 1, lr=lr, batch_loss=total_loss, grad_norm=float(norm_g), step_s=time.perf_counter() - t_step,
                   wait_prefetch_s=t_wait, enc_s_per_draw=pf.enc_s / max(1, pf.n), head_bwd_s=t_head, dp_sum_s=dp_s, skipped_all_unknown=skipped, infeasible=infeasible,
                   nll_by_species={k: float(np.mean(v)) for k, v in by_sp.items()}, time=time.time())
        log.write(json.dumps(row, allow_nan=False) + '\n')
        if (step + 1) % 16 == 0 or step < 3:
            print(json.dumps({k: row[k] for k in ('step', 'lr', 'batch_loss', 'grad_norm', 'step_s', 'wait_prefetch_s', 'enc_s_per_draw', 'head_bwd_s', 'dp_sum_s')}), flush=True)
        done = step + 1
        if args.mode == 'run' and (done % min(args.ckpt_every, 64) == 0 or done == end):     # GPU is shared (OOM kills happen): checkpoint every 64 steps
            st = dict(model=model.state_dict(), optimizer=opt.state_dict(), step=done, rng_torch=torch.get_rng_state(),
                      rng_numpy=np.random.get_state(), config=dict(CFG), draw_position=done * args.batch)
            M.atomic_save(st, ck)
            if done % args.keep_every == 0 or done == end:
                M.atomic_save(dict(model=model.state_dict(), step=done, config=st['config']), out / f'model_{done:06d}.pt')
    train_s = time.perf_counter() - t_run
    json.dump(dict(completed=True, steps=args.steps, reached=end, seconds=time.perf_counter() - began, train_seconds=train_s,
                   train_kbp_per_s=bases_seen / 1000 / train_s, peak_gpu_alloc_GB=torch.cuda.max_memory_allocated() / 2**30,
                   peak_gpu_reserved_GB=torch.cuda.max_memory_reserved() / 2**30), open(out / f'DONE_{end}.json', 'w'))
    pool.close()
    os._exit(0)


def DP_submit(pool, it, head_logits):
    cs, ce = it['core']
    return pool.submit(it['context'], head_logits.astype(np.float32), cs, ce, it['targets'], CFG['epsilon'])


if __name__ == '__main__':
    main()
