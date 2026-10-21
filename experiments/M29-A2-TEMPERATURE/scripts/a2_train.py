"""M29-A2 training (temperature-scaled CRF target; copy of a1_train with: E_tau = tau*(offset+head), bias-free role layer, held-out chromosomes excluded from the sampler, --stop-at).
Original A1 doc: M29-A1 training: fixed Carbon-A offset + zero-initialised residual C0 head, full-grammar CRF NLL (log Zall - log Zobs).
Resumable: checkpoint every --ckpt-every steps (model, AdamW, step, RNG states, schedule position).  Usage:
  a1_train.py run  --out DIR --steps N [--batch 8] [--workers 12] [--cpus 2-31] [--ckpt-every 128] [--resume]
  a1_train.py smoke --out DIR --steps N ...        (same code, writes to DIR, no resume)
"""
import argparse, json, math, os, queue, re, sys, threading, time
import numpy as np
import torch
import a2_common as C2
import a1_common as K
import a2_common as C2
import a1_model as MD
import a2_model as MD2
import a1_offset as O
import a1_dp as DP
import fit_full_grammar_objective_contrast as C
import prepare_full_grammar_fresh_draws as PP
K.repoint()
G, Q = C.G, C.Q
from src.m28.full_grammar_observed_partition import InfeasibleObservation   # noqa: E402

CFG = dict(tau=None, heldout=list(C2.HELDOUT), seed=0, sampling_seed=1, lr_peak=1e-4, lr_final=1e-5, warmup_frac=.05, weight_decay=.01, betas=[.9, .999], eps=1e-8,
           clip=1., epsilon=.1, max_exons=128, checkpoint_bp=256, encoder='frozen GENERanno-1.2b revision of the frozen pipeline, joint 24,576 bp BF16',
           init='C0 seed0, role layer zero and bias-free', offset='A_full Carbon-A fixed offset (parity gate)')


def lr_at(step, total, cfg=CFG):
    if cfg.get('const_lr'):
        w = max(1, total // 10)
        return cfg['lr_peak'] * min(1., (step + 1) / w)
    w = max(1, int(math.ceil(cfg['warmup_frac'] * total)))
    if step < w:
        return cfg['lr_peak'] * (step + 1) / w
    t = (step - w) / max(1, total - w)
    return cfg['lr_final'] + .5 * (cfg['lr_peak'] - cfg['lr_final']) * (1 + math.cos(math.pi * t))


def example(w, sequences, catalogs):
    """Copy of frozen C.example minus the comparison with a stored plan row (extra draws beyond the frozen 4,608 have none)."""
    chromosome = sequences[w['species']][w['seqid']]
    geometry = G.context_geometry(w['start'], w['end'], len(chromosome))
    a, b = geometry['genomic_context']
    raw = chromosome[w['start']:w['end']]
    core = raw if w['strand'] == '+' else Q.M.rc(raw)
    targets = C.P.P.build_targets(w, catalogs[w['species']], core)
    n = w['valid_bases']
    fine, region = (torch.as_tensor(targets[k])[:n] for k in ('seg', 'region'))
    counts = (int((fine >= 0).sum()), int(((fine < 0) & (region >= 0)).sum()), int(((fine < 0) & (region < 0)).sum()))
    offset = w['start'] - a if w['strand'] == '+' else b - w['end']
    context = re.sub('[^ACGT]', 'N', chromosome[a:b])
    if w['strand'] == '-':
        context = Q.M.rc(context)
    return chromosome, geometry, context, offset, targets, counts


def schedule(total_draws, windows, strata, catalogs, original):
    passes = int(math.ceil(total_draws / len(original)))
    strata2 = dict(strata)
    for (species, seqid, strand) in list(strata):
        if seqid in C2.HELDOUT:      # internal-validation chromosomes are never sampled: reuse the same-species/strand pool of a training chromosome
            strata2[(species, seqid, strand)] = strata[(species, C2.SUBSTITUTE[seqid], strand)]
    return PP.fresh_schedule(original, strata2, passes=passes, seed=CFG['sampling_seed'])[:total_draws]


class Prefetch(threading.Thread):
    """Frozen-encoder features do not depend on the head weights, so encode future draws while the DP workers run."""
    def __init__(self, draws, windows, sequences, catalogs, tracks, backbone, tokenizer, depth=16):
        super().__init__(daemon=True)
        self.q = queue.Queue(depth)
        self.args = (draws, windows, sequences, catalogs, tracks, backbone, tokenizer)
        self.err = None

    def run(self):
        draws, windows, sequences, catalogs, tracks, backbone, tokenizer = self.args
        try:
            for d in draws:
                w = windows[d['window_id']]
                chromosome, geometry, context, offset, targets, counts = example(w, sequences, catalogs)
                a, b = geometry['genomic_context']
                seqs = MD.oriented_windows(chromosome, geometry, w['strand'])
                feats = MD.encode_bf16(backbone, tokenizer, seqs)
                off = O.context_offset(tracks[w['seqid']], w['strand'], a, b, len(chromosome))
                item = dict(draw=d, w=w, geometry=geometry, context=context, core=(offset, offset + w['valid_bases']), counts=counts,
                            targets={k: np.asarray(targets[k]) for k in ('seg', 'region')}, seqs=seqs, feats=feats, offset=off)
                self.q.put(item)
        except BaseException as e:      # surface in the main thread
            self.err = e
            self.q.put(None)


def atomic_save(obj, path):
    tmp = str(path) + '.tmp'
    torch.save(obj, tmp)
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=('run', 'smoke'))
    ap.add_argument('--out', required=True)
    ap.add_argument('--steps', type=int, required=True)
    ap.add_argument('--batch', type=int, default=8)
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--cpus', default='')
    ap.add_argument('--ckpt-every', type=int, default=128)
    ap.add_argument('--keep-every', type=int, default=256)
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--lr-peak', type=float, default=None)
    ap.add_argument('--tau', type=float, required=True)
    ap.add_argument('--stop-at', type=int, default=None)
    ap.add_argument('--const-lr', action='store_true', help='pilot only: 10%% linear warmup then constant lr')
    args = ap.parse_args()
    out = K.Path(args.out) if hasattr(K, 'Path') else __import__('pathlib').Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cpus = None
    if args.cpus:
        cpus = []
        for part in args.cpus.split(','):
            lo, _, hi = part.partition('-')
            cpus += list(range(int(lo), int(hi or lo) + 1))
        os.sched_setaffinity(0, set(cpus))
    torch.set_num_threads(2)
    total_draws = args.steps * args.batch
    CFG['tau'] = args.tau
    if args.lr_peak is not None:
        CFG['lr_peak'] = args.lr_peak
    if args.const_lr:
        CFG['const_lr'] = True
    began = time.perf_counter()
    windows, strata, catalogs, original = PP.load_population()
    sequences, _ = C.P.read_train_DNA()
    sched = schedule(total_draws, windows, strata, catalogs, original)
    assert not any(windows[d['window_id']]['seqid'] in C2.HELDOUT for d in sched), 'held-out chromosomes must not be trained on'
    tracks = {s: O.load_tracks(s) for s in {w['seqid'] for w in windows.values()}}
    model = MD2.new_model(CFG['seed']).to('cuda').train()
    opt = torch.optim.AdamW(model.parameters(), lr=CFG['lr_peak'], weight_decay=CFG['weight_decay'], betas=tuple(CFG['betas']), eps=CFG['eps'])
    start = 0
    ck = out / 'ckpt_latest.pt'
    if args.resume and ck.exists():
        st = torch.load(ck, map_location='cpu', weights_only=False)
        model.load_state_dict(st['model']); opt.load_state_dict(st['optimizer'])
        torch.set_rng_state(st['rng_torch']); np.random.set_state(st['rng_numpy'])
        start = st['step']
        assert st['config'] == dict(CFG, steps=args.steps, batch=args.batch), 'config changed since checkpoint'
        print(json.dumps(dict(event='resumed', from_step=start, time=time.time())), flush=True)
    backbone, tokenizer = G.load_native()
    pool = DP.DPPool(args.workers, cpus)
    pf = Prefetch(sched[start * args.batch:], windows, sequences, catalogs, tracks, backbone, tokenizer)
    pf.start()
    log = (out / 'training.jsonl').open('a', buffering=1)
    print(json.dumps(dict(event='start', start_step=start, steps=args.steps, batch=args.batch, setup_s=time.perf_counter() - began)), flush=True)
    params = list(model.parameters())
    end = args.stop_at or args.steps
    for step in range(start, end):
        t_step = time.perf_counter()
        lr = lr_at(step, args.steps)
        for g in opt.param_groups:
            g['lr'] = lr
        opt.zero_grad(set_to_none=True)
        items, futs, hmeans = [], [], []
        t_wait = 0.
        for _ in range(args.batch):
            t0 = time.perf_counter()
            it = pf.q.get()
            t_wait += time.perf_counter() - t0
            if it is None:
                raise RuntimeError('prefetch failed') from pf.err
            with torch.no_grad():
                hl = MD.head_context_logits(model, it['feats'], it['seqs'], it['geometry'], it['w']['strand']).cpu().numpy()
            logits = (CFG['tau'] * (hl + it['offset'])).astype(np.float32)
            hmeans.append(hl.mean(0))
            cs, ce = it['core']
            futs.append(pool.submit(it['context'], logits, cs, ce, it['targets'], CFG['epsilon']) if it['counts'][0] + it['counts'][1] > 0 else None)
            items.append(it)
        total_loss, rows, dp_s = 0., [], 0.
        t_head = 0.
        for it, fut in zip(items, futs):
            if fut is None:      # no fine/coarse observation: loss and gradient are exactly 0 (frozen verification; tested) -> skip the DP
                rows.append(dict(window=it['w']['id'], species=it['w']['species'], raw_NLL=0., known=0, skipped_all_unknown=True)); continue
            try:
                r = fut.get()
            except InfeasibleObservation:
                rows.append(dict(window=it['w']['id'], infeasible=True)); continue
            n_known = max(1, it['counts'][0] + it['counts'][1])
            mult = 1. / n_known / args.batch
            t0 = time.perf_counter()
            hl = MD.head_context_logits(model, it['feats'], it['seqs'], it['geometry'], it['w']['strand'])
            g = torch.from_numpy(r['grad']).to('cuda') * (mult * CFG['tau'])
            hl.backward(g)
            torch.cuda.synchronize(); t_head += time.perf_counter() - t0
            total_loss += r['loss'] / n_known / args.batch
            dp_s += r['seconds']
            rows.append(dict(window=it['w']['id'], species=it['w']['species'], raw_NLL=r['loss'], known=n_known,
                             normalized_NLL=r['loss'] / n_known, dp_s=r['seconds'], fine=it['counts'][0], coarse=it['counts'][1],
                             unknown=it['counts'][2]))
        norm = torch.nn.utils.clip_grad_norm_(params, CFG['clip'], error_if_nonfinite=True)
        opt.step()
        if any(not bool(torch.isfinite(p).all()) for p in params):
            raise RuntimeError('nonfinite parameter')
        row = dict(step=step + 1, lr=lr, batch_loss=total_loss, grad_norm=float(norm), step_s=time.perf_counter() - t_step,
                   head_mean=[float(x) for x in np.mean(hmeans, 0)], wait_prefetch_s=t_wait, head_bwd_s=t_head, dp_sum_s=dp_s, windows=rows, time=time.time())
        log.write(json.dumps(row, allow_nan=False) + '\n')
        if (step + 1) % 16 == 0 or step < 3:
            print(json.dumps({k: row[k] for k in ('step', 'lr', 'batch_loss', 'grad_norm', 'step_s', 'wait_prefetch_s', 'head_bwd_s', 'dp_sum_s')}), flush=True)
        done = step + 1
        if args.mode == 'run' and (done % args.ckpt_every == 0 or done == end):
            st = dict(model=model.state_dict(), optimizer=opt.state_dict(), step=done, rng_torch=torch.get_rng_state(),
                      rng_numpy=np.random.get_state(), config=dict(CFG, steps=args.steps, batch=args.batch),
                      draw_position=done * args.batch)
            atomic_save(st, ck)
            if done % args.keep_every == 0 or done == end:
                atomic_save(dict(model=model.state_dict(), step=done, config=st['config']), out / f'model_{done:06d}.pt')
    json.dump(dict(completed=True, steps=args.steps, reached=end, seconds=time.perf_counter() - began), open(out / (f'DONE_{end}.json'), 'w'))
    pool.close()
    os._exit(0)


if __name__ == '__main__':
    main()
