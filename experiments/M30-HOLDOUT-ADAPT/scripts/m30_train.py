"""M30 partial-label adaptation of the Carbon-A-offset residual head (protocol 13 §3-4).  Copy of a2_train with: new-species labeler/windows/sampler.
  m30_train.py run --species S --k 0.25 --treatment U|F|FULL --out DIR --steps 1024 [--stop-at N] [--resume]
"""
import argparse, json, math, os, queue, threading, time
import numpy as np
import torch
import m30_common as M
import a1_common as K
import a1_model as MD
import a2_model as MD2
import a1_dp as DP
import a1_offset as O
import m30_labels as LB
K.repoint()
G, Q = MD.G, MD.Q
import re
from src.m28.full_grammar_observed_partition import InfeasibleObservation

CFG = dict(tau=0.03, seed=0, sampling_seed=1, lr_peak=1e-4, lr_final=1e-5, warmup_frac=.05, weight_decay=.01, betas=[.9, .999], eps=1e-8,
           clip=1., epsilon=.1, max_exons=128, init='C0 seed0, role layer zero and bias-free', offset='A_full Carbon-A fixed offset (parity gate)',
           uniform_fraction=.5)


def lr_at(step, total, cfg=CFG):
    w = max(1, int(math.ceil(cfg['warmup_frac'] * total)))
    if step < w:
        return cfg['lr_peak'] * (step + 1) / w
    t = (step - w) / max(1, total - w)
    return cfg['lr_final'] + .5 * (cfg['lr_peak'] - cfg['lr_final']) * (1 + math.cos(math.pi * t))


def make_schedule(labeler, total):
    """Protocol §4 sampler: 50% uniform over (cell,strand) of the training contigs, 50% labeled-gene-first.  Depends only on (species, k)."""
    rng = np.random.default_rng(CFG['sampling_seed'])
    cells = [(c, s, a, b) for c in labeler.contigs for s in '+-' for a, b in labeler.grid(c)]
    gene_cells = labeler.labeled_windows()
    genes = sorted(gene_cells)
    out = []
    for _ in range(total):
        if rng.random() < CFG['uniform_fraction'] or not genes:
            out.append(cells[int(rng.integers(len(cells)))])
        else:
            g = genes[int(rng.integers(len(genes)))]
            cl = gene_cells[g]
            out.append(cl[int(rng.integers(len(cl)))])
    return out


def example(labeler, treatment, cell):
    c, s, a, b = cell
    tg, win, core = labeler.targets(c, s, a, b, treatment)
    chromosome = labeler.seqs[c]
    geometry = G.context_geometry(a, b, len(chromosome))
    ca, cb = geometry['genomic_context']
    n = b - a
    seg, region = tg['seg'][:n], tg['region'][:n]
    fine = int((seg >= 0).sum()); coarse = int(((seg < 0) & (region >= 0)).sum())
    offset = a - ca if s == '+' else cb - b
    context = re.sub('[^ACGT]', 'N', chromosome[ca:cb])
    if s == '-':
        context = LB.MF.rc(context)
    return chromosome, geometry, context, offset, tg, (fine, coarse, n - fine - coarse), win


class Prefetch(threading.Thread):
    def __init__(self, cells, labeler, treatment, tracks, backbone, tokenizer, depth=16):
        super().__init__(daemon=True)
        self.q = queue.Queue(depth); self.err = None
        self.args = (cells, labeler, treatment, tracks, backbone, tokenizer)

    def run(self):
        cells, labeler, treatment, tracks, backbone, tokenizer = self.args
        try:
            for cell in cells:
                chromosome, geometry, context, offset, tg, counts, win = example(labeler, treatment, cell)
                ca, cb = geometry['genomic_context']
                seqs = MD.oriented_windows(chromosome, geometry, cell[1])
                feats = MD.encode_bf16(backbone, tokenizer, seqs) if counts[0] + counts[1] > 0 else None
                off = O.context_offset(tracks[cell[0]], cell[1], ca, cb, len(chromosome)) if feats is not None else None
                self.q.put(dict(cell=cell, w=win, geometry=geometry, context=context, core=(offset, offset + win['valid_bases']), counts=counts,
                                targets={k: np.asarray(tg[k]) for k in ('seg', 'region')}, seqs=seqs, feats=feats, offset=off))
        except BaseException as e:
            self.err = e
            self.q.put(None)


def atomic_save(obj, path):
    tmp = str(path) + '.tmp'
    torch.save(obj, tmp)
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--species', required=True)
    ap.add_argument('--k', type=float, required=True)
    ap.add_argument('--treatment', choices=('U', 'F', 'FULL'), required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--steps', type=int, default=1024)
    ap.add_argument('--stop-at', type=int, default=None)
    ap.add_argument('--batch', type=int, default=8)
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--cpus', default='')
    ap.add_argument('--ckpt-every', type=int, default=256)
    ap.add_argument('--resume', action='store_true')
    args = ap.parse_args()
    out = K.Path(args.out); out.mkdir(parents=True, exist_ok=True)
    cpus = None
    if args.cpus:
        cpus = []
        for part in args.cpus.split(','):
            lo, _, hi = part.partition('-'); cpus += list(range(int(lo), int(hi or lo) + 1))
        os.sched_setaffinity(0, set(cpus))
    torch.set_num_threads(2)
    began = time.perf_counter()
    labeler = LB.Labeler(args.species, M.SPECIES[args.species]['train'], args.k)
    cells = make_schedule(labeler, args.steps * args.batch)
    tracks = {c: {s: np.load(M.offset_path(c, s), mmap_mode='r') for s in '+-'} for c in labeler.contigs}
    cfg = dict(CFG, species=args.species, k=args.k, treatment=args.treatment, steps=args.steps, batch=args.batch, n_eligible=labeler.n_eligible, n_labeled=labeler.n_labeled)
    model = MD2.new_model(CFG['seed']).to('cuda').train()
    opt = torch.optim.AdamW(model.parameters(), lr=CFG['lr_peak'], weight_decay=CFG['weight_decay'], betas=tuple(CFG['betas']), eps=CFG['eps'])
    start, ck = 0, out / 'ckpt_latest.pt'
    if args.resume and ck.exists():
        st = torch.load(ck, map_location='cpu', weights_only=False)
        assert st['config'] == cfg, 'config changed since checkpoint'
        model.load_state_dict(st['model']); opt.load_state_dict(st['optimizer'])
        torch.set_rng_state(st['rng_torch']); np.random.set_state(st['rng_numpy'])
        start = st['step']
        print(json.dumps(dict(event='resumed', from_step=start, time=time.time())), flush=True)
    backbone, tokenizer = G.load_native()
    pool = DP.DPPool(args.workers, cpus)
    pf = Prefetch(cells[start * args.batch:], labeler, args.treatment, tracks, backbone, tokenizer)
    pf.start()
    log = (out / 'training.jsonl').open('a', buffering=1)
    print(json.dumps(dict(event='start', start_step=start, cfg=cfg, setup_s=time.perf_counter() - began)), flush=True)
    params = list(model.parameters())
    end = args.stop_at or args.steps
    tau = CFG['tau']
    for step in range(start, end):
        t_step = time.perf_counter()
        lr = lr_at(step, args.steps)
        for g in opt.param_groups: g['lr'] = lr
        opt.zero_grad(set_to_none=True)
        items, futs, hmeans = [], [], []
        for _ in range(args.batch):
            it = pf.q.get()
            if it is None:
                raise RuntimeError('prefetch failed') from pf.err
            if it['feats'] is not None:
                with torch.no_grad():
                    hl = MD.head_context_logits(model, it['feats'], it['seqs'], it['geometry'], it['cell'][1]).cpu().numpy()
                logits = (tau * (hl + it['offset'])).astype(np.float32)
                hmeans.append(hl.mean(0))
                futs.append(pool.submit(it['context'], logits, it['core'][0], it['core'][1], it['targets'], CFG['epsilon']))
            else:
                futs.append(None)
            items.append(it)
        total_loss, rows, dp_s = 0., [], 0.
        for it, fut in zip(items, futs):
            if fut is None:
                rows.append(dict(window=it['w']['id'], raw_NLL=0., known=0, skipped_no_observation=True)); continue
            try:
                r = fut.get()
            except InfeasibleObservation:      # label path not representable in the frozen grammar (e.g. >128 exons): skip window, logged
                rows.append(dict(window=it['w']['id'], raw_NLL=0., known=0, infeasible=True)); continue
            n_known = max(1, it['counts'][0] + it['counts'][1])
            hl = MD.head_context_logits(model, it['feats'], it['seqs'], it['geometry'], it['cell'][1])
            hl.backward(torch.from_numpy(r['grad']).to('cuda') * (tau / n_known / args.batch))
            total_loss += r['loss'] / n_known / args.batch
            dp_s += r['seconds']
            rows.append(dict(window=it['w']['id'], raw_NLL=r['loss'], known=n_known, normalized_NLL=r['loss'] / n_known, dp_s=r['seconds'],
                             fine=it['counts'][0], coarse=it['counts'][1], unknown=it['counts'][2]))
        torch.cuda.synchronize()
        norm = torch.nn.utils.clip_grad_norm_(params, CFG['clip'], error_if_nonfinite=True)
        opt.step()
        if any(not bool(torch.isfinite(p).all()) for p in params):
            raise RuntimeError('nonfinite parameter')
        row = dict(step=step + 1, lr=lr, batch_loss=total_loss, grad_norm=float(norm), step_s=time.perf_counter() - t_step,
                   head_mean=[float(x) for x in np.mean(hmeans, 0)] if hmeans else None, dp_sum_s=dp_s, windows=rows, time=time.time())
        log.write(json.dumps(row, allow_nan=False) + '\n')
        if (step + 1) % 16 == 0 or step < 3:
            print(json.dumps({k: row[k] for k in ('step', 'lr', 'batch_loss', 'grad_norm', 'step_s', 'dp_sum_s')}), flush=True)
        done = step + 1
        if done % args.ckpt_every == 0 or done == end:
            st = dict(model=model.state_dict(), optimizer=opt.state_dict(), step=done, rng_torch=torch.get_rng_state(), rng_numpy=np.random.get_state(),
                      config=cfg, draw_position=done * args.batch)
            atomic_save(st, ck)
            atomic_save(dict(model=model.state_dict(), step=done, config=cfg), out / f'model_{done:06d}.pt')
    json.dump(dict(completed=True, reached=end, seconds=time.perf_counter() - began), open(out / f'DONE_{end}.json', 'w'))
    pool.close()
    os._exit(0)


if __name__ == '__main__':
    main()
