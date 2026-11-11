"""M40-SEMISUP-MLM trainer = the M34 pc2 trainer (`M34-NTV2-SCALEUP/scripts/m34_train_lora.py`: PlantCAD2-Small + M35 LoRA, C0 head warm-started from
M32 runA, partition-function CRF loss with unknown positions marginalised, parallel window DP, 3,072-step cosine schedule) + a masked-language-model
auxiliary loss on the unlabeled positions of the new species' draws (protocol docs/review_20261016/34_semisup_MLM_protocol.md).

What is identical to M34 (imported, not copied): labelers / targets / example geometry (m34_train_common), encoder + LoRA + head (m34_backend, M35 pc2_lora),
the DP pool, the optimiser groups and learning-rate schedule, checkpoint format (m34_ckpt), the stage / selection / validation / DEV tools.
What is new, per optimiser step s with lambda(s) = m40_mlm.lam_at(...) (0 for the first 10 % of the steps, linear ramp over the next 10 %, then lambda_max):
  * lambda(s) = 0: every draw goes through the exact M34 code path (`m34_backend.window_pots`, clean input); no masking.
  * lambda(s) > 0: a new-species draw is encoded from its label-blind masked token ids (m40_mlm.draw_masks; both the no-grad DP pass and the gradient pass
    see the same corrupted input, so the DP gradient is exact for the encoded input); the loss of the draw is
        CRF_NLL / n_known  +  lambda(s) * sum_p w_p CE_p / sum_p w_p           (p = selected, unlabeled, owned bases; w = 1 or 1e-4 in repeats)
    divided by the batch size as in M34.  Legacy-species draws are unchanged (fully labeled, clean input).
  * treatment X (arm M40_k0): new-species draws have no CRF target at all; their only loss is the MLM term (no DP, no no-grad pass).
The MLM head is PlantCAD2's own pre-trained, tied, frozen RCPS LM head (m40_mlm.lm_logits on the backbone's embedding), so the MLM gradient reaches only
the LoRA parameters (the C0 head is not on that path).
Usage: m40_train.py run --arm M40_k25U_MLM --out DIR [--steps 3072] [--stop-at N] [--resume] [--diag-grad]
"""
import argparse
import json
import os
import queue
import threading
import time

import numpy as np

import m40_common as W
import torch

import m40_mlm as MLM
import m40_repeats as RP
T = W.T
import m34_backend as BK          # noqa: E402
import m34_ckpt as CK             # noqa: E402
import a1_dp as DP                # noqa: E402
from src.m28.full_grammar_observed_partition import InfeasibleObservation   # noqa: E402
from src.m28.labels import IGNORE                                             # noqa: E402

M, G, A, MD = T.M, T.G, T.A, T.MD
CFG = dict(seed=0, sampling_seed=1, lr_head_peak=1e-4, lr_lora_peak=3e-4, lr_final=1e-5, warmup_frac=.05, weight_decay=.01, betas=[.9, .999],
           eps=1e-8, clip=1., epsilon=.1, max_exons=128, checkpoint_bp=256, chunk=2)          # = m34_train_lora.CFG
CFG.update(W.MLM_CFG)


# ---- draw -> item (CPU) -------------------------------------------------------------------------------------------------------------------------------
def build_item(r, labs, arm, rep, cfg):
    sp = r['species']
    lab = labs[sp]
    treat = W.species_treatment(arm, sp)
    chromosome, geometry, context, offset, targets, counts, win = T.example(lab, 'U' if treat == 'X' else treat, r)
    if treat == 'X':                                         # no CRF label on any base of the draw
        n = win['valid_bases']
        targets = dict(seg=np.full_like(targets['seg'], IGNORE), region=np.full_like(targets['region'], IGNORE))
        counts = (0, 0, n)
    a, b = geometry['genomic_context']
    geo = BK.geometry(a, b)
    it = dict(draw=r, w=win, strand=r['strand'], geo=geo, context=context, core=(offset, offset + win['valid_bases']), counts=counts,
              targets={k: np.asarray(targets[k]) for k in ('seg', 'region')}, mlm=None, **BK.item_windows(chromosome, geo, r['strand']))
    if cfg['mlm_lambda_max'] > 0 and W.mlm_on_species(sp, cfg):
        unknown = (it['targets']['seg'] < 0) & (it['targets']['region'] < 0)
        lo, hi = r['start'] - a, r['end'] - a
        if cfg['mlm_scope'] == 'unlabeled':
            scope = MLM.scope_genomic(b - a, lo, hi, unknown, r['strand'])
        else:
            scope = np.ones(b - a, dtype=bool)
        rng = np.random.default_rng([cfg['mlm_seed'], int(r['draw'])])
        it['mlm'] = MLM.draw_masks(MLM.ids_of(chromosome[a:b]), scope, rep.mask(sp, r['seqid'], a, b), geo['local_windows'], MLM.owned_ranges(geo), rng,
                                   cfg['mlm_rate'], cfg['mlm_mask_frac'], cfg['mlm_rand_frac'], cfg['mlm_repeat_weight'])
    return it


class DataPrefetch(threading.Thread):
    def __init__(self, rows, labs, arm, rep, cfg, depth=16):
        super().__init__(daemon=True)
        self.q = queue.Queue(depth)
        self.args = (rows, labs, arm, rep, cfg)
        self.err = None

    def run(self):
        rows, labs, arm, rep, cfg = self.args
        try:
            for r in rows:
                self.q.put(build_item(r, labs, arm, rep, cfg))
        except BaseException as e:
            self.err = e
            self.q.put(None)


# ---- encoder forward with masked ids -------------------------------------------------------------------------------------------------------------
def encode_ids(enc, ids_list, grad, device='cuda', batch=8):
    """As m34_backend._pc2_encode, but from (possibly corrupted) token-id arrays instead of strings.  -> list of [len, 1536] BF16 (autograd iff grad)."""
    out = [None] * len(ids_list)
    by_len = {}
    for i, x in enumerate(ids_list):
        if not 0 < len(x) <= A.PC_WIDTH:
            raise ValueError('Nonempty window of at most PC_WIDTH bp required')
        by_len.setdefault(len(x), []).append(i)
    cuda = str(device).startswith('cuda')
    with (torch.enable_grad() if grad else torch.no_grad()):
        for n, idx in by_len.items():
            for j in range(0, len(idx), batch):
                part = idx[j:j + batch]
                ids = torch.stack([torch.from_numpy(np.asarray(ids_list[i], dtype=np.int64)) for i in part]).to(device)
                with torch.autocast('cuda' if cuda else 'cpu', dtype=torch.bfloat16, enabled=cuda):
                    h = enc(ids).last_hidden_state
                if h.shape != (len(part), n, A.HIDDEN):
                    raise ValueError(f'unexpected feature shape {tuple(h.shape)}')
                if not bool(torch.isfinite(h).all()):
                    raise ValueError('Finite PlantCAD2 features required')
                for k, i in enumerate(part):
                    out[i] = h[k].to(torch.bfloat16)
    return out


def make_fwd(head, ctx, device='cuda'):
    E, comp = MLM.lm_params(ctx.enc)

    def fwd(it, idx, grad, use_mask, need_pots):
        """-> (potentials list or None, mlm stats dict or None)."""
        m = it['mlm'] if use_mask else None
        if m is None:
            return BK.window_pots(head, ctx, it, idx, grad, device=device), None         # exact M34 path
        h = encode_ids(ctx.enc, [m['ids'][i] for i in idx], grad, device)
        pots = MD.window_potentials(head, A.oriented_tokens(h, it['strand']), [it['seqs'][i] for i in idx], device=device) if need_pots else None
        st = dict(num=h[0].new_zeros((), dtype=torch.float32), wsum=0., correct=0, n=0)
        if grad or not need_pots:
            for k, i in enumerate(idx):
                pos = torch.from_numpy(m['pos'][i]).to(device)
                s, ws, c = MLM.window_mlm(h[k], pos, torch.from_numpy(m['tgt'][i]).to(device), torch.from_numpy(m['w'][i]).to(device), E, comp)
                st['num'] = st['num'] + s
                st['wsum'] += ws; st['correct'] += c; st['n'] += len(m['pos'][i])
        return pots, st
    return fwd


# ---- one optimiser step's gradient accumulation (device-agnostic: unit-tested on CPU with a fake encoder and a fake DP) -----------------------
def accumulate(items, fwd, submit, lam, batch, ch, epsilon, device='cuda', diag_params=None):
    """Forward/backward of one step over `items`; gradients are ACCUMULATED into .grad (caller zeroes / clips / steps).
    submit(context, logits, cs, ce, targets, eps) -> future with .get() -> dict(grad, loss, seconds) or raising InfeasibleObservation."""
    use_mask = lam > 0
    futs, pots = [], []
    for it in items:
        if it['counts'][0] + it['counts'][1] == 0:
            futs.append(None); pots.append(None); continue
        n = len(it['seqs'])
        with torch.no_grad():
            p = []
            for a in range(0, n, 4):
                p += [x.detach() for x in fwd(it, list(range(a, min(a + 4, n))), False, use_mask, True)[0]]
        ctx_logits = G.owned_context_logits(p, it['geo'], it['w']['strand']).cpu().numpy()
        cs, ce = it['core']
        futs.append(submit(it['context'], ctx_logits.astype(np.float32), cs, ce, it['targets'], epsilon))
        pots.append(p)
    st = dict(total_loss=0., dp_s=0., bwd_s=0., skipped=0, infeasible=0, by_sp={}, mlm_num=0., mlm_wsum=0., mlm_correct=0, mlm_n=0, mlm_draws=0,
              mlm_loss=0., crf_draws=0)
    g_crf = [torch.zeros_like(q) for q in diag_params] if diag_params else None
    g_mlm = [torch.zeros_like(q) for q in diag_params] if diag_params else None
    for it, fut, p in zip(items, futs, pots):
        sp = it['draw']['species']
        g = None
        if fut is None:
            st['skipped'] += 1; st['by_sp'].setdefault(sp, []).append(0.)
        else:
            try:
                r = fut.get()
            except InfeasibleObservation:
                st['infeasible'] += 1; r = None
            if r is not None:
                n_known = max(1, it['counts'][0] + it['counts'][1])
                g = torch.from_numpy(r['grad']).to(device) * (1. / n_known / batch)
                st['total_loss'] += r['loss'] / n_known / batch
                st['dp_s'] += r['seconds']
                st['by_sp'].setdefault(sp, []).append(r['loss'] / n_known)
                st['crf_draws'] += 1
        m = it['mlm'] if use_mask else None
        mlm_on = m is not None and m['wsum'] > 0
        if g is None and not mlm_on:
            continue
        st['mlm_draws'] += int(mlm_on)
        t0 = time.perf_counter()
        n = len(it['seqs'])
        for a in range(0, n, ch):
            idx = list(range(a, min(a + ch, n)))
            if g is None and not any(len(m['pos'][i]) for i in idx):
                continue
            new, ms = fwd(it, idx, True, use_mask, g is not None)
            crf_term = mlm_term = None
            if g is not None:
                pieces = list(p)
                for i, x in zip(idx, new):
                    pieces[i] = x
                crf_term = (G.owned_context_logits(pieces, it['geo'], it['w']['strand']) * g).sum()
            if mlm_on and ms is not None and ms['n']:
                mlm_term = ms['num'] * (lam / m['wsum'] / batch)
                st['mlm_num'] += float(ms['num'].detach()); st['mlm_wsum'] += ms['wsum']; st['mlm_correct'] += ms['correct']; st['mlm_n'] += ms['n']
                st['mlm_loss'] += float(mlm_term.detach())
            terms = [x for x in (crf_term, mlm_term) if x is not None]
            if not terms:
                continue
            if diag_params:
                for x, acc in ((crf_term, g_crf), (mlm_term, g_mlm)):
                    if x is not None:
                        gs = torch.autograd.grad(x, diag_params, retain_graph=True, allow_unused=True)
                        for q, gg in zip(acc, gs):
                            if gg is not None:
                                q += gg
            sum(terms).backward()
            del new, terms, crf_term, mlm_term
        if str(device).startswith('cuda'):
            torch.cuda.synchronize()
        st['bwd_s'] += time.perf_counter() - t0
    if diag_params:
        st['diag_grad_norm_crf'] = float(torch.sqrt(sum((q.float() ** 2).sum() for q in g_crf)))
        st['diag_grad_norm_mlm'] = float(torch.sqrt(sum((q.float() ** 2).sum() for q in g_mlm)))
    return st


def lr_pair(step, total):
    return (M.lr_at(step, total, CFG['lr_head_peak'], CFG['lr_final'], CFG['warmup_frac']),
            M.lr_at(step, total, CFG['lr_lora_peak'], CFG['lr_final'], CFG['warmup_frac']))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=('run',))
    ap.add_argument('--arm', required=True, choices=sorted(W.ARMS40))
    ap.add_argument('--out', required=True)
    ap.add_argument('--schedule', default=None, help='explicit draw schedule file (smoke); default: m40_common.schedule_path(arm)')
    ap.add_argument('--init', default=str(T.INIT))
    ap.add_argument('--steps', type=int, default=3072)
    ap.add_argument('--stop-at', type=int, default=None)
    ap.add_argument('--batch', type=int, default=8)
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--cpus', default='')
    ap.add_argument('--ckpt-every', type=int, default=64)
    ap.add_argument('--keep-every', type=int, default=768)
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--mlm-start-frac', type=float, default=None, help='smoke only (goes into the config)')
    ap.add_argument('--mlm-ramp-frac', type=float, default=None, help='smoke only (goes into the config)')
    ap.add_argument('--mlm-lambda-max', type=float, default=None, help='smoke only (goes into the config)')
    ap.add_argument('--diag-grad', action='store_true', help='smoke only: separate CRF / MLM LoRA gradient norms per step (2 extra autograd calls per chunk)')
    args = ap.parse_args()
    out = M.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cpus = M.parse_cpus(args.cpus)
    if cpus:
        os.sched_setaffinity(0, set(cpus))
    torch.set_num_threads(2)
    init = torch.load(args.init, map_location='cpu', weights_only=False)
    BK.check_init(init)
    if args.mlm_start_frac is not None:
        CFG['mlm_start_frac'] = args.mlm_start_frac
    if args.mlm_ramp_frac is not None:
        CFG['mlm_ramp_frac'] = args.mlm_ramp_frac
    if args.mlm_lambda_max is not None:
        CFG['mlm_lambda_max'] = args.mlm_lambda_max
    began = time.perf_counter()
    spec40 = W.ARMS40[args.arm]
    sched_path = M.Path(args.schedule) if args.schedule else W.schedule_path(args.arm)
    rows = T.read_schedule(sched_path)
    assert len(rows) >= args.steps * args.batch, (len(rows), args.steps, args.batch)
    rows = rows[:args.steps * args.batch]
    forbidden = T.forbidden_contigs()
    assert not any(r['seqid'] in forbidden for r in rows), 'schedule contains DEV / held-out / validation contigs'
    spec = T.arm_spec(args.arm)
    assert {r['species'] for r in rows} <= set(spec['species'])
    CFG.update(BK.encoder_cfg(init['config']), lora_mode='on', init_file=str(args.init), steps=args.steps, batch=args.batch, init_step=init['step'], arm=args.arm,
               species=list(spec['species']), new_species_label_fraction=spec['new_k'], treatment=spec['treatment'], schedule_md5=W.md5_file(sched_path),
               schedule=sched_path.name, heldout=list(T.VAL_REGIONS_LEGACY), m40_control_arm=spec40['control'], diag_grad=bool(args.diag_grad))
    labs, _ = T.arm_labelers(args.arm)
    assert all(r['seqid'] in labs[r['species']].lengths for r in rows)
    rep = RP.RepeatIndex([k for k in spec['species'] if W.mlm_on_species(k, CFG)])
    label_load_s = time.perf_counter() - began
    head = BK.new_head(CFG, CFG['seed'])
    head.load_state_dict(init['model'])
    head = head.to('cuda').train()
    ctx, lora_params = BK.load(CFG, lora=True)
    CFG['head_parameters'] = sum(p.numel() for p in head.parameters())
    CFG['lora_parameters'] = sum(p.numel() for p in lora_params)
    groups = [dict(params=list(head.parameters()), name='head'), dict(params=lora_params, name='lora')]
    opt = torch.optim.AdamW(groups, lr=CFG['lr_head_peak'], weight_decay=CFG['weight_decay'], betas=tuple(CFG['betas']), eps=CFG['eps'])
    start = 0
    ck = out / 'ckpt_latest.pt'
    if args.resume and ck.exists():
        start = CK.restore(ck, head, lambda st: BK.copy_lora(ctx, st), opt, CFG)
        print(json.dumps(dict(event='resumed', from_step=start, time=time.time())), flush=True)
    pool = DP.DPPool(args.workers, cpus)
    pf = DataPrefetch(rows[start * args.batch:], labs, args.arm, rep, CFG)
    pf.start()
    log = (out / 'training.jsonl').open('a', buffering=1)
    print(json.dumps(dict(event='start', arm=args.arm, start_step=start, steps=args.steps, head_parameters=CFG['head_parameters'],
                          lora_parameters=CFG['lora_parameters'], setup_s=time.perf_counter() - began, label_load_s=label_load_s,
                          eligible={k: v.n_eligible for k, v in labs.items()}, labeled={k: v.n_labeled for k, v in labs.items()},
                          mlm={k: CFG[k] for k in W.MLM_CFG})), flush=True)
    params = list(head.parameters()) + lora_params
    end = args.stop_at or args.steps
    fwd = make_fwd(head, ctx)
    for step in range(start, end):
        t_step = time.perf_counter()
        lr_h, lr_l = lr_pair(step, args.steps)
        lam = MLM.lam_at(step, args.steps, CFG['mlm_lambda_max'], CFG['mlm_start_frac'], CFG['mlm_ramp_frac'])
        opt.param_groups[0]['lr'] = lr_h
        opt.param_groups[1]['lr'] = lr_l
        opt.zero_grad(set_to_none=True)
        items, t_wait = [], 0.
        for _ in range(args.batch):
            t0 = time.perf_counter()
            it = pf.q.get()
            t_wait += time.perf_counter() - t0
            if it is None:
                raise RuntimeError('prefetch failed') from pf.err
            items.append(it)
        st = accumulate(items, fwd, pool.submit, lam, args.batch, CFG['chunk'], CFG['epsilon'], 'cuda', lora_params if args.diag_grad else None)
        norm_g = torch.nn.utils.clip_grad_norm_(params, CFG['clip'], error_if_nonfinite=True)
        opt.step()
        if any(not bool(torch.isfinite(q).all()) for q in params):
            raise RuntimeError('nonfinite parameter')
        row = dict(step=step + 1, lr_head=lr_h, lr_lora=lr_l, mlm_lambda=lam, batch_loss=st['total_loss'], mlm_loss=st['mlm_loss'],
                   mlm_ce=st['mlm_num'] / st['mlm_wsum'] if st['mlm_wsum'] else None, mlm_acc=st['mlm_correct'] / st['mlm_n'] if st['mlm_n'] else None,
                   mlm_positions=st['mlm_n'], mlm_draws=st['mlm_draws'], crf_draws=st['crf_draws'], grad_norm=float(norm_g), step_s=time.perf_counter() - t_step,
                   wait_prefetch_s=t_wait, bwd_s=st['bwd_s'], dp_sum_s=st['dp_s'], skipped_all_unknown=st['skipped'], infeasible=st['infeasible'],
                   nll_by_species={k: float(np.mean(v)) for k, v in st['by_sp'].items()}, time=time.time(),
                   **({k: st[k] for k in ('diag_grad_norm_crf', 'diag_grad_norm_mlm')} if args.diag_grad else {}))
        log.write(json.dumps(row, allow_nan=False) + '\n')
        if (step + 1) % 16 == 0 or step < 3:
            print(json.dumps({k: row[k] for k in ('step', 'mlm_lambda', 'batch_loss', 'mlm_ce', 'mlm_acc', 'grad_norm', 'step_s', 'bwd_s', 'infeasible')}), flush=True)
        done = step + 1
        if done % args.ckpt_every == 0 or done == end:
            stt = CK.save_resumable(ck, head, BK.lora_state(ctx), opt, done, CFG, done * args.batch)
            if done % args.keep_every == 0 or done == end:
                CK.save_model(out / f'model_{done:06d}.pt', stt)
    json.dump(dict(completed=True, arm=args.arm, steps=args.steps, reached=end, seconds=time.perf_counter() - began,
                   peak_gpu_alloc_GB=torch.cuda.max_memory_allocated() / 2**30, peak_gpu_reserved_GB=torch.cuda.max_memory_reserved() / 2**30),
              open(out / f'DONE_{end}.json', 'w'))
    pool.close()
    os._exit(0)


if __name__ == '__main__':
    main()
