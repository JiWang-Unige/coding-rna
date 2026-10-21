"""Training-step DP backward with / without the reusable alpha arena (CPU, 1 thread): same real TRAIN windows as m29_g0b_dp_profile."""
import json, os, sys, time
import numpy as np, torch
import m29_g0b_dp_profile as P      # re-points hard-coded Baobab paths
import m29_common as M
import src.m28.full_grammar_boundary_compiled as BC
from m29_reverse_arena import _reverse_messages_arena
C, PP, coarse_channel_nll = P.C, P.PP, P.coarse_channel_nll
torch.set_num_threads(1)
plan, picks = P.load_plan_rows(int(sys.argv[1]))
windows, strata, catalogs, original = PP.load_population()
sequences, _ = C.P.read_train_DNA()
orig_rm = BC._reverse_messages
rows = []
def run(ordinal, fn):
    BC._reverse_messages = fn
    expected = plan[ordinal]; w = windows[expected['window_id']]
    chromosome, geometry, context, offset, targets = C.example(w, expected, sequences, catalogs)
    L, n = len(context), w['valid_bases']
    cuts = tuple(x for x in (offset, offset + n) if 0 < x < L)
    g = torch.Generator().manual_seed(ordinal)
    logits = (torch.randn(L, 15, generator=g) * 2.0).requires_grad_(True)
    t = time.perf_counter()
    res = coarse_channel_nll(context, logits, 128, offset, offset + n, targets, .1, cuts=cuts, checkpoint_bp=256)
    tf = time.perf_counter() - t; t = time.perf_counter()
    res['loss'].backward()
    return tf, time.perf_counter() - t, logits.grad.clone(), float(res['loss'].detach())
run(picks[0], orig_rm); run(picks[0], _reverse_messages_arena)     # JIT warm-up for both
for ordinal in picks:
    a = run(ordinal, orig_rm); b = run(ordinal, _reverse_messages_arena)
    rows.append(dict(draw=ordinal, fwd_s=a[0], bwd_original_s=a[1], bwd_arena_s=b[1], total_original_s=a[0] + a[1], total_arena_s=b[0] + b[1],
                     max_abs_grad_diff=float((a[2] - b[2]).abs().max()), loss_equal=a[3] == b[3], loadavg=os.getloadavg()[0]))
    print(json.dumps(rows[-1]), flush=True)
(M.M29 / 'G0b' / f'dp_arena_bench_{os.getpid()}.json').write_text(json.dumps(rows, indent=1))
