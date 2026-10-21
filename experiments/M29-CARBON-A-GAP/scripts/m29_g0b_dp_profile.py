"""G0b: where does a training step spend its time?  (a) CPU-only: FP64 full-grammar partition function (logZ_all, logZ_obs, backward)
on REAL TRAIN windows (real DNA context + real observation targets, random logits - the DP cost depends on DNA reachability and
allowed-role masks, not on logit values); (b) optional GPU part: encoder forward (joint24576 BF16) + head, then the same DP with
the real logits, optimizer step; GPU peak memory.   Uses fresh-draw plan rows 0..K-1 (first draws of both species), no gating.
Usage: m29_g0b_dp_profile.py cpu K   |  m29_g0b_dp_profile.py gpu K
"""
import json
import os
import sys
import time
from pathlib import Path
import torch
import m29_common as M   # puts pipeline/m28 on sys.path
sys.path.insert(0, str(M.ROOT / 'pipeline/m28/scripts'))
import fit_full_grammar_objective_contrast as C   # noqa: E402
import prepare_full_grammar_fresh_draws as PP     # noqa: E402
from src.m28.full_grammar_coarse_channel import coarse_channel_nll   # noqa: E402
from src.m28.full_grammar_boundary_compiled import compiled_boundary_full_grammar_log_partition   # noqa: E402

# Some frozen modules hard-code the Baobab ROOT; re-point every such module-level Path at this checkout (no frozen file is edited).
BAOBAB = '/home/users/j/jwang/coding-rna'
for _m in list(sys.modules.values()):
    for _k, _v in list(getattr(_m, '__dict__', {}).items()):
        if isinstance(_v, Path) and str(_v).startswith(BAOBAB):
            setattr(_m, _k, M.ROOT / Path(str(_v)).relative_to(BAOBAB))
OUT = M.M29 / 'G0b'
OUT.mkdir(parents=True, exist_ok=True)


def load_plan_rows(K):
    plan = [json.loads(x) for x in (M.OUTS / 'M28-FULL-GRAMMAR-FRESH-DRAWS-PREPARATION/training_plan.jsonl').read_text().splitlines()]
    # first K//2 of each species (fixed, before looking at any timing)
    sel = {}
    for i, r in enumerate(plan):
        sel.setdefault(r['species'], []).append(i)
    species = sorted(sel)
    picks = []
    for j in range(K):
        picks.append(sel[species[j % len(species)]][j // len(species)])
    return plan, picks


def main():
    mode, K = sys.argv[1], int(sys.argv[2])
    torch.set_num_threads(1)
    plan, picks = load_plan_rows(K)
    windows, strata, catalogs, original = PP.load_population()
    sequences, _ = C.P.read_train_DNA()
    rows = []
    cpu_only = mode == 'cpu'
    if not cpu_only:
        model = None
    for ordinal in picks:
        expected = plan[ordinal]
        w = windows[expected['window_id']]
        chromosome, geometry, context, offset, targets = C.example(w, expected, sequences, catalogs)
        L, n = len(context), w['valid_bases']
        cuts = tuple(x for x in (offset, offset + n) if 0 < x < L)
        g = torch.Generator().manual_seed(ordinal)
        logits = (torch.randn(L, 15, generator=g) * 2.0).requires_grad_(True)
        row = dict(draw=ordinal, species=w['species'], context_bp=L, core_bp=n, windows_in_context=len(geometry['genomic_windows']),
                   fine=expected['fine_positions'], coarse=expected['coarse_positions'], unknown=expected['fully_unknown_positions'],
                   loadavg=os.getloadavg()[0])
        c0, t0 = time.process_time(), time.perf_counter()
        with torch.no_grad():
            allm = compiled_boundary_full_grammar_log_partition(context, logits.detach(), 128, cuts=cuts, checkpoint_bp=256)
        row['logZ_all_forward_s'] = time.perf_counter() - t0
        t1 = time.perf_counter()
        res = coarse_channel_nll(context, logits, 128, offset, offset + n, targets, .1, cuts=cuts, checkpoint_bp=256)
        row['loss_forward_s(all+obs)'] = time.perf_counter() - t1
        t2 = time.perf_counter()
        res['loss'].backward()
        row['backward_s'] = time.perf_counter() - t2
        row['cpu_process_s_total'] = time.process_time() - c0
        row['wall_s_total'] = time.perf_counter() - t0
        st = res['statistics_all']; row['state_count'] = (st[0] if isinstance(st, list) else st)['state_count']
        row['peak_rss_gb'] = int(open('/proc/self/status').read().split('VmHWM:')[1].split()[0]) / 1e6
        rows.append(row)
        print(json.dumps(row), flush=True)
    (OUT / f'dp_profile_{mode}_{os.getpid()}.json').write_text(json.dumps(rows, indent=1))


if __name__ == '__main__':
    main()
