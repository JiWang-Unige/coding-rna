"""CPU dry check of the trainer data path without the encoder: the smoke stress draws (most-exons / longest-span genes of all 10 species) plus real T10 draws go through the
trainer's example builder and the parallel window DP with random 15-role potentials.  Counts InfeasibleObservation (independent of the potentials) and reports
loss / DP time per species.  usage: m34_cpu_dp_check.py <smoke_schedule.jsonl> [n_real=24] [workers=4] [cpus=25-28]"""
import json, sys, time
import numpy as np
import m34_train_common as T
import a1_dp as DP
from src.m28.full_grammar_observed_partition import InfeasibleObservation

def main():
    sched_path = sys.argv[1]
    n_real = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    cpus = T.M.parse_cpus(sys.argv[4] if len(sys.argv) > 4 else '25-28')
    arm = sys.argv[5] if len(sys.argv) > 5 else 'T10'
    rows = T.read_schedule(sched_path)
    special = [r for r in rows[:160] if 'note' in r]
    real = rows[160:160 + n_real]
    labs, treat = T.arm_labelers(arm)
    rng = np.random.default_rng(0)
    pool = DP.DPPool(workers, cpus)
    jobs = []
    t0 = time.time()
    for r in special + real:
        chromosome, geometry, context, offset, targets, counts, win = T.example(labs[r['species']], treat, r)
        lg = (rng.standard_normal((len(context), 15)) * 1.0).astype(np.float32)
        n = win['valid_bases']
        fut = pool.submit(context, lg, offset, offset + n, {k: np.asarray(targets[k]) for k in ('seg', 'region')}) if counts[0] + counts[1] > 0 else None
        jobs.append((r, counts, fut, len(context)))
    out, bad = [], 0
    for r, counts, fut, L in jobs:
        if fut is None:
            out.append(dict(species=r['species'], note=r.get('note'), skipped_all_unknown=True)); continue
        try:
            res = fut.get()
            out.append(dict(species=r['species'], note=r.get('note'), gene=r.get('gene'), ctx_bp=L, fine=counts[0], coarse=counts[1], nll_per_known=res['loss'] / max(1, counts[0] + counts[1]), dp_s=round(res['seconds'], 1)))
        except InfeasibleObservation as e:
            bad += 1
            out.append(dict(species=r['species'], note=r.get('note'), gene=r.get('gene'), infeasible=str(e)[:100]))
    pool.close()
    print(json.dumps(dict(arm=arm, windows=len(jobs), special=len(special), infeasible=bad, seconds=round(time.time() - t0), details=out)))



if __name__ == '__main__':
    main()
