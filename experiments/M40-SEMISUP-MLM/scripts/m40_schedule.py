"""M40 draw schedules (CPU only).  M40_k25U_MLM / M40_k50U_MLM reuse the M34 files schedule_T10_k25U.jsonl / schedule_T10_k50U.jsonl unchanged (md5 recorded);
M40_k0 gets its own file from the UNCHANGED M34 sampler law (m34_schedule.draw_schedule, sampling seed 1) at new-species label fraction 0: with no labeled
gene the gene-first half of the sampler falls back to uniform windows for the new species, so sampling uses no annotation of a new species.  The two legacy
species are sampled as in every M34 arm (k = 1).
Output: ~/coding-rna/outputs/M40-SEMISUP-MLM/schedules/{schedule_M40_k0.jsonl, schedule_summary.json}
usage: m40_schedule.py
"""
import json
from collections import Counter

import m40_common as W

S, C = W.S, W.C


def main():
    W.SCHEDULES40.mkdir(parents=True, exist_ok=True)
    summ = {}
    for arm, spec in W.ARMS40.items():
        if spec['schedule'] != 'M40':
            p = W.schedule_path(arm)
            summ[arm] = dict(arm=arm, schedule=str(p), md5=W.md5_file(p), reused_from=spec['schedule'])
            continue
        sspec = S.ARMS[arm]
        pools = {k: S.load_species(k, 1.0 if k in S.LEGACY else sspec['new_k']) for k in sspec['species']}
        assert all(not pools[k]['genes'] for k in W.NEW), 'k0 must have no labeled new-species gene'
        rows = S.draw_schedule(arm, sspec, pools)
        forb = W.T.forbidden_contigs()
        assert not any(r['seqid'] in forb for r in rows)
        p = W.schedule_path(arm)
        with open(p, 'w') as f:
            for r in rows:
                f.write(json.dumps(r) + '\n')
        summ[arm] = dict(arm=arm, schedule=str(p), md5=W.md5_file(p), draws=len(rows), unique_windows=len({r['id'] for r in rows}),
                         draws_per_species=dict(Counter(r['species'] for r in rows)),
                         gene_first_per_species=dict(Counter(r['species'] for r in rows if r['mode'] == 'G')), sampling_seed=S.SAMPLING_SEED)
        print(json.dumps({k: summ[arm][k] for k in ('arm', 'md5', 'unique_windows')}), flush=True)
    json.dump(summ, open(W.SCHEDULES40 / 'schedule_summary.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
