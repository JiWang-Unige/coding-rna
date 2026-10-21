"""Compact training curve: per 256-step block mean normalized NLL (by species), step time, grad norm, nonzero-loss share."""
import json, sys, numpy as np
rows = {}
for l in open(sys.argv[1]):
    r = json.loads(l); rows[r['step']] = r          # last record per step (resume duplicates)
steps = sorted(rows)
out = []
for a in range(0, len(steps), 256):
    b = [rows[s] for s in steps[a:a + 256]]
    sp = {}
    for r in b:
        for k, v in r['nll_by_species'].items():
            sp.setdefault(k, []).append(v)
    out.append(dict(steps=f'{b[0]["step"]}-{b[-1]["step"]}', mean_batch_loss=round(float(np.mean([r['batch_loss'] for r in b])), 5),
                    nll_by_species={k: round(float(np.mean(v)), 5) for k, v in sp.items()}, mean_step_s=round(float(np.mean([r['step_s'] for r in b])), 2),
                    median_grad_norm=round(float(np.median([r['grad_norm'] for r in b])), 3), lr_end=round(b[-1].get('lr', b[-1].get('lr_head')), 7)))
tot = sum(r['step_s'] for r in rows.values())
print(json.dumps(dict(blocks=out, steps=len(rows), sum_step_seconds=round(tot)), indent=1))
