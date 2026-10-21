"""Smoke draw schedule (CPU only): 20 steps x 8 draws that hit the stress cases named in protocol 19 §3.13 -- the genes with the most exons and the genes
with the longest CDS span (> one 24,576-bp window) of every species of the arm -- followed by the arm's real schedule (so the file still has 3,072 x 8
rows and the trainer's length assertion holds).  Usage: m34_smoke_schedule.py <arm> <out.jsonl>"""
import json
import sys

import m34_train_common as T
import m34_common as C

arm, outp = sys.argv[1], sys.argv[2]
spec = T.arm_spec(arm)
real = T.read_schedule(T.SCHEDULES / f'schedule_{arm}.jsonl')
tab = C.species_table()
special = []
for key in spec['species']:
    tr = set(tab[key]['train'])
    elig = [t for t in map(json.loads, open(C.CATALOG / key / 'transcripts.jsonl')) if t['eligible'] and t['seqid'] in tr]
    by_exons = sorted(elig, key=lambda t: (-t['n_exons'], t['id']))[:1]
    by_span = sorted(elig, key=lambda t: (-(t['span'][1] - t['span'][0]), t['id']))[:1]
    for t in by_exons + by_span:
        cells = C.gene_cells(t['span'], tab[key]['lengths'][t['seqid']])
        a, b = cells[0]
        special.append(dict(species=key, seqid=t['seqid'], strand=t['strand'], start=a, end=b, mode='G', gene=t['id'],
                            id=f"{key}|{t['seqid']}|{a}|{t['strand']}", note=f"n_exons={t['n_exons']} span={t['span'][1] - t['span'][0]}"))
rows = special[:160]
rows += real[:160 - len(rows)]
rows = [dict(r, draw=i, step=i // 8 + 1) for i, r in enumerate(rows)] + real[160:]
assert len(rows) == len(real)
with open(outp, 'w') as f:
    for r in rows:
        f.write(json.dumps(r) + '\n')
print(json.dumps(dict(arm=arm, special=len(special), total=len(rows))))
