"""M36 deterministic draw schedule of one (arm, species) adaptation run (protocol 22 §3: 'seed-first + uniform' half/half on the label-set A training chromosomes,
isomorphic to the M34 gene-prior sampler `m34_schedule.draw_schedule`).  CPU only, no model.
Law (one law for every arm, `default_rng(1)`): per draw, with probability 0.5 a window uniform over the species' (A training contig, strand, 24,576-bp grid window) pool,
else 'seed-first': a uniform seed locus, then a uniform grid window among the (<= 5) windows covering its CDS span (`m34_common.gene_cells`), on the seed's strand.
S1 and F1 use the same seed set and therefore the SAME draw sequence (as M34's U / F arms share one); S3 draws from its own seed set.  The validation pseudo-contig
and every evaluation-set B contig are not in the pool.
usage: m36_schedule.py make <arm> <species> [out.jsonl]      (reads the arm's seeds as m36_train does; default OUT/schedules/schedule_<arm>_<species>.jsonl)
"""
import json
import sys
from pathlib import Path

import numpy as np

import m36_common as W
import m34_common as C


def draw_schedule(species, seeds, lengths, contigs=None, steps=W.TRAIN['steps'], batch=W.TRAIN['batch'], seed=W.TRAIN['sampling_seed'],
                  uniform_fraction=W.TRAIN['uniform_fraction']):
    """seeds: [seed rows] (any order; sorted here); lengths: {contig: length} of the sampled contigs (training contigs only).  -> [draw rows]"""
    contigs = list(contigs if contigs is not None else lengths)
    if species in W.TARGETS:
        W.assert_not_B(species, contigs)
    seeds = sorted((r for r in seeds if r['seqid'] in contigs), key=lambda r: (r['seqid'], r['span'][0], r['span'][1], r['strand'], r['id']))
    nwin = [len(C.grid(lengths[c])) for c in contigs]
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(steps * batch):
        mode, gid = 'U', None
        if rng.random() < uniform_fraction or not seeds:
            tot = 2 * sum(nwin)
            u = int(rng.integers(tot))
            strand = '+' if u < tot // 2 else '-'
            u %= tot // 2
            for c, n in zip(contigs, nwin):
                if u < n:
                    a, b = C.grid(lengths[c])[u]
                    break
                u -= n
        else:
            mode = 'G'
            r = seeds[int(rng.integers(len(seeds)))]
            c, strand, gid = r['seqid'], r['strand'], r['id']
            cells = C.gene_cells(r['span'], lengths[c])
            a, b = cells[int(rng.integers(len(cells)))]
        rows.append(dict(draw=i, step=i // batch + 1, species=species, seqid=c, strand=strand, start=a, end=b, mode=mode,
                         id=f'{species}|{c}|{a}|{strand}', **({'gene': gid} if mode == 'G' else {})))
    return rows


def write_schedule(rows, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = str(path) + '.tmp'
    with open(tmp, 'w') as f:
        for r in rows:
            f.write(json.dumps(r) + '\n')
    Path(tmp).replace(path)


def schedule_path(arm, species):
    return W.OUT / 'schedules' / f'schedule_{W.arm_key(arm, species)}.jsonl'


def make(arm, species, out=None):
    """build (or reuse, if byte-identical) the schedule of an arm from its seeds"""
    import m36_train_inputs as TI
    inputs = TI.load_arm_inputs(arm, species)
    lens = {c: W.a_lengths(species)[c] for c in W.TARGETS[species]['train']}
    rows = draw_schedule(species, inputs['rows'], lens)
    path = Path(out) if out else schedule_path(arm, species)
    if path.exists():
        old = [json.loads(x) for x in open(path)]
        assert old == rows, f'{path} exists and differs from the schedule rebuilt from the current seeds'
    else:
        write_schedule(rows, path)
    return path, rows


if __name__ == '__main__':
    if sys.argv[1] == 'make':
        p, rows = make(sys.argv[2], sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else None)
        print(json.dumps(dict(schedule=str(p), draws=len(rows), seed_first=sum(r['mode'] == 'G' for r in rows) / len(rows),
                              unique_windows=len({r['id'] for r in rows}))))
