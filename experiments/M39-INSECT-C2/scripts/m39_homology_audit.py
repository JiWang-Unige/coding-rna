"""Protocol 29 §5.3 V9 (descriptive only, filters nothing): how close is each NEW species' proteome to the two fully labeled OLD species?  Adapted from
`m34_homology_audit.py`: queries = a fixed random sample (random.Random(0)) of 1,500 longest-per-gene proteins of every new species (data/m39_insect/busco/<key>.faa,
written by `m39_busco.py faa`); targets = all longest-per-gene proteins of Drosophila + Tribolium; MMseqs2 easy-search -s 4.0 -e 1e-5; best hit with query
coverage >= 0.7 -> fractions with best-hit identity >= 50 / 70 / 90 %, against both old species together and against each one.
usage: m39_homology_audit.py [threads]
"""
import json
import random
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import m39_common as C

MMSEQS = str(Path.home() / 'envs/mmseqs/bin/mmseqs')
FAA = C.DATA / 'busco'
OUTD = C.OUT / 'homology'


def read_faa(path):
    seqs, name = {}, None
    for line in open(path):
        if line[0] == '>':
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
    return {k: ''.join(v) for k, v in seqs.items()}


def main(threads=4):
    OUTD.mkdir(parents=True, exist_ok=True)
    old, new = C.old_keys(), C.new_keys()
    target, query = OUTD / 'old_species.faa', OUTD / 'new_species_sample.faa'
    with open(target, 'w') as g:
        for key in old:
            for pid, s in read_faa(FAA / f'{key}.faa').items():
                g.write(f'>{key}|{pid}\n{s}\n')
    nq = {}
    with open(query, 'w') as g:
        for key in new:
            seqs = read_faa(FAA / f'{key}.faa')
            ids = sorted(seqs)
            pick = sorted(random.Random(0).sample(ids, min(1500, len(ids))))
            nq[key] = len(pick)
            for pid in pick:
                g.write(f'>{key}|{pid}\n{seqs[pid]}\n')
    res = OUTD / 'best_hits.tsv'
    if not res.exists():
        subprocess.run([MMSEQS, 'easy-search', str(query), str(target), str(res), str(OUTD / 'tmp'), '--threads', str(threads), '-s', '4.0', '-e', '1e-5',
                        '--max-seqs', '300', '--format-output', 'query,target,pident,qcov,tcov,evalue,bits'], check=True)      # pident in percent
    groups = {'both_old': set(old), **{k: {k} for k in old}}
    best = defaultdict(lambda: defaultdict(float))
    for line in open(res):
        q, t, pid, qc, tc, ev, bits = line.rstrip('\n').split('\t')
        if float(qc) < 0.7:
            continue
        for g, sp in groups.items():
            if t.split('|')[0] in sp and float(pid) / 100 > best[q][g]:
                best[q][g] = float(pid) / 100
    out = {}
    for key in new:
        out[key] = {}
        for g in groups:
            vals = [best[q][g] for q in best if q.split('|')[0] == key]
            out[key][g] = dict(n_queries=nq[key], **{f'frac_best_identity_ge_{x}': sum(v >= x / 100 for v in vals) / nq[key] for x in (50, 70, 90)})
    json.dump(out, open(OUTD / 'homology_audit.json', 'w'), indent=1)
    for key in out:
        print(key, {g: {x: round(out[key][g][f'frac_best_identity_ge_{x}'], 3) for x in (50, 70, 90)} for g in ('both_old',)})


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 4)
