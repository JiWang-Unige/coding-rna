"""Protein-homology audit: how close is each PANEL proteome to the training proteomes of every nested arm (T2/T4/T6/T10)?  CPU only (MMseqs2).
Per panel protein: best hit among the arm's training proteins -> fractions with best-hit identity >= 50/70/90% at query coverage >= 70%.  Queries = a fixed random sample (seed 0) of 1,500 proteins per panel species (isoforms included; ratios are reported); targets = all training proteins.  Purely descriptive: it does not filter anything.
usage: m34_homology_audit.py [threads]
"""
import gzip
import json
import random
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import m34_common as C
from m34_schedule import LEGACY, NEW

MMSEQS = str(Path.home() / 'envs/mmseqs/bin/mmseqs')
OUTD = C.DATA / 'homology'
RAW = C.RAW


def faa(tag, path, out, sample=None):
    """Append a proteome with '<tag>|' prefixed ids; sample = n -> a fixed random subset (random.Random(0)) of n records (queries only)."""
    keep = None
    if sample:
        with gzip.open(path, 'rt') as f:
            n = sum(1 for line in f if line[0] == '>')
        keep = set(random.Random(0).sample(range(n), min(sample, n)))
    i = -1
    on = True
    with gzip.open(path, 'rt') as f, open(out, 'a') as g:
        for line in f:
            if line[0] == '>':
                i += 1
                on = keep is None or i in keep
                if on:
                    g.write(f'>{tag}|{line[1:].split()[0]}\n')
            elif on:
                g.write(line)


def main(threads=4):
    OUTD.mkdir(parents=True, exist_ok=True)
    rows = C.new_species()
    train_faa = OUTD / 'train_T10.faa'
    train_faa.unlink(missing_ok=True)
    for key in NEW:
        r = rows[key]
        faa(key, C.raw_dir(r) / f"{r['accession']}_{r['assembly_name']}_protein.faa.gz", train_faa)
    faa('arabidopsis_thaliana', RAW / 'legacy_proteins/GCF_000001735.4_TAIR10.1_protein.faa.gz', train_faa)
    faa('oryza_sativa', RAW / 'legacy_proteins/GCF_034140825.1_ASM3414082v1_protein.faa.gz', train_faa)
    panel = {}
    for f in sorted((RAW / 'panel_proteins').glob('*_protein.faa.gz')):
        panel[f.name.split('_protein')[0]] = f
    m30 = Path.home() / 'coding-rna/data/m30/raw'
    for f in m30.glob('*/*_protein.faa.gz'):
        panel[f.name.split('_protein')[0]] = f
    query = OUTD / 'panel_all.faa'
    query.unlink(missing_ok=True)
    for tag, f in panel.items():
        faa(tag, f, query, sample=1500)
    res = OUTD / 'best_hits.tsv'
    if not res.exists():
      subprocess.run([MMSEQS, 'easy-search', str(query), str(train_faa), str(res), str(OUTD / 'tmp'), '--threads', str(threads), '-s', '4.0', '-e', '1e-5',
                    '--max-seqs', '300', '--format-output', 'query,target,pident,qcov,tcov,evalue,bits'], check=True)            # mmseqs pident is in PERCENT (0-100)
    arms = {'T2': set(LEGACY), 'T4': set(LEGACY + NEW[:2]), 'T6': set(LEGACY + NEW[:4]), 'T10': set(LEGACY + NEW)}
    nq = defaultdict(int)
    with open(query) as f:
        for line in f:
            if line[0] == '>':
                nq[line[1:].split('|')[0]] += 1
    best = defaultdict(lambda: defaultdict(float))        # (panel, query) -> arm -> best pident among hits with qcov >= .7
    for line in open(res):
        q, t, pid, qc, tc, ev, bits = line.rstrip('\n').split('\t')
        if float(qc) < 0.7:
            continue
        ptag = q.split('|')[0]
        for arm, sp in arms.items():
            if t.split('|')[0] in sp:
                d = best[ptag, q]
                if float(pid) / 100 > d[arm]:
                    d[arm] = float(pid) / 100
    out = {}
    for tag in panel:
        out[tag] = {}
        for arm in arms:
            vals = [best[k][arm] for k in best if k[0] == tag]
            n = nq[tag]
            out[tag][arm] = dict(n_proteins=n, **{f'frac_best_identity_ge_{t}': sum(v >= t / 100 for v in vals) / n for t in (50, 70, 90)} if n else {})
    json.dump(out, open(OUTD / 'homology_audit.json', 'w'), indent=1)
    for tag in out:
        print(tag, {a: round(out[tag][a]['frac_best_identity_ge_70'], 3) for a in arms})


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 4)
