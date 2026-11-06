"""PANEL-P6: whole-nuclear-genome FASTA inputs of the baselines (CPU only, stdlib only, streaming).  Same rule as M38 `vip_make_genomes.py`.

For every core-10 species, from the frozen NCBI download of PANEL-PREP (data/panel/raw/<accession>_<name>/*_genomic.fna.gz, md5-checked in report 25 §2):
  <out>/<key>.nuclear.fa         every sequence except the assembly report's 'non-nuclear' unit (organelles); sequences, NCBI soft-masking and line width unchanged
                                 -> Vipsania (its masked-LM loss down-weights soft-masked bases, as in M38)
  <out>/<key>.nuclear.upper.fa   the same sequences upper-cased -> Helixer / Tiberius / ANNEVO / Carbon-A (the M28 DEV native runs and the M29/M30 Carbon-A runs were given
                                 upper-case DNA; Tiberius additionally runs with --no_softmasking)
The scored (main) sequences of the panel (prepared/<key>/lengths.json, 119 for the core 10) are checked to be identical to the ones our T2 run read (md5 of the upper-cased
sequence text vs prepared/<key>/genome.fa).  Output manifest <out>/genomes_manifest.json.
  split KEY DEST   one upper-case FASTA per main sequence (<contig>.fasta) + unplaced.fasta (all other nuclear sequences) for the per-file Carbon-A inference (M30 layout)
Usage (ont):  taskset -c 25-31 nice -n 10 python3.12 -I p6_genomes.py build [--only KEY]   |   ... p6_genomes.py split KEY DEST
"""
import argparse
import gzip
import hashlib
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_common as P          # noqa: E402

M38_MANIFEST = P.HOME / 'coding-rna' / 'outputs' / 'M38-VIPSANIA' / 'genomes' / 'genomes_manifest.json'


def opener(p):
    return gzip.open(p, 'rt') if str(p).endswith('.gz') else open(p)


def nonnuclear(report):
    """sequence names (RefSeq accession for GCF, GenBank accession for GCA assemblies, as in the FASTA headers) of the 'non-nuclear' assembly unit"""
    out = set()
    for line in open(report):
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip('\n').split('\t')
        if f[7] == 'non-nuclear':
            out.update(x for x in (f[4], f[6]) if x != 'na')
    return out


def seq_md5(path, wanted):
    """md5 of the concatenated upper-cased sequence text of the wanted sequence names"""
    res, cur, h = {}, None, None
    with opener(path) as f:
        for line in f:
            if line[0] == '>':
                if cur in wanted:
                    res[cur] = h.hexdigest()
                cur = line[1:].split()[0]
                h = hashlib.md5() if cur in wanted else None
            elif h is not None:
                h.update(line.strip().upper().encode())
    if cur in wanted:
        res[cur] = h.hexdigest()
    return res


def build(key):
    d = P.raw_dir(key)
    src = next(d.glob('*_genomic.fna.gz'))
    report = next(d.glob('*_assembly_report.txt'))
    drop = nonnuclear(report)
    main = P.lengths(key)
    out, up = P.GENOMES / f'{key}.nuclear.fa', P.GENOMES / f'{key}.nuclear.upper.fa'
    kept, dropped, lower = {}, {}, 0
    cur, keep, n = None, False, 0
    sha, shu = hashlib.sha256(), hashlib.sha256()
    with opener(src) as f, open(str(out) + '.tmp', 'w') as w, open(str(up) + '.tmp', 'w') as wu:
        for line in f:
            if line[0] == '>':
                if cur is not None:
                    (kept if keep else dropped)[cur] = n
                cur = line[1:].split()[0]
                keep, n = cur not in drop, 0
                if keep:
                    w.write(line); sha.update(line.encode())
                    wu.write(line); shu.update(line.encode())
            else:
                n += len(line.strip())
                if keep:
                    w.write(line); sha.update(line.encode())
                    lu = line.upper()
                    wu.write(lu); shu.update(lu.encode())
                    lower += sum(map(line.count, 'acgtn'))
    (kept if keep else dropped)[cur] = n
    os.replace(str(out) + '.tmp', out)
    os.replace(str(up) + '.tmp', up)
    missing = [c for c in main if c not in kept]
    assert not missing, (key, 'scored sequences missing from the nuclear genome', missing)
    assert all(kept[c] == main[c] for c in main), (key, 'length differs from prepared lengths.json')
    a = seq_md5(up, set(main))
    b = seq_md5(P.PREP / key / 'genome.fa', set(main))
    assert a == b and set(a) == set(main), (key, 'scored sequence differs from prepared/genome.fa')
    tot = sum(kept.values())
    rec = dict(key=key, source=str(src), report=str(report), fasta=str(out), fasta_upper=str(up), bytes=out.stat().st_size, bytes_upper=up.stat().st_size,
               sha256=sha.hexdigest(), sha256_upper=shu.hexdigest(), n_sequences=len(kept), total_bp=tot, main_sequences=len(main), main_bp=sum(main.values()),
               other_nuclear_sequences=len(kept) - len(main), other_nuclear_bp=tot - sum(main.values()), soft_masked_fraction=round(lower / max(tot, 1), 4),
               dropped_non_nuclear=dropped, scored_sequence_md5_equal_to_prepared=True)
    if M38_MANIFEST.exists():                                   # Dryas / Eschscholzia: the M38 whole genomes (Vipsania pre on set B) came from the same NCBI files
        m38 = json.load(open(M38_MANIFEST)).get(key)
        if m38:
            rec['same_as_m38_nuclear_fa'] = (m38['sha256'] == rec['sha256'])
    return rec


def split(key, dest):
    """per-file Carbon-A inputs (upper case): <dest>/<contig>.fasta for every main sequence, <dest>/unplaced.fasta for the rest (if any)"""
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    main = set(P.lengths(key))
    h, unpl = None, None
    with open(P.GENOMES / f'{key}.nuclear.upper.fa') as f:
        for line in f:
            if line[0] == '>':
                name = line[1:].split()[0]
                if h is not None and h is not unpl:
                    h.close()
                if name in main:
                    h = open(dest / f'{name}.fasta.tmp', 'w')
                else:
                    if unpl is None:
                        unpl = open(dest / 'unplaced.fasta.tmp', 'w')
                    h = unpl
            h.write(line)
    if h is not None and h is not unpl:
        h.close()
    if unpl is not None:
        unpl.close()
    for p in dest.glob('*.tmp'):
        os.replace(p, str(p)[:-4])
    got = sorted(p.name for p in dest.glob('*.fasta'))
    assert len([g for g in got if g != 'unplaced.fasta']) == len(main), (key, got)
    print(json.dumps(dict(key=key, dest=str(dest), files=len(got), unplaced=(dest / 'unplaced.fasta').exists())))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    b = sub.add_parser('build'); b.add_argument('--only', default=None)
    s = sub.add_parser('split'); s.add_argument('key'); s.add_argument('dest')
    a = ap.parse_args()
    if a.cmd == 'split':
        return split(a.key, a.dest)
    P.GENOMES.mkdir(parents=True, exist_ok=True)
    mp = P.GENOMES / 'genomes_manifest.json'
    man = json.load(open(mp)) if mp.exists() else {}
    for key in P.by_size(P.CORE10):
        if a.only and key != a.only:
            continue
        t0 = time.time()
        man[key] = build(key)
        man[key]['build_seconds'] = round(time.time() - t0, 1)
        json.dump(man, open(mp, 'w'), indent=1)
        r = man[key]
        print(f"{key}: {r['n_sequences']} sequences {r['total_bp']:,} bp (main {r['main_sequences']} / {r['main_bp']:,}; other {r['other_nuclear_sequences']} / "
              f"{r['other_nuclear_bp']:,}; dropped non-nuclear {len(r['dropped_non_nuclear'])}) soft-masked {r['soft_masked_fraction']} "
              f"same_as_m38={r.get('same_as_m38_nuclear_fa')} {r['build_seconds']} s", flush=True)


if __name__ == '__main__':
    main()
