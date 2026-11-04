"""M38-VIPSANIA: build the four whole-nuclear-genome FASTA files that Vipsania annotates (CPU only, stdlib only, streaming).

The M30/M34 working FASTAs of Dryas / Eschscholzia (`data/m1_screen/<sp>/genome.fa`) and of rice (`genome.fa`) are NOT whole genomes (they hold the chromosomes of the
experiment only, with a 20-Mb crop of the validation chromosome), so the whole assemblies are rebuilt from the original NCBI files:
  arabidopsis_thaliana      data/m1_screen/arabidopsis_thaliana/genome.fa  (GCF_000001735.4; drop the two organelles)
  oryza_sativa              data/m1_screen/oryza_sativa/genome.full.fa     (GCF_034140825.1; drop mitochondrion, plastid and the mitochondrial plasmid B1)
  dryas_octopetala          data/m30/raw/GCF_963921425.1_drDryOcto1.1/*_genomic.fna.gz   (assembly report: drop 'non-nuclear' unit)
  eschscholzia_californica  data/m30/raw/GCF_041429995.1_dmEscCali1.0.p/*_genomic.fna.gz (assembly report: all sequences are Primary Assembly scaffolds)
Sequences, soft-masking (lowercase, as distributed by NCBI) and line width are copied unchanged.  The sequences of the scored contigs (DEV / set B) are checked to be identical
(md5 of the UPPER-CASED sequence text: the M30 working FASTAs of Dryas / Eschscholzia are upper-case, NCBI's are soft-masked) to those the scoring reference
(data/m1_screen/<sp>/genome.fa) uses.  Output: <out>/<key>.nuclear.fa and <out>/genomes_manifest.json.
Usage (ont):  nice -n 19 python3.12 vip_make_genomes.py [--out DIR]
"""
import argparse
import gzip
import hashlib
import json
import os
import time
from pathlib import Path

DATA = Path(os.path.expanduser('~/coding-rna-reorg/data'))
GENOMES = {
    'arabidopsis_thaliana': dict(
        assembly='GCF_000001735.4', name='TAIR10.1', src=DATA / 'm1_screen/arabidopsis_thaliana/genome.fa', report=None,
        exclude=['NC_037304.1', 'NC_000932.1'], scored=['NC_003074.8'], role='DEV',
        score_src=DATA / 'm1_screen/arabidopsis_thaliana/genome.fa'),
    'oryza_sativa': dict(
        assembly='GCF_034140825.1', name='ASM3414082v1', src=DATA / 'm1_screen/oryza_sativa/genome.full.fa', report=None,
        exclude=['NC_011033.1', 'NC_001320.1', 'NC_001751.1'], scored=['NC_089041.1'], role='DEV',
        score_src=DATA / 'm1_screen/oryza_sativa/genome.fa'),
    'dryas_octopetala': dict(
        assembly='GCF_963921425.1', name='drDryOcto1.1', src=DATA / 'm30/raw/GCF_963921425.1_drDryOcto1.1/GCF_963921425.1_drDryOcto1.1_genomic.fna.gz',
        report=DATA / 'm30/raw/GCF_963921425.1_drDryOcto1.1/GCF_963921425.1_drDryOcto1.1_assembly_report.txt', exclude=[],
        scored=['NC_141367.1', 'NC_141369.1', 'NC_141371.1', 'NC_141373.1'], role='B',
        score_src=DATA / 'm1_screen/dryas_octopetala/genome.fa'),
    'eschscholzia_californica': dict(
        assembly='GCF_041429995.1', name='dmEscCali1.0.p', src=DATA / 'm30/raw/GCF_041429995.1_dmEscCali1.0.p/GCF_041429995.1_dmEscCali1.0.p_genomic.fna.gz',
        report=DATA / 'm30/raw/GCF_041429995.1_dmEscCali1.0.p/GCF_041429995.1_dmEscCali1.0.p_assembly_report.txt', exclude=[],
        scored=['NW_028238309.1', 'NW_028238461.1', 'NW_028238483.1'], role='B',
        score_src=DATA / 'm1_screen/eschscholzia_californica/genome.fa'),
}


def opener(p):
    return gzip.open(p, 'rt') if str(p).endswith('.gz') else open(p)


def report_nonnuclear(path):
    """RefSeq accessions of the assembly-report rows that belong to the 'non-nuclear' assembly unit"""
    out = set()
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip('\n').split('\t')
        if f[7] == 'non-nuclear':
            out.add(f[6])
    return out


def seq_md5(path, wanted):
    """md5 of the concatenated UPPER-CASED sequence text (newlines removed) of the wanted sequence names"""
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


def build(key, g, outdir):
    drop = set(g['exclude'])
    if g['report']:
        drop |= report_nonnuclear(g['report'])
    out = outdir / f'{key}.nuclear.fa'
    tmp = out.with_suffix('.fa.tmp')
    kept, dropped = {}, {}
    cur, keep, n, lower = None, False, 0, 0
    sha = hashlib.sha256()
    with opener(g['src']) as f, open(tmp, 'w') as w:
        for line in f:
            if line[0] == '>':
                if cur is not None:
                    (kept if keep else dropped)[cur] = n
                cur = line[1:].split()[0]
                keep, n = cur not in drop, 0
                if keep:
                    w.write(line); sha.update(line.encode())
            else:
                n += len(line.strip())
                if keep:
                    w.write(line); sha.update(line.encode())
                    lower += sum(map(line.count, 'acgtn'))
    (kept if keep else dropped)[cur] = n
    os.replace(tmp, out)
    for c in g['scored']:
        assert c in kept, (key, c, 'scored contig missing from the nuclear genome')
    a, b = seq_md5(out, set(g['scored'])), seq_md5(g['score_src'], set(g['scored']))
    assert a == b and set(a) == set(g['scored']), (key, a, b)
    return dict(key=key, assembly=g['assembly'], assembly_name=g['name'], source=str(g['src']), fasta=str(out), bytes=out.stat().st_size, sha256=sha.hexdigest(),
                n_sequences=len(kept), total_bp=sum(kept.values()), soft_masked_fraction=round(lower / max(sum(kept.values()), 1), 4), dropped=dropped, scored_role=g['role'], scored_contigs={c: kept[c] for c in g['scored']},
                scored_sequence_md5_equal_to_scoring_reference=True, scored_sequence_md5=a,
                largest=sorted(kept.items(), key=lambda kv: -kv[1])[:14])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=os.path.expanduser('~/coding-rna/outputs/M38-VIPSANIA/genomes'))
    ap.add_argument('--only', default=None)
    args = ap.parse_args()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    man = {}
    mp = outdir / 'genomes_manifest.json'
    if mp.exists():
        man = json.load(open(mp))
    for key, g in GENOMES.items():
        if args.only and key != args.only:
            continue
        t0 = time.time()
        man[key] = build(key, g, outdir)
        man[key]['build_seconds'] = round(time.time() - t0, 1)
        print(key, man[key]['n_sequences'], 'sequences', man[key]['total_bp'], 'bp', man[key]['fasta'], flush=True)
    json.dump(man, open(mp, 'w'), indent=1)


if __name__ == '__main__':
    main()
