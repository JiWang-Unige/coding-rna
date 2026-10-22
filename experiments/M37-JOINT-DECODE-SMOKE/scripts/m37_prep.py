"""M37 data preparation (no model output): 9311 chromosome 10 (CM000135.2) from the NCBI raw download ->
data/m37/oryza_sativa_9311/{genome.fa, reference.gff3, lengths.json, inputs/CM000135.2.fasta, source.json};
and the assembly-A chromosome (NC_089044.1, from the frozen data/m1_screen/oryza_sativa/genome.fa) as data/m37/align/NC_089044.1.fa for minimap2."""
import gzip, hashlib, json
from pathlib import Path
import m37_common as C
import manifest as MF

RAW = C.DATA / 'raw' / C.B_ASM
OUTD = C.DATA / C.B_SPECIES


def read_fasta_gz(path, wanted):
    seqs, name, buf = {}, None, []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line[0] == '>':
                if name == wanted:
                    seqs[name] = ''.join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif name == wanted:
                buf.append(line.strip())
    if name == wanted:
        seqs[name] = ''.join(buf).upper()
    return seqs[wanted]


def write_fa(path, name, s):
    with open(path, 'w') as f:
        f.write('>' + name + '\n')
        for i in range(0, len(s), 80):
            f.write(s[i:i + 80] + '\n')


def main():
    (OUTD / 'inputs').mkdir(parents=True, exist_ok=True)
    (C.DATA / 'align').mkdir(parents=True, exist_ok=True)
    s = read_fasta_gz(RAW / 'chr10.fna.gz', C.B_CONTIG)
    write_fa(OUTD / 'genome.fa', C.B_CONTIG, s)
    write_fa(OUTD / 'inputs' / f'{C.B_CONTIG}.fasta', C.B_CONTIG, s)
    json.dump({C.B_CONTIG: len(s)}, open(OUTD / 'lengths.json', 'w'), indent=1)
    n_feat = 0
    with open(OUTD / 'reference.gff3', 'w') as f:
        f.write('##gff-version 3\n')
        f.write(f'##sequence-region {C.B_CONTIG} 1 {len(s)}\n')
        for line in gzip.open(RAW / f'{C.B_ASM}_genomic.gff.gz', 'rt'):
            if line.startswith('#'):
                continue
            r = line.rstrip('\n').split('\t')
            if len(r) == 9 and r[0] == C.B_CONTIG:
                f.write(line); n_feat += 1
    # assembly A chromosome for the alignment (frozen file read only)
    M = C.register()
    a = MF.read_allowed(C.DATA.parent / 'm1_screen' / 'oryza_sativa' / 'genome.fa', {C.A_CONTIG: M.LENGTH[C.A_CONTIG]})[C.A_CONTIG].upper()
    write_fa(C.DATA / 'align' / f'{C.A_CONTIG}.fa', C.A_CONTIG, a)
    md5 = {p.name: hashlib.md5(p.read_bytes()).hexdigest() for p in (RAW / 'chr10.fna.gz', RAW / f'{C.B_ASM}_genomic.gff.gz')}
    ncbi = dict(l.split()[::-1] for l in open(RAW / 'md5checksums.txt') if l.strip())
    ncbi = {Path(k).name: v for k, v in ncbi.items()}
    src = dict(assembly='GCA_000004655.2', assembly_name='ASM465v1', cultivar='93-11 (Oryza sativa Indica Group)', submitter='Beijing Institute of Genomics, CAS',
               assembly_release='2011-02-25', annotation='submitted by BIG CAS, release 2015-03-23, 37,358 protein-coding genes (NCBI feature_count)',
               contig=C.B_CONTIG, chromosome=10, length=len(s), n_bases=s.count('N'), gff_features=n_feat,
               download_time=open(RAW / 'download_date.txt').read().strip(), local_md5=md5, ncbi_md5={k: ncbi.get(k) for k in md5},
               md5_match={k: ncbi.get(k) == v for k, v in md5.items()},
               assembly_A=dict(contig=C.A_CONTIG, length=len(a), n_bases=a.count('N'), source='data/m1_screen/oryza_sativa/genome.fa (GCF_034140825.1 AGIS-1.0)'))
    json.dump(src, open(OUTD / 'source.json', 'w'), indent=1)
    print(json.dumps(src, indent=1))


if __name__ == '__main__':
    main()
