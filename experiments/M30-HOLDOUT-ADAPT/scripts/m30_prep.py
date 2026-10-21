"""M30 data preparation (no model output): selected chromosomes -> data/m1_screen/<species>/{genome.fa, reference.gff3, lengths.json, splits.json}
and inputs/<contig>.fasta for Carbon-A.  The validation chromosome is replaced by its first 20 Mb as pseudo-contig '<seqid>__val20M';
GFF features of that contig keep only genes that lie fully inside the crop (coordinates unchanged: crop starts at 1)."""
import gzip, json, re, sys
from pathlib import Path
import m30_common as M


def read_fasta(path, wanted):
    seqs, name, buf = {}, None, []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line[0] == '>':
                if name in wanted:
                    seqs[name] = ''.join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif name in wanted:
                buf.append(line.strip())
    if name in wanted:
        seqs[name] = ''.join(buf).upper()
    return seqs


def attrs(s):
    return dict(x.split('=', 1) for x in s.split(';') if '=' in x)


def main():
    for species, sp in M.SPECIES.items():
        raw = M.RAW / sp['asm']
        out = M.DATA / species
        (out / 'inputs').mkdir(parents=True, exist_ok=True)
        seqs = read_fasta(raw / (sp['asm'] + '_genomic.fna.gz'), set(sp['chroms']))
        assert set(seqs) == set(sp['chroms']), set(sp['chroms']) - set(seqs)
        contigs = {}
        for c in sp['chroms']:
            if c == sp['val_chrom']:
                contigs[sp['val_contig']] = seqs[c][:M.VAL_BP]
            elif c in sp['A'] or c in sp['B']:
                contigs[c] = seqs[c]
        with open(out / 'genome.fa', 'w') as f:
            for k, s in contigs.items():
                f.write('>' + k + '\n')
                for i in range(0, len(s), 80):
                    f.write(s[i:i + 80] + '\n')
                with open(out / 'inputs' / (k + '.fasta'), 'w') as g:
                    g.write('>' + k + '\n')
                    for i in range(0, len(s), 80):
                        g.write(s[i:i + 80] + '\n')
        lens = {k: len(s) for k, s in contigs.items()}
        json.dump(lens, open(out / 'lengths.json', 'w'), indent=1)
        # GFF
        lines = [l.rstrip('\n').split('\t') for l in gzip.open(raw / (sp['asm'] + '_genomic.gff.gz'), 'rt') if not l.startswith('#')]
        valc, crop = sp['val_chrom'], M.VAL_BP
        keep_ids_outside = set()
        inside_ids = set()
        parent = {}
        for r in lines:
            if len(r) < 9:
                continue
            a = attrs(r[8])
            if 'ID' in a and 'Parent' in a:
                parent[a['ID']] = a['Parent'].split(',')[0]
        def root(i):
            seen = 0
            while i in parent and seen < 10:
                i = parent[i]; seen += 1
            return i
        # genes of the validation chromosome fully inside the crop
        gene_ok = {}
        for r in lines:
            if len(r) >= 9 and r[0] == valc and r[2] in ('gene', 'pseudogene'):
                gene_ok[attrs(r[8])['ID']] = int(r[4]) <= crop
        with open(out / 'reference.gff3', 'w') as f:
            f.write('##gff-version 3\n')
            for k, n in lens.items():
                f.write(f'##sequence-region {k} 1 {n}\n')
            n_feat = 0
            for r in lines:
                if len(r) < 9 or r[0] not in sp['chroms']:
                    continue
                if r[0] == valc:
                    if r[2] == 'region':
                        continue
                    a = attrs(r[8])
                    rid = root(a.get('ID', '')) if 'Parent' not in a else root(a['Parent'].split(',')[0])
                    if 'Parent' not in a and a.get('ID') in gene_ok:
                        ok = gene_ok[a['ID']]
                    elif rid in gene_ok:
                        ok = gene_ok[rid]
                    else:
                        ok = int(r[4]) <= crop
                    if not ok:
                        continue
                    r = [sp['val_contig']] + r[1:]
                elif r[0] not in lens:
                    continue
                f.write('\t'.join(r) + '\n'); n_feat += 1
        json.dump(dict(species=species, assembly=sp['asm'], A=sp['A'], B=sp['B'], validation_chromosome=valc, validation_contig=sp['val_contig'],
                       train=sp['train'], lengths=lens), open(out / 'splits.json', 'w'), indent=1)
        print(species, {k: v for k, v in lens.items()}, 'gff features', n_feat, flush=True)


if __name__ == '__main__':
    main()
