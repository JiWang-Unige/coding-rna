"""M34 data preparation for the NEW training species (no model output): raw NCBI files -> data/m34_train/prepared/<key>/
   genome.fa, reference.gff3, lengths.json, splits.json   (same layout as data/m1_screen/<species>/ used by M28/M30/M31).
Rules (fixed before any model output):
  * main nuclear chromosomes only: assembled-molecule + Chromosome + Primary Assembly + length >= 5 Mb (organelles/unplaced/unlocalized dropped);
  * internal-validation chromosome = the lower-median-length chromosome; only its first 8 Mb is kept, as pseudo-contig '<chrom>__val8M'
    (GFF features of genes fully inside the crop, coordinates unchanged); the remainder of that chromosome is discarded;
  * all other main chromosomes are TRAINING chromosomes.
Also writes the panel-exclusion audit to prepared/panel_exclusion_audit.json.
usage: m34_prep.py [species_key ...]
"""
import gzip
import json
import sys
from datetime import datetime

import m34_common as C
import panel_exclusion as PX


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


def write_fasta(path, contigs):
    with open(path, 'w') as f:
        for k, s in contigs.items():
            f.write('>' + k + '\n')
            for i in range(0, len(s), 80):
                f.write(s[i:i + 80] + '\n')


def prepare(key, row):
    raw = C.raw_dir(row)
    stem = f"{row['accession']}_{row['assembly_name']}"
    rep = C.parse_assembly_report(raw / (stem + '_assembly_report.txt'))
    chroms = C.main_chromosomes(rep)
    held = C.choose_heldout(chroms)
    hname = held['refseq']
    crop = min(C.VAL_BP, held['length'])
    vcontig = f'{hname}__val{crop // 1_000_000}M'
    train = sorted(r['refseq'] for r in chroms if r['refseq'] != hname)
    wanted = set(train) | {hname}
    seqs = read_fasta(raw / (stem + '_genomic.fna.gz'), wanted)
    assert set(seqs) == wanted, wanted - set(seqs)
    for r in chroms:
        assert len(seqs[r['refseq']]) == r['length'], (r['refseq'], len(seqs[r['refseq']]), r['length'])
    contigs = {c: seqs[c] for c in train}
    contigs[vcontig] = seqs[hname][:crop]
    out = C.PREP / key
    out.mkdir(parents=True, exist_ok=True)
    write_fasta(out / 'genome.fa', contigs)
    lens = {k: len(s) for k, s in contigs.items()}
    json.dump(lens, open(out / 'lengths.json', 'w'), indent=1)
    # ---- GFF (single pass; the same parent-closure logic as M30 prep) ----
    lines = [l.rstrip('\n').split('\t') for l in gzip.open(raw / (stem + '_genomic.gff.gz'), 'rt') if not l.startswith('#')]
    parent = {}
    for r in lines:
        if len(r) >= 9:
            a = attrs(r[8])
            if 'ID' in a and 'Parent' in a:
                parent[a['ID']] = a['Parent'].split(',')[0]

    def root(i):
        n = 0
        while i in parent and n < 10:
            i = parent[i]
            n += 1
        return i
    gene_ok = {r2: None for r2 in ()}
    for r in lines:
        if len(r) >= 9 and r[0] == hname and r[2] in ('gene', 'pseudogene'):
            gene_ok[attrs(r[8])['ID']] = int(r[4]) <= crop
    n_feat, n_gene = 0, 0
    with open(out / 'reference.gff3', 'w') as f:
        f.write('##gff-version 3\n')
        for k, n in lens.items():
            f.write(f'##sequence-region {k} 1 {n}\n')
        for r in lines:
            if len(r) < 9 or r[0] not in wanted:
                continue
            if r[0] == hname:
                if r[2] == 'region':
                    continue
                a = attrs(r[8])
                rid = root(a['ID']) if 'Parent' not in a and 'ID' in a else root(a['Parent'].split(',')[0]) if 'Parent' in a else ''
                if 'Parent' not in a and a.get('ID') in gene_ok:
                    ok = gene_ok[a['ID']]
                elif rid in gene_ok:
                    ok = gene_ok[rid]
                else:
                    ok = int(r[4]) <= crop
                if not ok:
                    continue
                r = [vcontig] + r[1:]
            elif r[0] not in lens:
                continue
            n_gene += r[2] == 'gene'
            f.write('\t'.join(r) + '\n')
            n_feat += 1
    splits = dict(species=key, organism=row['organism'], accession=row['accession'], assembly=stem,
                  rule='main nuclear chromosomes (assembled-molecule, Chromosome, Primary Assembly, >=5 Mb); validation chromosome = lower-median length; '
                       'validation crop = first 8 Mb (genes fully inside); remainder of that chromosome discarded',
                  main_chromosomes={r['refseq']: r['length'] for r in chroms}, heldout_chromosome=hname, heldout_length=held['length'],
                  validation_crop=dict(chromosome=hname, start0=0, end=crop), validation_contig=vcontig,
                  discarded_bp=held['length'] - crop, train=train, lengths=lens,
                  train_bp=sum(lens[c] for c in train), dropped_sequences=len(rep) - len(chroms),
                  gff_features=n_feat, gff_genes=n_gene, prepared=datetime.now().astimezone().isoformat(timespec='seconds'))
    json.dump(splits, open(out / 'splits.json', 'w'), indent=1)
    print(key, 'train', len(train), 'contigs', splits['train_bp'], 'bp; val', vcontig, 'held-out chrom', hname, held['length'], 'features', n_feat, 'genes', n_gene, flush=True)
    return splits


def main():
    rows = C.new_species()
    keys = sys.argv[1:] or list(rows)
    C.PREP.mkdir(parents=True, exist_ok=True)
    for k in keys:
        prepare(k, rows[k])
    allrows = C.new_species()
    training = [(r['organism'].split()[0], r['accession'], C.FAMILY[k]) for k, r in allrows.items()]
    training += [(v['organism'].split()[0], v['accession'], v['family']) for v in C.LEGACY.values()]
    audit = PX.audit(training)
    audit['training_species'] = [dict(genus=g, accession=a, family=f) for g, a, f in training]
    json.dump(audit, open(C.PREP / 'panel_exclusion_audit.json', 'w'), indent=1)
    assert not audit['violations'], audit['violations']
    print('panel audit: no genus/accession violation; family overlaps:', sorted({(x['training_genus'], x['panel_genus'], x['panel_tier']) for x in audit['family_overlaps']}))


if __name__ == '__main__':
    main()
