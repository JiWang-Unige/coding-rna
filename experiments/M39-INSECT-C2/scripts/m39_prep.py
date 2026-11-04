"""M39 data preparation (no model output): raw NCBI files -> data/m39_insect/prepared/<key>/ genome.fa, reference.gff3, lengths.json, splits.json
(same layout as data/m1_screen/<species>/ and the M34 prepared species).  Copy of `m34_prep.py` with the protocol 29 §3.2 rules:
  * main nuclear chromosomes only: assembled-molecule + Chromosome / Linkage Group + Primary Assembly + length >= 5 Mb;
  * NEW species: validation chromosome = the lower-median-length chromosome; only its first L bp are kept, as pseudo-contig '<chrom>__val<L/1e6>M'
    ('<chrom>__valAll' when L is the whole chromosome), genes wholly inside, coordinates unchanged; L = the smallest whole-Mb length >= 8 Mb whose crop holds
    >= 500 CDS-assessable primary references of the frozen ruler (counted on the whole chromosome, a reference counts when its gene ends inside the crop),
    capped at the whole chromosome; the remainder of that chromosome is discarded;
  * OLD species: selection chromosome S = the lower-median-length chromosome, descriptive held-out D = the next shorter one, both whole;
  * all other main chromosomes are TRAINING chromosomes;
  * trans-spliced transcripts (GFF attribute exception=trans-splicing; only Drosophila mod(mdg4), 23 of its 31 mRNAs) and their child features are dropped: the frozen
    parser rejects a transcript with exons on both strands, and the M28 grammar cannot represent one (protocol 29 §4.12).
usage: m39_prep.py [species_key ...]
"""
import gzip
import json
import shutil
import sys
from datetime import datetime

import m39_common as C
import common_ruler as CR                     # noqa: E402  (frozen; path set by a1_common)


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


def drop_trans_spliced(lines):
    """Remove trans-spliced transcripts (non-gene features carrying exception=trans-splicing) and every feature whose parent is one of them; the gene line
    stays (Drosophila mod(mdg4) keeps its 8 cis-spliced isoforms).  -> (kept lines, removed transcript ids)."""
    def trans(r):
        return len(r) >= 9 and r[2] != 'gene' and 'exception=trans-splicing' in r[8]
    bad = {attrs(r[8]).get('ID') for r in lines if trans(r) and r[2] not in ('exon', 'CDS')}
    bad.discard(None)
    keep = [r for r in lines if not (trans(r) or (len(r) >= 9 and attrs(r[8]).get('Parent', '').split(',')[0] in bad))]
    return keep, sorted(bad)


def reference_gene_ends(out, hname, seq, lines):
    """Frozen-ruler CDS-assessable primary references of the whole held-out chromosome -> sorted list of their genes' GFF end coordinates."""
    tmp = out / '_heldout_full'
    tmp.mkdir(parents=True, exist_ok=True)
    gff = tmp / 'reference.gff3'
    with open(gff, 'w') as f:
        f.write('##gff-version 3\n')
        for r in lines:
            if len(r) >= 9 and r[0] == hname and r[2] != 'region':
                f.write('\t'.join(r) + '\n')
    length = {hname: len(seq)}
    coding = CR.S.E.parse_annotation(gff, length, protein_coding_only=True)
    with open(gff) as handle:
        _, _, raw_cds, _ = CR.P.raw_features(handle, hname)
    refs, _ = CR.select_reference(coding, raw_cds, {hname: seq})
    ends = sorted(coding['genes'][t['gene_id']]['end'] for t in refs)
    shutil.rmtree(tmp)
    return ends


def crop_length(ends, chrom_len):
    L = C.VAL_MIN_BP
    while L < chrom_len and sum(e <= L for e in ends) < C.VAL_TARGET_REFS:
        L += C.VAL_STEP_BP
    return min(L, chrom_len)


def prepare(key, row):
    raw = C.raw_dir(row)
    stem = row['ftp_stem']
    rep = C.parse_assembly_report(raw / (stem + '_assembly_report.txt'))
    chroms = C.main_chromosomes(rep)
    held = C.choose_heldout(chroms)
    hname = held['refseq']
    old = row['role'] == 'old'
    dname = C.choose_d(chroms)['refseq'] if old else None
    train = sorted(r['refseq'] for r in chroms if r['refseq'] not in (hname, dname))
    wanted = {r['refseq'] for r in chroms}
    seqs = read_fasta(raw / (stem + '_genomic.fna.gz'), wanted)
    assert set(seqs) == wanted, wanted - set(seqs)
    for r in chroms:
        assert len(seqs[r['refseq']]) == r['length'], (r['refseq'], len(seqs[r['refseq']]), r['length'])
    lines = [l.rstrip('\n').split('\t') for l in gzip.open(raw / (stem + '_genomic.gff.gz'), 'rt') if not l.startswith('#')]
    lines, trans_spliced = drop_trans_spliced(lines)
    out = C.PREP / key
    out.mkdir(parents=True, exist_ok=True)
    crop_info, val_refs_est = None, None
    if old:
        contigs = {c: seqs[c] for c in train + [hname, dname]}
        vcontig, kept = hname, set(contigs)
    else:
        ends = reference_gene_ends(out, hname, seqs[hname], lines)
        crop = crop_length(ends, held['length'])
        val_refs_est = sum(e <= crop for e in ends)
        vcontig = f'{hname}__val{crop // 1_000_000}M' if crop < held['length'] else f'{hname}__valAll'
        contigs = {c: seqs[c] for c in train}
        contigs[vcontig] = seqs[hname][:crop]
        kept = set(train) | {hname}
        crop_info = dict(chromosome=hname, start0=0, end=crop, whole_chromosome=crop == held['length'],
                         references_whole_chromosome=len(ends), references_in_crop_estimate=val_refs_est)
    write_fasta(out / 'genome.fa', contigs)
    lens = {k: len(s) for k, s in contigs.items()}
    json.dump(lens, open(out / 'lengths.json', 'w'), indent=1)
    # ---- GFF (single pass; the M34 / M30 parent-closure crop logic for the new-species validation chromosome) ----
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
    crop_end = crop_info['end'] if crop_info else None
    gene_ok = {}
    if not old:
        for r in lines:
            if len(r) >= 9 and r[0] == hname and r[2] in ('gene', 'pseudogene'):
                gene_ok[attrs(r[8])['ID']] = int(r[4]) <= crop_end
    n_feat, n_gene = 0, 0
    with open(out / 'reference.gff3', 'w') as f:
        f.write('##gff-version 3\n')
        for k, n in lens.items():
            f.write(f'##sequence-region {k} 1 {n}\n')
        for r in lines:
            if len(r) < 9 or r[0] not in kept or r[2] == 'region':
                continue
            if not old and r[0] == hname:
                a = attrs(r[8])
                rid = root(a['ID']) if 'Parent' not in a and 'ID' in a else root(a['Parent'].split(',')[0]) if 'Parent' in a else ''
                if 'Parent' not in a and a.get('ID') in gene_ok:
                    ok = gene_ok[a['ID']]
                elif rid in gene_ok:
                    ok = gene_ok[rid]
                else:
                    ok = int(r[4]) <= crop_end
                if not ok:
                    continue
                r = [vcontig] + r[1:]
            n_gene += r[2] == 'gene'
            f.write('\t'.join(r) + '\n')
            n_feat += 1
    rule = ('main nuclear chromosomes (assembled-molecule, Chromosome or Linkage Group, Primary Assembly, >=5 Mb); '
            + ('OLD species: selection chromosome S = lower-median length, descriptive held-out D = next shorter, both whole' if old else
               'NEW species: validation chromosome = lower-median length; crop = smallest whole-Mb length >= 8 Mb with >= 500 CDS-assessable primary references '
               '(genes wholly inside), capped at the whole chromosome; remainder discarded'))
    order = C.length_order(chroms)
    splits = dict(species=key, role=row['role'], organism=row['organism'], accession=row['accession'], assembly=stem, rule=rule,
                  main_chromosomes={r['refseq']: r['length'] for r in order}, main_chromosome_names={r['refseq']: r['name'] for r in order},
                  heldout_chromosome=hname, heldout_length=held['length'], validation_contig=vcontig, validation_crop=crop_info,
                  heldout_d=[dname] if old else [], discarded_bp=0 if old else held['length'] - crop_info['end'], train=train, lengths=lens,
                  train_bp=sum(lens[c] for c in train), main_bp=sum(r['length'] for r in chroms), assembly_bp=sum(r['length'] for r in rep if r['refseq'] != 'na'),   # RefSeq rows only (Apis cerana lists GenBank copies too)
                 
                  dropped_sequences=len(rep) - len(chroms), dropped_trans_spliced_transcripts=trans_spliced, gff_features=n_feat, gff_genes=n_gene,
                  prepared=datetime.now().astimezone().isoformat(timespec='seconds'))
    json.dump(splits, open(out / 'splits.json', 'w'), indent=1)
    print(json.dumps(dict(species=key, role=row['role'], train=len(train), train_bp=splits['train_bp'], val=vcontig, heldout=hname, d=dname,
                          crop=crop_info, features=n_feat, genes=n_gene)), flush=True)
    return splits


def main():
    rows = C.roster()
    keys = sys.argv[1:] or list(rows)
    C.PREP.mkdir(parents=True, exist_ok=True)
    allrows = {r['key']: r for r in C.species_rows()}
    for k in keys:
        prepare(k, allrows[k])


if __name__ == '__main__':
    main()
