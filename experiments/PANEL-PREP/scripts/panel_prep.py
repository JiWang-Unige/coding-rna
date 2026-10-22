"""Independent-panel preparation (no model output): raw NCBI files -> data/panel/prepared/<key>/{genome.fa, reference.gff3, lengths.json, splits.json}
(the m1_screen / M34 `prepared` layout).  The whole prepared genome is EVALUATION data: there is no train/val split and no crop.
Sequence rule (fixed before any model output, identical to M34/M30; see panel_common.main_sequences):
  assembled-molecule + Chromosome + Primary Assembly + >= 5 Mb; if the assembly has no assembled-molecule, Primary Assembly sequences >= 10 Mb.
Organelles, unplaced/unlocalized scaffolds and short sequences are dropped (counted, with their annotated genes, in splits.json).
`stats`: reference statistics with the frozen ruler through the unchanged M34 Labeler (primary transcripts, CDS-assessable references, exon strata).
usage: panel_prep.py prep  [key ...]      (CPU, a few minutes per species; Dipteris/Chenopodium need ~6 GB RAM)
       panel_prep.py stats [key ...]      (CPU; writes prepared/<key>/reference_stats.json and results summary)
"""
import gzip
import json
import sys
from collections import Counter
from datetime import datetime

import panel_common as P
import m34_common as C
from m34_prep import read_fasta, write_fasta, attrs


def ncbi_report(raw):
    p = raw / 'ncbi_dataset_report.json'
    if not p.exists():
        return {}
    rep = json.load(open(p)).get('reports', [{}])[0]
    ai, an = rep.get('assembly_info', {}), rep.get('annotation_info', {})
    busco = an.get('busco', {})
    return dict(organism=rep.get('organism', {}).get('organism_name'), source_database=rep.get('source_database'),
                assembly_level=ai.get('assembly_level'), assembly_release_date=ai.get('release_date'), assembly_status=ai.get('assembly_status'),
                submitter=ai.get('submitter'), sequencing_tech=ai.get('sequencing_tech'),
                annotation_name=an.get('name'), annotation_release_date=an.get('release_date'), annotation_provider=an.get('provider'),
                annotation_method=an.get('method'), annotation_pipeline=an.get('pipeline'), annotation_software_version=an.get('software_version'),
                gene_counts=an.get('stats', {}).get('gene_counts'), busco=dict(lineage=busco.get('busco_lineage'), complete=busco.get('complete'), busco_ver=busco.get('busco_ver')) if busco else None,
                total_sequence_length=rep.get('assembly_stats', {}).get('total_sequence_length'))


def local_md5(raw):
    out = {}
    p = raw / 'local_md5.txt'
    if p.exists():
        for line in open(p):
            f = line.split()
            out[f[1]] = dict(md5=f[0], status=f[-1])
    return out


def prepare(key, row):
    raw = P.raw_dir(row)
    stem = P.stem(row)
    rep = C.parse_assembly_report(raw / (stem + '_assembly_report.txt'))
    main, rule = P.main_sequences(rep)
    names = {P.seq_name(r): r for r in main}
    assert len(names) == len(main) and 'na' not in names, 'unresolved sequence ids'
    seqs = read_fasta(raw / (stem + '_genomic.fna.gz'), set(names))
    assert set(seqs) == set(names), set(names) - set(seqs)
    for n, r in names.items():
        assert len(seqs[n]) == r['length'], (n, len(seqs[n]), r['length'])
    out = P.PREP / key
    out.mkdir(parents=True, exist_ok=True)
    contigs = {n: seqs[n] for n in sorted(names)}
    write_fasta(out / 'genome.fa', contigs)
    lens = {k: len(s) for k, s in contigs.items()}
    json.dump(lens, open(out / 'lengths.json', 'w'), indent=1)
    del seqs, contigs
    # ---- GFF: every feature of the kept sequences, unchanged (no crop) ----
    n_feat, n_gene, dropped_genes, types, biotypes = 0, 0, Counter(), Counter(), Counter()
    with gzip.open(raw / (stem + '_genomic.gff.gz'), 'rt') as g, open(out / 'reference.gff3', 'w') as f:
        f.write('##gff-version 3\n')
        for k, n in lens.items():
            f.write(f'##sequence-region {k} 1 {n}\n')
        for line in g:
            if line[0] == '#':
                continue
            r = line.rstrip('\n').split('\t')
            if len(r) < 9:
                continue
            if r[0] not in lens:
                if r[2] == 'gene':
                    dropped_genes[r[0]] += 1
                continue
            if r[2] == 'region':
                continue
            types[r[2]] += 1
            if r[2] == 'gene':
                n_gene += 1
                biotypes[attrs(r[8]).get('gene_biotype', 'NA')] += 1
            f.write('\t'.join(r) + '\n')
            n_feat += 1
    by_role = Counter()
    for r in rep:
        if P.seq_name(r) not in lens:
            by_role[f"{r['role']}/{r['location']}/{r['unit']}"] += 1
    dropped_bp = sum(r['length'] for r in rep if P.seq_name(r) not in lens)
    md5 = local_md5(raw)
    splits = dict(species=key, organism=row['organism'], accession=row['accession'], assembly=stem, tier=row['tier'], source=row['source'],
                  rule=rule, purpose='independent evaluation panel: all kept sequences are evaluation-only (no training, no validation crop)',
                  main_sequences={P.seq_name(r): dict(length=r['length'], name=r['name'], genbank=r['genbank'], refseq=r['refseq'], role=r['role']) for r in main},
                  test=sorted(lens), train=[], val=[], lengths=lens, test_bp=sum(lens.values()),
                  dropped_sequences=len(rep) - len(main), dropped_bp=dropped_bp, dropped_by_role=dict(by_role),
                  dropped_annotated_genes=sum(dropped_genes.values()), dropped_annotated_genes_by_sequence=dict(dropped_genes.most_common(20)),
                  gff_features=n_feat, gff_genes=n_gene, gff_feature_types=dict(types), gene_biotypes=dict(biotypes),
                  ncbi=ncbi_report(raw), md5={k: v for k, v in md5.items()},
                  annotation_frozen_by=dict(gff_md5=md5.get(stem + '_genomic.gff.gz', {}).get('md5'), fna_md5=md5.get(stem + '_genomic.fna.gz', {}).get('md5')),
                  prepared=datetime.now().astimezone().isoformat(timespec='seconds'))
    json.dump(splits, open(out / 'splits.json', 'w'), indent=1)
    print(key, rule, '| kept', len(lens), 'seqs', splits['test_bp'], 'bp | dropped', splits['dropped_sequences'], 'seqs', dropped_bp, 'bp,',
          splits['dropped_annotated_genes'], 'genes | kept genes', n_gene, dict(biotypes).get('protein_coding'), 'pc', flush=True)
    return splits


def _patched_species_table(orig):
    """Expose the prepared panel genomes to the unchanged M34 Labeler: every kept sequence is a 'val' contig (evaluation), train = []."""
    def table():
        t = orig()
        for key, row in P.species().items():
            p = P.PREP / key / 'splits.json'
            if p.exists():
                s = json.load(open(p))
                t[key] = dict(key=key, organism=row['organism'], accession=row['accession'], order=100 + row['order'], dir=P.PREP / key,
                              train=[], val=list(s['test']), lengths={k: int(v) for k, v in s['lengths'].items()}, val_crop=None, family=row['family'])
        return t
    return table


def stats(key):
    import a1_common as K
    K.repoint()
    C.species_table = _patched_species_table(C.species_table)
    import m34_labels as L
    sp = C.species_table()[key]
    contigs = list(sp['val'])
    lab = L.Labeler(key, contigs, 1.0)
    st = L.stats(lab, {c: 'val' for c in contigs})['val']
    cat = [t for t in lab.catalog.values() if t['assessable']]
    cds = sorted(t['cds_bp'] for t in cat)
    exons = Counter(min(t['n_exons'], 20) for t in cat)
    per_contig = Counter(t['seqid'] for t in cat)
    out = dict(species=key, organism=sp['organism'], accession=sp['accession'], contigs=len(contigs), bp=st['bp'],
               coding_genes_parsed=lab.n_coding_genes, primary_coding_genes=st['primary_coding_genes'], cds_assessable_references=st['cds_assessable_references'],
               assessable_per_Mb=1e6 * st['cds_assessable_references'] / st['bp'],
               complete_supported_nonoverlapping=st['eligible_positive'], exon_hist=st['exon_hist'], exon_count_hist={str(k): v for k, v in sorted(exons.items())},
               median_exons=st['median_exons'], median_cds_bp=(cds[len(cds) // 2] if cds else None),
               cds_bp_hist={'<300': sum(x < 300 for x in cds), '300-1000': sum(300 <= x < 1000 for x in cds), '1000-3000': sum(1000 <= x < 3000 for x in cds), '>=3000': sum(x >= 3000 for x in cds)},
               spans_longer_than_window=st['spans_longer_than_window'], nonstandard_splice_fraction=st['nonstandard_splice_fraction'],
               unsupported_reasons=st['unsupported_reasons'], assessable_per_contig=dict(per_contig),
               label_policy='CDS_assessable_primary_R4 (frozen ruler common_ruler.select_reference; M34 Labeler unchanged)',
               computed=datetime.now().astimezone().isoformat(timespec='seconds'))
    json.dump(out, open(P.PREP / key / 'reference_stats.json', 'w'), indent=1)
    print(json.dumps({k: out[k] for k in ('species', 'bp', 'coding_genes_parsed', 'cds_assessable_references', 'complete_supported_nonoverlapping', 'exon_hist', 'nonstandard_splice_fraction')}), flush=True)
    return out


def main():
    cmd, keys = sys.argv[1], sys.argv[2:]
    rows = P.species()
    keys = keys or list(rows)
    P.PREP.mkdir(parents=True, exist_ok=True)
    for k in keys:
        if cmd == 'prep':
            prepare(k, rows[k])
        elif cmd == 'stats':
            stats(k)
        else:
            raise SystemExit('usage: panel_prep.py prep|stats [key ...]')
    print('PANEL_PREP_DONE', cmd, flush=True)


if __name__ == '__main__':
    main()
