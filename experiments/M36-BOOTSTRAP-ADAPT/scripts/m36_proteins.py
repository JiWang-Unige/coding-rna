"""Build the cross-species protein query set for the S2 arm: one RefSeq protein per protein-coding gene (the isoform with the longest CDS) of each of
the 10 M34 training species.  Legacy species (Arabidopsis / rice) contribute only genes on their M34 training + internal-validation chromosomes
(never the DEV / split_test chromosomes); organelle contigs (GFF region attribute genome=chloroplast|mitochondrion|plastid) are dropped; proteins
whose CDS carries partial=true are dropped.  The target species' own proteins are never used (protocol 22 §2).
Output: data/m36/proteins/queries_10sp.faa  (headers '<species_key>|<protein_id>')  +  queries_stats.json
usage: m36_proteins.py
"""
import gzip
import json
import re
import sys
import time
from collections import defaultdict

import m36_common as W


def opener(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def attrs(s):
    return dict(x.split('=', 1) for x in s.rstrip(';').split(';') if '=' in x)


def select_proteins(gff, keep_contigs=None):
    """-> ({protein_id: (gene_id, cds_bp)} longest-CDS isoform per gene, stats)"""
    organelle = set()
    rna2gene, cds_len, partial, contig_of = {}, defaultdict(int), set(), {}
    with opener(gff) as f:
        for line in f:
            if line[0] == '#':
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9:
                continue
            typ = p[2]
            if typ == 'region':
                a = attrs(p[8])
                if a.get('genome') in ('chloroplast', 'mitochondrion', 'plastid', 'apicoplast'):
                    organelle.add(p[0])
                continue
            if typ == 'CDS':
                a = attrs(p[8])
                pid = a.get('protein_id')
                if not pid:
                    continue
                parent = a.get('Parent', '').split(',')[0]
                cds_len[pid] += int(p[4]) - int(p[3]) + 1
                contig_of[pid] = p[0]
                if a.get('partial') == 'true':
                    partial.add(pid)
                rna2gene.setdefault(pid, parent)
            elif typ in ('mRNA', 'transcript'):
                a = attrs(p[8])
                if 'ID' in a:
                    rna2gene[a['ID']] = a.get('Parent', a['ID'])
    per_gene = {}
    n_all, n_contig, n_org, n_partial = 0, 0, 0, 0
    for pid, bp in cds_len.items():
        n_all += 1
        c = contig_of[pid]
        if c in organelle:
            n_org += 1
            continue
        if keep_contigs is not None and c not in keep_contigs:
            n_contig += 1
            continue
        if pid in partial:
            n_partial += 1
            continue
        parent = rna2gene.get(pid, pid)                 # rna-XM_... (or the gene itself when the CDS hangs directly under the gene)
        gene = rna2gene.get(parent, parent)
        cur = per_gene.get(gene)
        if cur is None or bp > cur[1] or (bp == cur[1] and pid < cur[0]):
            per_gene[gene] = (pid, bp)
    chosen = {pid: (gene, bp) for gene, (pid, bp) in per_gene.items()}
    return chosen, dict(proteins_in_gff=n_all, dropped_other_contigs=n_contig, dropped_organelle=n_org, dropped_partial=n_partial,
                        genes=len(per_gene), selected=len(chosen), organelle_contigs=sorted(organelle))


def read_faa(path, wanted):
    out, name, buf = {}, None, []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line[0] == '>':
                if name in wanted:
                    out[name] = ''.join(buf)
                name, buf = line[1:].split()[0], []
            elif name in wanted:
                buf.append(line.strip())
    if name in wanted:
        out[name] = ''.join(buf)
    return out


def main():
    W.pin_and_nice()
    t0 = time.time()
    outd = W.DATA / 'proteins'
    outd.mkdir(parents=True, exist_ok=True)
    stats = {}
    n_written = 0
    with open(outd / 'queries_10sp.faa.tmp', 'w') as out:
        for src in W.protein_sources():
            chosen, st = select_proteins(src['gff'], src['keep_contigs'])
            seqs = read_faa(src['faa'], set(chosen))
            missing = len(chosen) - len(seqs)
            n_sp = 0
            for pid in sorted(seqs):
                s = re.sub('[^A-Z]', 'X', seqs[pid].upper()).rstrip('*')
                if len(s) < 30:
                    continue
                out.write(f">{src['key']}|{pid}\n")
                for i in range(0, len(s), 80):
                    out.write(s[i:i + 80] + '\n')
                n_sp += 1
            n_written += n_sp
            st.update(written=n_sp, missing_in_faa=missing, residues=sum(len(seqs[p]) for p in seqs), organism=src['organism'],
                      faa=str(src['faa']), gff=str(src['gff']), keep_contigs=sorted(src['keep_contigs']) if src['keep_contigs'] else 'all nuclear')
            stats[src['key']] = st
            print(src['key'], json.dumps({k: st[k] for k in ('proteins_in_gff', 'genes', 'written', 'dropped_other_contigs', 'dropped_organelle', 'dropped_partial')}), flush=True)
    (outd / 'queries_10sp.faa.tmp').rename(outd / 'queries_10sp.faa')
    json.dump(dict(total_written=n_written, rule='one protein per gene (longest CDS); legacy species restricted to M34 train+val chromosomes; organelle and partial CDS dropped',
                   species=stats, seconds=time.time() - t0), open(outd / 'queries_stats.json', 'w'), indent=1)
    print('total', n_written, f'{time.time() - t0:.0f}s')


if __name__ == '__main__':
    main()
