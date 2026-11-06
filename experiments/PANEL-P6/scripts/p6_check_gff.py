"""PANEL-P6: sanity check of one native prediction GFF against its input FASTA (stdlib; same checks as the M28 DEV runner `run_native_baseline.check_output`, for a
whole genome): 9 columns, seqid in the FASTA, 1 <= start <= end <= length, CDS strand / phase / Parent valid.  Prints a JSON summary; exit 1 on a violation.
usage: python3.12 -I p6_check_gff.py GFF FASTA [MAIN_LENGTHS_JSON]
"""
import json
import sys
from collections import Counter


def fasta_lengths(path):
    out, cur, n = {}, None, 0
    with open(path) as f:
        for line in f:
            if line[0] == '>':
                if cur is not None:
                    out[cur] = n
                cur, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if cur is not None:
        out[cur] = n
    return out


def main():
    gff, fasta = sys.argv[1], sys.argv[2]
    main_seqs = set(json.load(open(sys.argv[3]))) if len(sys.argv) > 3 else None
    L = fasta_lengths(fasta)
    counts, genes_main, seqs = Counter(), 0, Counter()
    with open(gff) as h:
        for n, line in enumerate(h, 1):
            if not line.strip() or line[0] == '#':
                continue
            r = line.rstrip('\n').split('\t')
            if len(r) != 9:
                sys.exit(f'{gff}:{n}: {len(r)} columns')
            if r[0] not in L:
                sys.exit(f'{gff}:{n}: seqid {r[0]} not in the FASTA')
            s, e = int(r[3]), int(r[4])
            if not 1 <= s <= e <= L[r[0]]:
                sys.exit(f'{gff}:{n}: coordinates outside the sequence')
            if r[2] == 'CDS' and (r[6] not in '+-' or r[7] not in ('0', '1', '2') or not any(a.startswith('Parent=') and a[7:] for a in r[8].split(';'))):
                sys.exit(f'{gff}:{n}: CDS strand / phase / Parent invalid')
            counts[r[2]] += 1
            if r[2] == 'gene':
                seqs[r[0]] += 1
                genes_main += main_seqs is not None and r[0] in main_seqs
    print(json.dumps(dict(gff=gff, features=dict(counts), sequences_with_genes=len(seqs), fasta_sequences=len(L), genes=counts['gene'],
                          genes_on_main_sequences=genes_main if main_seqs is not None else None)))


if __name__ == '__main__':
    main()
