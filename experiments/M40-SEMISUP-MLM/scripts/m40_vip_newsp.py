"""Vipsania external 0 % point on the 8 new-species internal-validation pseudo-contigs (protocol 32 section 6, decision m40-vipsania-newsp, user 2026-11-11).
  fasta  OUT.fa      write the 8 '<chrom>__val8M' pseudo-contigs (exactly the sequences the M34 validation reads) as one FASTA, headers = pseudo-contig names
  regions            print the comma-separated pseudo-contig list (= the own_species regions of every M34 / M40 arm)
  concat RUN OUT.gff concatenate RUN/<contig>/prediction/full_grammar.gff3 of the 8 contigs (CPU check of the convert + score path on an arm's own predictions)
Vipsania annotates the 8 Mb crops only (pretrained-only mode chunks the sequence into independent 200-kb blocks, so a crop changes nothing but the ends)."""
import sys

import m40_common as W

C, T = W.C, W.T


def contigs():
    return [T.val_contig(k) for k in W.NEW]


def main():
    cmd = sys.argv[1]
    if cmd == 'regions':
        print(','.join(contigs()))
    elif cmd == 'fasta':
        tab = C.species_table()
        with open(sys.argv[2], 'w') as f:
            for k in W.NEW:
                c = T.val_contig(k)
                s = T.MF.read_allowed(C.PREP / k / 'genome.fa', {c: int(tab[k]['lengths'][c])})[c]
                f.write(f'>{c}\n')
                for i in range(0, len(s), 80):
                    f.write(s[i:i + 80] + '\n')
    elif cmd == 'concat':
        with open(sys.argv[3], 'w') as out:
            out.write('##gff-version 3\n')
            for c in contigs():
                for line in open(f'{sys.argv[2]}/{c}/prediction/full_grammar.gff3'):
                    if not line.startswith('#'):
                        out.write(line)
    else:
        raise SystemExit('unknown command')


if __name__ == '__main__':
    main()
