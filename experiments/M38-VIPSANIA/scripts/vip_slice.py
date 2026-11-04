"""M38-VIPSANIA: cut a slice (1-based inclusive) out of a contig of a FASTA file, keeping the soft-masking, for the CPU / GPU smoke tests.
Usage: vip_slice.py GENOME.fa CONTIG START1 END1 OUT.fa      (header: >CONTIG__START1_END1 ; sidecar OUT.fa.json = contig / start1 / end1)"""
import json
import sys


def main():
    fa, contig, s1, e1, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    parts, cur = [], None
    with open(fa) as f:
        for line in f:
            if line[0] == '>':
                cur = line[1:].split()[0]
                if parts and cur != contig:
                    break
            elif cur == contig:
                parts.append(line.strip())
    seq = ''.join(parts)[s1 - 1:e1]
    assert len(seq) == e1 - s1 + 1, (contig, s1, e1, len(seq))
    name = f'{contig}__{s1}_{e1}'
    with open(out, 'w') as w:
        w.write(f'>{name}\n')
        for i in range(0, len(seq), 80):
            w.write(seq[i:i + 80] + '\n')
    json.dump(dict(name=name, contig=contig, start1=s1, end1=e1, length=len(seq), lowercase_fraction=round(sum(c.islower() for c in seq) / len(seq), 4)), open(out + '.json', 'w'))
    print(name, len(seq))


if __name__ == '__main__':
    main()
