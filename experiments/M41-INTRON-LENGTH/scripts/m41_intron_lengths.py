"""Training-chromosome intron lengths (bp) of the eligible catalog transcripts of the two T2 species -> JSON for `m41_prior.py fit`.
Only the four M31/M34 TRAIN chromosomes are read (DEV, validation and never-used chromosomes are not).  Run on ont (CPU, seconds):
  nice -n 19 python3 m41_intron_lengths.py ~/coding-rna-reorg/data/m34_train/catalog > ../results/train_intron_lengths.json
"""
import collections
import json
import sys

TRAIN = {'arabidopsis_thaliana': ['NC_003070.9', 'NC_003076.8'], 'oryza_sativa': ['NC_089043.1', 'NC_089046.1']}


def main():
    base = sys.argv[1].rstrip('/') + '/'
    out = {}
    for sp, cs in TRAIN.items():
        lengths, n, bins = [], 0, collections.Counter()
        for line in open(base + sp + '/transcripts.jsonl'):
            t = json.loads(line)
            if t['seqid'] not in cs or not t.get('eligible'):
                continue
            c = sorted(tuple(x) for x in t['CDS'])
            n += 1
            bins[min(len(c), 6)] += 1
            lengths += [c[i + 1][0] - c[i][1] for i in range(len(c) - 1)]
        out[sp] = dict(n_tx=n, exon_bins=dict(bins), introns=lengths)
    json.dump(out, sys.stdout)


if __name__ == '__main__':
    main()
