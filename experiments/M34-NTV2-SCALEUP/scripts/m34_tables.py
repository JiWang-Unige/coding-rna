"""Collect catalog summaries / schedule summaries / download manifest into results/ (tsv + json) for the report.  CPU only.
usage: m34_tables.py <results_dir>"""
import json
import sys
from pathlib import Path

import m34_common as C
from m34_schedule import ARMS, LEGACY, NEW

res = Path(sys.argv[1])
res.mkdir(parents=True, exist_ok=True)
rows = C.new_species()
order = LEGACY + NEW
lines = ['species\taccession\ttrain_contigs\ttrain_Mb\tval_contig\tval_Mb\tprimary_coding_genes_train\tcds_assessable_train\teligible_train\teligible_per_Mb\tcds_assessable_val\teligible_val\texon1\texon2_5\texon6p\tmedian_exons\tspans_gt_window\tnonstd_splice_frac\twindows_train\twindows_val']
tot = dict(eligible=0, ass=0, bp=0, val=0)
for k in order:
    s = json.load(open(C.CATALOG / k / 'summary.json'))
    t, v = s['stats']['train'], s['stats']['val']
    n_train = sum(1 for c in s['split_of'].values() if c == 'train')
    vc = [c for c, x in s['split_of'].items() if x == 'val'][0]
    lines.append('\t'.join(str(x) for x in (k, s['accession'], n_train, round(t['bp'] / 1e6, 1), vc, round(v['bp'] / 1e6, 1), t['primary_coding_genes'], t['cds_assessable_references'],
        t['eligible_positive'], round(t['eligible_per_Mb'], 1), v['cds_assessable_references'], v['eligible_positive'], t['exon_hist']['1'], t['exon_hist']['2-5'], t['exon_hist']['>=6'],
        t['median_exons'], t['spans_longer_than_window'], round(t['nonstandard_splice_fraction'], 4), s['windows']['train'], s['windows']['val'])))
open(res / 'species_supervision.tsv', 'w').write('\n'.join(lines) + '\n')
p = C.OUT / 'schedules' / 'schedule_summary.json'
if p.exists():
    sch = json.load(open(p))
    L = ['arm\tspecies\tnew_k\ttreatment\tdraws\tunique_windows\tlabeled_genes_total\tlabeled_genes_exposed_unique\tgene_exposures_total\tgene_first_fraction\twindows_with_positive_fraction']
    for a, s in sch.items():
        L.append('\t'.join(str(x) for x in (a, len(s['species']), s['new_species_label_fraction'], s['treatment'], s['draws'], s['unique_windows'], sum(s['labeled_genes'].values()),
                 s['labeled_genes_exposed_unique'], s['labeled_gene_exposures_total'], round(s['gene_first_fraction'], 3), round(s['windows_with_positive_fraction'], 3))))
    open(res / 'arms_exposure.tsv', 'w').write('\n'.join(L) + '\n')
    json.dump(sch, open(res / 'schedule_summary.json', 'w'), indent=1)
man = ['species\tfile\tmd5\tncbi_md5_ok\tdownload_date']
for k, r in rows.items():
    d = C.raw_dir(r)
    date = open(d / 'download_date.txt').readline().strip()
    for line in open(d / 'local_md5.txt'):
        f = line.split()
        man.append('\t'.join((k, f[1], f[0], f[3], date)))
open(res / 'download_manifest.tsv', 'w').write('\n'.join(man) + '\n')
print(open(res / 'species_supervision.tsv').read())
