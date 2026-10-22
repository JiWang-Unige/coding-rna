"""PANEL-PREP shared setup: data preparation of the independent plant panel (docs/review_20261016/10 §2.2, 12 species) for Fig. 4 / Table 1.
CPU + network only; no model, no inference, no scoring.  Re-uses the M34 preparation code unchanged (m34_common / m34_prep / m34_labels / panel_exclusion).

Run on ont:  cd ~/coding-rna-reorg/experiments/PANEL-PREP/scripts && taskset -c 25-31 nice -n 10 ~/envs/generanno/bin/python <script>.py
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
M34 = HERE.parents[1] / 'M34-NTV2-SCALEUP' / 'scripts'
if str(M34) not in sys.path:
    sys.path.insert(0, str(M34))
import m34_common as C            # noqa: E402  (sets up a1_common -> pipeline/m28 paths; ROOT = ~/coding-rna-reorg)
import panel_exclusion as PX      # noqa: E402

ROOT = C.ROOT
DATA = ROOT / 'data' / 'panel'                 # -> ~/coding-rna/data/panel
RAW = DATA / 'raw'                             # raw/<accession>_<assembly_name>/ : NCBI files, md5, dataset report, download date
PREP = DATA / 'prepared'                       # prepared/<key>/ : genome.fa, reference.gff3, lengths.json, splits.json (m1_screen / M34 layout)
HOMOLOGY = DATA / 'homology'
MIN_CHROM_BP = C.MIN_CHROM_BP                  # 5 Mb: assembled-molecule chromosomes shorter than this are dropped (M34 rule)
MIN_SCAFFOLD_BP = 10_000_000                   # fallback when the assembly has NO assembled-molecule: sequences >= 10 Mb (M30 rule for Eschscholzia)


def species():
    rows = json.load(open(HERE / 'panel_species.json'))
    return {r['key']: r for r in sorted(rows, key=lambda r: r['order'])}


def stem(row):
    return f"{row['accession']}_{row['assembly_name']}"


def raw_dir(row):
    return RAW / stem(row)


def seq_name(r):
    """Sequence id as it appears in the NCBI *_genomic.fna.gz: RefSeq accession for GCF assemblies, GenBank accession for GCA ones."""
    return r['refseq'] if r['refseq'] != 'na' else r['genbank']


def main_sequences(report_rows):
    """Objective nuclear-sequence rule, fixed before any model output and identical to M34/M30:
    (a) if the assembly has assembled-molecule chromosomes: M34 rule (assembled-molecule + Chromosome + Primary Assembly + >= 5 Mb);
    (b) otherwise (scaffold-level, e.g. Eschscholzia): every Primary Assembly sequence >= 10 Mb (M30 rule).
    Returns (rows, rule_name)."""
    chroms = C.main_chromosomes(report_rows)
    if chroms:
        return chroms, 'assembled-molecule+Chromosome+Primary Assembly+>=5Mb (M34)'
    big = [r for r in report_rows if r['unit'] == 'Primary Assembly' and r['length'] >= MIN_SCAFFOLD_BP]
    return big, 'no assembled-molecule: Primary Assembly sequences >=10Mb (M30)'


def exposure_sets():
    """Training-side species that the panel must be audited against: M34's 10 training species (8 new + 2 legacy), the M30 adaptation species
    (Dryas, Eschscholzia: their label-set-A chromosomes carried labels in M30), and the project-sealed Setaria."""
    rows = C.new_species()
    train = [dict(key=k, genus=r['organism'].split()[0], species=r['organism'], accession=r['accession'], family=C.FAMILY[k], role='M34 training (new)') for k, r in rows.items()]
    train += [dict(key=k, genus=v['organism'].split()[0], species=v['organism'], accession=v['accession'], family=v['family'], role='M34 training (legacy: M28/M31)') for k, v in C.LEGACY.items()]
    m30 = [dict(key='dryas_octopetala', genus='Dryas', species='Dryas octopetala', accession='GCF_963921425.1', family='Rosaceae', role='M30 adaptation (label set A)'),
           dict(key='eschscholzia_californica', genus='Eschscholzia', species='Eschscholzia californica', accession='GCF_041429995.1', family='Papaveraceae', role='M30 adaptation (label set A)')]
    sealed = [dict(key='setaria', genus='Setaria', species='Setaria (project sealed)', accession=None, family='Poaceae', role='project-sealed (never read)')]
    return dict(m34_training=train, m30_adaptation=m30, sealed=sealed)
