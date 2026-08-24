# Dossier: Helixer: ab initio prediction of primary eukaryotic gene models combining deep learning and a hidden Markov model

- slug: `helixer-2025` · type: sota · added: 2026-06-09
- Paper: https://www.nature.com/articles/s41592-025-02939-1
- PDF: `refs/pdfs/helixer-2025.pdf`
- Repo: `refs/repos/helixer-2025/` cloned @ `f7eb4dd`
- GitHub: https://github.com/usadellab/Helixer
- Weights: https://zenodo.org/records/17850139
- Why relevant: required open HMM-postprocessing baseline.

## Dataset source

- Fungi: RefSeq 2022-03-04.
- Plants: train/val from Phytozome13 2021-06-07; test from RefSeq 2022-07-14 excluding train/val species.
- Vertebrates and invertebrates: RefSeq 2022-05-06.
- Mammal model: RefSeq 2025-03-13, following Tiberius mammal species selection/partition.
- `docs/model_overview.md` lists training and validation species per lineage.

## Metric implementation

- Main paper reports phase F1 and feature-level exon/intron/intron-chain/transcript F1.
- Repo comparison docs use gffcompare for annotation-level precision/recall/F1 and `scripts/accs_genic_intergenic.py` for base-wise genic class / coding phase / subgenic class metrics.
- Label classes include intergenic, UTR, CDS, intron, plus coding phase labels.

## Split scheme

- Species-level split, not random sequence split.
- Validation species use random subsets of subsequences for several lineages; exact species are in `docs/model_overview.md`.
- Final benchmark contract must freeze model version and species list because Zenodo has older and current model records.

## Weights / license

- Repo license: GPL-3.0.
- Weights: Zenodo `10.5281/zenodo.17850139`, v2, 2025-12-07, CC BY 4.0; includes fungi, land_plant, vertebrate, invertebrate and mammal `.h5` models.
- Older model list in repo also points to Zenodo `10836346`; do not mix weight records in one benchmark.

## Reported values

- Main Table 1 phase F1, HelixerPost group means: fungi 0.9540, plants 0.8099, vertebrates 0.8829, invertebrates 0.8562.
- Main Table 2 feature F1:
  - fungi exon/intron/intron-chain/transcript: 0.6678 / 0.6061 / 0.4431 / 0.5386
  - plants: 0.7143 / 0.7232 / 0.4338 / 0.4618
  - vertebrates: 0.7405 / 0.6912 / 0.1740 / 0.1977
  - invertebrates: 0.6608 / 0.6416 / 0.2939 / 0.3066

## Reproducibility notes

- Reproducibility: `moderate`.
- Code, weights, and commands exist, but full benchmark reproduction spans paper, supplement, Zenodo, and repo notebooks/scripts.

## Relevance to our project

- Main baseline for two-stage base-wise prediction plus HMM postprocessing.
- Important source for HMM states, transition penalties, minimum intron lengths, and base-wise class/phase evaluation.
