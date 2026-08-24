# Dossier: Highly accurate ab initio gene annotation with ANNEVO

- slug: `annevo-2026` · type: sota · added: 2026-06-09
- Paper: https://www.nature.com/articles/s41592-026-03036-7
- Supplement/PDF: `refs/pdfs/annevo-2026.pdf`
- Repo: `refs/repos/annevo-2026/` cloned @ `37bdd9a`
- GitHub: https://github.com/xjtu-omics/ANNEVO
- Why relevant: strongest full-gene SOTA-anchor candidate for cross-species eukaryotic ab initio annotation.

## Dataset source

- RefSeq / Ensembl-derived eukaryotic gene annotations; article reports broad benchmark across 566 phylogenetically diverse species.
- Current repo performance table uses 12 model species across Mammalia, Insecta, Fungi, Magnoliopsida, Aves, and Actinopteri.
- Exact frozen train/test species lists are in supplementary tables and still need expansion before final benchmark locking.

## Metric implementation

- Paper reports nucleotide-level mean F1 / recall / precision; repo current release reports exon recall/precision, locus recall/precision, and BUSCO.
- Repo evaluation notes use `gffcompare -r ${path_to_ref} ${path_to_pred} --no-exon-merge --strict-match`.
- Repo notes exclude non-coding transcripts and UTRs and remove invalid start/stop codons or erroneous intron lengths.
- For exon-level evaluation, the longest transcript per gene is used as ground truth; for locus-level evaluation, a predicted transcript is correct if it strictly matches any transcript from the same reference gene.

## Split scheme

- Repo notes state ANNEVO was intentionally trained without including the listed model species; none of the test species appeared in training or validation.
- Other tools in the ANNEVO repo comparison may have overlap with training/validation species and are marked as such in repo tables.
- Before claim, expand supplementary species tables and freeze exact train/validation/test species.

## Weights / license

- Weights are in repo under `saved_model/*.pt`: Mammalia, Insecta, Aves, Actinopteri, Magnoliopsida, Fungi.
- License: ANNEVO Non-Commercial License; academic/non-profit use allowed; commercial use prohibited without separate license; not OSI-approved.

## Reported values

- Paper: nucleotide-level mean F1 = 0.92, recall = 0.922, precision = 0.919 across 12 model species.
- Paper: 43 mammalian species comparison reports mean improvements over Tiberius of NT(CDS)-F1 +5.9%, gene-F1 +1.0%, BUSCO +3.5%.
- Current repo release table: ANNEVO average exon recall 91.4, exon precision 90.2, locus recall 76.3, locus precision 74.3, BUSCO 97.8; Tiberius 89.8/88.7/74.0/68.8/96.3; Helixer 86.1/75.3/50.2/47.0/92.5.

## Reproducibility notes

- Reproducibility: `moderate`.
- Setup and weights are public, but paper vs current GitHub release metrics/runtime have drifted.
- Treat paper values as published anchor candidates; treat current repo values as a separately pinned reproduction target.

## Relevance to our project

- Candidate published SOTA anchor.
- Main ideas to inspect: long/distal modeling, clade/evolution MoE routing, resolution restoration, Viterbi decoding, and strict transcript/locus evaluation.
