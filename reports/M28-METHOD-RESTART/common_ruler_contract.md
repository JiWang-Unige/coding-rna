# M28 R3 common evaluation ruler — frozen before rescore

2026-09-30. Reuse M26 CDS-assessable semantics, not a new outcome-dependent rule.
DEV only: Arabidopsis NC_003074.8 (23,459,830 bp), rice NC_089041.1
(29,936,421 bp). No test/Setaria sequence or annotation inspection.

Select protein-coding transcripts with no partial flag on their own CDS rows,
first transcription-order phase 0, ACGT-only complete ATG-to-stop CDS,
modulo-three length, no internal in-frame stops and no overlapping CDS pieces.
Ignore parent gene/transcript partial tags for this CDS-specific endpoint.
Per gene select longest assessable CDS, tie transcript ID. Phase consistency is
a diagnostic stop, not a reason to silently drop a difficult reference.
Nonstandard splice, split terminal codon, overlapping genes and long spans remain
in the evaluation denominator, even if M28 training cannot support them.

Use unchanged cached M25R A/B, ANNEVO, Helixer and Tiberius predictions. Replay
the historical 4151/2299 parent-filtered chain counts exactly first. New primary
metric is unique exact CDS-chain precision/recall/F1, with separate per-species
and micro-pooled counts. Secondary metrics: unstranded CDS-base union PRF,
alternative-assessable-isoform exact matches (not primary TP), and predicted
span/CDS coverage plus wholly-background chain density outside all annotated
gene-feature spans on either strand. Background disagreement is not biological
false-gene proof. No threshold-penalized metric is presented as ordinary accuracy.

One 2-CPU/8-GiB/10-minute job; no GPU, new baseline runs, parameter search or fitting.
Export the reference ledger and compact result JSON. Measure the primary-choice
difference from R2 longest-all policy explicitly before formal fitting.
