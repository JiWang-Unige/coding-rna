# M24 direct structure diagnostic

## Scope audit

- `arabidopsis_thaliana`: 18,585,056 bp; 3,161 complete primary transcripts; 17,131 CDS intervals; 1,019 partial genes excluded from exact metrics.
- `oryza_sativa`: 28,605,474 bp; 1,981 complete primary transcripts; 9,397 CDS intervals; 6 partial genes excluded from exact metrics.

## Structural metrics

Candidate rows are coordinate-only diagnostics. Their strand/phase/transcript/gene fields are `not_applicable`.

| method | species | exact CDS F1 | pseudo/exact CDS-chain F1 | gbF1 | intergenic FPR | gene-count ratio |
|---|---|---:|---:|---:|---:|---:|
| M19_s0 | arabidopsis_thaliana | 0.0691 | 0.0123 | 0.9222 | 0.0133 | 1.114 |
| M19_s0 | oryza_sativa | 0.0531 | 0.0107 | 0.7226 | 0.0151 | 1.740 |
| M19_s1 | arabidopsis_thaliana | 0.1498 | 0.0095 | 0.9358 | 0.0080 | 0.912 |
| M19_s1 | oryza_sativa | 0.0576 | 0.0082 | 0.8038 | 0.0118 | 1.213 |
| M8_3C_s0 | arabidopsis_thaliana | 0.0029 | 0.0000 | 0.8123 | 0.0750 | 0.875 |
| M8_3C_s0 | oryza_sativa | 0.0025 | 0.0000 | 0.6508 | 0.0449 | 1.878 |
| M8_3C_s2 | arabidopsis_thaliana | 0.0117 | 0.0000 | 0.8020 | 0.0745 | 0.812 |
| M8_3C_s2 | oryza_sativa | 0.0082 | 0.0004 | 0.6552 | 0.0434 | 1.694 |
| M8_3C_s4 | arabidopsis_thaliana | 0.0207 | 0.0000 | 0.8003 | 0.0465 | 0.881 |
| M8_3C_s4 | oryza_sativa | 0.0198 | 0.0004 | 0.6290 | 0.0278 | 1.614 |
| ANNEVO | arabidopsis_thaliana | 0.8614 | 0.7479 | 0.9467 | 0.0161 | 0.893 |
| ANNEVO | oryza_sativa | 0.8882 | 0.7324 | 0.9037 | 0.0182 | 1.004 |
| Helixer | arabidopsis_thaliana | 0.8117 | 0.6339 | 0.9520 | 0.0341 | 0.972 |
| Helixer | oryza_sativa | 0.8121 | 0.5850 | 0.8840 | 0.0289 | 1.114 |
| Tiberius | arabidopsis_thaliana | 0.8547 | 0.7254 | 0.9571 | 0.0129 | 0.914 |
| Tiberius | oryza_sativa | 0.8631 | 0.6660 | 0.8842 | 0.0128 | 1.039 |

## SegmentNT released feature cache

Thresholds were selected independently on each species' validation seqid and applied once to its test seqid. The cache uses independent 6,000-bp tiles.

| species | view | element | test AUCPR | test F1 | prevalence |
|---|---|---|---:|---:|---:|
| arabidopsis_thaliana | primary_transcript | exon | 0.6569 | 0.6092 | 0.294158 |
| arabidopsis_thaliana | primary_transcript | intron | 0.2500 | 0.3369 | 0.128167 |
| arabidopsis_thaliana | primary_transcript | splice_donor | 0.0421 | 0.1275 | 0.001578 |
| arabidopsis_thaliana | primary_transcript | splice_acceptor | 0.0443 | 0.1358 | 0.001578 |
| arabidopsis_thaliana | all_isoform_union | exon | 0.6615 | 0.6130 | 0.309117 |
| arabidopsis_thaliana | all_isoform_union | intron | 0.2584 | 0.3429 | 0.134902 |
| arabidopsis_thaliana | all_isoform_union | splice_donor | 0.0425 | 0.1296 | 0.001668 |
| arabidopsis_thaliana | all_isoform_union | splice_acceptor | 0.0448 | 0.1376 | 0.001679 |
| oryza_sativa | primary_transcript | exon | 0.5866 | 0.5490 | 0.126588 |
| oryza_sativa | primary_transcript | intron | 0.3243 | 0.3865 | 0.137100 |
| oryza_sativa | primary_transcript | splice_donor | 0.0314 | 0.1053 | 0.000566 |
| oryza_sativa | primary_transcript | splice_acceptor | 0.0388 | 0.1252 | 0.000566 |
| oryza_sativa | all_isoform_union | exon | 0.5875 | 0.5487 | 0.128900 |
| oryza_sativa | all_isoform_union | intron | 0.3265 | 0.3870 | 0.140923 |
| oryza_sativa | all_isoform_union | splice_donor | 0.0319 | 0.1068 | 0.000591 |
| oryza_sativa | all_isoform_union | splice_acceptor | 0.0390 | 0.1263 | 0.000600 |

## Interpretation boundary

This report measures saved artifacts only. Candidate `+` strand and phase `0` placeholders are not ranked as structural predictions. A weak SegmentNT row applies only to the existing 6-kb tiled cache, not to a longer-context extraction or the checkpoint in general.
