# Current research state — 2026-09-01

## Executive conclusion

The project has established a useful but limited positive result: the CDS-specialized `GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview` backbone can be adapted into a low-FPR coarse coding-region detector on the Arabidopsis/rice development panel. The strongest frozen M19 seed reaches gene-body F1 `0.8815`, intergenic FPR `0.0065`, and predicted-gene-count ratio `0.830` after validation-only calibration.

That result is not yet a direct gene annotator. M24 showed that M19's saved coordinate candidates have low exact CDS-interval F1 (`0.0531–0.1498`) and near-zero coordinate pseudo-chain F1 (`0.0082–0.0123`), while released callers on the identical held-out seqids reach exact CDS-interval F1 `0.8117–0.8882` and exact CDS-chain F1 `0.5850–0.7479`. Strand and phase in M19 are placeholders and cannot be scored as model outputs.

M25R was the first numerically valid attempt to turn GENERanno 1.2B features into complete, strand-aware and phase-aware CDS models through explicit region, boundary and phase heads plus a constrained decoder. It is a valid development-stage **no-go under the frozen admission contract**: all `5,625` checkpoint/threshold tuples failed the predicted-gene-count constraint (`0.80–1.20`), with observed ratios only `0.0953–0.3305`. Setaria inference was therefore not run. This does not prove that the GENERanno representation is intrinsically inadequate; it shows that the present combination of representation, structural supervision and decoder under-recovers complete models.

The project is not yet ready for a Nature Communications claim. The immediate task is a read-only development-set error decomposition using the existing M25R checkpoints. No new training, Setaria annotation access, SegmentNT promotion, extra seeds or baseline reruns are justified until that diagnostic identifies where reference genes disappear.

## Frozen scientific scope

- Product: raw nuclear genome FASTA to primary protein-coding gene annotation.
- Current organisms in structural development: *Arabidopsis thaliana* and *Oryza sativa*.
- Frozen blind organism: *Setaria viridis* assembly `GCF_005286985.2`.
- Included sequence scope: predefined primary nuclear chromosomes only.
- Excluded: chloroplast, mitochondrion and unplaced/unlocalized scaffolds.
- Blind-test rule: do not read or use Setaria annotation until model, checkpoint, decoder and full/ablation predictions are frozen.
- GENERanno provenance: generic pretraining and CDS post-training species/accession overlap remain `overlap_unknown`; no clean no-overlap claim is allowed.

## Research progression

| Phase | Main question | Decisive result | Interpretation |
|---|---|---|---|
| M9–M11 | Can a DNA foundation model produce coherent coding-region spans with strict FP control? | M11 validation-only calibration: gene-body F1 `0.8178`, FPR `0.0087`, count ratio `1.003`. | Coarse region detection and operating-point control are feasible. |
| M12–M17 | Does the NT-v2 route transfer as a universal fixed caller, and how does it compare with released tools? | Arabidopsis-to-rice M12A: F1 `0.6556`, FPR `0.0311`, count `1.755`; close-plant M13 and animal M14 also fail. M17 released callers remain much stronger. | Current NT-v2 model is not a broad fixed cross-species caller. |
| M15–M19 | Does GENERanno CDS specialization provide a stronger representation? | M18 1.2B: F1 `0.8494`, FPR `0.0071`, count `0.864`; M19 calibrated seeds remain low-FPR. GENERanno 0.5B base is markedly worse. | CDS specialization carries useful coarse coding signal; the evidence is adaptation evidence, not a clean held-out claim. |
| M20–M23 | Can calibration, CRF or a scalar objective close the remaining gap? | M20 oracle replay finds little unused calibration gain; M21 CRF worsens FPR; M22 Tversky fails; M23 reproduces the older NT-v2 result. | More threshold tuning, CRF retuning and scalar-loss swaps are low-value directions. |
| M24 | Are existing outputs already structurally exact, and does cached SegmentNT provide boundary signal? | M19 exact interval/chain are far below ANNEVO, Helixer and Tiberius. SegmentNT 6-kb cache has exon AUCPR `0.5866–0.6569` but donor/acceptor AUCPR only `0.0314–0.0443`. | The project needs explicit structural prediction. The cache result does not reject SegmentNT generally. |
| M25 | Can explicit structural heads produce full models? | Four of 1,536 training windows had empty structural masks; the old boundary-loss reduction produced non-finite losses. | Implementation-invalid; old checkpoints cannot support scientific inference. |
| M25R | Does the minimally repaired, frozen structural experiment pass development admission? | Finite training and complete grid, but `0/5,625` admissible because every tuple under-predicts complete genes. | Valid no-go for this exact pipeline; bottleneck location remains unresolved. |

## Decisive quantitative evidence

| Result | Exact or coarse metric | FPR | Gene-count ratio | Status |
|---|---:|---:|---:|---|
| M19 GENERanno 1.2B calibrated s0 | gene-body F1 `0.8421` | `0.0083` | `1.083` | positive coarse adaptation |
| M19 GENERanno 1.2B calibrated s1 | gene-body F1 `0.8815` | `0.0065` | `0.830` | strongest coarse adaptation |
| M23 NT-v2 clean-provenance s0 | gene-body F1 `0.8427` | `0.01673` | `0.867` | valid but not frontier-moving |
| M24 M19 candidates | exact CDS interval F1 `0.0531–0.1498`; pseudo-chain F1 `0.0082–0.0123` | `0.0080–0.0151` | `0.912–1.740` | not complete structural output |
| M24 ANNEVO | exact interval `0.8614–0.8882`; chain `0.7324–0.7479` | `0.0161–0.0182` | `0.893–1.004` | same-scope structural baseline |
| M24 Helixer | exact interval `0.8117–0.8121`; chain `0.5850–0.6339` | `0.0289–0.0341` | `0.972–1.114` | same-scope structural baseline |
| M24 Tiberius | exact interval `0.8547–0.8631`; chain `0.6660–0.7254` | `0.0128–0.0129` | `0.914–1.039` | same-scope structural baseline |
| M25R best-ranked development tuple | exact interval `0.1204`; chain `0.3250` | `0.01247` | `0.3253` | rejected by count constraint |

## M25R terminal result

- Slurm job: `12116383`; state `COMPLETED`; exit `0:0`; elapsed `1-13:03:39`.
- Experiment status: `STOP_M25_BRANCH`.
- Training loss: `0.8205 -> 0.6669 -> 0.6177` across three epochs.
- Checkpoints: all three epoch checkpoints exist on Baobab and contain finite tensors; they are intentionally not versioned in Git.
- Validation grid: `3` epochs × `3` region thresholds × `5^4` boundary thresholds = `5,625` complete finite rows.
- All `5,625` rows pass intergenic FPR `<=0.020`.
- All `5,625` rows fail gene-count ratio `0.80–1.20`.
- Exact CDS interval F1 range: `0.03509–0.12057`.
- Exact CDS-chain F1 range: `0.12733–0.32499`.
- Intergenic FPR range: `0.00303–0.01266`.
- Predicted-gene-count ratio range: `0.09535–0.33054`.
- Recorded structural-valid fraction is always `1.0`, but the current implementation makes this non-independent and therefore uninformative; it must be replaced by per-transcript validation in the diagnostic.
- Best-ranked tuple: epoch `1`, region `0.4`, start `0.5`, stop `0.5`, donor `0.1`, acceptor `0.1`; interval F1 `0.12041`, chain F1 `0.32499`, FPR `0.01247`, count ratio `0.32527`.
- Later epochs reduce training loss but worsen complete-gene recovery; the best-ranked count ratio drops from `0.3253` at epoch 1 to `0.1829` at epoch 2 and `0.1034` at epoch 3.

Because no validation tuple was admissible, the following do **not** exist and must not be reported as zero:

- selected checkpoint;
- frozen validation decoder parameters;
- Setaria raw scores or full/ablation GFF3 predictions;
- embargo release marker;
- Setaria exact interval, chain, gene, strand or phase metrics;
- full-versus-unchanged-input ablation differences.

The Setaria annotation embargo remains intact.

## What is established

1. GENERanno 1.2B CDS specialization is empirically more useful than the tested 0.5B generic base checkpoint under this project's adaptation setup.
2. Validation-only calibration can control coarse intergenic FPR, but it does not create exact gene structures.
3. Released callers produce much more exact CDS intervals and chains on the same Arabidopsis/rice held-out ranges.
4. M25R training is numerically valid; its terminal result is not attributable to the M25 empty-mask bug.
5. The current M25R pipeline strongly under-recovers complete genes even while FPR is acceptable.
6. The blind Setaria test has not been consumed.

## What remains unresolved

- **Representation:** whether the backbone scores contain recoverable start/stop/donor/acceptor and strand/phase information before decoding.
- **Structural supervision:** whether sparse boundary and phase targets are learned well enough, and why performance degrades with later epochs.
- **Decoder:** where candidates are lost during region formation, event snapping, transition constraints, complete-model filtering and uniqueness filtering.
- **Transfer:** whether a development-passing model transfers to Setaria; this cannot be answered before development admission.
- **SegmentNT:** whether longer-context extraction or direct plant adaptation provides stronger structural emissions than the existing independent 6-kb cache.
- **Provenance:** the species/accession composition of GENERanno generic pretraining and CDS post-training.
- **Publication relevance:** whether a future method provides a reproducible advantage over Tiberius, Helixer and ANNEVO rather than merely approaching their accuracy.

## Current decision

The GENERanno structural branch is **no-go for new training under the present evidence**, but not rejected as a backbone. The only approved next analysis is `M25R-DEV-REDECODE-ERROR-DECOMPOSITION`: re-decode existing epoch checkpoints on Arabidopsis/rice development data, reproduce frozen aggregates, account for every reference gene and emitted prediction through decoder stages, validate structures independently, and then stop for review. The diagnostic must not read Setaria annotation or change any model parameter.

The decision tree and publication route are defined in `docs/29_end_to_end_technical_roadmap_2026-09-01.md`.

## Primary evidence paths

- `reports/M24-DIRECT-STRUCTURE-DIAGNOSTIC/report.md`
- `reports/M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0/terminal_summary.md`
- `reports/M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0/validation_grid_summary.json`
- `docs/25_direct_annotator_execution_plan.md`
- `docs/26_m25_generanno_structural_heads_execution_plan.md`
- `docs/27_m25r_implementation_repair_plan.md`
- `reports/M19-COMPARABILITY-EVIDENCE/comparison_tables.md`
- `reports/M20-SOTA-ERROR-ANALYSIS/report.md`
