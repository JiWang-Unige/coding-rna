# M25R implementation-repair execution plan

Status: pending final ChatGPT Pro approval. Do not submit before approval.

## Decision

Job 12094731 ended `COMPLETED` with experiment status `STOP_M25_BRANCH`, but it is not a scientific no-go. Four of the 1,536 frozen training windows have an all-false `structural_mask`; taking the mean of the empty boundary-loss tensor produced non-finite batch losses. The run is implementation-invalid. Its three checkpoints remain archived and must not be used for Setaria inference.

The only next experiment is:

`M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0`

M25R repeats the same registered scientific experiment after a minimal implementation repair. It is not a new architecture or a tuning experiment.

## Allowed code changes

1. In the structural-head loss, keep the region loss active for every window. If `structural_mask.any()` is false, set only the boundary-loss contribution to `boundary_logits.sum() * 0.0`. For non-empty masks, preserve the existing calculation exactly.
2. Before `backward()`, require the total loss to be finite and raise an error otherwise. The error must identify epoch, batch index, and the region/structural/boundary/phase valid counts available from that batch. Do not skip a non-finite batch.
3. Add one regression test containing valid region supervision and an all-false structural mask. Assert finite total loss, zero boundary contribution, computable region gradients, and finite parameters after one optimizer step.
4. Save a validation-grid diagnostic record without changing selection. Each tuple records checkpoint/epoch, thresholds, exact CDS interval F1, exact CDS-chain F1, intergenic FPR, gene-count ratio, structural-valid fraction, each admission-constraint result, admissibility, and frozen ranking/tie-break fields.

No other trainer, decoder, evaluator, data, model, loss, threshold, sampling, or gate change is allowed.

## Frozen contract

- Seed 0; one fit; exactly three epochs.
- GENERanno 1.2B CDS annotator, official 6-mer tokenizer, 6,144-bp windows.
- LoRA r=8, alpha=16, dropout=0.05 on q/k/v/o projections.
- Existing region, boundary and phase heads; existing loss coefficients and learning rates.
- The identical 1,536 training windows and 768 validation windows, batch size 1, bf16.
- Historical Arabidopsis/rice split and current primary-nuclear-chromosome allowlists.
- Frozen validation threshold grid, admission constraints, ranking and tie-break.
- Frozen Setaria FASTA and chromosome allowlist; annotation embargo remains active.
- Existing full and unchanged-input ablation definitions and all registered Setaria success/stop gates.
- No Setaria inference before an admissible development-validation tuple is selected and frozen.
- No use of old M25 checkpoints for Setaria prediction.

## Pre-submit verification

Run only checks that can block submission:

1. The regression test above passes.
2. A one-step GPU smoke includes an all-false structural-mask window and produces finite component and total losses, finite gradients, and finite updated parameters. It must not perform Setaria inference.
3. The resolved M25R config differs from M25 only in experiment/output identifiers and the approved audit fields; the scientific contract is unchanged.
4. The M25R output and Slurm log paths are new and empty.
5. Read the current remote `cluster_config.yaml` and live Slurm state immediately before submission; select the smallest eligible resource satisfying the unchanged run.

## M25R terminal decisions

### Implementation-invalid stop

Stop without scientific interpretation if any batch or epoch loss is non-finite, a checkpoint is missing/non-finite, the validation grid is incomplete without an explicit error, or the Setaria embargo is violated.

### Valid stop before Setaria

If finite training completes but no validation tuple simultaneously satisfies intergenic FPR <= 0.020, gene-count ratio 0.80-1.20, and structural-valid fraction >= 0.99, write `STOP_M25_BRANCH`, retain tuple-level rejection diagnostics, and do not access Setaria annotation.

### Setaria discovery gate

Only after a validation tuple is frozen, run the unchanged full and ablation predictions and evaluate the frozen nine primary chromosomes. `PASSED_DISCOVERY_GATE` requires all of:

- strand-aware exact CDS interval F1 >= 0.80;
- exact CDS-chain/coding-transcript F1 >= 0.55;
- exact coding-gene F1 >= 0.50;
- matched-gene strand accuracy >= 0.98;
- exact-matched CDS phase accuracy >= 0.90;
- structurally valid complete-transcript fraction >= 0.99;
- intergenic FPR <= 0.020;
- predicted gene-count ratio 0.80-1.20;
- full minus unchanged-input ablation interval F1 >= +0.10;
- full minus unchanged-input ablation chain F1 >= +0.10;
- full minus ablation FPR <= +0.005.

Any failed registered gate stops the branch without changing loss, decoder, thresholds, or sampling and without adapting to Setaria labels.

## Result package

At terminal state, report Slurm status, logs, finite-loss evidence, three checkpoints, selected checkpoint if any, resolved config, allowlist, complete validation-grid diagnostics, frozen decoder parameters, embargo marker, full/ablation GFF3, and all available metrics. Diagnose errors as representation, structural supervision, decoder/threshold, or species transfer only to the extent supported by raw emissions, decoded structures, full/ablation differences, and development-to-Setaria changes.

Helixer, Tiberius, ANNEVO, SegmentNT, extra seeds, alternative backbones, and decoder experiments remain blocked until M25R reaches a valid terminal scientific result.
