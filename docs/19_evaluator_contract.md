# Evaluator Contract / 评估器与可比性合约

> v4.1 central ledger. Migrated on 2026-06-14 from legacy `docs/11_evaluator_contract.md`; the legacy file is preserved for compatibility, but this file is now the preferred read-first contract for evaluator/metric/split/claim comparability.

## 0. Contract Status
- Status: active_migrated
- Owner skill: framework-upgrade migration; future changes by `/benchmark-roadmap`, `/reproduce-baselines`, `/code-review-gate`, `/result-log`
- Last verified: 2026-06-14 migration from existing M1 evaluator contract
- Applies to goal: cross-species ab initio protein-coding gene annotation, dual co-primary intergenic specificity + gene-body F1 contract in `ACTIVE_GOAL.json`

## 1. Migration Notes
- Source of truth before v4.1: `docs/11_evaluator_contract.md`.
- v4.1 source of truth from now on: `docs/19_evaluator_contract.md`.
- The migrated contract below intentionally preserves historical wording, thresholds, commands, and open decisions without reinterpretation.
- Any future metric/split/claim change should update this file and register evidence in `docs/15_evidence_register.md`.

## 2. Migrated Legacy Contract

# M1 Evaluator Contract

> Frozen during M1 setup after `BASE-TIBERIUS-MINISMOKE-EVALFIX`.
> This contract defines the project-level gene-body/FPR evaluator consumed by `scripts/validate_goal.py`.

## Scope

- Profiles: `smoke`, `screen`, `full`, `scale`.
- Claim policy: `smoke` and `screen` never claim SOTA. `full` and `scale` can only claim after frozen `sota_benchmark`, strict exceedance, comparability review, and human gate.
- Current implementation: `scripts/eval_gene_body_mask.py`.

## Gene-Body Mask

Input:
- reference GFF/GTF
- prediction GFF/GTF
- genome FASTA

Feature set used for spans:
- `CDS`
- `exon`
- `intron`
- `start_codon`
- `stop_codon`

Grouping:
- Use `transcript_id` when present.
- Fall back to `gene_id`.
- Fall back to a line-local synthetic ID only for malformed records.

For each `(seqid, strand, group_id)`, gene-body span is `min(start)..max(end)` over the selected features, using 0-based half-open intervals internally. Spans are then merged per `seqid` before computing base-mask lengths and overlaps. This intentionally treats protein-coding gene body as the transcript span implied by CDS/intron/start/stop evidence, not as exon-only bases.

## Primary Metrics

Definitions:
- `reference_gene_body_bases`: total merged reference gene-body mask length.
- `predicted_gene_body_bases`: total merged prediction gene-body mask length.
- `gene_body_overlap_bases`: base overlap between reference and prediction masks.
- `gene_body_precision = overlap / predicted_gene_body_bases`.
- `gene_body_recall = overlap / reference_gene_body_bases`.
- `gene_body_F1_unconstrained`: harmonic mean of precision and recall.
- `reference_intergenic_bases = genome_bases - reference_gene_body_bases`.
- `predicted_intergenic_false_positive_bases = predicted_gene_body_bases - overlap`.
- `intergenic_FPR = predicted_intergenic_false_positive_bases / reference_intergenic_bases`.

`constrained_gene_body_F1`:
- equals `gene_body_F1_unconstrained` when `intergenic_FPR` passes the active profile threshold;
- otherwise equals `0.0`.

Sensitivity fields are always emitted:
- `constrained_gene_body_F1_at_0.005`
- `constrained_gene_body_F1_at_0.01`
- `constrained_gene_body_F1_at_0.02`
- matching `intergenic_guardrail_pass_at_*` fields

## Profile Thresholds

`ACTIVE_GOAL.json` is authoritative.

Current thresholds:
- `smoke`: `intergenic_FPR <= 0.02`
- `screen`: `intergenic_FPR <= 0.02`
- `full`: `intergenic_FPR <= 0.01`
- `scale`: `intergenic_FPR <= 0.01`

Rationale:
- `smoke/screen` use the roadmap sensitivity upper bound while M1 evaluator/baseline anchors are being established.
- `full/scale` retain the stricter claim guardrail.

## Count Guardrails

`predicted_gene_count_ratio_vs_reference` is:

```text
unique predicted gene_id count / unique reference gene_id count
```

It is not divided by reference transcript count. Multi-isoform references can otherwise create false underprediction warnings for single-transcript gene callers such as Tiberius.

Transcript count is reported separately as:

```text
predicted_transcript_count_ratio_vs_reference
```

Current active guardrail:
- `predicted_gene_count_ratio_vs_reference <= 1.25`

No lower-bound guardrail is active yet. Add one only after M1 screen evidence shows underprediction is systematic rather than an annotation/transcript multiplicity artifact.

## Required Metrics JSON Fields

Required for `validate_goal.py`:
- `primary_metric`
- `constrained_gene_body_F1`
- `intergenic_FPR`
- `predicted_gene_count_ratio_vs_reference`
- `semantic_success`

Required for result interpretation:
- `gene_body_F1_unconstrained`
- `gene_body_precision`
- `gene_body_recall`
- `reference_gene_body_bases`
- `predicted_gene_body_bases`
- `gene_body_overlap_bases`
- `reference_intergenic_bases`
- `predicted_intergenic_false_positive_bases`
- `constrained_gene_body_F1_at_0.005`
- `constrained_gene_body_F1_at_0.01`
- `constrained_gene_body_F1_at_0.02`
- `intergenic_guardrail_pass_at_0.005`
- `intergenic_guardrail_pass_at_0.01`
- `intergenic_guardrail_pass_at_0.02`

## Current Limitations

- This contract is a screen/freeze evaluator for gene-body/FPR comparability. Full gene/transcript/locus-level SOTA claims still require SOTA-native metric cross-checks such as gffcompare where appropriate.
- The exact full-eval `sota_benchmark` remains pending M2 ANNEVO-compatible reproduction.

## Semantic-Gate Disposition (M1-AGGREGATION-GATE-AUDIT, 2026-06-10)

`scripts/validate_goal.py` distinguishes two kinds of `constrained_gene_body_F1 == 0.0`:

- **completed_poor** — a finite, valid-but-poor baseline whose constrained primary was intentionally hard-zeroed by the `intergenic_FPR` guardrail. Evidence required: the evaluator's own `semantic_success` flag is true, OR `gene_body_F1_unconstrained > 0` (configured under the goal contract's `semantic_success` block). Result: `status=not_yet`, `disposition=completed_poor`, exit 1 — iteration continues; NOT a stop. Not anchor- or claim-eligible.
- **failed_run** — genuinely degenerate output with no underlying signal, OR a real infrastructure failure (run-status OOM/TIMEOUT/missing/STALE, or non-finite/NaN metric). Result: `status=failed_run`, exit 3 — `/pursue` STOPS and notifies. This tripwire is unchanged.

The gate also emits `per_species_summary` and a heterogeneity warning when `per_species` is present.

## Aggregation Policy (required before any `screen_anchor` freeze)

Multi-species screen results MUST report all three views (`scripts/aggregate_gene_body_metrics.py` emits them):

1. **base-weighted** — base-count-weighted aggregate; the authoritative gate metric.
2. **macro** — unweighted per-species mean (`macro_*` fields); makes a single weak species visible when base-weighting is dominated by a large genome.
3. **per-species** — `per_species` rows with per-species constrained/unconstrained F1 + FPR.

Do NOT freeze `screen_anchor` from a mixed aggregate while any species fails the FPR guardrail. Example (PILOT-M1): base-weighted constrained F1 = 0.0 (hard-zeroed), macro = 0.4925, S. cerevisiae = 0.985, D. melanogaster = 0.0.

## Profile-Aware Guardrail Resolution (implemented in this audit)

`check()` in `validate_goal.py` now resolves, per rule:
- `threshold_by_profile`: `{profile: thr}` overrides base `threshold` for the active `--profile`.
- `profiles`: `[...]`; a rule with this key applies ONLY to the listed profiles and is SKIPPED (not failed) elsewhere (e.g. `nucleotide_gene_body_F1_drop_vs_anchor` is skipped for smoke). Skipped rules are marked `"skipped": true` in the output.

(Prior to this audit, profile resolution was documented but not actually implemented in the file; its regression test was red. Now green.)

## Gene-Body Span Mode (M1-SPAN-HARMONIZE-CDS, 2026-06-10)

`scripts/eval_gene_body_mask.py --span-mode {transcript,cds}` controls which features define the per-transcript min..max gene-body extent:

- `transcript` (default; backward-compatible): span over `CDS, exon, intron, start_codon, stop_codon` → full mRNA span INCLUDING UTR (via exon). Fair ONLY among tools that all emit UTR. **NOT apples-to-apples across tools**: Tiberius is a CDS-only caller (no UTR), so a transcript-span metric structurally penalizes it.
- `cds` (REQUIRED for cross-tool screen): span over `CDS, start_codon, stop_codon` only → coding-region span, the one layer EVERY gene caller emits (RefSeq, Tiberius, Helixer, ANNEVO). Canonical fair common denominator for protein-coding annotation; removes the tool-dependent UTR confound.

**Why CDS, not transcript**: comparing Helixer (GFF3, emits UTR) vs Tiberius (GTF, CDS-only) under `transcript` gave S. cerevisiae FPR 0.6544 (Helixer) vs 0.0164 (Tiberius) on the SAME reference — a pure schema artifact (Helixer predicted/reference gene-body ratio 1.224 on gene-dense yeast). Under `cds`, the same artifacts give S.cer FPR 0.0333 (Helixer) / 0.0186 (Tiberius), ratio ~0.99 — comparable. CDS-only is not "coarser"; it is the biologically canonical and tool-fair layer for a protein-coding annotation task. Per-feature-class (CDS/exon/intron/transcript) multi-resolution reporting is a future extension; CDS is the primary cross-tool layer.

**Cross-tool screen rule**: all screen/anchor baseline evaluations MUST use `--span-mode cds`. `transcript` mode is retained for within-tool diagnostics only.

## Open contract decision: screen FPR guardrail (surfaced by M1-SPAN-HARMONIZE-CDS)

Under the FAIR `cds` span, BOTH published SOTA tools' base-weighted aggregate FPR on the two pilot species narrowly EXCEEDS the current screen guardrail `0.02` (Tiberius 0.0225, Helixer 0.0228), hard-zeroing constrained F1 for both. Both PASS at `0.025` and `0.03`. This indicates the `0.02` screen guardrail is too strict to register any SOTA baseline on these gene-dense pilot species (small intergenic denominator inflates FPR). Pending user decision (a `/revise-goal`-style guardrail change): raise the screen FPR guardrail to ~0.025, and/or choose anchor species with more typical intergenic fractions, and/or report sensitivity-based anchor. `screen_anchor` stays BLOCKED until resolved. Unconstrained CDS F1 (informative): Tiberius base-weighted 0.8608 / macro 0.9150; Helixer base-weighted 0.9213 / macro 0.9494.

