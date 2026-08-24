# Pivot Decisions

> 由 /pivot append。每个 pivot 一段。

每个 entry 用 # Pivot Decision: <exp_id> 开头。模板见 /pivot SKILL.md。

---

# Pivot Decision: BASE-TIBERIUS-MINISMOKE

## Inputs consumed

- `/tri-review`: `docs/07_tri_review.md#tri-review-base-tiberius-minismoke`
- `/result-log`: `docs/06_results_log.md#result-base-tiberius-minismoke`
- Metrics: `outputs/BASE-TIBERIUS-MINISMOKE/metrics/metrics.json`
- Resource profile: screen-style mini-smoke; not claim eligible.

## Current evidence summary

The Tiberius bundled mini-smoke is a successful baseline infrastructure reproduction. It passed the official repo thresholds: CDS exact F1 `0.8594 >= 0.75` and transcript-chain exact F1 `0.3124 >= 0.28`. The project-level active primary became `0.0` only because the provisional evaluator hard-zeroed `gene_body_F1_unconstrained = 0.9196` after `intergenic_FPR = 0.0187` exceeded the draft guardrail `<= 0.01`.

The false-positive intergenic bases are approximately `495,567` bp over roughly `26,499,138` reference-intergenic bp. This is a strict threshold miss, not a sign of uncontrolled gene inflation: gene-body precision is `0.9654`, recall is `0.8779`, and predicted gene count ratio is `0.382x`.

## SOTA gap

| Metric | Current | SOTA / anchor | Gap (abs) | Gap (rel %) | Severity |
|---|---:|---:|---:|---:|---|
| `constrained_gene_body_F1` | `0.0` | `screen_anchor` pending | N/A | N/A | Not interpretable; anchor and evaluator not frozen |
| `gene_body_F1_unconstrained` | `0.9196` | `screen_anchor` pending | N/A | N/A | Smoke-only support metric |
| `intergenic_FPR` | `0.0187` | draft guardrail `0.01` | `+0.0087` | `+87.0%` over guardrail | Guardrail sensitivity / metric-contract issue |

## Sanity check

- [x] At least two independent CLI reviewers succeeded: 3/3.
- [x] Any reviewer raised comparability blocker: yes, all three.
- [ ] Metric implementation matches SOTA: no. The M1 frozen evaluator is not yet implemented, and the provisional gene-body mask is asymmetric.
- [x] Artifacts are present and parseable.
- [x] Official Tiberius mini-smoke thresholds pass.
- [x] No evidence of leakage or prediction explosion in this smoke result.
- [ ] Seed variance: not applicable / unknown for bundled single smoke.

Because metric implementation/comparability is not yet resolved, the decision must solve that before further anchor-setting runs.

## Tri-review summary

| Reviewer | Judgment | Next action proposed | Main concern | Confidence |
|---|---|---|---|---|
| A · Claude | `comparability-blocker` | Revise metric/guardrail, make gene-body mask symmetric, re-evaluate mini-smoke, then continue M1. | Provisional evaluator mixes reference CDS/intron-derived spans with prediction gene/transcript spans; `intergenic_FPR <= 0.01` hard zeroing is too strict for smoke. | High |
| B · Codex | `comparability-blocker` | Metric-contract revision plus sanity check; report unconstrained F1 and FPR thresholds separately; do not use hard-zeroed primary for unfrozen smoke evaluator. | Draft/provisional hard guardrail caused `validate_goal.py` to classify a successful baseline smoke as `failed_run`. | High |
| C · Antigravity | `comparability-blocker` | Relax/remove hard-zeroing for non-claim runs, implement frozen M1 evaluator, then continue baseline roadmap. | `0.01` is exceptionally strict for mini-smoke and evaluator span derivation is mismatched. | High |

Consensus: Tiberius mini-smoke succeeded; current `constrained_gene_body_F1=0.0` is an artificial metric-contract failure, not architecture failure.

Disagreement: whether to immediately relax screen guardrail to `0.02` or first disable hard-zeroing for non-claim profiles until M1 evaluator is frozen. This is a contract implementation detail, not a strategic disagreement.

Quorum / degraded review status: 3/3, no degradation.

## Reviewer-proposed directions

| # | From reviewer | Direction | major_axis | mechanism_delta | Orthogonal to others? | Into this round's cohort? |
|---:|---|---|---|---|---|---|
| 1 | A · Claude | Symmetric gene-body mask derivation and mini-smoke re-evaluation. | metric_contract | Same transcript-collapsing/span rule for reference and prediction. | Yes | Yes |
| 2 | B · Codex | Separate `gene_body_F1_unconstrained`, `intergenic_FPR`, and pass/fail at `0.005/0.01/0.02`; avoid hard-zero failed_run for unfrozen smoke evaluator. | metric_contract | Profile-aware guardrail instead of destructive score collapse. | Complementary to #1 | Yes |
| 3 | C · Antigravity | Revise goal/guardrail for smoke/screen and implement frozen M1 evaluator before setting anchors. | benchmark_contract | Frozen evaluator cross-checked against SOTA scripts plus profile-aware guardrail. | Yes | Yes |

## Is tuning justified?

Premature. This was a baseline reproduction and metric-contract failure, not a model-training result. Hyperparameter tuning or architecture replacement would be unsupported.

## Architecture hypothesis status

Unknown / not tested. The result supports that Tiberius inference and structured output pipeline are usable, but it does not evaluate our architecture paths.

## DECISION

- [ ] Continue current architecture as-is
- [ ] Tune current architecture
- [ ] Scale data / training
- [ ] Replace component
- [ ] Change backbone
- [ ] Change objective / loss
- [x] Comparability audit first
- [ ] Sanity check first
- [ ] Abandon this route
- [ ] Return to literature

## Why this decision

The run has no training or architecture failure to fix. It passed Tiberius' own mini-smoke metrics, and all three reviewers agree the blocker is the project evaluator/guardrail contract. Continuing M1 without fixing this would propagate an unfrozen, asymmetric, hard-zeroing metric into `screen_anchor`. Tuning is irrelevant; replacing architecture is premature; abandoning Tiberius reproduction is contradicted by the official smoke pass.

## Best next moves

| Priority | Move | Expected mechanism | Goes to which EXP / Track |
|---:|---|---|---|
| 1 | Freeze profile-aware metric contract for smoke/screen: always report unconstrained F1 and FPR separately; treat hard-zeroing as claim/full gate or report `pass@0.005/0.01/0.02` rather than semantic failed_run while evaluator is provisional. | Prevent guardrail threshold sensitivity from destroying otherwise valid smoke/screen evidence. | M1 evaluator prerequisite |
| 2 | Make provisional gene-body span derivation symmetric and re-run only evaluation on existing Tiberius mini-smoke artifacts. | Distinguish real extra spans from reference/prediction feature-schema mismatch. | `BASE-TIBERIUS-MINISMOKE-EVALFIX` or same B0 audit note |
| 3 | Implement frozen M1 evaluator cross-checking Tiberius, ANNEVO, and Helixer metric scripts before screen_anchor runs. | Ensure all screen baselines compare under one fair contract. | M1 `$reproduce-baselines` |

## Parallel cohort this round

No training cohort. The next work is a metric/comparability audit and evaluator sanity task before further M1 baseline screen runs.

| Slot | EXP ID (new) | Direction | major_axis | mechanism_delta | Track | Resource profile |
|---|---|---|---|---|---|---|
| primary | `BASE-TIBERIUS-MINISMOKE-EVALFIX` | Re-evaluate existing mini-smoke with symmetric masks and non-destructive guardrail reporting. | metric_contract | symmetric span derivation + threshold sensitivity report | baseline/M1 prerequisite | smoke |

## TODO update

- [x] Update `docs/08_pivot_decisions.md`.
- [ ] Update `docs/05_todo.md` to put metric-contract audit before further M1 baseline screen runs.
- [ ] Update `docs/04_experiment_iterations.md` with tri-review consensus and pivot decision.

---

# Pivot Decision: BASE-TIBERIUS-MINISMOKE-EVALFIX

## Inputs consumed

- `/tri-review`: `docs/07_tri_review.md#tri-review-base-tiberius-minismoke-evalfix`
- `/result-log`: `docs/06_results_log.md#result-base-tiberius-minismoke-evalfix`
- Metrics: `outputs/BASE-TIBERIUS-MINISMOKE-EVALFIX/metrics/metrics.json`
- Resource profile: smoke; not claim eligible.

## Current evidence summary

The profile-aware metric-contract fix worked as intended. Under smoke profile, `intergenic_FPR <= 0.02`, the existing Tiberius mini-smoke artifacts now yield `constrained_gene_body_F1 = 0.9196` with `run_ok=true`, `semantic_ok=true`, and guardrails passing. The same artifacts still fail stricter `0.005/0.01` sensitivity thresholds, so full/scale claim strictness remains intact.

## SOTA gap

Not evaluated. This is smoke-only and cannot set `screen_anchor` or compare to published SOTA.

## Sanity check

- [x] At least two independent CLI reviewers succeeded: 3/3.
- [x] Metric implementation change is deterministic and profile-scoped.
- [x] Full/scale strictness is preserved at `intergenic_FPR <= 0.01`.
- [x] Artifacts and metrics are parseable and finite.
- [x] No training/inference rerun was performed.
- [ ] `screen_anchor` established: no, still pending M1.

## Tri-review summary

| Reviewer | Judgment | Next action proposed | Main concern | Confidence |
|---|---|---|---|---|
| A · Claude | `continue-current-route` | Run unified screen baselines under revised evaluator, starting with Tiberius-like. | Apparent low count ratio (`0.328x`) should be checked. It was later resolved as transcript-count ratio, not gene-count ratio. | Medium |
| B · Codex | `continue-current-route` | Freeze M1 evaluator and run unified screen baselines, starting with Tiberius-like. | Do not use this smoke value as `screen_anchor` or SOTA evidence. | High |
| C · Antigravity | `continue-current-route` | Establish and freeze unified `screen_anchor` baselines. | None blocking. | High |

Consensus: continue M1. The threshold adjustment is acceptable for smoke/screen and does not weaken full/scale claim criteria.

Disagreement: none blocking. Reviewer A's low-count warning was advisory and is resolved by the evaluator contract: gene count ratio uses unique gene IDs; transcript multiplicity is reported separately.

Quorum / degraded review status: 3/3, no degradation.

## Reviewer-proposed directions

| # | From reviewer | Direction | major_axis | mechanism_delta | Orthogonal to others? | Into this round's cohort? |
|---:|---|---|---|---|---|---|
| 1 | A · Claude | Check apparent low predicted count ratio. | metric_contract | Distinguish gene-count ratio from transcript-count ratio. | Yes | Resolved in evaluator contract |
| 2 | B · Codex | Freeze M1 evaluator and start Tiberius-like unified screen baseline. | benchmark_contract | Move from bundled mini-smoke to common screen protocol. | Yes | Primary |
| 3 | C · Antigravity | Establish unified `screen_anchor` baselines. | benchmark_contract | Max over Tiberius-like / Helixer-like / ANNEVO-light under one protocol. | Complements #2 | Primary sequence |

## Is tuning justified?

No. This was a metric-contract sanity fix, not a model-training result.

## Architecture hypothesis status

Unknown / not tested.

## DECISION

- [x] Continue current route as-is: proceed to M1 evaluator freeze and unified baseline screen.
- [ ] Tune current architecture
- [ ] Scale data / training
- [ ] Replace component
- [ ] Change backbone
- [ ] Change objective / loss
- [ ] Comparability audit first
- [ ] Sanity check first
- [ ] Abandon this route
- [ ] Return to literature

## Why this decision

The blocker identified in the previous pivot is resolved for smoke/screen without weakening full/scale claim standards. The next milestone is not another metric debate; it is M1: freeze the evaluator under the revised contract and run comparable small-budget baselines to establish `screen_anchor`.

## Best next moves

| Priority | Move | Expected mechanism | Goes to which EXP / Track |
|---:|---|---|---|
| 1 | Freeze M1 evaluator contract and implementation. | Ensure all baselines use the same gene-body/FPR/sensitivity metrics. | M1 prerequisite |
| 2 | Run Tiberius-like screen baseline under unified protocol. | First component of `screen_anchor`. | `BASE-TIBERIUS-SCREEN-M1` |
| 3 | Run Helixer-like and ANNEVO-light/available screen baselines under the same protocol. | Complete `screen_anchor = max(...)`. | `BASE-HELIXER-SCREEN-M1`, `BASE-ANNEVO-SCREEN-M1` |

## Parallel cohort this round

No architecture cohort yet. Next work is M1 baseline/evaluator work. Use cluster-aware submission for any real screen inference/training.

## TODO update

- [x] Update `docs/08_pivot_decisions.md`.
- [x] Update `docs/05_todo.md` with M1 next-step and transcript multiplicity TODO.
- [ ] Update `docs/04_experiment_iterations.md` with tri-review and pivot links.

---

## Goal Revision: 2026-06-10 profile-aware intergenic FPR guardrail

### Trigger

User approved threshold adjustment after `BASE-TIBERIUS-MINISMOKE` tri-review/pivot found that `intergenic_FPR=0.0187` should not zero smoke/screen evidence while the evaluator is being frozen.

### Diff applied

- `ACTIVE_GOAL.json guardrails[intergenic_FPR].threshold`: kept base/full value `0.01`.
- Added `threshold_by_profile`: `smoke=0.02`, `screen=0.02`, `full=0.01`, `scale=0.01`.
- Added profile scoping for `nucleotide_gene_body_F1_drop_vs_anchor`: skipped for `smoke`.
- `scripts/validate_goal.py`: added support for `threshold_by_profile` and profile-scoped guardrails.
- `docs/03_benchmark_roadmap.md`, `docs/00_active_goal.md`, `CLAUDE.md`, and generated `AGENTS.md`: updated to reflect profile-aware rule.

### Comparability check

| Dimension | Verdict | Notes |
|---|---|---|
| Dataset version | N/A | This revision changes smoke/screen gating only, not a benchmark value. |
| Split scheme | N/A | No split or anchor value changed. |
| Metric implementation | PASS for gating semantics | Same metric is reported; sensitivity at `0.005/0.01/0.02` is explicit. |
| Preprocessing | N/A | No preprocessing change. |
| External weights | N/A | No model or weight change. |
| Test-time inference | N/A | No inference change. |
| Claim strictness | PASS | Full/scale claim candidates still use `intergenic_FPR <= 0.01` and cannot claim without frozen `sota_benchmark` + human gate. |

### Review evidence

- `BASE-TIBERIUS-MINISMOKE` tri-review: 3/3 judged the zeroed primary as metric-contract/comparability blocker, not model failure.
- `BASE-TIBERIUS-MINISMOKE-EVALFIX` tri-review: 3/3 judged the revision acceptable and recommended continuing M1.

### Effect

Smoke/screen runs can now surface valid progress instead of being classified as failed solely by `0.01` FPR sensitivity. Full/scale SOTA claim criteria remain strict.

---

# Pivot Decision: BASE-TIBERIUS-PILOT-M1

## Inputs consumed

- `/tri-review`: `docs/07_tri_review.md#tri-review-base-tiberius-pilot-m1`
- `/result-log`: `docs/06_results_log.md#result-base-tiberius-pilot-m1`
- Metrics: `outputs/BASE-TIBERIUS-PILOT-M1/metrics/metrics.json`
- Validation: `outputs/BASE-TIBERIUS-PILOT-M1/metrics/validate_goal.json`
- Resource profile: screen; not claim eligible.

## Current evidence summary

`BASE-TIBERIUS-PILOT-M1` is a completed Tiberius inference/evaluation pilot that was classified as `failed_run` by the project semantic gate. The failure is not an infrastructure failure: both species produced GTF predictions, metrics are parseable and finite, and Slurm exit `3:0` reflects the validator result being propagated from the sbatch script.

The biological/evaluation result is mixed. S. cerevisiae is strong (`constrained_gene_body_F1=0.9850`, `intergenic_FPR=0.0164`), while D. melanogaster fails the screen FPR guardrail (`intergenic_FPR=0.0295`) and has low recall (`0.5172`). The aggregate base-weighted `intergenic_FPR=0.0287` exceeds the screen threshold `0.02`, so aggregate `constrained_gene_body_F1` is hard-zeroed.

## SOTA gap

| Metric | Current | SOTA / anchor | Gap (abs) | Gap (rel %) | Severity |
|---|---:|---:|---:|---:|---|
| `constrained_gene_body_F1` | `0.0` | `screen_anchor` pending | N/A | N/A | Not anchor-eligible |
| `gene_body_F1_unconstrained` | `0.7087` | `screen_anchor` pending | N/A | N/A | Useful negative-control support metric |
| `intergenic_FPR` | `0.0287` | screen threshold `0.02` | `+0.0087` | `+43.5%` over threshold | Guardrail fail |

## Sanity check

- [x] At least two independent CLI reviewers succeeded: 2/3.
- [x] Degraded review status recorded: `DEGRADED_REVIEW`.
- [x] Metrics file exists and is parseable.
- [x] Values are finite; no NaN/Inf.
- [x] Prediction artifacts exist for both species.
- [x] Slurm failure is validator propagation, not OOM/timeout.
- [x] Any reviewer raised comparability / semantic-gate blocker: yes, Antigravity.
- [ ] Metric implementation matches published SOTA: unknown; this is the project M1 screen evaluator, not official paper metrics.
- [ ] `screen_anchor` established: no.

## Tri-review summary

| Reviewer | Judgment | Next action proposed | Main concern | Confidence |
|---|---|---|---|---|
| A · Claude | failed-after-retry; not counted | N/A | Raw text failed the required structured marker check, so it cannot contribute to quorum. | N/A |
| B · Codex | `continue-current-route` | Run Helixer two-species smoke/screen under the same evaluator; keep aggregation/per-species gate as TODO. | Tiberius pilot is a valid negative control but cannot update anchor; official comparability and preprocessing remain unresolved. | Medium |
| C · Antigravity | `comparability-blocker` | Fix evaluator/aggregation or validator semantics so valid poor finite baselines are completed-but-poor, then run Helixer. | Current `validate_goal.py` converts hard-zeroed poor baseline results into `failed_run`, which can break autonomous iteration. | High |

Consensus: no tuning, no anchor update, no claim, no abandon. The run is useful evidence that Tiberius current multi-clade release is heterogeneous under this FP-controlled screen objective.

Disagreement: whether Helixer should run immediately or after a lightweight aggregation/semantic-gate fix. Pivot resolves this by making the fix/audit the primary action and scheduling Helixer next.

Quorum / degraded review status: 2/3 `DEGRADED_REVIEW`; maximum confidence Medium.

## Reviewer-proposed directions

| # | From reviewer | Direction | major_axis | mechanism_delta | Orthogonal to others? | Into this round's cohort? |
|---:|---|---|---|---|---|---|
| 1 | B · Codex | Run Helixer smoke/screen on the same two species with the same evaluator. | benchmark_baseline | Add broad-lineage Helixer evidence to the baseline matrix. | Yes | Next after sanity fix |
| 2 | B · Codex | Add species-level gate or macro/per-species reporting to evaluator contract. | metric_contract | Prevent base-weighted aggregate from hiding single-species failure. | Complements #1 | Yes, primary |
| 3 | B · Codex | D. melanogaster sanity audit: feature mapping, contig filtering, strand/span conversion, softmasking comparability. | data_contract | Check whether the failed insect result is metric/data conversion rather than model behavior. | Yes | Optional |
| 4 | C · Antigravity | Change validator/result semantics so finite poor baselines are not classified as pipeline `failed_run`. | metric_contract | Distinguish semantic infrastructure failure from a valid baseline that fails the FPR-constrained objective. | Complements #2 | Yes, primary |
| 5 | C · Antigravity | Prepare ANNEVO smoke after Helixer to complete M1 baseline portfolio. | benchmark_baseline | Add ANNEVO current-release comparison while keeping paper/current metrics separate. | Yes | Later |

## Is tuning justified?

Premature / no. This is a frozen external baseline reproduction; the gap is to an unfrozen screen objective and the run failed a guardrail. Hyperparameter tuning would violate the workflow's architecture-first and comparability-first discipline.

## Architecture hypothesis status

Unknown for our architecture. For the baseline inventory, Tiberius current release is operational but weakened as a broad-eukaryote screen anchor candidate under this two-species pilot because performance is highly species-dependent.

## DECISION

- [ ] Continue current architecture as-is
- [ ] Tune current architecture
- [ ] Scale data / training
- [ ] Replace component
- [ ] Change backbone
- [ ] Change objective / loss
- [ ] Comparability audit first
- [x] Sanity check first -> lightweight M1 aggregation / semantic-gate audit
- [ ] Abandon this route
- [ ] Return to literature

## Why this decision

The key blocker is not the Tiberius command path; that path is now proven on two species. The blocker is how the harness should represent finite poor baseline performance under a hard FPR constraint. Treating every hard-zero constrained F1 as `failed_run` is appropriate for degenerate model outputs, but too coarse for baseline reproduction: a method can be reproducibly poor and still provide a valid screen-anchor comparison point. We need this distinction before `$pursue` or additional baseline runs repeatedly stop on valid negative controls.

This is not a reason to abandon Tiberius reproduction: it is already useful as a negative control and mechanism baseline. It is also not a reason to update `screen_anchor`: the aggregate fails the screen guardrail and the run is pilot-only.

## Best next moves

| Priority | Move | Expected mechanism | Goes to which EXP / Track |
|---:|---|---|---|
| 1 | Add/record M1 aggregation and semantic-gate policy: finite poor baseline with parseable metrics should be closed as completed-poor/not-anchor, not confused with OOM/missing-output failure. | Prevent false failed-run stops while preserving guardrail strictness. | `M1-AGGREGATION-GATE-AUDIT` |
| 2 | Add explicit per-species/macro reporting requirement before any `screen_anchor` update. | Make species heterogeneity visible instead of relying only on base-weighted aggregate. | evaluator contract / M1 |
| 3 | Run Helixer two-species smoke/screen with existing SIF and fungi/invertebrate weights. | Add a broad-lineage baseline comparator under the same project evaluator. | `BASE-HELIXER-SAC-DMEL-SMOKE-M1` |
| 4 | Keep ANNEVO setup next, but do not mix paper and current-release metrics. | Complete M1 baseline portfolio without corrupting the published SOTA anchor. | `BASE-ANNEVO-SMOKE-M1` |

## Parallel cohort this round

No model-training cohort. This is still B0/M1 baseline reproduction.

| Slot | EXP ID (new) | Direction | major_axis | mechanism_delta | Track | Resource profile |
|---|---|---|---|---|---|---|
| primary | `M1-AGGREGATION-GATE-AUDIT` | Lightweight evaluator/validator semantic audit for poor finite baselines. | metric_contract | completed-poor vs failed-run distinction; per-species/macro reporting | B0/M1 | local/no GPU |
| next | `BASE-HELIXER-SAC-DMEL-SMOKE-M1` | Helixer two-species smoke under same evaluator. | benchmark_baseline | broad-lineage baseline comparator | B0/M1 | smoke/screen |

## TODO update

- [x] Update `docs/08_pivot_decisions.md`.
- [x] Update `docs/04_experiment_iterations.md` with tri-review and pivot links.
- [x] Update `docs/05_todo.md`: close this failed run as reviewed/pivoted, add aggregation-gate audit, then Helixer smoke.

---

## Goal Revision: 2026-06-10 CDS-span harmonization + unconstrained screen anchor (user option-2)

### Trigger
`BASE-HELIXER-SAC-DMEL-SMOKE-M1` + `M1-SPAN-HARMONIZE-CDS` showed (a) the transcript-span gene-body metric is not apples-to-apples across tools (Helixer GFF3 incl UTR vs Tiberius CDS-only → S.cer FPR 0.654 vs 0.016), and (b) under the fair CDS span, both published SOTA tools' aggregate FPR (~0.0225) narrowly exceeds the 0.02 screen guardrail on gene-dense pilot species. User chose option-2: use unconstrained CDS F1 as the screen direction-selection anchor, with FPR advisory for screen.

### Diff applied (ACTIVE_GOAL.json)
- `screen_anchor.metric`: `constrained_gene_body_F1` → `gene_body_F1_unconstrained`; `value`: `0.0` → `0.9213` (PROVISIONAL = max(Tiberius 0.8608, Helixer 0.9213) base-weighted CDS-span unconstrained gene-body F1 on pilot species; macro Helixer 0.9494). Marked provisional pending ANNEVO-light + frozen anchor species.
- `guardrails[intergenic_FPR]`: added `profiles: [full, scale]` → HARD claim guardrail at ≤0.01 for full/scale ONLY; ADVISORY (reported in metrics + per_species_summary, not gating) for smoke/screen.
- `scripts/eval_gene_body_mask.py`: added `--span-mode {transcript,cds}`; cross-tool screen MUST use `cds` (docs/11 contract).

### Comparability check
| Dimension | Verdict | Notes |
|---|---|---|
| Dataset version | N/A | No dataset change; re-eval of existing predictions. |
| Split scheme | N/A | Unchanged (pilot species). |
| Metric implementation | PASS | CDS-only span makes ref/Tiberius/Helixer apples-to-apples; same construction for all. |
| Preprocessing | N/A | None. |
| External weights | N/A | None. |
| Test-time inference | N/A | None. |
| Claim strictness | PASS | full/scale claim guardrail UNCHANGED at intergenic_FPR ≤0.01; screen never claims. Only the screen direction-selection bar changed. |

### Status
- status stays `draft` (sota_benchmark still pending M2). screen_anchor now functional/provisional for Track A direction-selection.
- Not yet tri-reviewed; user-directed (option-2). Optional tri-review on the combined contract change (gate audit + CDS span + screen anchor) before any Track A promotion or screen_anchor finalization.

### Next
- Run ANNEVO-light under `--span-mode cds` (may raise the provisional anchor); re-derive screen_anchor on the frozen anchor species set (pilot yeast/fly are gene-dense outliers).

---

# Pivot Decision: M1-CONTRACT-REVIEW

## Inputs consumed
- `/tri-review`: `docs/07_tri_review.md#tri-review-m1-contract-review` (2/3 DEGRADED; Claude + Codex both `comparability-blocker`, High)
- Reproduction record: `docs/12_baseline_reproduction.md`
- Findings: `docs/10_findings.md` (screen_anchor semantics)
- Resource profile: baseline/M1, non-claim.

## Current evidence summary
The session's four M1 changes (completed_poor gate; CDS-span harmonization; screen_anchor→unconstrained CDS + FPR advisory; three pretrained-inference baselines) are individually sound fixes — CDS-span is the highest-value (removed a real UTR-vs-CDS unfairness, Helixer S.cer FPR 0.654→0.033). BUT both independent reviewers and the user converge on one blocker: `screen_anchor=0.9213` was the MAX of three PRETRAINED-inference baselines, i.e. a full-data ceiling, not the same-budget reference the two-tier anchor system requires. Gating Track A against it would unfairly kill every from-scratch small-sample candidate.

## DECISION
- [x] **Comparability blocker first** → build the true same-budget `screen_anchor` before any Track A architecture screening.

## Why this decision
The project's own anchor discipline (CLAUDE.md §4.5 #2) forbids judging screen candidates against a full-data SOTA value. The current anchor violated exactly that. This is not tuning, not abandon, not a model result — it is a benchmark-contract repair. Both reviewers High-confidence agree; the user identified it. Track A cannot start with a dishonest ruler.

## User-selected scope (the open disagreement, resolved)
- Reviewer A: 1 standard backbone + FLOOR. Reviewer B: ≥2 family refs. **User chose Reviewer B scope: TWO family references** (Tiberius-like + Helixer-like; ANNEVO-light if tractable else deferred), plus a cheap FLOOR baseline.

## Best next moves (the M1-SAMEBUDGET-SCREEN-ANCHOR sub-stage)
| Priority | Move | Mechanism |
|---:|---|---|
| 1 | Freeze ONE unified small-sample screen protocol (train species/fraction, val/test split, window/step, epochs, patience, seeds≥3, CDS-span metric, preprocessing) shared by anchor refs AND future Track A candidates. | One harness, only architecture varies → fair. |
| 2 | Implement + train (random-init) `Tiberius-like` (CNN+biLSTM + structure-aware labels / CRF-or-semiCRF-style constraint) and `Helixer-like` (per-base gene-body segmentation backbone) light references under that protocol. ANNEVO-light if its training path is tractable, else record `ANNEVO-light deferred` (do NOT mix pretrained-ANNEVO numbers into the screen gate). | Same-budget references. |
| 3 | Add a cheap FLOOR baseline (ORF/GC/majority-class) → bracket floor < screen_anchor < pretrained_ceiling(0.9213). | Guards against a too-weak anchor. |
| 4 | `screen_anchor :=` seed-averaged max of the same-budget refs (CDS unconstrained gene-body F1). Report FPR / pred-ref ratio / macro as advisory diagnostics. | Honest screen bar. |
| 5 | Tighten `completed_poor` exemption: the evaluator's `semantic_success` flag is a constant (always true) so it is NOT real evidence; require `gene_body_F1_unconstrained` ≥ a non-trivial floor (~0.05) AND a sane predicted/reference count ratio. Add a regression test. | Keep the failed_run tripwire strict. |
| 6 | Keep a soft-warn FPR threshold (advisory, e.g. > 2× reference FPR) + recalibrate guardrails on real screen species (not gene-dense yeast/fly). | Catch over-prediction. |
| 7 | `/revise-goal` (human gate) → set `status: draft → active` once a real same-budget `screen_anchor` exists. | Enables success judging. |

## Track A / Track B
- Both PAUSED until the same-budget `screen_anchor` exists. status stays `draft`.

## Parallel cohort this round
- No model-training cohort yet. Next work is the M1-SAMEBUDGET-SCREEN-ANCHOR harness + reference training (B0/M1).

## Contract changes applied at pivot
- `ACTIVE_GOAL.json`: 0.9213 moved to a non-gating `pretrained_ceiling` field; `screen_anchor.value` reset to 0.0 (pending) with the same-budget definition.
- Confidence: Medium (2/3 DEGRADED_REVIEW; both reviewers independently High on the central conclusion).

## TODO
- [x] Update `docs/07`, `docs/08`, `docs/12`, `ACTIVE_GOAL.json`.
- [ ] Build M1-SAMEBUDGET-SCREEN-ANCHOR (next phase; user to direct start — implementing 2 light architectures is a substantial `/implement` effort).
- [ ] Tighten completed_poor exemption + add test (do with the anchor harness).

---

## Mid-iteration note 2026-06-10: m1-screen-anchor-bracket

### Mini-retrospective on m1-screen-anchor-bracket
- Relevance: directly-attacks-current-gap (this IS the same-budget screen_anchor the pivot required).
- Does it change our hypothesis? no — it operationalizes the agreed correction; quantifies same-budget 0.56 vs pretrained 0.92.
- Conflicts with an abandoned route (docs/09)? no.
- Recommendation: continue-as-is — Track A is now unblocked against an honest bar (beat 0.5579).
- Urgency: low — proceed to Track A when ready; re-derive on frozen typical-intergenic species before heavy reliance.
- advisory only; not a pivot, no docs/03 change.

---

## Goal Revision 2026-06-10: M1-GOAL-REVISION (R1-R5)

### Trigger
User /revise-goal: "adjustments from this session not yet landed in the contract." Audit of ACTIVE_GOAL.json vs session decisions found 4 stale/missing items; tri-review added a 5th (required).

### tri-review
- docs/07#tri-review-m1-goal-revision — 2/3 DEGRADED (Claude+Codex both continue-current-route; antigravity absent). Unanimous: revision sound, claim strictness preserved, two-tier design intact. Both required the R5 promotion-discipline addition; Claude required R3 to mark the inert guardrail disabled.

### Comparability check
- Dataset/split/metric/preprocessing/weights/test-time: all UNCHANGED. Claim guardrails (intergenic_FPR<=0.01, predicted_gene_count_ratio<=1.25) remain HARD at full/scale. sota_benchmark unchanged (0.0, M2). status stays draft. Claim strictness PRESERVED — only the never-claim screen profile is relaxed.

### Diff applied (ACTIVE_GOAL.json)
- R1 `_comment`: clarified screen_anchor established (0.5576); draft now blocks only full-claim (sota_benchmark, M2).
- R2 guardrails[predicted_gene_count_ratio_vs_reference]: + profiles=[full,scale] (advisory/reported for smoke/screen; hard claim guardrail). Rationale: same-budget anchor refs themselves fail 1.25 (tiberius 1.8-4.1) -> hard-gating screen is unfair; mirrors intergenic_FPR.
- R3 guardrails[nucleotide_gene_body_F1_drop_vs_anchor]: DISABLED via sentinel profile + _status=inert_pending_evaluator (evaluator emits placeholder 0.0 -> false assurance). Re-enable after evaluator computes true drop vs screen_anchor.
- R4 + top-level `_metrics_note`: documents the two-metric design (screen=unconstrained CDS F1; claim=constrained).
- R5 + `track_a_promotion` (tri-review REQUIRED): seed-wise gene-count reporting; beats-anchor != promotable; severe fragmentation blocks Track B promotion unless plan includes a structural-decoder/coherence fix.

### Effect
- screen/Track-A candidates judged on base-F1 vs screen_anchor (fair, same terms as the anchor refs); fragmentation reported + governed by promotion discipline, not a premature hard screen gate. Verified: SCREENREF runs -> not_yet/completed_poor with all guardrails advisory-skipped at screen. 7/7 tests pass. status remains draft.

---

# Pivot Decision: TA-DECODER-M3
## Inputs: docs/07#tri-review-ta-decoder-m3 (2/3 DEGRADED, 1-1 split), docs/06 result. Profile screen, non-claim.
## DECISION (autonomous part — common ground)
- [x] Record CONSTR as a VALIDATED same-budget screen result: seed-mean 0.5791 > gate 0.5676 (beats screen_anchor 0.5576) AND gene_count_ratio 2.74->1.12 (now < 1.25 claim guardrail). M3 primary_progress_gate MET. But it is constrained-Viterbi POST-PROCESSING, NOT a learned structured decoder.
- [x] Highest-SOTA-ROI next direction (BOTH reviewers agree): vectorize CRF/semi-CRF to fairly test the CORE learned-structure bet (the project's primary architecture hypothesis), which tractability — not science — left untested. Add CONSTR per-seed paired delta vs softmax.
## DEFERRED TO USER (decision-autonomy exception: 1-1 reviewer tie/no-leader + Track-B resource + new long sub-iteration)
The contested choice between:
  (A) [Codex] promote CONSTR to Track B NOW (scale data/seeds/CI) + parallel vectorized learned-decoder batch; vs
  (B) [Claude] HOLD CONSTR Track-B; vectorize CRF/semi-CRF FIRST for a fair learned-structure test, then decide.
Both keep the CONSTR screen win on record and both schedule the vectorized learned-decoder batch; they differ only on whether CONSTR scales now or after the fair learned-decoder test.
## Anti-tuning / R5
CONSTR meets R5 (beat anchor + fixed fragmentation). No tuning. The learned-structure axis (CRF/semi-CRF) is the architecture direction to actually validate next.
## docs updated: 04, 05, 06, 07, 08, exp-log.


## Pivot UPDATE: TA-DECODER-M3 — tie broken (Reviewer C = B)
3/3 tri-review, majority **B** (Codex + Antigravity): promote CONSTR to Track B NOW (scale data/seeds/CI; labeled post-processing structured-inference candidate, NOT learned-decoder success) + in PARALLEL run a vectorized CRF/semi-CRF learned-decoder batch. Claude's A (vectorize-first) is the minority but its core caution is RETAINED in the plan: the vectorized learned-decoder batch is mandatory and CONSTR must NOT be framed as the learned-structure success.
RESOLVED next direction = B. LAUNCH still pending user go-ahead (decision-autonomy exception: Track B scale-up + new vectorization sub-iteration = new long sub-iteration / possible >24h compute).


---

# Pivot Decision: TA-DECODER-VEC-M3
## Inputs: docs/07#tri-review-ta-decoder-vec-m3 (2/3 DEGRADED, 2-0 consensus), docs/06 result. Screen, non-claim.
## DECISION (consensus, autonomous): promote-learned-decoder-to-Track-B
CRF-vec (vectorized learned linear-chain CRF) is the Track A winner: seed-mean base-w gene_body_F1_unconstrained (CDS) 0.6186 > gate 0.5676 > anchor 0.5576, AND > CONSTR 0.5791, with best coherence (ratio 0.88). LEARNED structure > post-processing > per-base baseline — the project's core architecture bet is validated at same-budget screen. semi-CRF deferred (vectorize segment DP later).
## Track B plan (job #1 = seed variance, per both reviewers)
≥5-8 seeds + mean±CI + paired test vs CONSTR (current spread 0.081 > edge 0.040; s2 loses to CONSTR); scale data/epochs to test if the CRF advantage GROWS (scalability bet) or is a small-sample artifact. Keep CONSTR as the in-Track-B baseline. Anti-tuning: this is architecture (decoder), not tuning.
## ⏸ LAUNCH pending user go-ahead (decision-autonomy exception: Track B scale-up = new long sub-iteration / likely >24h compute).
## Anti-tuning/R5: CRF-vec meets R5 (beat anchor + CONSTR + coherent). docs updated: 04,05,06,07,08,exp-log.

## Goal Revision 2026-06-11 — REVISE-INTERGENIC-PRIMARY-M1 (FOUNDATIONAL ruler change)
### Trigger (user, foundational)
"intergenic 稳定性升为主指标 + intergenic 用 full-transcript（含 UTR）补集"。User insight: the old numbers (anchor 0.5576, CONSTR 0.5791, CRF-vec 0.6186) were all measured by the CDS-only-complement ruler that wrongly counts UTR as intergenic; change the ruler BEFORE changing architecture (SegmentNT) or it is wasted compute ("先换架构再改尺子 = 白跑").
### Comparability (revise-goal Step 2) — this is a RE-EVALUATION of EXISTING predictions, no retrain/GPU
Deterministic re-scoring of frozen prediction GFFs on the SAME held-out test subsets. Strong reproducibility. FLOOR re-done on the identical test subset (was whole-genome) to fix the comparability blocker all 3 reviewers raised.
### Code changes landed (eval ruler)
- scripts/eval_gene_body_mask.py: intergenic = genome - FULL-transcript span (incl UTR via exon), DECOUPLED from gene-body-F1 span_mode. Adds intergenic_specificity (=1-FPR, the new primary), reference_full_transcript_bases, intergenic_FPR_cds_complement_diag. eval tests 2 pass.
- scripts/aggregate_gene_body_metrics.py: emits intergenic_specificity (base-weighted) + macro_intergenic_specificity + per-species; primary_metric -> intergenic_specificity.
### Recompute (NEW ruler, 3 seeds, base-weighted; identical test subsets) — RANKING FLIPS
FLOOR(ORF) spec 0.8805 / F1 0.3735 (BLOCKED by F1 floor) | tiberius_like(ANCHOR) 0.8710 / 0.5576 | CONSTR 0.8369 / 0.5791 | helixer_like 0.7954 / 0.5579 (frag 99.5) | CRF-vec 0.7138 / 0.6186. CRF-vec went from OLD-ruler WINNER (0.6186 F1) to WORST on the new primary (0.7138 spec, highest FPR 0.2862) — structured decoders raise recall by spilling into intergenic DNA.
### tri-review 3/3 (docs/07#REVISE-INTERGENIC-PRIMARY-M1): approve-with-modifications. UTR redefinition CORRECT (code-verified). Required mods: anti-gaming floor (0.50 too loose), macro as a GATE, anchor provisional (low-UTR outliers), reframe CRF-vec as FP-aware-objective need not "decoders out".
### DECISION (user human gate)
User delegated the primary-metric structure to the agent ("选最合理、贴合发表的；CLAUDE §0 当初没考虑这么细，及时修正"), and chose anchor=PROVISIONAL. Agent chose **DUAL CO-PRIMARY (Pareto)**: AXIS-1 headline/ranking = intergenic_specificity (honors user's elevation); AXIS-2 = gene-level F1 (keeps the SOTA-comparable, publishable claim + §0 north star + protects the structured-decoder bet). Promotable iff specificity STRICTLY > anchor AND gene_body_F1 >= floor (0.5276 screen / 0.5576 promotion) AND macro_specificity >= 0.7978. Both M1 gaming modes blocked (FLOOR by F1, CRF-vec by specificity).
### ACTIVE_GOAL.json diff applied
- primary_metric: constrained_gene_body_F1 -> intergenic_specificity.
- success_criteria: now TWO (intergenic_specificity AND constrained_gene_body_F1), both pending M2 (Pareto claim).
- screen_anchor: gene_body_F1_unconstrained 0.5576 -> intergenic_specificity 0.8710 (macro 0.8278), status=provisional, anchor_gene_body_F1=0.5576 retained.
- guardrails +gene_body_F1_unconstrained>=0.5276 [screen,full,scale] +macro_intergenic_specificity>=0.7978 [screen,full,scale]; intergenic_FPR now full-transcript-based; nucleotide_drop still inert (superseded by absolute floor).
- pretrained_ceiling: metric -> intergenic_specificity, value null pending recompute (0.9213 was old-ruler gene_body_F1).
- track_a_promotion: Pareto rewrite. architecture_hypothesis: reframe decoders + add foundation+FP-aware path. status stays draft (sota_benchmark pending M2).
### CONSEQUENCE — CRF-vec Track-B promotion (TA-DECODER-VEC-M3) INVALIDATED/PAUSED
The prior pivot's promote-CRF-vec-to-Track-B is invalidated under the corrected ruler (CRF-vec is now the worst candidate on the primary). It had NOT launched (was ⏸ pending user go-ahead) — caught BEFORE any compute was spent, exactly as the user warned. CRF-vec retained as an ablation (recall↑/specificity↓), NOT scaled. Next direction: foundation-probe (SegmentNT/GENERanno) → cut intergenic FP without losing recall → then semi-CRF + FP-aware objective.


## Pivot Decision: FP-SEGMENTNT-PROBE-M1 (2026-06-11)
### Inputs: docs/07#FP-SEGMENTNT-PROBE-M1 (3/3 quorum, all iterate/change-objective), docs/06 result. Track A screen, NON-CLAIM, validate=not_yet.
### Result recap: frozen SegmentNT features -> anchor-matched conv+biLSTM head. AXIS-2 gene_body_F1 0.6888 >> anchor 0.5576 (+0.13, PASS); AXIS-1 intergenic_specificity 0.8416 < anchor 0.8710, macro 0.7543 < gate 0.7978 (FAIL, yeast/fungus over-prediction). Not Pareto-dominant. High spec seed variance (s1=0.897 > anchor).
### DECISION (3/3 consensus, autonomous): ITERATE-PROBE via change-objective-or-loss. NOT abandon (features add real recall), NOT promote (fails AXIS-1), NOT scale.
### Next round (NEW goal; parallel ≤3 orthogonal, all reuse outputs/FP-SEGMENTNT-FEATCACHE, same same-budget screen protocol, NEW ruler):
  - Direction A (MAIN, change-objective-or-loss): FP-aware / specificity-targeted objective on the frozen-feature probe — asymmetric intergenic-FP cost / precision-biased (focal or boundary precision re-balance) added to class-weighted CE, penalizing predicted-genic bases that fall in true intergenic (full-transcript complement). Convert the +0.16 F1 margin into specificity; target the yeast over-prediction directly. mechanism_delta = loss_design.
  - Direction B (ORTHOGONAL, data_view/training_signal): FUSE raw-DNA one-hot ⊕ frozen SegmentNT logits into the same conv+biLSTM head (gated fusion) — anchor(raw-DNA) already has spec 0.871, foundation has recall; combine to get both. mechanism_delta = input fusion.
  - Optional per-clade calibration (threshold/prior/temperature) layered on either; report per-species + gene_count_ratio as HARD diagnostic.
  - STATISTICS: run ≥5 seeds + mean±CI + paired test vs anchor on intergenic_specificity (current spec spread 0.084 >> edge 0.040; AXIS-1 verdict is variance-fragile — one seed already beats the anchor).
### DEFER (record, do not pursue this round):
  - semi-CRF / structured decoder: M1 evidence (CRF-vec spec 0.7138 << anchor) shows it HURTS specificity; on spillover-prone emissions it only makes coherent wrong genes. Revisit AFTER FP-aware emissions control the spillover.
  - unfreeze / fine-tune SegmentNT: expensive, breaks the clean frozen ablation, small-cross-species overfit risk -> Track B only, after a frozen route wins.
  - GENERanno (different foundation model): parallel cheap literature/probe branch (its cross-clade base/CDS signal may transfer better to fungi), NOT the main path. = replace-component.
### Anti-tuning: gap to AXIS-1 anchor (0.8710-0.8416=0.029) < tuning_gap_threshold 0.05, BUT the chosen move is loss_design (architecture axis), not lr/batch/dropout tuning -> compliant. Direction B is a structural input change.
### Pre-claim HARD guard (Claude): before ANY full/scale claim, verify the test species/clade are NOT in SegmentNT's pretraining corpus (it saw vertebrate genomes) — else the foundation-feature advantage is leakage-contaminated. Pilot yeast+fly are NOT in a human/vertebrate pretraining test set, so the screen comparison is clean; the guard is for future vertebrate held-out evals.
### Docs updated: 04,05,06,07,08,10,00.


## Pivot Decision: TA-FOUNDATION-DECODER-M4 (2026-06-11)
### Inputs: docs/07#TA-FOUNDATION-DECODER-M4 (3/3 quorum), docs/06 result. Track A screen, NON-CLAIM.
### Result: FPLOSS (FP-aware specificity-targeted loss on frozen SegmentNT features) PARETO-beats the same-budget anchor on the dual co-primary: intergenic_specificity 0.9303 > 0.8710 (all 5 seeds > anchor mean), gene_body_F1 0.6157 > 0.5576, macro 0.8431 > 0.7978. FIRST candidate to strictly exceed the anchor on the new ruler -> MAIN architecture bet (foundation features + FP-aware objective) VALIDATED at screen. FUSION 0.8615 (just below anchor) no; CRF 0.8298 (high variance ±0.119, one seed 0.59) no but best gene_count coherence 0.90.
### DECISION (autonomous, reconciling 3/3 split 2 promote / 1 iterate): ITERATE one cheap screen round to FIX FRAGMENTATION before Track-B promotion. FPLOSS is the validated lead; promote-as-is WITHHELD.
### Rationale: ALL 3 reviewers flagged FPLOSS gene_count_ratio 2.25 = fragmentation; the full/scale HARD guardrail predicted_gene_count_ratio<=1.25 would BLOCK it -> promoting a 2.25-fragmented candidate sends up something that fails at the promotion ruler. base-weighted spec + base gbF1 are blind to fragmentation; only gene_count_ratio catches it. AGENT CORRECTION (reviewers missed): the CRF candidate ALREADY = FP-loss + CRF decoder (their recommended 'synthesis'), and it TRADED specificity for coherence (spec 0.830, coherent 0.90, +variance). So the learned CRF is not a free coherence fix.
### Next round (NEW goal; cheap Track A screen, reuse FEATCACHE, same protocol):
  - Direction A (MAIN): FPLOSS (FP-aware loss winner) + CHEAP constrained-decode POST-PROCESSING (reuse src/screen_anchor/decoders.py constrained_decode: merge small intergenic gaps / drop tiny CDS) to fix the 2.25 fragmentation WITHOUT the learned-CRF's specificity cost + variance. Target: keep spec ~0.93, pull gene_count_ratio 2.25 -> <=1.25.
  - Direction B (stabilize the structured decoder): diagnose the CRF collapsed seed (0.593); regularize / warm-start emissions from FPLOSS / stronger FP-aware weight inside CRF, to recover spec while keeping coherence 0.90.
  - Direction C (validity, REQUIRED): rerun the ANCHOR to 5 seeds (currently 3, one collapse 0.773) -> valid paired test vs FPLOSS on intergenic_specificity. Report transcript-span precision/recall + gene-length/exon-count distributions (fragmentation diagnostics, per B).
  - >=5 seeds + CI. Confirm FP-loss lambda was NOT tuned on test (it was hardcoded 1.0 — confirmed; document).
### THEN promote the coherence-fixed winner to Track B (scale-up = new long sub-iteration -> USER GO-AHEAD required per decision-autonomy >24h-compute exception). Track-B job#1: scale data/epochs/seeds + CI + the richer Tiberius-style multi-class output (CDS/intron/intergenic/phase/splice) — multi-class gives the CRF transitions real biological meaning (reviewer A: likely the ceiling-approach step).
### DEFER: chasing ceiling 0.9917 (different pretrained+full-data regime, not screen-comparable); unfreeze/fine-tune SegmentNT (Track B); GENERanno (parallel probe). Pre-claim guard: verify test clade not in SegmentNT pretraining.
### FUSION: dropped as a standalone mechanism (fails spec+macro, worst fragmentation 3.40); its data_view idea may return only fused with FP-loss later.
### Anti-tuning: gap FPLOSS->anchor is FPLOSS ABOVE anchor (it won); the iterate is a structural/coherence fix (decoder/post-proc axis), not lr/batch tuning -> compliant.
### Docs updated: 04,05,06,07,08,10,00.


## Pivot Decision: TA-COHERENCE-FIX-M5 (2026-06-11)
### Inputs: docs/07#TA-COHERENCE-FIX-M5 (3/3 quorum), docs/06 result. Track A screen, NON-CLAIM.
### Result: FP-FRAGFIX-CONSTR (FPLOSS + deterministic constrained post-proc) PARETO-beats the 5-seed anchor, PAIRED-SIGNIFICANT: intergenic_specificity 0.9272 vs anchor 0.8436 (paired +0.0836±0.037, all 5 seeds positive, t~5.0 p<0.01), gbF1 0.6581 > anchor 0.5768, macro 0.8555 > gate, gene_count 2.25 -> 1.28 (de-fragmented, 0.03 above full/scale guardrail 1.25). Clean net constrained contribution (vs M4 FPLOSS): spec kept (0.930->0.927), F1 up (0.616->0.658), fragmentation fixed.
### DECISION (autonomous, 3/3 consensus): PROMOTE-READY pending ONE cheap constrained-param sweep. CONSTR is the validated winner; the only blemish (gene_count 1.28 > 1.25) is cleared by a deterministic, no-retrain param sweep BEFORE spending Track-B compute.
### Next (NEW goal -> ③; the cheap sweep is its first step):
  - STEP 0 (cheap, deterministic, gate): save raw pre-constrained per-base predictions (small code add) OR re-run the 5 CONSTR seeds; sweep constrained_decode (max_fill_gap/min_cds_len) on TRAIN/VAL, require gene_count_ratio<=1.25 AND intergenic_specificity>=anchor AND gbF1>=floor. (spec & gene_count are COUPLED — verify, since merging more into gene-body can raise intergenic FP.) Confirm params chosen on train/val NOT test (current defaults 30/20 are non-test — OK).
  - THEN promote to Track B (= USER GO-AHEAD, new long sub-iteration, >24h-compute exception): scale data/epochs/seeds + CI; add richer Tiberius-style multi-class output (CDS/intron/intergenic/phase/splice) — gives the structure real meaning; UNFREEZE/fine-tune SegmentNT as a SEPARATE staged axis (NOT mixed into the first scale run — attribution). Keep CONSTR (deterministic post-proc) as the coherence layer.
### Optional /revise-goal (human-gated, does NOT change promotion): update screen_anchor 0.8710(3-seed) -> 0.8436(5-seed, more representative + higher variance); keep both recorded. CONSTR beats both via the paired test.
### Recorded: constrained params = defaults (min_cds_len 30 / max_fill_gap 20), NOT test-tuned (no leakage). FP-loss lambda 1.0 hardcoded, not test-tuned. Clean attribution baseline = M4 FPLOSS (anchor comparison mixes SegmentNT pretraining dividend). CRFSTAB deferred (CONSTR superseded the need this round). FUSION dropped (M4).
### Docs updated: 04,05,06,07,08,10,00.


## Pivot Decision: TA-FRAGFIX-SWEEP-M6 (2026-06-11)
### Inputs: docs/07#TA-FRAGFIX-SWEEP-M6 (3/3), docs/06. Track A screen, NON-CLAIM, STEP-0 promote-gate.
### Result: VAL-chosen constrained params (no test leakage) clear FP-FRAGFIX-CONSTR's gene_count. Initial max-spec-s.t.<=1.25 rule picked mfg=20/mcl=90 -> under-prediction (0.55-0.70 seeds); 2/3 tri-review flagged the one-sided rule. ZERO-GPU re-selection with a TWO-SIDED band [1.0,1.25] -> ADOPTED mfg=20/mcl=60: TEST spec 0.9218 (all 5 > anchor 0.8710 & 0.8436), gbF1 0.6439 (>floor, recovered), macro 0.8331 (>gate), gene_count 1.037 (~1.0, <=1.25). ALL 4 GATES PASS.
### DECISION (autonomous, 3/3 consensus): FP-FRAGFIX-CONSTR (frozen SegmentNT feats + FP-aware loss + constrained post-proc mfg=20/mcl=60) is PROMOTE-READY. STEP-0 gate cleared with the reviewer-corrected two-sided selection rule.
### Next = ③ Track-B promotion (NEW goal, USER GO-AHEAD required, new long sub-iteration / >24h):
  - scale data/epochs/seeds (>=8 seed) + report intergenic_specificity & gene_count_ratio with CI + TWO-SIDED gene_count band (not one-sided).
  - richer Tiberius-style multi-class output (CDS/intron/intergenic/phase/splice) — gives the structure real meaning.
  - staged UNFREEZE / fine-tune SegmentNT as a SEPARATE axis (not mixed with the first scale run — attribution).
  - keep the deterministic constrained post-proc (mfg=20/mcl=60) as the coherence layer; re-select the band per Track-B protocol on val.
### Optional /revise-goal (human-gated, doesn't change conclusion): screen_anchor 0.8710(3-seed)->0.8436(5-seed).
### Recorded: no leakage (VAL-chosen, TEST applied once); FP-loss lambda 1.0 non-test; constrained deterministic; gene_count seed variance (0.6-1.5) to manage in Track B with CI + band. Docs: 04,05,06,07,08,10,00.

---

## Retrospective 2026-06-11

### Scope
- Iterations covered: ITER-B0-008 .. ITER-FP-004 (the architecture-relevant arc since the same-budget anchor was set; B0-001..007 were metric-contract/baseline-reproduction, summarized but not re-litigated).
- Trigger: every-5 (≥7 completed iterations since last retrospective; first retrospective of the project).
- Focus: empty (full retrospective).

### Inline timeline summary (evidence base)
- docs/04: B0-001..007 = metric contract + baseline trio (Tiberius/Helixer/ANNEVO); B0-008 = same-budget screen_anchor 0.5579 (OLD ruler); B0-009 CONSTR / B0-010 CRF-vec (OLD ruler F1 winners); RULER CHANGE (REVISE-INTERGENIC-PRIMARY-M1) flipped primary→intergenic_specificity and INVALIDATED CRF-vec before launch; FP-001..004 = foundation-probe → FP-aware loss → coherence-fix → promote-gate.
- docs/08: every pivot was tri-reviewed; FP-series consensus = foundation-features + FP-aware objective is the validated lead.
- docs/09: EMPTY — no formally abandoned route (CRF-vec retained as ablation, not abandoned).

### Are we doing marginal tuning?

| Verdict | Evidence |
|---|---|
| **partially (mostly NO)** | FP-001 input=frozen SegmentNT features (data_view) → FP-002 FP-aware loss (loss_design) → FP-003 + deterministic constrained post-proc (decoder/post-proc axis): each adds a genuine mechanism, not lr/batch/dropout. BUT FP-004 (TA-FRAGFIX-SWEEP-M6) was a constrained-decode PARAM sweep (max_fill_gap × min_cds_len) — deterministic param SELECTION, not a new axis. It is defensible as a ONE-TIME promote-gate to clear the gene_count≤1.25 guardrail (VAL-chosen, no leakage), but it traded headline spec down (0.9272→0.9218) for coherence. ⚠️ The NEXT post-proc param tweak would cross into marginal tuning — Track B must move on a structural axis (multi-class output / staged unfreeze), not more mfg/mcl fiddling. |

### Gap trajectory (NEW ruler: intergenic_specificity; anchor 0.8710 (3-seed) / 0.8436 (5-seed); ceiling 0.9917)

| ITER | Track | Candidate | spec | gene_body_F1 | gene_count | Δspec vs anchor |
|---|---|---|---:|---:|---:|---:|
| FP-001 | A-screen | SegmentNT-probe | 0.8416 | 0.6888 | — | −0.029 (BELOW) |
| FP-002 | A-screen | FP-SEGNT-FPLOSS | 0.9303 | 0.6157 | 2.25 | +0.059 |
| FP-003 | A-screen | FP-FRAGFIX-CONSTR | 0.9272 | 0.6581 | 1.28 | +0.056 |
| FP-004 | A-screen | FRAGFIX mfg20/mcl60 | 0.9218 | 0.6439 | 1.037 | +0.051 |

- Trend: **big jump then deliberate plateau.** The decisive move was FP-001→FP-002 (+0.089 spec from adding the FP-aware objective). FP-003/004 held spec roughly flat while pulling gene_count 2.25→1.037 — i.e. we stopped chasing the spec axis and paid down the coherence guardrail. This is correct sequencing, NOT stagnation.
- Gap-to-ceiling: 0.150 (FP-001) → 0.061 (FP-002) → 0.070 (FP-004). The slight widening at FP-004 is the spec-for-coherence trade, not regression. Half-life: one architecture move (FP-002) closed ~60% of the FP-001 gap-to-ceiling. Indeterminate beyond that until a NEW axis is tried.

### Repeated failure pattern

| Pattern | Affected ITERs | Evidence | Root cause |
|---|---|---|---|
| Learned structured decoders HURT specificity | B0-010 CRF-vec (0.7138, worst), FP-002 FP-SEGNT-CRF (0.8298, high-var) | structured decoders raise recall by coherent over-prediction into intergenic DNA | without FP-aware emissions, a coherence head only makes *coherent wrong genes*. Already in docs/10. The deterministic constrained post-proc (not a learned CRF) is what worked. |
| Per-base FP-aware loss FRAGMENTS genes | FP-002 FPLOSS 2.25, FUSION 3.40 | gene_count_ratio explosion | FP-aware per-base objective with no coherence term shatters gene bodies; fixed only by the post-proc layer (FP-003/004). |

### Early signal we skipped (the one that matters)

| Signal | First appeared | Why it matters now | Suggested re-examination |
|---|---|---|---|
| **screen_anchor + ceiling are PROVISIONAL: yeast+fly are low-UTR gene-dense OUTLIERS where the full-transcript UTR ruler-fix is ~no-op; and both are likely IN the baselines' training/pretraining corpora** | tri-review 3/3 at REVISE-INTERGENIC-PRIMARY-M1; restated in docs/05 "Next up" | **The entire promote-ready conclusion rests on 2 in-corpus, low-UTR species.** The project's north star is *cross-species intergenic stability on held-out clades* — which has NEVER been evaluated. spec 0.9218 may not survive a UTR-rich / held-out-clade species. | Before any HARD Track-B gating: re-derive anchor + ceiling on a UTR-rich / held-out-clade set (deterministic once those prediction GFFs exist). docs/05 lists this TODO but it has been deferred across FP-001..004. |
| pretrained_ceiling recompute | docs/05 line 59 still open, but ACTIVE_GOAL already has 0.9917 | minor staleness | close the docs/05 TODO (already done in contract). |

### Abandoned route worth reconsidering?
- docs/09 EMPTY → nothing to revive. CRF-vec is an ablation-on-record (recall↑/spec↓), correctly NOT scaled under the new ruler. No re-entry needed.

### Subagent / scout fan-out gaps
- Pending integration queue empty; no literature re-scan since /sota-inventory (pre-ruler-change). A cheap read-only lit sweep for *FP-aware / specificity-targeted gene-annotation objectives* + *UTR-rich cross-clade benchmarks* could inform Track B but is NOT blocking.
- FP-series were run-and-evaluate (short) → no submit-and-handoff wait windows to fill → no missed scout tasks.

### Recommendation (advisory only)
- [x] **escalate to user decision** (primary): the Track-B promotion go-ahead is already a user gate; that stands. The retrospective does NOT change the promote-ready verdict for *direction-selection*.
- [x] **run focused ablation to isolate root cause** (secondary, STRONG): make "re-derive anchor + ceiling on a UTR-rich / held-out-clade species set" the FIRST step of Track B (or a cheap pre-Track-B gate), NOT a later afterthought. Rationale: every spec number to date is on yeast+fly; the north-star claim (cross-species intergenic stability) is literally untested. This is the highest-value de-risking move before spending >24h Track-B compute.
- [ ] continue current path unchanged / pivot axis / revisit abandoned route / return to literature — not selected.

### Advisory boundary
This retrospective is advisory only. It does NOT overwrite docs/03, cancel jobs, override promotion, or write docs/09. If the user wants the held-out-clade re-anchoring to become a binding pre-Track-B requirement (a benchmark-contract change), it must go: /tri-review (with this entry as context) → /pivot → user confirmation.


## Pivot Decision: REANCHOR-HELDOUT-M7 (2026-06-12)
### Inputs: docs/07#tri-review-reanchor-heldout-m7 (2/3 DEGRADED, 2-0 consensus), docs/06 result. Track A screen, NON-CLAIM, retrospective-derived re-anchor gate.
### DECISION (autonomous within screen scope; 2/2 consensus): RE-ANCHOR GATE PASSED -> ③ Track-B is GREEN-LIT pending USER GO-AHEAD, with a MANDATORY job#1 redirection.
- Held-out re-anchor gate did its job: the candidate's intergenic_specificity advantage TRANSFERS cross-clade (0.9604 vs held-out anchor 0.8054, +0.155, all 5 seeds; LARGER margin than yeast+fly +0.078; near ANNEVO ceiling 0.9824). The retrospective worry ("spec numbers only on low-UTR in-corpus outliers") is REFUTED for the specificity axis. Methodology clean (leakage/fairness/VAL-band PASS).
### CORRECTION (both reviewers, ADOPTED): the M7 verdict was OVER-STATED. The candidate is **Pareto-ADMISSIBLE** (passes the R6 screen promotion contract: spec strictly>anchor AND gbF1>=floor AND macro>=gate AND gene_count<=1.25), NOT "Pareto-beat BOTH co-primary axes" — it DOMINATES AXIS-1 (spec) but on AXIS-2 gbF1 0.6664 < anchor 0.7099 (it loses the publishable axis vs the raw-DNA anchor). docs/06/00/10 wording corrected.
### ③ Track-B (USER GO-AHEAD required — >24h compute / new long sub-iteration):
  - **job#1 (MANDATORY redirection, both reviewers): gbF1 RECOVERY, not more spec.** richer strand/phase/splice-aware MULTI-CLASS structured output (semi-CRF / segment-level + FP-aware objective) on frozen SegmentNT features; target constrained_gene_body_F1 climbing toward ANNEVO ceiling 0.8976 WHILE intergenic_specificity stays >=~0.95. The gbF1->ceiling gap 0.231 is ARCHITECTURAL (>>0.05 anti-tuning threshold) — multi-class is the lever, tuning will not close it.
  - **mandatory eval upgrades (parallel, non-blocking):** (a) add Gallus gene-sparse MACROCHROMOSOME stratum (the untested, hardest spec regime — current chicken subset is gene-dense microchromosomes); report 3 strata (Arabidopsis / Gallus-micro / Gallus-macro). (b) deterministically AUDIT SegmentNT pretraining species membership for arabidopsis+gallus -> pre-claim leakage gate (held-out novelty is feature-level-contaminated; matters for full/scale claim, NOT for this non-claim screen).
  - pass conditions for the gated entry: spec still >> same-budget anchor; gbF1 no longer < raw-DNA anchor (or clear recovery trend); no macrochromosome specificity collapse; gene_count<=1.25 + under-prediction not worsening.
  - staged UNFREEZE/fine-tune SegmentNT = SEPARATE later axis (attribution); keep deterministic constrained post-proc as coherence layer; >=5-8 seeds + CI.
### Optional /revise-goal (human-gated): record held-out anchor (spec 0.8054 / macro 0.7804 / gbF1 0.7099) + ANNEVO ceiling 0.9824 alongside the yeast+fly anchor (do not replace; both are valid same-budget references on different species sets).
### Anti-tuning: gbF1 gap 0.231 >> 0.05 -> tuning_allowed=false on the gbF1 axis -> Track-B job#1 MUST be a structural (multi-class output) change, not lr/batch. Compliant.
### Docs updated: 04,05,06,07,08,10,00.


## Pivot Decision: TB-GBF1-MULTICLASS-M8 (2026-06-12)
### Inputs: docs/07#tri-review-tb-gbf1-multiclass-m8 (2/3 DEGRADED, 2-0 consensus), docs/06 result. Track B, NON-CLAIM.
### DECISION (autonomous, 2/2 consensus): DROP multi-class (M8 bet refuted); NEXT AXIS = STAGED UNFREEZE / fine-tune SegmentNT, entered via a BOUNDED screen-profile PREFLIGHT (not a direct >24h scale-up).
- M8 result: multi-class structured output did NOT recover gbF1 on CLEAN held-out plants {arabidopsis,rice}: mc gbF1 0.7189 NOT > 3c 0.7392 (worse + gcount 0.66 under-prediction). The gbF1->ANNEVO-ceiling gap (~0.16) is NOT closed by richer decoder labels -> structural; the most likely cause (both reviewers) is the FROZEN features (ANNEVO ceiling 0.8976 is end-to-end-trained; frozen-head caps ~0.74).
- Clean POSITIVE (new honest headline): 3c-candidate (frozen SegmentNT + FP-aware + constrained) PARETO-beats the raw-DNA same-budget anchor on CLEAN plants on BOTH co-primary axes (spec 0.9663 vs 0.9045 +0.062; gbF1 0.7392 vs 0.6960 +0.043) — leakage-free (SegmentNT backbone excludes plants), replacing M7's chicken-contaminated +0.155 with an honest clean dual-axis win. The foundation-feature route is VALIDATED clean; the open problem is gbF1 headroom to SOTA.
### NEXT = ④ STAGED-UNFREEZE PREFLIGHT (bounded screen; the bounded screen compute is within autonomy, but the IMPLEMENTATION is non-trivial -> user go-ahead recommended; the full >24h unfreeze scale-up AFTER is a hard user gate):
  - Mechanism: unfreeze the TOP N layers of SegmentNT + the existing 3c FP-aware constrained head; low LR; backprop into the backbone. Tests directly whether frozen features are the gbF1 cap.
  - **IMPLEMENTATION REALITY (key)**: SegmentNT is JAX/Haiku; the current head is torch with SEPARATE frozen jax feature extraction. Unfreezing needs an in-process trainable path: either (a) a JAX/Haiku head + jax fine-tune of SegmentNT, or (b) a torch port of SegmentNT (e.g. HF AutoModel if available). This is a substantial new implementation (the M7/M8 jax-extract / torch-head split was deliberate because they don't coexist). So ④ is a real new goal, not a config flag.
  - Bounded-screen success (both reviewers): gbF1 directionally > frozen 3c 0.7392 on clean plants, intergenic_specificity not collapsing, gene_count sane (avoid mc-style under-call). NOT set near ANNEVO. If no directional gbF1 gain -> backbone-only domain-adaptation (masked/self-sup or pseudo-label) or a different foundation model.
  - Leakage discipline (codex): clean species/chrom split; no test labels in early-stop/decode-tuning; stay raw-DNA ab-initio; same 3-class collapse ruler; no test-truth gene_count calibration. Pre-claim: keep evaluating on segmentation-clean species (plants) — chicken/fly stay contaminated.
### NOT pursued: multi-class scaling (refuted); accept-frozen-ceiling (gap too big -> can't reach north star); evidence/RNA (breaks ab-initio purity). chicken-macrochromosome stratum DEFERRED (the binding finding is the clean-species gbF1 negative; chicken is contaminated so its robustness doesn't change the route pivot).
### Anti-tuning: gbF1 gap 0.16 >> 0.05 -> structural; unfreeze is an architecture axis (not lr/batch) -> compliant.
### Optional /revise-goal (human-gated): record the clean 3c-candidate dual-axis result (spec 0.966 / gbF1 0.739 on clean plants) + the frozen-feature gbF1 ceiling finding.
### Docs updated: 04,06,07,08,10,00.

## failed_run: M9-UNFREEZE-BACKBONE CK3 (2026-06-13)

- **事件**: 9-job 3-arm batch (M9-UNFREEZE-L{0,2,4}-s{0,1,2}, jobs 8575441-49) 全部 **TIMEOUT** (Elapsed 11:50:13 撞 shared-gpu 12h 上限被 SIGKILL @ 2026-06-13T09:55-09:58, gpu034/048/053 等)。0 metrics 产出。
- **幽灵态**: `#SBATCH --time=11:50:00` 被杀时 trap EXIT 未触发 → STATUS 卡 RUNNING。已确定性重置为 FAILED_TIMEOUT。
- **根因 (确定)**: (1) 算力超预算——8 epoch × 2物种(arab+rice) × sample 0.3 的 backbone fine-tune 在 12h 内跑不完 (smoke arab-only/0.1/2ep 已 1h56m)。(2) trainer `train_unfreeze_backbone.py:228-230` best 权重仅存内存、从不 torch.save → TIMEOUT 全白跑、无法 resume。
- **正面信号 (CK2 smoke, job 8565782)**: arabidopsis 2 epoch 已 gbF1 0.808 ≥ frozen-full 0.805, spec 0.959 — unfreeze 抬 gbF1 的假设有强信号, **只是预算配置错**, 不是方向错。
- **裁决**: failed_run → STOP, 不静默续跑。修复=减 epoch (4 足够) + 加磁盘 checkpoint 兜底; 重跑预算选择上交主人 (大算力重跑 = 需确认节点)。NON-CLAIM, screen profile 不变。

## Pivot: TB-UNFREEZE-BACKBONE-M9 (CK3, 2026-06-14)

- **Inputs**: docs/06 M9-CK3 result (3 jobs 8667188/89/90 COMPLETED single-species arabidopsis), validate_goal=not_yet (run_ok+semantic_ok true, guardrails pass, primary_progress_gate blocked by constrained_gene_body_F1=0). Screen, NON-CLAIM.
- **Scientific verdict**: unfreeze backbone VALIDATED — monotonic DUAL-axis gain with depth: gbF1 0.8284(L0)->0.8544(L2)->0.8759(L4), spec 0.9656->0.9669->0.9754. L4 gbF1 0.876 nears ANNEVO ceiling 0.898. Confirms M8 diagnosis (frozen features capped gbF1); fine-tuning emissions lifts BOTH axes (no trade-off).
- **Why not_yet (key)**: constrained_gene_body_F1 = gbF1 if intergenic_FPR<=threshold else 0 (eval line 237; screen threshold 0.02). L4 FPR=0.0246 JUST exceeds 0.02 -> constrained=0 -> goal success_criteria (constrained>0) unmet. L4 spec 0.9754 needs only 0.98 (FPR 0.0246->0.02, gap 0.0046).
- **DECISION (autonomous, reversible single-direction)**: continue-current-route -> DEEPER unfreeze. Since spec rises monotonically with unfreeze depth, L6/L8/full-unfreeze should push FPR below 0.02, unlocking constrained F1 + passing the gate. This is the clear next lever (not tuning — it's the same structural axis extended).
- **Track B scope (needs USER gate — >24h compute)**: promote unfreeze to Track B = {deeper unfreeze L6/L8/full} x {multi-seed >=3} x {cross-species: add rice + held-out clade} for CI + cross-species generalization (the north-star). This is the scale-up that tests whether the dual-axis advantage holds with seeds/species and whether deeper unfreeze breaks the 0.98 spec / 0.02 FPR barrier.
- **Caveats carried forward**: (1) gcount drifts down with depth (1.03->0.90->0.82) — watch under-prediction at deeper unfreeze. (2) single-species/1-seed — CI needed. (3) transcript_count_ratio 0.445 (predicts ~half the transcripts) — isoform under-prediction, investigate at scale.
- **tri-review**: formal 3-CLI tri-review deferred to the Track B PROMOTION gate (claim-adjacent); this screen-level continue-direction pivot is reversible/non-claim/single-axis (SINGLE_REVIEW_CONTINUATION-class).

## Goal Revision PROPOSAL (advisory, 2026-06-14): constrained_gene_body_F1 硬门→软化

- **触发**: 指标审计(docs/10 2026-06-14) — constrained 硬门(eval:237 `gbF1 if FPR<=thr else 0`)是 R6 前单-primary 遗留, R6 后 FPR/unconstrained-gbF1 已各自独立 guardrail → 硬门功能冗余 + 悬崖效应误杀 M9 L4(FPR 0.0246 仅超 0.02 门 0.0046, 真实 gbF1 0.876 被归零)。
- **可比性复核(6维)**: dataset/split/preprocessing/weights/test-time 均不变(同一 eval 管道、同一 ruler); 仅 constrained_gbF1 的**聚合函数**变 → 改的是指标定义本身, 影响所有历史 constrained 值的可比性。anchor_gene_body_F1=0.5576 是 UNCONSTRAINED gbF1, **不受影响**(软化只动 constrained)。
- **方案 A (软截断/平滑惩罚)**: constrained_gbF1 = gbF1 × min(1, threshold/FPR) — FPR<=thr 时=gbF1(连续), 超标时按比例衰减不归零。M9 L4 → 0.876×(0.02/0.0246)=0.712。保留排序+平滑+不作弊(纯打分变换, 不用 reference 裁预测)。
- **方案 B (推荐, 弃冗余门)**: 承认硬门冗余, success_criteria AXIS-2 改用 `gene_body_F1_unconstrained`(已是 guardrail floor 0.5276), constrained 降级为 diagnostic; FPR 由独立 guardrail 守(screen advisory 0.02 / full hard 0.01)。最符合 R6 dual-co-primary 两轴独立连续哲学。M9 L4: gbF1 0.876>>0.5276 ✓ + spec 0.975>anchor 0.871 ✓ → PASS(FPR 0.0246 仅 screen advisory 超标, 不阻塞 screen; full 仍须 hard 0.01)。
- **对 M9 影响(关键)**: 两方案下 M9 L4 均过门 → M9-DEEP(L6/L8/L12)从"破门必需"降级为"验证更深双轴趋势/为 full hard FPR 0.01 做准备"(仍有价值, 不 kill)。
- **状态**: PROPOSAL only. 须 tri-review 复核 + 用户确认才落盘(改 ACTIVE_GOAL success_criteria + eval_gene_body_mask.py + 重算受影响历史记录)。未确认不改。

## Goal Revision LANDED (2026-06-14): REVISE-CONSTRAINED-SOFTGATE
- tri-review 3/3 APPROVE-B + 用户确认 → 落盘。改动 (最小面, 零代码改动, 利用已有 profile-scoping):
  - ACTIVE_GOAL.json success_criteria: `constrained_gene_body_F1` rule 加 `profiles:[full,scale]` → screen 时 check() 自动 skip; AXIS-1 intergenic_specificity 保持全 profile。
  - screen AXIS-2 质量由已存在的 HARD guardrails 守 (gene_body_F1_unconstrained>=0.5276 + macro_intergenic_specificity>=0.7978)。
  - full/claim AXIS-2 保持 constrained + SOTA-comparable (不动); full/scale intergenic_FPR<=0.01 HARD 保留。
  - eval_gene_body_mask.py 不动 (constrained 仍产出作 full gate + screen diagnostic); validate_goal.py 不动 (check 已支持 profiles)。
  - 备份: ACTIVE_GOAL.json.bak-20260614。_revision_log 已记。
- **验证生效**: M9-UNFREEZE-L4-s0 重判 screen → status `not_yet`→**`progress`** (primary_progress_gate.ok=True: spec 0.975>0 ✓ + constrained skipped[full,scale]; guardrails.ok=True)。
- **M1-M9 新口径重判**: screen 历史值不再被 constrained 硬零误杀; AXIS-2 screen 质量看 unconstrained gbF1 (已记录值不变, 仅 gate 判读变)。M9 L4 screen PASS。

## Pivot Decision: TB-UNFREEZE-BACKBONE-M9-DEEP (2026-06-14)

### Inputs consumed
- `$result-log`: docs/06 `Result: TB-UNFREEZE-BACKBONE-M9-DEEP`
- `$tri-review`: docs/07 `Tri-Review: TB-UNFREEZE-BACKBONE-M9-DEEP`, 3/3 quorum
- validate_goal: `outputs/M9-UNFREEZE-L12-s0/metrics/validate_goal.json` -> `status=progress`, screen non-claim
- Resource profile: screen / Track-B preflight, NON-CLAIM

### Current evidence summary
M9-DEEP directly answered the M9 CK3 blocker: L4 had gbF1 0.8759 but FPR 0.0246; deeper unfreeze made all L6/L8/L12 pass FPR<=0.02. Best arm L12 reached intergenic_specificity 0.9810, FPR 0.0190, gbF1 0.9035, constrained_gbF1 0.9035, gene_count_ratio 0.792 on arabidopsis seed0. This is the strongest NT-v2 route evidence so far, but it is not a claim because it is single species / single seed and full/scale `sota_benchmark` is still draft.

### SOTA gap
| Metric | Current | Reference | Gap (abs) | Severity |
|---|---:|---:|---:|---|
| intergenic_specificity vs screen_anchor | 0.9810 | 0.8710 | +0.1100 | strong screen win |
| gene_body_F1 vs arabidopsis ANNEVO ceiling reference | 0.9035 | ~0.8980 | +0.0055 | promising but local/non-claim |
| intergenic_FPR vs full/scale hard guardrail | 0.0190 | 0.0100 | +0.0090 | claim blocker remains |
| gene_count_ratio vs full/scale upper guardrail | 0.792 | <=1.25 | pass | under-prediction watch |

### Sanity check
- [x] Three independent CLI reviewers succeeded (3/3).
- [x] No reviewer raised a semantic failure.
- [x] No reviewer raised immediate leakage proof, but all require NT-v2 pretraining species coverage audit before claim.
- [x] Metrics finite; loss decreases; no OOM/NaN/Traceback.
- [x] Resource profile does not support claim.
- [x] `validate_goal.py` profile-aware regression fixed and tests pass.

### Tri-review summary (record ALL reviewer conclusions)
| Reviewer | Judgment | Next action proposed | Main concern | Confidence |
|---|---|---|---|---|
| A · Claude | `continue-current-route` | M9-L12 multi-seed + clean plants as primary; GENERanno LoRA + 3-class head as parallel second direction | screen only, single seed, FPR still >0.01, NT-v2 species audit, gcount under-prediction | Medium-High |
| B · Codex | `scale-to-track-b` | Promote M9-L12 to multi-seed/cross-species Track B; start GENERanno LoRA in parallel if budget allows | benchmark not frozen, single species/seed, FPR 0.019 > 0.01, backbone coverage audit, gcount 0.792 | Medium-High |
| C · Antigravity | `scale-to-track-b` | M9-L12 primary Track B; GENERanno LoRA parallel hedge because native FPR already satisfies full/scale specificity | FPR hard gate, SOTA benchmark, NT-v2 contamination audit, gene_count under-prediction | High |

Consensus: primary next run = M9-L12 multi-seed + clean plants `{arabidopsis,rice}`. Parallel challenger = GENERanno LoRA + 3-class intron-aware head if GPU budget allows.
Disagreement: label only (`continue-current-route` vs `scale-to-track-b`); no material disagreement in proposed next action.
Quorum / degraded review status: 3/3.

### Reviewer-proposed directions
| # | From reviewer | Direction | major_axis | mechanism_delta | Orthogonal to others? | Into cohort? |
|---:|---|---|---|---|---|---|
| 1 | A · Claude | M9-L12 multi-seed + clean plants | training_signal / scale-validation | scale validated L12 unfreeze across seeds/species | primary | yes |
| 2 | A · Claude | GENERanno LoRA + 3-class head | backbone | replace NT-v2 with CDS-annotator-pretrained backbone and add intron-aware head | yes | yes, if budget |
| 3 | B · Codex | L12 + stronger FP-aware objective later | objective | target full/scale FPR<=0.01 after stability evidence | yes but downstream | not this round |
| 4 | B · Codex | L12 + segment/gene-count regularization later | decoder/regularization | counter L12 gene_count under-prediction | yes but downstream | not this round |
| 5 | C · Antigravity | Enhanced structured decoder as ablation | decoder | lightweight CRF/Viterbi after emissions improve | yes but lower priority | not this round |

### Is tuning justified?
No as a primary decision. This is a structural route-validation point, not an LR/dropout tuning point. Limited objective/postproc calibration may be justified only after M9-L12 multi-seed/cross-species shows the route is stable and the remaining blocker is specifically full/scale FPR<=0.01.

### Architecture hypothesis status
Supported. Deeper NT-v2 unfreeze improves emissions enough to raise gbF1 and reduce intergenic spillover simultaneously. GENERanno evidence suggests an orthogonal backbone with even stronger specificity but missing intron/coherence.

### DECISION
- [x] **Scale data / training: M9-L12 multi-seed + clean plants `{arabidopsis,rice}` as the single primary direction.**
- [x] **Parallel cohort: GENERanno LoRA + 3-class intron-aware head as a challenger backbone direction if GPU budget permits.**
- [ ] Tune current architecture
- [ ] Abandon route
- [ ] Claim SOTA

### Why this decision
M9-L12 is already in the project's native 3-class/intron-aware pipeline and has the strongest gbF1 evidence. It solved the immediate L4 blocker and is the lowest-risk path to a claim candidate. GENERanno should not replace it yet because native gbF1/coherence are poor, but its FPR 0.004-0.005 is uniquely strong; LoRA+3-class directly tests whether that specificity can be combined with intron-aware coherence. Pure tuning is premature because the evidence gap is still dominated by species/seed/full-scale validation and claim-contract blockers.

### Parallel cohort this round
| Slot | EXP ID (new) | Direction | major_axis | mechanism_delta | Track | Resource profile |
|---|---|---|---|---|---|---|
| primary | `M10-M9L12-CLEANPLANTS` | M9-L12 multi-seed on clean plants | training_signal / scale-validation | same L12 architecture, seeds + arabidopsis/rice validation | B | full/screen-to-full preclaim |
| parallel-1 | `M10-GENERANNO-LORA-3C` | GENERanno LoRA + 3-class head | backbone | CDS-annotator backbone with intron-aware output | A/B challenger | screen |

Shared-code conflict: likely yes if both alter shared trainer/data code; if implemented simultaneously, use careful exp_id isolation or `$workspace-matrix`. If run sequentially, no worktree required.

### Required blockers before claim
- Freeze ANNEVO-compatible `sota_benchmark` under the same evaluator/ruler.
- Audit NT-v2 and GENERanno pretraining species/sequence overlap for arabidopsis/rice claim cleanliness.
- Show multi-seed and multi-species stability, including macro specificity.
- Clear full/scale FPR<=0.01 or explicitly revise the claim contract through human gate.

### TODO update
- [ ] update docs/05 TODO status for tri-review/pivot completed and next M10 cohort.
- [ ] update docs/04 M9-DEEP tri-review/pivot fields.
- [ ] run `$note-gate` / evidence register capture for the route decision.

## Pivot Decision: M10-COMBINED-M9L12-GENERANNO (2026-06-15)

### Inputs consumed
- `$result-log`: docs/06 `Result: M10-M9L12-CLEANPLANTS` and `Result: M10-GENERANNO-LORA-3C-SMOKE`
- `$tri-review`: docs/07 `Tri-Review: M10-COMBINED-M9L12-GENERANNO`, 3/3 quorum
- validate_goal: `outputs/M10-M9L12-CLEANPLANTS-s{0,1,2}/metrics/validate_goal.json` -> all `status=progress`, screen non-claim
- Resource profile: screen/smoke, NON-CLAIM

### Current evidence summary
M10 mainline validates M9-L12 across clean plants and seeds: mean intergenic_specificity `0.9826`, FPR `0.0174`, macro_specificity `0.9801`, gbF1 `0.8398`, gene_count_ratio `0.897`. It is far above the same-budget screen anchor, but formal claim remains blocked by full/scale FPR `<=0.01`. The most visible failure mode is arabidopsis: unconstrained gbF1 is high (`0.896-0.903`), but constrained_gbF1 is zero for all seeds because FPR `0.022-0.0285` exceeds the 0.02 sensitivity threshold. GENERanno LoRA smoke is engineering-positive but metric-negative and should not be submitted as the prepared screen run.

### SOTA gap
| Metric | Current | Reference | Gap (abs) | Severity |
|---|---:|---:|---:|---|
| intergenic_specificity vs screen_anchor | 0.9826 | 0.8710 | +0.1116 | strong screen win |
| gene_body_F1 vs screen floor | 0.8398 | 0.5276 | +0.3122 | strong pass |
| intergenic_FPR vs full/scale hard guardrail | 0.0174 | 0.0100 | +0.0074 | claim blocker |
| macro_specificity vs gate | 0.9801 | 0.7978 | +0.1823 | strong pass |
| gene_count_ratio vs full/scale upper guardrail | 0.897 | <=1.25 | pass | coherent |

### Sanity check
- [x] Three independent CLI reviewers succeeded (3/3).
- [x] No reviewer raised semantic failure for M10 mainline.
- [x] Metrics finite and parseable; all STATUS files COMPLETED; validate rerun with STATUS file path.
- [x] No OOM/NaN/Traceback in logs.
- [x] Resource profile does not support claim.
- [x] Main blocker is known and measurable: FPR/constrained operating point.

### Tri-review summary
| Reviewer | Judgment | Next action proposed | Main concern | Confidence |
|---|---|---|---|---|
| A · Claude | `continue-current-route` | `M11-CONSTRAINED-FIX-FP-CALIBRATE`: diagnose constrained postproc, adjust decode params/soft constraints, increase FP pressure | arabidopsis constrained zero; FPR>0.01; SOTA benchmark draft; LoRA route damages specificity | Medium-High |
| B · Codex | `continue-current-route` | `M11-L12-SPEC-CALIBRATION`: validation-only operating point sweep, stronger FP loss, two-stage decoder calibration | FPR must drop 0.0174->0.01; avoid test-set tuning; LoRA overpredicts | Medium |
| C · Antigravity | `continue-current-route` | Decode/calibration sweep first; if insufficient, retrain M9-L12 with higher `fp_lambda` or asymmetric focal loss | FPR guardrail; rice gbF1 may drop; LoRA cannot learn introns | Medium-High |

Consensus: continue M9-L12 as primary; next experiment must target specificity/FPR calibration. Park GENERanno LoRA until redesigned.
Disagreement: only implementation order (post-hoc decode first vs retrain with stronger loss). This can be encoded as a staged M11: cheap decode/diagnostic first, then retraining if needed.
Quorum / degraded review status: 3/3.

### Reviewer-proposed directions
| # | From reviewer | Direction | major_axis | mechanism_delta | Orthogonal to others? | Into this round's cohort? |
|---:|---|---|---|---|---|---|
| 1 | A · Claude | constrained-postproc diagnosis/fix | decoder | repair hard postproc and arabidopsis constrained zero | complements FP loss | yes |
| 2 | A · Claude | stronger FP-aware loss | objective | `fp_lambda` increase / hard-negative weighting | complements decoder | yes |
| 3 | B · Codex | validation-only operating point sweep | decoder/calibration | threshold/margin/min length/gap fill selected on val only | complements retrain | yes |
| 4 | B · Codex | two-stage decoder calibration | decoder | preserve recall emissions but calibrate constrained output | overlaps 1/3 | maybe |
| 5 | C · Antigravity | post-hoc intergenic emission penalty | decoder/calibration | bias against genic transitions to hit FPR<=0.01 | overlaps 1/3 | yes |
| 6 | C · Antigravity | retrain with higher FP pressure if post-hoc fails | objective | higher `fp_lambda` / asymmetric focal loss | complements decode | yes |

### Is tuning justified?
Yes, but only targeted near-claim calibration is justified. This is not generic LR/dropout tuning: the architecture is validated, and the remaining gap is a concrete operating-point/FPR blocker. Test labels must not be used for parameter selection.

### Architecture hypothesis status
Supported. NT-v2 L12 unfreeze + intron-aware head remains the mainline. GENERanno LoRA as currently configured is weakened/parked, not abandoned forever.

### DECISION
- [x] **Tune current architecture: targeted M9-L12 specificity calibration / constrained-FPR repair before full-scale promotion.**
- [ ] Change backbone now.
- [ ] Submit prepared GENERanno LoRA screen now.
- [ ] Claim SOTA.
- [ ] Abandon M9-L12.

### Why this decision
M10-M9L12 already solves the large architecture gap: it is stable across clean plants/seeds, has strong gbF1, and coherent gene counts. Replacing the backbone now would discard the best evidence. Direct full/scale would likely reproduce the FPR blocker expensively. The right next move is a bounded M11 that uses validation-only calibration and/or stronger FP objective to push FPR toward `<=0.01` while monitoring gbF1, especially rice recall. GENERanno LoRA is not ready because the smoke damaged the very specificity advantage that justified it.

### Best next architecture moves
| Priority | Move | Expected mechanism | Goes to which EXP / Track |
|---:|---|---|---|
| 1 | Add raw-pred / validation decode sweep for M9-L12 | Separate emissions quality from postproc operating point; choose thresholds/`min_cds_len`/`max_fill_gap` without test leakage | `M11-L12-SPEC-CALIBRATION` |
| 2 | Constrained decoder/intergenic penalty | Reduce arabidopsis genic spillover and avoid constrained zero | `M11-L12-SPEC-CALIBRATION` |
| 3 | Stronger FP-aware objective if decode sweep insufficient | Train emissions to be more conservative in intergenic / UTR-adjacent hard negatives | follow-up M11 retrain arm |
| 4 | Redesigned GENERanno LoRA only after mainline calibration | Preserve native specificity while learning intron continuity | parked challenger, no immediate submit |

### Parallel cohort this round
- **Primary direction (single)**: `M11-L12-SPEC-CALIBRATION`, M9-L12 targeted specificity calibration.
- **Parallel cohort**: none by default. GENERanno LoRA is parked; do not spend another GPU slot until a redesigned schedule is written and reviewed.

| Slot | EXP ID (new) | Direction | major_axis | mechanism_delta | Track | Resource profile |
|---|---|---|---|---|---|---|
| primary | `M11-L12-SPEC-CALIBRATION` | M9-L12 FPR/constrained repair | decoder/objective | val-only decode sweep + optional stronger FP loss | B preflight | screen |

Shared-code conflict: no parallel shared-code conflict if GENERanno is parked. M11 likely needs trainer support for raw prediction/logit saving or configurable calibration; run `$implement` + `$code-review-gate` before submission.

### Required blockers before claim
- Freeze ANNEVO-compatible `sota_benchmark`.
- Audit NT-v2/GENERanno pretraining species or sequence overlap.
- Hit full/scale FPR `<=0.01` without using test labels for calibration.
- Confirm constrained_gbF1 is not zeroed on arabidopsis and macro specificity remains strong.

### TODO update
- [x] update docs/05_todo.md with M10 done + M11 next.
- [x] update docs/07_tri_review.md with all reviewer conclusions.
- [x] update docs/08_pivot_decisions.md with single primary decision.
- [x] update docs/11_master_plan.md navigation to M11 calibration.
- [ ] if launching M11, run `$implement` then `$code-review-gate` before `sbatch`.

## Pivot Decision: M11-L12-SPEC-CALIBRATION (2026-06-16)

### Inputs consumed
- `$result-log`: docs/06 `Result: M11-L12-SPEC-CALIBRATION`
- `$tri-review`: docs/07 `Tri-Review: M11-L12-SPEC-CALIBRATION`, 2/3 `DEGRADED_REVIEW` quorum (Codex + Antigravity succeeded; Claude failed marker heuristic)
- validate_goal: `outputs/M11-L12-SPEC-CALIBRATION-s{0,1,2}/metrics/validate_goal.json` -> all `status=progress`, screen non-claim
- Resource profile: screen / Track-B preflight, NON-CLAIM

### Current evidence summary
M11 answers the M10 blocker directly. Validation-only calibration over saved raw emissions reduces aggregate FPR from M10 mean `0.0174` to M11 mean `0.0087`; all seeds clear `FPR<=0.01`, constrained gbF1@0.01 equals unconstrained gbF1, and gene_count_ratio remains coherent (`1.003` mean). This indicates the immediate FPR problem was primarily an operating-point/decode calibration issue, not a failure of NT-v2 L12 emissions or the FP-aware objective.

### SOTA gap
| Metric | Current | Reference | Gap (abs) | Severity |
|---|---:|---:|---:|---|
| intergenic_specificity vs screen_anchor | 0.9913 | 0.8710 | +0.1203 | strong screen win |
| gene_body_F1 vs screen floor | 0.8178 | 0.5276 | +0.2902 | strong pass |
| intergenic_FPR vs full/scale hard guardrail | 0.0087 | <=0.0100 | pass | blocker cleared at screen scale |
| macro_specificity vs gate | 0.9909 | 0.7978 | +0.1931 | strong pass |
| gene_count_ratio vs full/scale upper guardrail | 1.003 | <=1.25 | pass | coherent |
| published SOTA benchmark | unknown | not frozen | unknown | claim blocker |

### Sanity check
- [x] At least two independent CLI reviewers succeeded: 2/3 `DEGRADED_REVIEW`.
- [x] No successful reviewer raised a semantic failure.
- [x] Metrics finite and parseable; all STATUS files COMPLETED.
- [x] Loss trend matches M10 and is sane.
- [x] Seed variance is small for screen: specificity `[0.9908,0.9921]`, gbF1 `[0.8086,0.8277]`.
- [x] Calibration is validation-only; no test-label parameter selection is recorded.
- [x] Screen resource profile does not support claim.
- [ ] Published SOTA benchmark and pretraining-overlap audit remain unresolved for claim.

### Tri-review summary
| Reviewer | Judgment | Next action proposed | Main concern | Confidence |
|---|---|---|---|---|
| A · Claude | failed for quorum | not counted; two non-empty attempts lacked required `Overall judgment` marker | unusable as independent structured review | N/A |
| B · Codex | `scale-to-track-b` | Freeze calibration rule, derive/reproduce ANNEVO-compatible SOTA benchmark, evaluate calibrated M9-L12 at full/scale with per-species FPR sensitivity | benchmark unfrozen; screen cannot claim; arabidopsis seed2 FPR `0.0111`; no checkpoint; NT-v2 overlap audit | Medium-High |
| C · Antigravity | `scale-to-track-b` | Promote current M9-L12 + validation-only calibration to full/scale testing; do not change objective now | benchmark draft; screen cannot claim; per-species FPR variance; small validation split | High |

Consensus: successful reviewers agree to scale/promote the calibrated M9-L12 route and defer stronger FP objective.
Disagreement: none material among successful reviewers. Quorum is degraded, so aggregate confidence is Medium.
Quorum / degraded review status: 2/3 `DEGRADED_REVIEW`.

### Reviewer-proposed directions
| # | From reviewer | Direction | major_axis | mechanism_delta | Orthogonal to others? | Into this round's cohort? |
|---:|---|---|---|---|---|---|
| 1 | B · Codex | full/scale calibrated M9-L12 comparability preflight | scale/comparability | freeze VAL-only calibration and evaluate full/scale with sensitivity bands | primary | yes |
| 2 | B · Codex | ANNEVO-compatible benchmark freeze | benchmark | reproduce/derive published SOTA under active evaluator before claim | prerequisite | yes |
| 3 | B · Codex | checkpoint/reproducibility hardening | engineering | save full/scale model checkpoints, not only raw-score artifacts | prerequisite | yes |
| 4 | C · Antigravity | increase validation chromosome/diversity | data_view/calibration | stabilize per-species operating point and reduce arabidopsis edge-case risk | complements primary | yes |
| 5 | C · Antigravity | stronger FP objective only if full/scale FPR fails | objective | asymmetric focal / hard-negative FP loss as fallback | orthogonal fallback | no immediate GPU |

### Is tuning justified?
Yes only for narrow validation-only calibration and full-scale confirmation. Generic hyperparameter tuning or stronger FP objective is not justified now because M11 already clears the direct FPR blocker with the existing emissions.

### Architecture hypothesis status
Supported. NT-v2 L12 unfreeze + 3-class intron-aware head + FP-aware loss remains the mainline. M11 strengthens the view that the residual failure was decode/calibration, not backbone/objective failure.

### DECISION
- [x] **Scale data / training: promote calibrated M9-L12 to Track-B/full-scale comparability preparation.**
- [ ] Tune current architecture broadly.
- [ ] Change objective / loss now.
- [ ] Change backbone now.
- [ ] Revive GENERanno LoRA now.
- [ ] Claim SOTA.

### Why this decision
The direct blocker named by M10 reviewers was FPR `<=0.01`. M11 cleared it across all seeds while retaining high gbF1 and coherent gene count. Launching stronger FP objective now would spend GPU on a fallback whose trigger condition has not occurred. The next limiting factor is no longer screen architecture quality; it is claim-readiness: frozen SOTA benchmark, full/scale split/protocol, pretraining-overlap audit, checkpointed reproducibility, and per-species FPR robustness.

### Best next architecture moves
| Priority | Move | Expected mechanism | Goes to which EXP / Track |
|---:|---|---|---|
| 1 | Freeze full/scale M9-L12 calibration protocol | Prevent test-set tuned decode parameters; keep M11 mechanism claim-safe | `M12-M9L12-FULLSCALE-CALIBRATED` preflight |
| 2 | ANNEVO-compatible benchmark freeze / comparability audit | Convert `ACTIVE_GOAL.status=draft` into a claim-capable benchmark contract | M12 prerequisite / `$reproduce-baselines` or `$revise-goal` |
| 3 | Full/scale calibrated M9-L12 with checkpoints and per-species sensitivity | Test whether screen FPR `<=0.01` and gbF1 hold at claim-like scale | `M12-M9L12-FULLSCALE-CALIBRATED` |
| 4 | Stronger FP objective fallback | If full/scale FPR fails again, train emissions more conservatively | follow-up only |

### Parallel cohort this round
- **Primary direction (single)**: calibrated M9-L12 full/scale/comparability preflight.
- **Parallel cohort**: none by default. Do not resume GENERanno LoRA or stronger FP-objective GPU work until the primary full/scale/comparability blockers are explicit.

| Slot | EXP ID (new) | Direction | major_axis | mechanism_delta | Track | Resource profile |
|---|---|---|---|---|---|---|
| primary | `M12-M9L12-FULLSCALE-CALIBRATED` | calibrated M9-L12 full/scale preflight | scale/comparability | frozen VAL-only calibration + checkpointed full/scale eval | Track B / full | full/scale |

Shared-code conflict: no, if M12 only extends config/sbatch/checkpointing/calibration protocol. If benchmark reproduction needs external baseline code, keep it as a separate prerequisite output namespace.

### Required blockers before claim
- Freeze ANNEVO-compatible published SOTA benchmark under the active evaluator.
- Audit NT-v2 pretraining/species overlap and document whether arabidopsis/rice or other test species are admissible for claim.
- Use a frozen validation-only calibration recipe before full/scale test inference.
- Save model checkpoints for full/scale runs.
- Report aggregate and per-species FPR sensitivity at `0.005/0.01/0.02`; pay attention to arabidopsis-like edge cases.

### TODO update
- [x] update docs/07_tri_review.md with M11 degraded-quorum review.
- [x] update docs/08_pivot_decisions.md with single primary decision.
- [x] update docs/04_experiment_iterations.md M11 tri-review/pivot status.
- [x] update docs/05_todo.md and docs/11_master_plan.md navigation.
- [x] route decision to docs/15 evidence register.

## Retrospective Review · 2026-06-17

### Scope
- Iterations covered: baseline reproduction M1 through `M11-L12-SPEC-CALIBRATION`.
- Trigger: manual / user-suspicion. User raised three strategic concerns: current work is over-focused on model performance instead of publishable story; lack of intuitive comparisons to Tiberius/Helixer/ANNEVO; current train/val/test use the same species pool, so it does not prove a fixed cross-species model generalizes.
- Focus: publication-oriented route correction before launching `M12-M9L12-FULLSCALE-CALIBRATED`.

### File timeline summary
- `docs/03_benchmark_roadmap.md`: original roadmap already required published SOTA reproduction, held-out species/clade evaluation, screen/full separation, and full/scale claim only after `sota_benchmark` freeze. The current next step should therefore close comparability and cross-species evidence, not only optimize M9.
- `docs/04_experiment_iterations.md`: M1 established evaluator/baselines; M4-M6 found SegmentNT+FP-aware/constrained route; M8 refuted multi-class but validated clean 3-class; M9-M11 focused on NT-v2 unfreeze then calibration.
- `docs/05_todo.md`: at review time, open pending items still centered on `M12-M9L12-FULLSCALE-CALIBRATED`; after this retrospective+council+tri-review, the pending queue should instead center on `M12-PREREQ-AUDIT`, `M12A-FIXEDMODEL-CROSSSPECIES`, `M12B-SAMEPANEL-BASELINES`, and bounded `M12C-GENERANNO-FAIR-CHALLENGER`. A stale/open `FP-SEGMENTNT-FEATCACHE-M7` row remains and should be reconciled before relying on old featcache status.
- `docs/06_results_log.md`: baseline runners for Tiberius/Helixer/ANNEVO exist on pilot species, but published/full comparability is still not frozen. M11 gives strong M9 metrics but still non-claim screen evidence.
- `docs/07_tri_review.md`: M9-DEEP reviewers recommended M9-L12 primary plus GENERanno LoRA challenger. M10 reviewers parked GENERanno after a smoke, not after a fair same-treatment screen. M11 reviewers shifted risk to benchmark/comparability/full-scale robustness.
- `docs/08_pivot_decisions.md`: latest binding pivot says scale calibrated M9-L12. This retrospective does not override it, but advises changing M12 design so scale-up answers publication questions.
- `docs/09_decisions_log.md`: no abandoned routes recorded; therefore GENERanno is parked, not abandoned.

### Are we doing marginal tuning?

| Verdict | Evidence |
|---|---|
| partially | M9-L12 unfreeze was a structural backbone change; M10 scaled it across clean plants/seeds; M11 validation-only calibration was a targeted operating-point fix. These are not LR/dropout tuning. However, the last two pivots narrowed into the same M9-L12 family and optimized the same arabidopsis/rice setting while benchmark comparability, fixed-model cross-species evaluation, and GENERanno as a peer pretrained backbone remained unresolved. |

### Gap trajectory

| ITER | Track | Path | Primary metric snapshot | Reference | Gap / status | Delta vs prev |
|---|---|---|---:|---:|---|---|
| M1 pretrained trio | baseline | Tiberius/Helixer/ANNEVO inference on yeast/fly | Helixer CDS F1 `0.9213`; ANNEVO `0.9197`; Tiberius `0.8608` | pretrained ceiling, non-gating | runner/evaluator validated, not same-budget anchor | established intuitive baseline harness |
| M1 same-budget anchor | screen | random-init refs | screen_anchor `0.5579` | floor `0.3735`, ceiling `0.9213` | fair Track-A ruler built | fair anchor corrected |
| M8 clean plants 3c | full/screen-like | frozen SegmentNT/3-class | spec about `0.966`, gbF1 about `0.739` | anchor gbF1 floor `0.5276` | clean Pareto-positive but gbF1 ceiling | validated clean held-out plant split |
| M9-DEEP L12 | screen | NT-v2 unfreeze, arabidopsis seed0 | spec `0.9810`, FPR `0.0190`, gbF1 `0.9035` | screen_anchor spec `0.8710` | strong single-species route signal | large jump |
| M10-M9L12 | screen/preflight | clean plants `{arabidopsis,rice}`, 3 seeds | spec `0.9826`, FPR `0.0174`, gbF1 `0.8398` | full FPR guardrail `<=0.01` | strong but claim-blocked by FPR | generalized within plant species pool, gbF1 lower than arab-only |
| M11 calibration | screen/preflight | M9-L12 + VAL-only decode calibration | spec `0.9913`, FPR `0.0087`, gbF1 `0.8178` | full FPR guardrail `<=0.01` | aggregate FPR blocker cleared at screen scale | specificity improved, gbF1 traded down modestly |

- Fitted trend: improving on the internal FP-controlled metric, but publication readiness is flat because `sota_benchmark`, same-output comparisons to Tiberius/Helixer/ANNEVO, and fixed-model cross-species validation remain open.
- Gap half-life: indeterminate. The internal metric improved quickly; the claim gap cannot be computed until the published benchmark is frozen.

### Repeated failure pattern

| Pattern | Affected ITERs | Evidence | Possible root cause |
|---|---|---|---|
| Publication evidence deferred behind architecture iteration | M9-M11 | `ACTIVE_GOAL.status=draft`; M11 pivot still lists ANNEVO-compatible benchmark freeze and overlap audit as blockers | Research loop optimized local metrics faster than it built the claim/comparison matrix |
| Same-species-pool evaluation masquerades as cross-species progress | M10-M11 | Train/val/test split is seqid/chromosome-aware within arabidopsis/rice, but train and validation include the same species as test | The current split validates within-species chromosome generalization, not a fixed model applied to unseen species/clades |
| GENERanno not given same-quality treatment as M9 | M10 | Only bounded smoke: spec `0.9491`, FPR `0.0509`, gbF1 `0.7525`, gene_count `4.43`; earlier native GENERanno specificity was strong but binary CDS fragmentation was structural | The LoRA schedule/data volume may be inadequate; one smoke cannot answer whether pretrained CDS models can become coherent with a 3-class head |
| Baseline comparison exists technically but not as a paper-facing result | M1, M11 | Tiberius/Helixer/ANNEVO runner/evaluator ledger exists, but no full/paper-like side-by-side figure using a fixed split and fixed model inference | The project has infrastructure evidence but not a compact benchmark story a reviewer can interpret |

### Early signal we skipped

| Signal | Where it first appeared | Why it matters now | Suggested re-examination |
|---|---|---|---|
| M9 tri-review recommended GENERanno as a parallel challenger | `docs/07`, `TB-UNFREEZE-BACKBONE-M9-DEEP` | User's current question exactly asks whether all pretrained models can do well, or whether M9 did something special | Run a fair GENERanno 3-class screen/preflight under the same fixed split/protocol before declaring M9-specific novelty |
| Baseline trio had practical tool differences and strong pretrained ceiling | M1 results / `docs/20` | A paper needs intuitive comparison to tools biologists know, not only our internal screen anchor | Produce side-by-side Tiberius/Helixer/ANNEVO/M9/GENERanno predictions on a fixed evaluation panel with common output metrics and qualitative gene-count/FP examples |
| Clean plant split was chosen for contamination avoidance, not final cross-clade claim | M8-M11 notes | A model trained and validated on arabidopsis/rice does not prove fixed-model cross-species utility | Freeze a leave-one-species/clade-out protocol: train/calibrate on one set, evaluate without further tuning on held-out species |
| Publication docs are still skeletal | `docs/12`, `docs/14` | The work has no contribution menu / figure plan / validation matrix despite strong screen results | Switch next stage to a publication-validation evidence map, even if still inside Discovery-Iteration |

### Abandoned route worth reconsidering?

| Route | Original abandon reason | New evidence | Reconsider? | Re-entry criteria check |
|---|---|---|---|---|
| GENERanno LoRA + 3-class head | Not abandoned; parked after metric-negative smoke | User explicitly asks to test whether pretrained models broadly work; native GENERanno specificity remains a key contrast; docs/09 has no abandon entry | Yes, but as fair challenger/control, not as replacement mainline | No docs/09 re-entry gate applies because the route was never abandoned |
| Stronger FP objective for M9 | Deferred because M11 calibration cleared FPR | No new FPR failure after M11 | No immediate | Reopen only if fixed-model full/held-out FPR exceeds `0.01` |
| Multi-class SegmentNT route | M8 refuted gbF1 recovery | No new evidence | No | Would need a new mechanism beyond label enrichment alone |

### Subagent / scout fan-out gaps
- Missing read-only publication audit: map possible paper claims to current evidence and missing figure/table experiments.
- Missing same-panel baseline comparison: Tiberius/Helixer/ANNEVO outputs should be evaluated on the same panel as M9/M11 and future GENERanno, with gene-count, FPR sensitivity, gbF1, locus/exon/gffcompare-style metrics where feasible.
- Missing pretraining-overlap audit for NT-v2 and GENERanno in the actual claim species.
- Missing fixed-model generalization protocol: train/calibrate once, evaluate on unseen species/clade without per-test-species tuning.

### Recommendation (advisory only)

- [x] **run focused ablation / validation package to isolate root cause and recover the publication story.**

Primary recommendation: replace the next single-run `M12-M9L12-FULLSCALE-CALIBRATED` with a small parallel **M12 publication-alignment portfolio**:

1. `M12A-FIXEDMODEL-CROSSSPECIES`: fixed calibrated M9-L12 protocol, train/calibrate on a defined species set, then evaluate on held-out species/clades without further tuning. This directly addresses whether the model is useful as a reusable gene caller.
2. `M12B-SAMEPANEL-BASELINES`: run/evaluate Tiberius, Helixer, ANNEVO, M9-L12, and later GENERanno on the same evaluation panel with common metrics plus tool-native support metrics. This answers the paper-facing "is it more practical?" question.
3. `M12C-GENERANNO-FAIR-CHALLENGER`: give GENERanno a fair 3-class/intron-aware training schedule comparable to M9-L12, not only an 8-window smoke. This tests whether the result is a generic pretrained-model effect or something specific about the M9/NT-v2 recipe.

Secondary recommendation: start a lightweight publication plan / validation matrix update before any new GPU submission so each run maps to a figure/table and reviewer objection.

### Advisory boundary
This retrospective does not cancel the M11 pivot and does not modify `docs/03/04/06/09`. It advises that the next tri-review/council/pivot consume this evidence before submitting further full/scale GPU work.

## Interim Human-Gated Decision: M12-PUBLICATION-PREFLIGHT-TWOSEED (2026-06-17)

### Status
- This started as an interim human-gated decision before seed2 completed. Seed2 has now completed and confirmed the same negative direction; the entry remains interim because no formal `$pivot` has consumed the completed M12 package yet.
- User explicitly approved proceeding from two M12A seeds because seeds 0/1 agree and seed2 is unlikely to change the direction.

### Evidence consumed
- `docs/06_results_log.md` entry `M12-PUBLICATION-PREFLIGHT-TWOSEED`.
- M12A seeds 0/1 metrics: `outputs/M12A-FIXEDMODEL-CROSSSPECIES-A2R-s{0,1}/metrics/metrics.json`.
- M12B same-panel external baselines: `outputs/M12B-SAMEPANEL-BASELINES-{ANNEVO,TIBERIUS,HELIXER}/metrics/metrics.json`.
- M12C GENERanno smoke: `outputs/M12C-GENERANNO-{1P2B-CDS,0P5B-BASE}-SMOKE/metrics/metrics.json`.

### Current evidence summary
M12A fixed Arabidopsis->rice transfer fails in both available seeds: mean `gbF1=0.6608`, `specificity=0.9661`, `FPR=0.0339`, constrained_gbF1@0.01 `0.0`, and `gene_count_ratio=1.832`. Same-panel external tools are much stronger, especially Tiberius (`gbF1=0.9252`, `specificity=0.9927`, `FPR=0.0073`) although it under-calls genes. GENERanno 1.2B CDS-preview has useful signal but fragments badly; GENERanno 0.5B base collapses.

### Interim decision
- **Primary interim decision**: stop treating M9-L12 micro-optimization as the mainline.
- **Next mode**: publication-validation reframing. Define a claim that survives M12: likely low-data/domain adaptation plus calibrated FP control, not a fixed universal cross-species gene caller.
- **Formal closure still needed**: once seed2 finishes or the user requests, run `$tri-review` + `$pivot` over the M12 package and update this entry or append a formal pivot.

### Why not continue current route
The gap to same-panel Tiberius/ANNEVO/Helixer is too large for ordinary tuning: M12A fixed-model `gbF1~0.66` and `FPR~0.034` versus external baselines `gbF1~0.92` and Tiberius `FPR~0.007`. The failure mode is structural generalization/calibration transfer, not a small optimization issue.

### 2026-06-17 closure update
M12A seed2 completed after this interim decision and confirmed the same direction. Three-seed fixed Arabidopsis->rice mean is `gbF1=0.6556`, `specificity=0.9689`, `FPR=0.0311`, constrained_gbF1@0.01 `0.0`, and `gene_count_ratio=1.755`. A follow-up formal tri-review of the user's M13 distance/generalization proposal reached 3/3 quorum for `run-sanity-check-first`: do zero-GPU M12A-vs-M11 failure-mode analysis first, then only if warranted run a bounded single-seed close-plant distance scan; animals are diagnostic/negative controls unless overlap-clean.

---

# Pivot Decision: M13/M14/M16 Combined Generalization Diagnostics (2026-06-18)

## Inputs consumed
- `$result-log`: docs/06 entries for `M13-DISTANCE-GENERALIZATION-SCAN-s0`, `M14-ANIMAL-DISTANCE-NEGCTRL-s0`, `M15-GENERANNO-LORA-PANEL-SCREEN`, and `M16-MULTISPECIES-TRAIN-DISTANCE-SCAN-s0`.
- `$tri-review`: docs/07 `Tri-Review: M13/M14/M16 Combined Generalization Diagnostics (+ M15 GENERanno Context)`, 3/3 quorum.
- Resource profile: screen / publication-alignment diagnostics, NON-CLAIM.

## Current evidence summary
M13 and M14 close the fixed single-species route negatively: Arabidopsis-only train/calibration fails close A. lyrata, rice, and animal negative controls. M16 shows that adding rice to training/calibration partially improves animal FPR/gene-count behavior, but the mixed model still fails the aggregate guardrails (`FPR=0.0197`, gene_count_ratio `1.326`) and remains far from same-panel external gene callers in gbF1. M15 shows GENERanno 1.2B CDS-preview is a real signal source and better than 0.5B base, but the current LoRA+3-class schedule is not guardrail-valid.

## SOTA gap

| Metric | Current | Reference/SOTA-like comparator | Gap (abs) | Gap (rel %) | Severity |
|---|---:|---:|---:|---:|---|
| M16 TEST gbF1 | 0.5615 | M12 same-panel Tiberius `0.9252` | -0.3637 | -39.3% | large / structural |
| M16 TEST FPR | 0.0197 | full/scale hard guardrail `0.0100` | +0.0097 | +96.8% over guardrail | blocker |
| M16 gene_count_ratio | 1.326 | guardrail `1.25` | +0.076 | +6.1% over guardrail | blocker |
| M15 1.2B GENERanno FPR | 0.0258 | guardrail `0.0100` | +0.0158 | +158% over guardrail | blocker |

## Sanity check
- [x] At least two independent CLI reviewers succeeded: 3/3.
- [x] Any reviewer raised comparability blocker: yes, all reviewers converged on it operationally.
- [x] Any reviewer raised leakage/reproducibility blocker: no semantic failure; overlap/claim comparability remains a blocker.
- [x] Metric implementation matches the current project evaluator contract for screen diagnostics.
- [x] Loss showed expected pattern for M13/M14/M16.
- [x] Seed variance is incomplete for M13/M14/M16, but screen directional failures are large enough for non-claim pivot.

## Tri-review summary (record ALL reviewer conclusions — drop none)

| Reviewer | Judgment | Next action proposed | Main concern | Confidence |
|---|---|---|---|---|
| A · Claude | `return-to-literature` | Run same-panel Tiberius/Helixer/ANNEVO on the M16-style panel; audit NT-v2 overlap before more architecture work | Cannot distinguish task difficulty from architecture failure without same-panel baseline/comparability; current fixed-model route structurally fails | Medium |
| B · Codex | `comparability-blocker` | Freeze same-evaluator/same-panel Tiberius/Helixer/ANNEVO benchmark and pretraining-overlap audit, then choose broader/adaptive NT-v2 or GENERanno redesign | Further M9 gains may remain internal screen metrics; claim benchmark is draft | High |
| C · Antigravity | `comparability-blocker` | Execute baseline/comparability run and freeze `ACTIVE_GOAL.sota_benchmark`; then return to decoder/objective + broader/adaptive panel | Screen-only diagnostics and draft benchmark make architecture redesign premature | High |

Consensus: no generic tuning, no Track-B promotion, no fixed single-species generalization claim. Close comparability/SOTA-freeze blocker first.
Disagreement: label only (`return-to-literature` vs `comparability-blocker`); operational next action is the same.
Quorum / degraded review status: 3/3.

## Reviewer-proposed directions (ordered A→B→C, verbatim before convergence)

| # | From reviewer | Direction | major_axis | mechanism_delta | Orthogonal to others? | Into this round's cohort? |
|---:|---|---|---|---|---|---|
| 1 | A · Claude | Same-panel Tiberius/Helixer/ANNEVO on M16 train/test panel | comparability | freeze external baseline behavior under identical evaluator before architecture claims | yes | yes |
| 2 | A · Claude | Broader phylodiverse training + species-conditioned/per-clade head | data_view/head_arch | make head/clade adapters robust to species distribution shift | yes | design only after blocker |
| 3 | B · Codex | Freeze ANNEVO-compatible same-panel benchmark + overlap audit | comparability | make SOTA target and pretraining caveats claim-grade | overlaps #1 | yes |
| 4 | B · Codex | GENERanno 1.2B specificity-preserving continuity objective | backbone/objective | preserve CDS-preview specificity while learning intron/gene-body coherence | yes | later challenger |
| 5 | C · Antigravity | Segment/semi-CRF decoder + FP/species-balanced objective | decoder/objective | enforce gene-structure coherence and cross-domain FP control | yes | later architecture |

## Is tuning justified?
- ❌ No. The gap to external gene callers and the cross-species FPR/gene-count failures are too large and too systematic for lr/dropout/batch-size tuning.

## Architecture hypothesis status
weakened. NT-v2 L12 + 3-class calibrated head remains useful internally and species-diverse training helps, but the fixed universal model hypothesis is not supported by current evidence.

## DECISION (choose exactly one)

- [ ] Continue current architecture as-is
- [ ] Tune current architecture
- [ ] Scale data / training (Track B with current architecture)
- [ ] Replace component
- [ ] Change backbone
- [ ] Change objective / loss
- [x] Comparability audit first
- [ ] Sanity check first
- [ ] Abandon this route
- [ ] Return to literature

## Why this decision (not another)
The main blocker is now paper-facing, not runtime or metric extraction. M13/M14/M16 show that the current fixed-model route is not claim-ready, but they do not yet tell us whether Tiberius/Helixer/ANNEVO also degrade on the same A. lyrata/animal diagnostic panels under the same evaluator. Without that, changing backbone/objective would be underdetermined. `return-to-literature` is too broad because we already have local baseline runners and a concrete audit to run; `comparability audit first` is the smallest decisive action.

## Best next architecture moves (if applicable)

| Priority | Move | Expected mechanism | Goes to which EXP / Track |
|---:|---|---|---|
| 1 | Same-panel external baseline/comparability audit | Establish whether cross-species failure is model-specific or task/panel-specific | `M17-SAMEPANEL-GENERALIZATION-BASELINES` |
| 2 | Pretraining-overlap audit for NT-v2/GENERanno/baselines | Prevent contaminated animal/plant diagnostics from steering claim decisions | `M17-PRETRAINING-OVERLAP-AUDIT` |
| 3 | Broader/adaptive NT-v2 multi-species route | Test whether multi-species positive signal in M16 scales with clade diversity/adapters | post-M17 candidate |
| 4 | GENERanno 1.2B specificity-preserving schedule | Test whether CDS-preview specificity can be retained while adding intron/coherence | post-M17 challenger |

## Parallel cohort this round (primary + orthogonal directions)
- **Primary direction (single)**: comparability audit first.
- **Parallel cohort**: run baseline panel and overlap audit in parallel; no shared training-code conflict expected.

| Slot | EXP ID (new) | Direction | major_axis | mechanism_delta | Track | Resource profile |
|---|---|---|---|---|---|---|
| primary | `M17-SAMEPANEL-GENERALIZATION-BASELINES` | Tiberius/Helixer/ANNEVO on M13/M16 diagnostic species where runnable | comparability | same evaluator / same species panel / released weights | baseline | screen |
| parallel-1 | `M17-PRETRAINING-OVERLAP-AUDIT` | NT-v2, GENERanno, and baseline overlap/provenance audit for A. lyrata/rice/gallus/drosophila | comparability | claim caveat ledger and clean/diagnostic species status | local audit | local |

Shared-code conflict? no.

## TODO update
- [x] update docs/07_tri_review.md.
- [x] update docs/08_pivot_decisions.md.
- [ ] update docs/05_todo.md and docs/11_master_plan.md with M17 cohort.
- [ ] run `$note-gate` to persist the combined review/pivot evidence.
- [ ] implement/submit `M17-SAMEPANEL-GENERALIZATION-BASELINES` if existing M12B baseline runners can be reused safely.

---

# Pivot Decision: M17+M18 Combined Evidence (2026-06-19)

## Inputs consumed
- `$result-log`: docs/06 entries for `M17-SAMEPANEL-GENERALIZATION-BASELINES`, `M18-MULTICLADE-TRAIN-DIAGNOSTIC-s0`, `M18-GENERANNO-0P5B-SPEC-OBJECTIVE-s0`, and `M18-GENERANNO-1P2B-SPEC-OBJECTIVE-s0`.
- `$tri-review`: docs/07 `Tri-Review: M17+M18 Combined Evidence`, 2/3 `DEGRADED_REVIEW`.
- Resource profile: screen / publication-alignment diagnostics, NON-CLAIM.

## Current evidence summary
M17 shows released callers remain much stronger than our broad fixed NT-v2 diagnostics on gbF1, though each released caller has practical tradeoffs. M18 broad NT-v2 confirms that adding Drosophila to the supervised panel helps nearby plant transfer but does not produce a broad fixed caller; gallus remains a severe emission/coherence failure and test-label oracle calibration cannot rescue it. M18 0.5B GENERanno base is a clean negative under the same stronger objective. M18 1.2B GENERanno CDS-preview is the route-changing positive: aggregate clean-plant FPR `0.0071`, specificity `0.9929`, gbF1/constrained gbF1@0.01 `0.8494`, macro specificity `0.9943`, gene_count_ratio `0.864`.

## SOTA gap

| Metric | Current | Reference/SOTA-like comparator | Gap (abs) | Gap (rel %) | Severity |
|---|---:|---:|---:|---:|---|
| M18 1.2B gbF1 clean plants | 0.8494 | released clean-plant baselines `0.922-0.927` | -0.073 to -0.078 | ~-8% | large enough to require structure/calibration, not generic tuning |
| M18 1.2B FPR clean plants | 0.0071 | full/scale hard guardrail `0.0100` | +0.0029 margin | 29% below guardrail | strong |
| M18 1.2B gene_count_ratio | 0.864 | guardrail `<=1.25` | pass | n/a | coherent but Arabidopsis under-call risk |
| M18 NT-v2 broad gbF1 | 0.6170 | M17 released callers `0.879-0.912` | -0.262 to -0.295 | ~-30% | broad fixed route not competitive |

## Sanity check
- [x] At least two independent CLI reviewers succeeded: 2/3 `DEGRADED_REVIEW`.
- [x] Any reviewer raised comparability blocker: yes, both reviewers say claim is blocked by benchmark/provenance/screen status.
- [x] Any reviewer raised leakage/reproducibility blocker: no semantic failure; GENERanno provenance/overlap is a claim blocker.
- [x] Metric implementation matches the current project evaluator contract for screen diagnostics.
- [x] Loss showed expected pattern for M18 1.2B (`0.8206 -> 0.6194`).
- [x] Seed variance is not known for M18 1.2B, acceptable for screen but not for claim.

## Tri-review summary (record ALL reviewer conclusions — drop none)

| Reviewer | Judgment | Next action proposed | Main concern | Confidence |
|---|---|---|---|---|
| A · Claude | `change-backbone` | Primary `M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS`: 2 seeds, raw scores, VAL-only calibration; parallel provenance audit; NT-v2 adaptive design only as backup | M18 1.2B is route-changing but still below released clean-plant gbF1 and blocked by single seed/provenance/draft benchmark | Medium-High |
| B · Codex | `scale-to-track-b` | `M19-GENERANNO-1P2B-RAWCAL-STRUCT-PREFLIGHT`: raw-score calibration, plain vs constrained/segment-aware decode, FPR sensitivity | gbF1 gap remains >0.05, so generic tuning is forbidden; need calibration/structure and claim blockers resolved | High |
| C · Antigravity | failed | none | agy OAuth login required; no valid independent review produced | n/a |

Consensus: switch the primary next work from broad fixed NT-v2 to GENERanno 1.2B CDS-preview preflight; no claim; no generic tuning.
Disagreement: label only (`change-backbone` vs `scale-to-track-b`). Operational next step converges.
Quorum / degraded review status: 2/3 `DEGRADED_REVIEW`.

## Reviewer-proposed directions (ordered A→B→C, verbatim before convergence)

| # | From reviewer | Direction | major_axis | mechanism_delta | Orthogonal to others? | Into this round's cohort? |
|---:|---|---|---|---|---|---|
| 1 | A · Claude | `M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS` | backbone/calibration | GENERanno 1.2B CDS-preview + stronger FP objective + raw-score VAL-only calibration + 2 seeds | primary | yes |
| 2 | A · Claude | `M19-GENERANNO-PROVENANCE-AUDIT` | comparability | audit GENERanno pretraining/species overlap and model provenance | yes | yes |
| 3 | A · Claude | `M19-NTV2-ADAPTIVE-DESIGN` | architecture design | clade-aware/adaptive NT-v2 fallback design, no GPU until needed | yes | yes, local only |
| 4 | B · Codex | `M19-GENERANNO-1P2B-RAWCAL-STRUCT-PREFLIGHT` | calibration/decoder | save raw scores, run VAL-only calibration, compare plain vs constrained/segment-aware decode | overlaps #1 | folded into primary after raw-score availability |
| 5 | B · Codex | Freeze clean-plant / published-SOTA comparability contract | comparability | prepare next Track-B claim contract and baseline ruler | yes | yes, folded into audit/docs |
| 6 | C · Antigravity | none | n/a | reviewer failed due OAuth | n/a | no |

## Is tuning justified?
- ❌ Generic tuning is not justified. The gbF1 gap to released clean-plant baselines remains around `0.07-0.08`, and broad fixed NT-v2 failures are structural. Validation-only calibration is allowed because it is a protocol/operating-point check, not LR/dropout tuning.

## Architecture hypothesis status
partially supported but redirected. The "foundation pretrained signals + FP-aware objective + structured decoder/calibration" hypothesis is supported, but the strongest foundation signal is now GENERanno 1.2B CDS-preview rather than NT-v2 broad fixed unfreeze. Fixed universal NT-v2 is weakened.

## DECISION (choose exactly one)

- [ ] Continue current architecture as-is
- [ ] Tune current architecture
- [ ] Scale data / training as Track B claim run
- [ ] Replace component
- [x] Change backbone: primary next challenger becomes GENERanno 1.2B CDS-preview + our 3-class/FP-aware route
- [ ] Change objective / loss
- [ ] Comparability audit first
- [ ] Sanity check first
- [ ] Abandon this route
- [ ] Return to literature

## Why this decision (not another)
The evidence no longer supports broad fixed NT-v2 as the main publication route: M18 broad training still fails gallus badly, and oracle calibration rules out a simple threshold rescue. By contrast, M18 1.2B GENERanno is the first pretrained-CDS route to pass the FPR/gene-count screen guardrails while maintaining strong gbF1. This is a backbone-level change, not a small parameter adjustment. It is also not a full Track-B claim scale-up yet because the result is single-seed, screen-only, lower-gbF1 than released clean-plant baselines, and blocked by GENERanno provenance plus draft SOTA benchmark.

## Best next architecture moves

| Priority | Move | Expected mechanism | Goes to which EXP / Track |
|---:|---|---|---|
| 1 | Re-run GENERanno 1.2B with raw-score saving and 2 seeds | confirm M18 stability and enable no-leak VAL-only operating-point calibration | `M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS` |
| 2 | VAL-only calibration over raw scores | test whether gbF1/FPR sensitivity improves without test leakage | same primary EXP post-processing |
| 3 | GENERanno provenance/overlap audit | clear or label claim blocker before further paper framing | `M19-GENERANNO-PROVENANCE-AUDIT` |
| 4 | Structured decoder/segment-aware preflight after raw scores exist | improve rice gbF1/coherence without losing specificity | follow-up or folded diagnostic if cheap |

## Parallel cohort this round (primary + orthogonal directions)
- **Primary direction (single)**: change backbone to GENERanno 1.2B CDS-preview route; run raw-score calibration preflight with 2 seeds.
- **Parallel cohort**: one local provenance/comparability audit and one local NT-v2 adaptive fallback design. No shared-code conflict expected.

| Slot | EXP ID (new) | Direction | major_axis | mechanism_delta | Track | Resource profile |
|---|---|---|---|---|---|---|
| primary | `M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS` | 1.2B CDS-preview LoRA, stronger FP objective, save raw scores, 2 seeds, VAL-only calibration | backbone/calibration | promote route-changing M18 signal into reproducible calibrated preflight | screen / Track-B preflight | GPU screen array |
| parallel-1 | `M19-GENERANNO-PROVENANCE-AUDIT` | audit GENERanno training data/species overlap and released model provenance | comparability | decide whether clean-plant claim is possible or only diagnostic | local | CPU/local |
| parallel-2 | `M19-NTV2-ADAPTIVE-DESIGN` | write fallback design for clade-aware/adaptive NT-v2, no GPU yet | architecture_design | preserve an NT-v2 route only if GENERanno hits provenance/ceiling blocker | local | CPU/local |

Shared-code conflict? no. Primary uses existing trainer with recently added `--save-raw-scores`; local audits/design docs do not modify shared training code.

## TODO update
- [ ] update docs/05_todo.md with M19 cohort.
- [ ] update docs/11_master_plan.md current step after M19 submission.
- [ ] run `$note-gate` / evidence register entry for this tri-review+pivot.
- [ ] create M19 config/sbatch/code-review gate and submit 2-seed GPU array.

# Pivot Decision: M19-GENERANNO-COMBINED-DECISION (2026-06-21)

## Inputs consumed
- `$result-log`: docs/06 `Result: M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s{0,1}`.
- `$tri-review`: docs/07 `Tri-Review: M19-GENERANNO-COMBINED-DECISION`, 3/3 quorum.
- Same-evaluator report: `reports/M19-COMPARABILITY-EVIDENCE/comparison_tables.md`.
- Provenance dossier: `refs/dossiers/m19_generanno_provenance_audit.md`.
- Resource profile: screen / Track-B-preflight, NON-CLAIM.

## Current evidence summary
M19 confirms that the M18 GENERanno 1.2B CDS-preview signal is stable rather than a seed accident. Both M19 seeds are aggregate FPR-valid and gene-count sane after validation-only calibration: s0 `gbF1=0.8421`, FPR `0.0083`, gene_count_ratio `1.083`; s1 `gbF1=0.8815`, FPR `0.0065`, gene_count_ratio `0.830`. The 1.2B CDS-preview route remains clearly stronger than the 0.5B base control (`gbF1=0.6561`, FPR `0.0967`, gene_count_ratio `1.617`).

The same-evaluator table also prevents overclaiming. Released clean-plant callers still have higher gbF1: Tiberius `0.9252`, ANNEVO `0.9269`, Helixer `0.9220`. M19's real practical advantage is the specificity/gene-count balance: it passes aggregate FPR<=0.01 with a more reasonable gene_count_ratio than Tiberius (`0.628`, under-call), ANNEVO (`0.726`, under-call), and Helixer (`FPR=0.0216`). However, GENERanno Arabidopsis/rice provenance is still `overlap_unknown`, so this panel cannot support a clean held-out SOTA claim.

## SOTA / utility gap

| Metric | Current best M19 | Comparator | Gap / margin | Severity |
|---|---:|---:|---:|---|
| gbF1 | 0.8815 | ANNEVO 0.9269 / Tiberius 0.9252 / Helixer 0.9220 | -0.0405 to -0.0454 | large enough to require structural improvement, not generic tuning |
| FPR | 0.0065 | hard guardrail 0.0100 / Tiberius 0.0073 | +0.0035 guardrail margin; slightly better than Tiberius | strong |
| gene_count_ratio | 0.830 | guardrail <=1.25 / Tiberius 0.628 | pass; less under-called than Tiberius | practical advantage |
| provenance | `overlap_unknown` | claim requires clean/defensible split | hard blocker | cannot claim |

## Sanity check
- [x] At least two independent CLI reviewers succeeded: 3/3 quorum.
- [x] Any reviewer raised comparability blocker: yes, all three.
- [x] Any reviewer raised leakage/provenance blocker: yes, all three; GENERanno Arabidopsis/rice overlap remains unknown.
- [x] Metric implementation matches the current same-evaluator screen contract.
- [x] Training semantics are valid: Slurm completed, loss decreased, raw scores exist, metrics finite.
- [x] Seed stability is adequate for screen route promotion, not for claim.

## Tri-review summary (record ALL reviewer conclusions — drop none)

| Reviewer | Judgment | Next action proposed | Main concern | Confidence |
|---|---|---|---|---|
| A · Claude | `claim-panel-first` + `replace-component` parallel | `M20-CLAIM-PANEL-FREEZE` locally, in parallel with `M20-STRUCTURED-DECODER-IMPL` smoke/design; then evaluate on the claim-clean panel | M19 is stable and has practical utility, but provenance is a hard claim blocker and gbF1 still trails released callers | High |
| B · Codex | `claim-panel-first` | `M20-CLAIM-CLEAN-PANEL-FREEZE`: freeze a claim-clean held-out panel and same-panel external comparability before scaling GENERanno | M19 clears stability gates but not claim gates; Arabidopsis/rice overlap unknown makes scale GPU premature | High |
| C · Antigravity | `freeze-as-adaptation-evidence` | `M20-CLEAN-STRUCTURED-DECODER`: freeze GENERanno as adaptation evidence unless clean panel is proven; move mainline GPU toward clean structured decoder/error analysis | Provenance is fatal for clean claim; gbF1 gap requires segment-level structured modeling, not more generic tuning | High |

Consensus: M19 is stable, useful adaptation/comparability evidence; current Arabidopsis/rice M19 must not be claimed or scaled as a claim run.
Disagreement: operational ordering. Two reviewers put the primary label on claim-clean panel freeze; one reviewer would immediately freeze GENERanno as adaptation evidence and move the GPU mainline to a clean structured-decoder route.
Quorum status: 3/3.

## Reviewer-proposed directions (ordered A→B→C, verbatim before convergence)

| # | From reviewer | Direction | major_axis | mechanism_delta | Orthogonal to others? | Into this round's cohort? |
|---:|---|---|---|---|---|---|
| 1 | A · Claude | `M20-CLAIM-PANEL-FREEZE` | comparability/data contract | freeze claim-clean held-out species panel before claim-grade GPU | primary blocker | yes, primary |
| 2 | A · Claude | `M20-STRUCTURED-DECODER-IMPL` | decoder/head_arch | segment/semi-CRF or constrained Viterbi-style decoder smoke/design on top of strong low-FPR emissions | yes | yes, parallel implementation/smoke only |
| 3 | A · Claude | GENERanno provenance escalation | comparability/provenance | contact/search for training/fine-tuning manifest or exclusion evidence | yes | yes, folded into primary |
| 4 | B · Codex | `M20-CLAIM-CLEAN-PANEL-FREEZE` | comparability/data contract | freeze a claim-clean held-out panel and re-freeze same-panel external comparability | overlaps #1 | yes, primary |
| 5 | B · Codex | `M20-STRUCTURED-DECODER-DESIGN` | decoder/head_arch | design segment/structured decoder, IO, calibration and guardrails before large GPU | overlaps #2 | yes, parallel |
| 6 | B · Codex | `M20-PROVENANCE-AUDIT-ESCALATION` | provenance | systematic manifest search / author-contact template / alternative panel list | overlaps #3 | yes, folded into primary |
| 7 | C · Antigravity | `M20-CLEAN-STRUCTURED-DECODER` | decoder/head_arch | shift mainline to clean architecture/structured decoder route, inheriting FP-aware and calibration lessons | yes but should wait for panel for claim | yes, parallel implementation; GPU only after panel |
| 8 | C · Antigravity | `CLAIM-PANEL-FIRST` | data contract | data mining for 1-2 provably uncontaminated held-out species | overlaps #1 | yes, primary |
| 9 | C · Antigravity | `SOTA-ERROR-ANALYSIS` | evaluation/error_analysis | analyze where Tiberius/ANNEVO gain gbF1 and what decoder should target | yes | yes, local parallel |

## Is tuning justified?
- ❌ Generic tuning is not justified. The best M19 gbF1 still trails released callers by `0.040-0.045`, and the claim blocker is logical/provenance-based. LR/dropout/LoRA rank sweeps would not make current Arabidopsis/rice evidence claimable.
- ✅ Structural development is justified: a segment/structured decoder is a mechanism-level move that targets gene coherence/gbF1 while preserving the low-FPR emissions M19 established.
- ✅ Local comparability/provenance work is mandatory: without a claim-clean panel or manifest, improved metrics remain adaptation evidence only.

## Architecture hypothesis status
Supported but incomplete. GENERanno 1.2B CDS-preview plus our 3-class FP-aware adaptation supplies strong low-FPR emissions, validating the "pretrained backbone + FP-aware objective" half of the hypothesis. The remaining gap to released callers is gene-level sensitivity/coherence, which supports the planned "structured decoder / segment-level objective" half. The route is claim-blocked by provenance, not semantically failed.

## DECISION (choose exactly one)

- [ ] Continue current architecture as-is
- [ ] Tune current architecture
- [ ] Scale data / training as Track B claim run
- [ ] Replace component
- [ ] Change backbone
- [ ] Change objective / loss
- [x] Comparability audit first: freeze claim-clean panel/provenance before claim-grade GPU
- [ ] Sanity check first
- [ ] Abandon this route
- [ ] Return to literature

## Why this decision (not another)
M19 has done its job: it proves the 1.2B GENERanno adaptation route is stable and useful, and it gives a strong specificity/gene-count story against released callers. But it does not give a defensible paper claim. The clean-plant panel is `overlap_unknown` for GENERanno, and the best M19 gbF1 is still below the released-caller frontier. Therefore a claim-grade GPU scale-up on the current panel would be scientifically unspendable even if metrics improved.

The next primary action must be a claim/comparability gate, not another performance chase. In parallel, it is reasonable to implement and smoke a structured decoder because that is the only plausible mechanism-level route for closing the gbF1 gap while preserving M19's specificity advantage. That implementation should not become a claim-scale GPU run until the claim-clean panel/provenance decision is closed.

## Best next moves

| Priority | Move | Expected mechanism | Goes to which EXP / Track |
|---:|---|---|---|
| 1 | Freeze claim-clean held-out panel and provenance status | convert M19 from screen evidence into a defensible claim path or explicitly downgrade it to adaptation evidence | `M20-CLAIM-CLEAN-PANEL-FREEZE` / local |
| 2 | Rebuild same-evaluator external-model comparison contract for the frozen panel | prevent another internal metric win that cannot be compared to Tiberius/ANNEVO/Helixer | same EXP or `M20-SAMEPANEL-BASELINE-CONTRACT` / local+software |
| 3 | Implement structured/segment decoder smoke | target gbF1/gene coherence without generic tuning; carry M19's FP-aware objective/calibration forward | `M20-STRUCTURED-DECODER-IMPL` / implement+smoke |
| 4 | Run SOTA error analysis | identify exactly where released callers gain gbF1 and where under/over-call tradeoffs occur | `M20-SOTA-ERROR-ANALYSIS` / local |

## Parallel cohort this round (primary + orthogonal directions)
- **Primary direction (single)**: comparability/provenance first, via claim-clean panel freeze. This is the claim gate.
- **Parallel cohort**: structured-decoder implementation/smoke and SOTA error analysis can proceed concurrently because they do not require the same claim panel to be final. No claim-grade GPU scale until the primary gate closes.

| Slot | EXP ID (new) | Direction | major_axis | mechanism_delta | Track | Resource profile |
|---|---|---|---|---|---|---|
| primary | `M20-CLAIM-CLEAN-PANEL-FREEZE` | freeze claim-clean species/accession panel and provenance status for GENERanno/current baselines | comparability/data_contract | convert or block M19 route as claim evidence before more scale GPU | claim preflight | local / CPU |
| parallel-1 | `M20-STRUCTURED-DECODER-IMPL` | implement and smoke segment/structured decoder route, preserving M19 FP-aware objective and calibration discipline | decoder/head_arch | structural gbF1/coherence improvement rather than generic tuning | screen prep | local + bounded GPU smoke only after code-review |
| parallel-2 | `M20-SOTA-ERROR-ANALYSIS` | analyze Tiberius/ANNEVO/Helixer vs M19 errors on existing clean-plant outputs | evaluation/error_analysis | define what the structured decoder must fix and how to report utility | local analysis | CPU/local |

Shared-code conflict? low. Primary is data/provenance/report work; parallel-2 is read-only analysis; parallel-1 may touch trainer/decoder code and must pass `$implement` + `$code-review-gate` before any smoke/full job.

## TODO update
- [x] update docs/05_todo.md with M20 cohort.
- [x] update docs/11_master_plan.md current step to M20 claim-clean panel + structured-decoder prep.
- [x] run `$note-gate` / evidence register entry for M19 tri-review+pivot.
- [x] do not submit any claim-grade M19/GENERanno GPU until `M20-CLAIM-CLEAN-PANEL-FREEZE` closes.

---

# Pivot Decision: M20-GENERANNO-COMBINED-DECISION (2026-06-21)

## Inputs consumed
- `$tri-review`: docs/07 `Tri-Review: M20-GENERANNO-COMBINED-DECISION`, 3/3 quorum.
- Result/evidence logs: docs/06 `Result: M20-STRUCTURED-DECODER-IMPL-SMOKE3`, `refs/dossiers/M20-CLAIM-CLEAN-PANEL-FREEZE.md`, `reports/M20-SOTA-ERROR-ANALYSIS/report.md`.
- Resource profile: local + smoke/screen prep, NON-CLAIM.

## Current evidence summary
M20 resolves the immediate M19 blockers. The provenance gate is negative for clean held-out GENERanno claims on the current Arabidopsis/rice panel: public sources do not provide a complete species/accession exclusion manifest. The same-evaluator table shows M19 GENERanno s1 is strong on FPR (`0.0065`) and gene-count behavior (`0.830`) but trails the hard-FPR released comparator Tiberius in gbF1 (`0.8815` vs `0.9252`). The remaining technical gap is recall/gene recovery; the CRF path is smoke-proven but not screen-proven.

## SOTA / utility gap

| Metric | Current best M19 | Comparator | Gap / margin | Severity |
|---|---:|---:|---:|---|
| gbF1 | 0.8815 | Tiberius 0.9252 / ANNEVO 0.9269 | -0.0437 to -0.0454 | requires structural improvement |
| FPR | 0.0065 | hard guardrail 0.0100 / Tiberius 0.0073 | passes; slightly better than Tiberius | strong |
| gene_count_ratio | 0.830 | Tiberius 0.628 / guardrail <=1.25 | less under-called than Tiberius | practical advantage |
| provenance | `overlap_unknown` | clean claim requires exclusion manifest or controlled provenance | hard blocker | cannot claim |

## Sanity check
- [x] At least two independent CLI reviewers succeeded: 3/3.
- [x] Any reviewer raised comparability blocker: yes, provenance remains a claim blocker.
- [x] Any reviewer raised leakage/reproducibility blocker: no run leakage signal; yes unresolved pretraining overlap.
- [x] Metric implementation matches same-evaluator CDS-span contract.
- [x] CRF smoke semantic success passed; metrics are finite.
- [x] Seed variance is known for non-CRF M19; unknown for CRF because only smoke exists.

## Tri-review summary (record ALL reviewer conclusions — drop none)

| Reviewer | Judgment | Next action proposed | Main concern | Confidence |
|---|---|---|---|---|
| A · Claude | `continue-current-route` | `M21-GENERANNO-1P2B-CRF-SCREEN` non-claim, then provenance/clean-backbone branch depending on result | provenance blocks claim; CRF may or may not improve recall without FPR cost | Medium |
| B · Codex | `replace-component` | replace/strengthen decoder/head via real CRF screen; compare against M19 s1 and Tiberius | CRF smoke proves only code path; gbF1 gap and provenance remain | Medium |
| C · Antigravity | `comparability-blocker` | full Track A CRF screen as mechanism test, followed by clean-provenance backbone if successful | GENERanno provenance is fatal for clean claim; FPR may break under CRF | High |

Consensus: no claim/full-scale GENERanno on current panel; run a real non-claim CRF screen as the next mechanism experiment; keep provenance as a separate claim prerequisite.
Disagreement: reviewer labels differ, but the concrete next action converges.
Quorum / degraded review status: 3/3.

## Reviewer-proposed directions (ordered A→B→C, verbatim before convergence)

| # | From reviewer | Direction | major_axis | mechanism_delta | Orthogonal to others? | Into this round's cohort? |
|---:|---|---|---|---|---|---|
| 1 | A · Claude | `M21-GENERANNO-1P2B-CRF-SCREEN` | decoder/head_arch | real CRF screen on top of M19 low-FPR emissions | primary | yes |
| 2 | A · Claude | transition regularization / FP-aware CRF loss | decoder/objective | constrain transition matrix and encode FP pressure into structured loss | extension | later if CRF screen promising |
| 3 | A · Claude | controlled-provenance backbone fallback | backbone/provenance | move effective decoder recipe to clean backbone if GENERanno claim stays blocked | orthogonal | scout only |
| 4 | B · Codex | CRF decoder Track A screen | decoder/head_arch | replace per-token decode with CRF Viterbi screen | overlaps #1 | yes |
| 5 | B · Codex | semi-CRF / segment-level decoder | decoder | model CDS segment length/boundaries directly | extension | later |
| 6 | B · Codex | recall-aware FP-constrained objective | objective | add recall/gene-count recovery without breaking FPR | extension | later |
| 7 | C · Antigravity | full Track A CRF mechanism screen | decoder/head_arch | quantify gbF1 lift and FPR cost of CRF | overlaps #1 | yes |
| 8 | C · Antigravity | clean-provenance backbone search/training | backbone/provenance | solve final claim route after mechanism validation | orthogonal | scout only |

## Is tuning justified?
- ❌ Generic tuning is not justified. The gap is recall/gene recovery and provenance, not LR/dropout noise.
- ✅ Structural component replacement is justified: CRF decoder directly targets the identified error mode.
- ❌ Track B/full claim is forbidden until provenance and benchmark contract are resolved.

## Architecture hypothesis status
Supported but incomplete. GENERanno 1.2B plus FP-aware 3-class adaptation supplies strong low-FPR emissions. The missing piece is structured gene recovery. CRF is the next minimal component replacement to test that hypothesis; final claim still needs provenance-clean data/model evidence.

## DECISION (choose exactly one)

- [ ] Continue current architecture as-is
- [ ] Tune current architecture
- [ ] Scale data / training as Track B claim run
- [x] Replace component: decoder/head via real CRF screen
- [ ] Change backbone
- [ ] Change objective / loss
- [ ] Comparability audit first
- [ ] Sanity check first
- [ ] Abandon this route
- [ ] Return to literature

## Why this decision (not another)
Choosing `comparability audit first` again would loop: M20 already performed the audit and found the blocker real. Choosing scale/claim would be invalid because GENERanno remains `overlap_unknown`. Choosing generic tuning would ignore the error analysis: the gap is recall/gene recovery under hard FPR. The only actionable non-claim step with a direct mechanism is to replace the local decoder/head behavior with CRF and measure whether it closes gbF1 without sacrificing FPR.

## Best next moves

| Priority | Move | Expected mechanism | Goes to which EXP / Track |
|---:|---|---|---|
| 1 | Run real CRF decoder screen | improve recall/gene recovery from strong low-FPR GENERanno emissions | `M21-GENERANNO-1P2B-CRF-SCREEN` / screen NON-CLAIM |
| 2 | Compare against M19 s1 and fixed baselines | quantify gbF1 lift, FPR cost, and gene_count behavior | same M21 result-log |
| 3 | Scout clean-provenance backbone/panel options | prevent successful CRF mechanism from becoming unclaimable | local provenance/backbone scout |
| 4 | If CRF fails, escalate to semi-CRF/objective or clean-backbone route | avoid repeated local tuning | next pivot |

## Parallel cohort this round (primary + orthogonal directions)
- **Primary direction (single)**: `M21-GENERANNO-1P2B-CRF-SCREEN`, non-claim screen.
- **Parallel cohort**: local provenance/clean-backbone scout can proceed while M21 runs; no full/scale claim GPU.

| Slot | EXP ID (new) | Direction | major_axis | mechanism_delta | Track | Resource profile |
|---|---|---|---|---|---|---|
| primary | `M21-GENERANNO-1P2B-CRF-SCREEN` | real CRF screen on Arabidopsis+rice, compare to M19 s1/Tiberius | decoder/head_arch | replace independent argmax with learned transition + Viterbi | screen | GPU, NON-CLAIM |
| parallel-1 | `M21-CLEAN-BACKBONE-PROVENANCE-SCOUT` | local scout for controlled-provenance backbone/panel options | provenance/backbone | identify claimable route if CRF mechanism works | local | CPU/local |

Shared-code conflict? no for local scout; yes if later adding semi-CRF/objective changes. M21 screen touches trainer/config/sbatch and needs fresh `$code-review-gate`.

## TODO update
- [ ] update docs/05_todo.md with M21 tasks.
- [ ] update docs/11_master_plan.md current step to M21 CRF screen + provenance scout.
- [ ] update docs/15 evidence register.
- [ ] implement/code-review/smart-sbatch M21 before any GPU submission.

---

# Pivot Decision: M21-GENERANNO-1P2B-CRF-SCREEN (2026-06-22)

## Inputs consumed
- `$result-log`: docs/06 `Result: M21-GENERANNO-1P2B-CRF-SCREEN-s{0,1}`.
- `$tri-review`: docs/07 `Tri-Review: M21-GENERANNO-1P2B-CRF-SCREEN`, 3/3 quorum.
- Resource profile: screen / NON-CLAIM.

## Current evidence summary
M21 was a clean model-quality test of the M20 CRF hypothesis after separating runtime failures from valid seeds. Seed0 completed with gbF1 `0.8544`, FPR `0.0273`, and gene_count_ratio `0.956`. Seed1 rescue completed with gbF1 `0.8744`, FPR `0.0192`, and gene_count_ratio `0.690`. The original shared seed1 timeout is not model evidence, and the duplicate fast-validation rescue was cancelled after seed1 rescue produced metrics.

The best M21 CRF seed is worse than M19 non-CRF seed1 on both axes that matter for the paper story: gbF1 `0.8744 < 0.8815`, and FPR `0.0192 > 0.0065`. CRF therefore destroys the low-FPR property that made M19 useful while failing to close the released-caller gbF1 gap.

## SOTA gap
| Metric | Current best M21 CRF | M19 non-CRF best | Released comparator | Gap / margin | Severity |
|---|---:|---:|---:|---:|---|
| gbF1 | 0.8744 | 0.8815 | ANNEVO 0.9269 | -0.0525 vs ANNEVO; -0.0071 vs M19 | High |
| intergenic_FPR | 0.0192 | 0.0065 | Tiberius 0.0073 | +0.0119 vs Tiberius; +0.0128 vs M19 | High |
| gene_count_ratio | 0.690 | 0.830 | Helixer 0.820 / Tiberius 0.628 | under-calls more than M19, less than Tiberius | Medium |

## Sanity check
- [x] At least two independent CLI reviewers succeeded: 3/3 succeeded.
- [x] Any reviewer raised comparability blocker: yes, GENERanno provenance remains `overlap_unknown`.
- [x] Any reviewer raised leakage/reproducibility blocker: no new leakage; yes unresolved pretraining overlap.
- [x] Metric implementation matches same-evaluator CDS-span contract.
- [x] Loss/validation trace is plausible; no NaN/OOM in valid seeds.
- [x] Seed variance reasonable enough to decide: both valid seeds fail hard FPR<=0.01 and neither beats M19 best.

## Tri-review summary (record ALL reviewer conclusions — drop none)
| Reviewer | Judgment | Next action proposed | Main concern | Confidence |
|---|---|---|---|---|
| A · Claude | `replace-component` | abandon CRF as current decoder component; return to non-CRF low-FPR route and improve recall through emissions/foundation features or objective | CRF worsened FPR by 3-4x and did not improve gbF1 | High |
| B · Codex | `replace-component` | drop CRF as primary; use FPR-controlled M19-like path and only revisit structure with lighter inference-time constraints | decoder/component failure plus unresolved provenance, not tuning | Medium-High |
| C · Antigravity | `abandon-route` | abandon CRF decoder direction; if continuing GENERanno, change objective/loss on non-CRF head and move lessons to clean-provenance backbone | CRF breaks FPR red line and cannot claim under unresolved provenance | High |

Consensus: M21 is semantically successful but refutes the GENERanno+CRF decoder bet.
Disagreement: only label strength; A/B say replace component, C says abandon CRF route.
Quorum / degraded review status: 3/3.

## Reviewer-proposed directions (ordered A→B→C, verbatim before convergence)
| # | From reviewer | Direction | major_axis | mechanism_delta | Orthogonal to others? | Into this round's cohort? |
|---:|---|---|---|---|---|---|
| 1 | A · Claude | return to non-CRF low-FPR route | decoder/head_arch | remove trained CRF; preserve M19-style independent head/calibration | primary | yes, as route reset |
| 2 | A · Claude | improve emissions/foundation features | backbone/data_view | add or transfer stronger pretrained signal before any structured decoder | yes | later design |
| 3 | A · Claude | deepen FP-aware objective on non-CRF head | objective | improve recall/FPR without CRF transition noise | yes | candidate next |
| 4 | B · Codex | lighter inference-time constraints | decoder/postproc | test rule/constrained postprocess without trained CRF | partly | only if clearly distinct |
| 5 | B · Codex | claim evidence chain | provenance/comparability | continue exclusion-list/split audit for claimable route | yes | ongoing/local |
| 6 | C · Antigravity | abandon CRF decoder direction | route_decision | stop CRF tuning/scale | primary | yes |
| 7 | C · Antigravity | change objective/loss on non-CRF head | objective | FP-aware loss instead of CRF | yes | candidate next |
| 8 | C · Antigravity | change backbone to clean-provenance route | backbone/provenance | transfer lessons to NT-v2 or another clean backbone | yes | candidate next |

## Is tuning justified?
- ❌ CRF tuning is not justified. Both counted seeds miss hard FPR by a wide margin and the best seed is worse than M19 non-CRF on gbF1.
- ✅ A different structural axis remains justified: objective/emission/backbone changes that preserve M19's low-FPR behavior.

## Architecture hypothesis status
Falsified for this local component: standard trained CRF on GENERanno 1.2B LoRA does not recover genes while preserving low FPR. The broader M19 hypothesis remains partially supported: GENERanno 1.2B CDS-preview plus non-CRF FP-aware adaptation is stable adaptation/comparability evidence, but not claim-clean.

## DECISION (choose exactly one)

- [ ] Continue current architecture as-is
- [ ] Tune current architecture
- [ ] Scale data / training as Track B claim run
- [ ] Replace component: decoder/head via CRF
- [ ] Change backbone
- [ ] Change objective / loss
- [ ] Comparability audit first
- [ ] Sanity check first
- [x] Abandon this route: `GENERanno-1.2B + trained CRF decoder`
- [ ] Return to literature

## Why this decision (not another)
The CRF screen answered its specific mechanism question and the answer is negative. Continuing with transition regularization, temperature, or LR tuning would chase a component that failed both key comparisons: lower gbF1 than M19 best and much higher FPR than the hard claim threshold. Choosing `change-backbone` immediately would overreach because the M19 non-CRF route still has useful adaptation evidence; choosing `comparability audit first` would repeat the M20 gate already closed. The precise decision is to abandon the CRF decoder route while preserving M19 as non-claim evidence.

## Best next moves
| Priority | Move | Expected mechanism | Goes to which EXP / Track |
|---:|---|---|---|
| 1 | Freeze CRF as negative/abandoned route | prevent repeated GPU spending on CRF transition tuning | `DEC-001` / decisions log |
| 2 | Design next non-CRF claim path | preserve M19 low-FPR behavior while improving gbF1 or moving to clean backbone | `M22` planning / local |
| 3 | If GPU resumes, use a distinct axis: non-CRF objective/emission or clean-provenance backbone | avoid re-running the failed CRF hypothesis | next screen |

## Parallel cohort this round
- **Primary direction (single)**: abandon GENERanno+CRF decoder route.
- **Parallel cohort**: no new GPU launched from this pivot. Local claim/provenance and M22 design may proceed; any future GPU must be a distinct route.

| Slot | EXP ID (new) | Direction | major_axis | mechanism_delta | Track | Resource profile |
|---|---|---|---|---|---|---|
| primary | `DEC-001-GENERANNO-CRF-DECODER` | abandon CRF decoder route | route_decision | stop trained CRF decoder tuning/scale | decision | local |
| parallel-1 | `M22-NONCRF-CLAIM-ROUTE-DESIGN` | choose non-CRF objective/emission or clean-backbone route | design/provenance | preserve low FPR without CRF | planning | local |

Shared-code conflict? no.

## TODO update
- [x] update docs/05_todo.md.
- [x] update docs/08_pivot_decisions.md.
- [x] run `$note-gate` / evidence register for M21 result + decision.
- [x] update docs/11_master_plan.md state.
- [x] write `docs/09_decisions_log.md` entry for abandoned CRF route.

---

## Mid-iteration note 2026-06-23: m22-noncrf-claim-route

### Mini-retrospective on m22-noncrf-claim-route
- Relevance: directly-attacks-current-gap.
- Does it change our hypothesis? refines: M21 falsified trained CRF, but did not falsify M19's non-CRF low-FPR adaptation evidence. M22 should compare non-CRF objective/emission improvement against clean-provenance backbone transfer before any new GPU.
- Conflicts with an abandoned route (docs/09)? no, provided M22 does not use trained CRF/HMM-style transition tuning from `DEC-001`.
- Recommendation: fold-into-next-batch.
- Urgency: high — this is the immediate M22 local design gate before a new GPU direction.
- Routed idea: `wiki/ideas/m22-noncrf-claim-route.md`.
