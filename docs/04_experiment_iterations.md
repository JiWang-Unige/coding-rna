# Experiment Iterations

> 由 `/goal` 在每轮迭代结束时维护。每条 iteration 一段。
> Track A 是小样本并行筛架构；Track B 是从 Track A 晋升候选后的 scale-up / full validation。

---

## ITER-B0-001: Tiberius bundled mini-smoke reproduction

- Date (UTC): 2026-06-10
- Linked command summary: `$reproduce-baselines` first baseline smoke for Tiberius.
- Experiment ID(s): `BASE-TIBERIUS-MINISMOKE`
- Path tested: baseline reproduction, not a technical roadmap path.
- Milestone: M1 prerequisite.
- Track: baseline.
- Execution mode: run-and-evaluate via `srun`.
- Resource profile: screen-style mini-smoke.
- Claim eligibility: cannot claim.

### Data readiness

- Dataset used: Tiberius repo-bundled `test_data/Panthera_pardus/inp.tar.gz`.
- Downloaded this iteration? no; already archived inside `refs/repos/tiberius-2024`.
- Path / version / split source: `outputs/BASE-TIBERIUS-MINISMOKE/data/inp/`; bundled workflow data, not a claim split.

### Sbatch / run status

- Job id(s): `8527962` completed; setup failures `8527907`, `8527908`; pending private job `8527944` cancelled before running.
- Output dir(s): `outputs/BASE-TIBERIUS-MINISMOKE/`
- Log path(s): `outputs/BASE-TIBERIUS-MINISMOKE/logs/srun_tiberius_predict_4090_bound.log`
- Status: completed.

### Result summary

- Primary metric: `constrained_gene_body_F1 = 0.0` for the provisional project metric because intergenic FPR guardrail failed.
- Supporting metrics: CDS exact F1 `0.8594`; transcript-chain exact F1 `0.3124`; unconstrained gene-body F1 `0.9196`; intergenic FPR `0.0187`.
- Gates: screen_anchor not established; SOTA claim not eligible; review decision not triggered for this B0 smoke.
- Semantic success: pass.

### Tri-review consensus

- Quorum: 3/3 independent CLI reviewers.
- Consensus judgment: `comparability-blocker`.
- Reviewers agreed Tiberius mini-smoke reproduction is semantically successful by official thresholds, while project primary `constrained_gene_body_F1 = 0.0` is an artificial result of provisional hard-zero guardrail failure.
- Main blocker: metric/evaluator comparability. Reference gene-body spans were derived from CDS/intron/start/stop features, while prediction spans used explicit gene/transcript features; the evaluator is not yet the frozen M1 implementation.
- Architecture implication: no negative conclusion about Tiberius-style architecture.

### Pivot decision

- Decision: `Comparability audit first`.
- Do not tune, replace architecture, set `screen_anchor`, or claim from this mini-smoke.
- Next action: re-evaluate existing mini-smoke artifacts with symmetric gene-body masks and non-destructive `intergenic_FPR` sensitivity reporting at `0.005/0.01/0.02`, then continue M1 frozen evaluator + unified baseline screen.

### Links

- result-log: `docs/06_results_log.md#result-base-tiberius-minismoke`
- tri-review: `docs/07_tri_review.md#tri-review-base-tiberius-minismoke`
- pivot: `docs/08_pivot_decisions.md#pivot-decision-base-tiberius-minismoke`
- metrics: `outputs/BASE-TIBERIUS-MINISMOKE/metrics/metrics.json`

---

## ITER-B0-002: Tiberius mini-smoke metric-contract evalfix

- Date (UTC): 2026-06-10
- Linked command summary: post-pivot follow-up to adjust profile-aware `intergenic_FPR` threshold and re-evaluate existing Tiberius mini-smoke artifacts.
- Experiment ID(s): `BASE-TIBERIUS-MINISMOKE-EVALFIX`
- Path tested: metric contract / evaluator sanity, not a technical roadmap architecture path.
- Milestone: M1 prerequisite.
- Track: baseline.
- Execution mode: local evaluation only; no inference/training rerun.
- Resource profile: smoke.
- Claim eligibility: cannot claim.

### Data readiness

- Dataset used: existing Tiberius repo-bundled `test_data/Panthera_pardus/inp.tar.gz` extracted under `outputs/BASE-TIBERIUS-MINISMOKE/data/inp/`.
- Downloaded this iteration? no.
- Path / version / split source: same B0 bundled workflow data, not a claim split.

### Sbatch / run status

- Job id(s): none; local deterministic evaluator over existing artifacts.
- Output dir(s): `outputs/BASE-TIBERIUS-MINISMOKE-EVALFIX/`
- Log path(s): none.
- Status: completed.

### Result summary

- Primary metric: `constrained_gene_body_F1 = 0.9196` under smoke threshold `intergenic_FPR <= 0.02`.
- Sensitivity: fails `0.005` and `0.01`, passes `0.02`.
- Supporting metrics: unconstrained gene-body F1 `0.9196`; intergenic FPR `0.0187`; precision `0.9654`; recall `0.8779`; predicted gene count ratio versus reference genes `0.9872`.
- Gates: smoke cannot claim; screen_anchor not established.
- Semantic success: pass; `validate_goal.py` returns `progress` because smoke/profile and draft goal with placeholder anchors cannot claim success.

### Tri-review consensus

- Quorum: 3/3 independent CLI reviewers.
- Consensus judgment: `continue-current-route`.
- Profile-aware `intergenic_FPR` threshold is acceptable for smoke/screen while preserving full/scale strictness.
- No blocker remains before continuing M1; the earlier low count warning was transcript multiplicity (`154/470` transcripts), not gene-count underprediction (`154/156` genes).

### Pivot decision

- Decision: `Continue current route`.
- Next action: freeze M1 evaluator and run unified Tiberius-like / Helixer-like / ANNEVO-light screen baselines to establish `screen_anchor`.

### Links

- result-log: `docs/06_results_log.md#result-base-tiberius-minismoke-evalfix`
- tri-review: `docs/07_tri_review.md#tri-review-base-tiberius-minismoke-evalfix`
- pivot: `docs/08_pivot_decisions.md#pivot-decision-base-tiberius-minismoke-evalfix`
- metrics: `outputs/BASE-TIBERIUS-MINISMOKE-EVALFIX/metrics/metrics.json`

---

## ITER-<N>: <short description>

- Date (UTC):
- Linked `/goal` command summary:
- Experiment ID(s):
- Path tested: <docs/03 §7.3 Path N>
- Milestone: <M1 / M2 / M3 / M4 / M5>
- Track: <Track A screen / Track B scale-up / baseline / generalization>
- Execution mode: <run-and-evaluate / submit-and-stop>
- Resource profile: <smoke / screen / full / scale>
- Claim eligibility: <cannot claim / claim candidate / robust claim support>

### Track A screen setting (if applicable)

| exp_id | Path | Architecture change | sample_fraction | epochs | patience | seed | config | output_dir |
|---|---|---|---:|---:|---:|---:|---|---|

Promotion rule to Track B:

### Track B scale-up setting (if applicable)

- Promoted from Track A exp_id:
- Promotion reason:
- Scale-up change:
- Success criterion:
- Fallback if fails:

### Hypothesis being tested

### Architecture changes

- What changed:
- Why this is structural rather than hyperparameter tuning:
- Which SOTA weakness it attacks:

### Data readiness

- Dataset(s) used:
- Downloaded this iteration? yes / no
- Path / version / hash / split source:

### Sbatch / run status

- Job id(s):
- Output dir(s):
- Log path(s):
- Status: submitted / running / completed / failed / waiting-for-job
- Resume instruction if submit-and-stop:

### Result summary (if run-and-evaluate or completed job)

- Primary metric: <value> (SOTA: <value>, gap: <abs>)
- Gates: primary <pass/fail>, sota_claim <pass/fail>, review_decision <triggered/not>
- Semantic success: ✅/❌

### Tri-review consensus

### Pivot decision

### Links

- result-log: docs/06_results_log.md#<exp_id>
- tri-review: docs/07_tri_review.md#<ref>
- pivot: docs/08_pivot_decisions.md#<ref>
- decisions-log if abandoned: docs/09_decisions_log.md#<ref>

---

## ITER-B0-003: Tiberius current multi-clade pilot on two RefSeq species

- Date (UTC): 2026-06-10
- Linked command summary: `$reproduce-baselines` M1 pilot for current Tiberius fungi/insecta models under the frozen gene-body evaluator.
- Experiment ID(s): `BASE-TIBERIUS-PILOT-M1`
- Path tested: baseline reproduction / screen-anchor candidate input.
- Milestone: M1 screen baseline reproduction.
- Track: baseline.
- Execution mode: submit-and-handoff, then result-processing.
- Resource profile: screen.
- Claim eligibility: cannot claim.

### Hypothesis being tested

Current Tiberius multi-clade release should provide an operational structured neural gene-caller baseline on small cross-clade pilot species, producing finite gene-body/intergenic metrics that can inform M1 screen-anchor design.

### Architecture changes

- What changed: none; this is a published/current-release baseline reproduction.
- Why structural rather than hyperparameter tuning: N/A.
- Which SOTA weakness it attacks: establishes whether Tiberius-like structured HMM decoding is a usable unified screen baseline outside the bundled mammal mini-smoke.

### Data readiness

- Dataset(s) used: RefSeq `GCF_000146045.2_R64` S. cerevisiae and `GCF_000001215.4_Release_6_plus_ISO1_MT` D. melanogaster.
- Downloaded this iteration? yes, via `scripts/download_refseq_accessions.py`.
- Data check: pass for both selected rows in `data/m1_screen/check_data_report.json`.
- Caveat: species are pilot/reproduction rows, not held-out generalization evidence; D. melanogaster and S. cerevisiae appear in Helixer training lists and ANNEVO current-release tables.

### Sbatch / run status

- Job id(s): `8528176`
- Partition / GPU: `private-teodoro-gpu`, RTX3090 24GB.
- Output dir(s): `outputs/BASE-TIBERIUS-PILOT-M1/`
- Log path(s): `outputs/BASE-TIBERIUS-PILOT-M1/logs/BASE-TIBERIUS-PILOT-M1_8528176.{out,err}`
- Status: failed by validator, not by missing artifacts or OOM.
- Slurm elapsed: `00:18:42`; exit `3:0`.

### Result summary

- Primary metric: `constrained_gene_body_F1 = 0.0`.
- Aggregate supporting metrics: unconstrained gene-body F1 `0.7087`, precision `0.9743`, recall `0.5569`, intergenic FPR `0.0287`.
- Per-species: S. cerevisiae constrained F1 `0.9850`; D. melanogaster constrained F1 `0.0` due to FPR `0.0295`.
- Gates: primary fail; screen_anchor update blocked; SOTA claim not eligible.
- Semantic success: fail under `validate_goal.py` because primary is exactly zero after guardrail hard-zero.

### Tri-review consensus

- 2/3 independent CLI reviewers succeeded; `DEGRADED_REVIEW`.
- Codex: continue M1 baseline reproduction, run Helixer on the same two species, and carry per-species/macro aggregation as an evaluator-contract TODO.
- Antigravity: comparability/semantic-gate blocker; finite poor baselines should not be treated as pipeline `failed_run` when metrics and predictions are valid.
- Claude CLI returned text but failed the required structured-output marker check after retry; not counted toward quorum.

### Pivot decision

- Decision: sanity check first.
- Primary next action: `M1-AGGREGATION-GATE-AUDIT`, a lightweight local audit to distinguish completed-poor finite baseline results from infrastructure/semantic failed runs and to require per-species/macro reporting before any `screen_anchor` update.
- Next baseline after audit: `BASE-HELIXER-SAC-DMEL-SMOKE-M1`.
- `screen_anchor` update remains blocked.

### Links

- result-log: `docs/06_results_log.md#result-base-tiberius-pilot-m1`
- tri-review: `docs/07_tri_review.md#tri-review-base-tiberius-pilot-m1`
- pivot: `docs/08_pivot_decisions.md#pivot-decision-base-tiberius-pilot-m1`
- metrics: `outputs/BASE-TIBERIUS-PILOT-M1/metrics/metrics.json`
- validate: `outputs/BASE-TIBERIUS-PILOT-M1/metrics/validate_goal.json`
- sbatch: `scripts/run_BASE-TIBERIUS-PILOT-M1.sbatch`

---

## ITER-B0-004: M1 aggregation / semantic-gate audit

- Date (UTC): 2026-06-10
- Linked command summary: `$reproduce-baselines` M1; pivot-designated "sanity check first" audit from `BASE-TIBERIUS-PILOT-M1`.
- Experiment ID(s): `M1-AGGREGATION-GATE-AUDIT`
- Path tested: metric/semantic-gate contract + multi-species aggregation, not a roadmap architecture path.
- Milestone: M1 prerequisite.
- Track: baseline.
- Execution mode: local audit (no cluster job), then result-processing.
- Resource profile: local/no-GPU.
- Claim eligibility: cannot claim.

### Hypothesis being tested

The semantic gate should distinguish a finite, valid-but-poor baseline (guardrail hard-zeroed primary) from a true infrastructure/degenerate failure, and the evaluator should surface per-species + macro views so a base-weighted aggregate cannot hide single-species failure before a `screen_anchor` freeze.

### Architecture changes

- What changed: deterministic gate logic + aggregation reporting, not a model.
- `validate_goal.py`: config-driven `completed_poor` vs `failed_run` distinction; `disposition` + `per_species_summary` + heterogeneity warning; implemented the previously-missing `threshold_by_profile`/`profiles` resolution in `check()`.
- `aggregate_gene_body_metrics.py`: added macro fields.
- `ACTIVE_GOAL.json`: `semantic_success` config block (status stays draft).
- Restored 7 M1 scripts lost from `scripts/` (from md5-verified `scripts.backup-20260610-102213/`).

### Data readiness

- Dataset(s) used: existing `BASE-TIBERIUS-PILOT-M1` per-species metrics (`saccharomyces_cerevisiae`, `drosophila_melanogaster`); no new download.

### Sbatch / run status

- Job id(s): none; local deterministic audit on `olympus`.
- Output dir(s): `outputs/M1-AGGREGATION-GATE-AUDIT/`
- Status: completed.

### Result summary

- Gate verdict on PILOT artifacts (corrected): `not_yet` + `disposition=completed_poor`, exit 1 (was `failed_run` exit 3).
- Aggregation: base-weighted constrained F1 0.0 (hard-zeroed) vs macro 0.4925; per-species S.cer 0.985 / D.mel 0.0.
- Tests: 5/5 pass (`tests/test_validate_goal_profiles.py` incl. 3 new audit tests, `tests/test_eval_gene_body_mask.py`).
- Semantic success: pass.

### Tri-review consensus

- Pending: deterministic gate/contract change backed by passing unit tests; tri-review/pivot not yet run for this audit (next gate before proceeding, at user's discretion).

### Pivot decision

- Pending tri-review. Provisional next: `BASE-HELIXER-SAC-DMEL-SMOKE-M1` via `ssh baobab` + `/smart-sbatch`. `screen_anchor` remains blocked; never freeze from a mixed aggregate while a species fails the guardrail.

### Links

- result-log: `docs/06_results_log.md#result-m1-aggregation-gate-audit`
- metrics: `outputs/M1-AGGREGATION-GATE-AUDIT/metrics/metrics.json`
- validate: `outputs/M1-AGGREGATION-GATE-AUDIT/metrics/validate_goal.json`
- contract: `docs/11_evaluator_contract.md`

---

## ITER-B0-005: Helixer two-species smoke (broad-lineage baseline comparator)

- Date (UTC): 2026-06-10
- Linked command summary: `$reproduce-baselines` M1; pivot direction "Helixer two-species smoke under same evaluator".
- Experiment ID(s): `BASE-HELIXER-SAC-DMEL-SMOKE-M1`
- Path tested: baseline reproduction / screen_anchor candidate input.
- Milestone: M1 screen baseline reproduction.
- Track: baseline.
- Execution mode: remote_ssh submit-and-handoff (ssh baobab -> sbatch), then result-processing.
- Resource profile: smoke.
- Claim eligibility: cannot claim.

### Hypothesis being tested

Helixer current lineage models (fungi / invertebrate) should provide an operational broad-lineage gene-caller baseline on the two pilot species under the frozen gene-body evaluator, adding a second reference architecture toward `screen_anchor = max(...)`.

### Architecture changes

- None; published baseline reproduction. Helixer fungi(S.cer) + invertebrate(D.mel) via `--model-filepath` (offline-safe) + `--subsequence-length` (21384 / 213840).

### Data readiness

- Same RefSeq S.cerevisiae + D.melanogaster as PILOT; pilot/runner species, not held-out generalization.

### Sbatch / run status

- Job id(s): `8530017` (FAILED@30s, missing --subsequence-length) -> `8530344` (COMPLETED).
- Partition / GPU: `private-teodoro-gpu` gpu035, RTX3090 24GB. Helixer D.mel 1.10 h.
- Output dir: `outputs/BASE-HELIXER-SAC-DMEL-SMOKE-M1/`
- Status: completed; validator exit 1 (`not_yet`, `disposition=completed_poor`).

### Result summary

- base-weighted: constrained 0.0 (hard-zeroed), unconstrained 0.9009, FPR 0.0877; macro: unconstrained 0.8944, FPR 0.3525.
- Per-species: S.cer unconstrained 0.8864 / FPR 0.6544; D.mel unconstrained 0.9025 / FPR 0.0506.
- Runner path validated end-to-end; gate `completed_poor` verified in production.
- HEADLINE: cross-tool incomparability — Helixer GFF3 includes UTR vs Tiberius CDS-only span; S.cer FPR 0.654 vs 0.016 on same reference. Harmonize gene-body span before screen_anchor freeze.

### Tri-review consensus

- Pending / optional (smoke runner validation + comparability finding; non-claim). Recommended to tri-review the span-harmonization decision before freezing screen_anchor.

### Pivot decision

- Provisional next: harmonize gene-body span (CDS-only or dual-report) and re-evaluate existing Tiberius + Helixer predictions locally; then unified screen + ANNEVO-light -> screen_anchor. screen_anchor remains blocked.

### Links

- result-log: `docs/06_results_log.md#result-base-helixer-sac-dmel-smoke-m1`
- metrics: `outputs/BASE-HELIXER-SAC-DMEL-SMOKE-M1/metrics/metrics.json`
- validate: `outputs/BASE-HELIXER-SAC-DMEL-SMOKE-M1/metrics/validate_goal.json`
- sbatch: `scripts/run_BASE-HELIXER-SAC-DMEL-SMOKE-M1.sbatch`

---

## ITER-B0-006: CDS-span harmonization + unconstrained screen anchor

- Date (UTC): 2026-06-10
- Linked command summary: `$reproduce-baselines` M1 evaluator-contract harmonization + goal revision (user option-2).
- Experiment ID(s): `M1-SPAN-HARMONIZE-CDS`
- Path tested: evaluator comparability contract (gene-body span) + screen-anchor metric definition.
- Milestone: M1 prerequisite (screen_anchor enablement).
- Track: baseline.
- Execution mode: local re-eval (no GPU) + goal revision.
- Resource profile: local.
- Claim eligibility: cannot claim.

### Hypothesis being tested

Cross-tool gene-body comparison must use a tool-fair span. CDS-only span (every caller emits CDS; UTR is tool-dependent) should remove the Helixer-vs-Tiberius incomparability surfaced by the Helixer smoke.

### Architecture changes

- `scripts/eval_gene_body_mask.py`: added `--span-mode {transcript,cds}`; cross-tool screen uses `cds`. Backward-compatible default `transcript`; 5/5 tests pass.
- Goal revision (ACTIVE_GOAL.json): `screen_anchor` -> `gene_body_F1_unconstrained` (CDS), value 0.9213 provisional; `intergenic_FPR` guardrail scoped to `[full,scale]` (advisory for smoke/screen). See docs/08 Goal Revision.

### Data readiness

- Re-evaluated existing Tiberius PILOT + Helixer smoke predictions (both species); no new data.

### Sbatch / run status

- Job id(s): none; local deterministic re-eval.
- Output dir: `outputs/M1-SPAN-HARMONIZE-CDS/`
- Status: completed.

### Result summary

- Comparability FIXED under CDS span: Helixer S.cer FPR 0.654->0.033, pred/ref ratio 1.224->0.99; Tiberius/Helixer comparable (Helixer slightly stronger).
- Finding: 0.02 screen FPR too strict for SOTA tools on gene-dense pilot species (both ~0.0225, pass@0.025). Resolved by user option-2 (unconstrained CDS F1 screen anchor + FPR advisory for screen).
- screen_anchor now SET (provisional): gene_body_F1_unconstrained 0.9213 (CDS). Gate verified.

### Tri-review consensus

- Pending/optional (deterministic, test-backed, user-directed contract change). Recommended before any Track A promotion / screen_anchor finalization.

### Pivot decision

- Next: ANNEVO-light under cds span (may raise anchor); re-derive on frozen anchor species; then Track A portfolio vs screen_anchor.

### Links

- result-log: `docs/06_results_log.md#result-m1-span-harmonize-cds`
- goal revision: `docs/08_pivot_decisions.md` (Goal Revision 2026-06-10 CDS-span)
- contract: `docs/11_evaluator_contract.md`
- metrics: `outputs/M1-SPAN-HARMONIZE-CDS/metrics/`

---

## ITER-B0-007: ANNEVO two-species smoke — completes the baseline trio

- Date (UTC): 2026-06-10
- Linked command summary: `$reproduce-baselines` M1, third full gene-caller (ANNEVO-light).
- Experiment ID(s): `SETUP-ANNEVO-ENV`, `BASE-ANNEVO-SAC-DMEL-SMOKE-M1`
- Path tested: baseline reproduction / screen_anchor completion.
- Milestone: M1 (reproduce-baselines core trio complete).
- Track: baseline.
- Execution mode: remote_ssh submit-and-handoff (observed to completion).
- Resource profile: screen (shared-gpu).
- Claim eligibility: cannot claim.

### Hypothesis being tested

ANNEVO (provisional published-SOTA anchor candidate) should run as the third operational baseline under the fair CDS evaluator, completing the screen_anchor = max(Tiberius, Helixer, ANNEVO) computation.

### Architecture changes

- None (published baseline). Dedicated `annevo` conda env (Py3.10/torch2.1/cu12.1, sanctioned exception); annotation.py one-step (predict+decode) Fungi(S.cer)/Insecta(D.mel); eval --span-mode cds.

### Data readiness

- Same RefSeq S.cer + D.mel pilot species; no new data.

### Sbatch / run status

- Job id(s): 8536914 (private PENDING→cancel), 8537041/8537118 (FAILED: set-u activate + AF_UNIX TMPDIR), 8537422 (COMPLETED, shared-gpu gpu021, 18:09).
- Output dir: `outputs/BASE-ANNEVO-SAC-DMEL-SMOKE-M1/`
- Status: completed; gate not_yet/completed_poor.

### Result summary

- ANNEVO CDS: base-w F1 0.9197 / macro 0.9429; S.cer FPR 0.0072 (lowest of 3), prec 0.9971.
- Trio (base-w CDS unconstrained F1): Tiberius 0.8608 / Helixer 0.9213 / ANNEVO 0.9197 → screen_anchor = max = Helixer 0.9213 (ANNEVO ties, doesn't raise).
- reproduce-baselines core trio COMPLETE under the fair CDS evaluator.

### Tri-review consensus

- Pending/optional. Recommended: tri-review the combined M1 contract (gate audit + CDS span + unconstrained anchor + 3-tool screen_anchor) before Track A promotion.

### Pivot decision

- Next: Track A architecture portfolio vs screen_anchor 0.9213 (primary track foundation_probe→semi-CRF decoder), OR re-derive screen_anchor on frozen typical-intergenic anchor species first. SegmentNT/GENERanno probes deferred.

### Links

- result-log: `docs/06_results_log.md#result-base-annevo-sac-dmel-smoke-m1`
- metrics: `outputs/BASE-ANNEVO-SAC-DMEL-SMOKE-M1/metrics/`
- sbatch: `scripts/run_BASE-ANNEVO-SAC-DMEL-SMOKE-M1.sbatch`; env: `scripts/setup_annevo_env.sh`

---

## ITER-B0-008: M1-SAMEBUDGET-SCREEN-ANCHOR — true same-budget screen_anchor

- Date (UTC): 2026-06-10
- Linked command summary: `$implement` from pivot M1-CONTRACT-REVIEW (user scope=2 family refs).
- Experiment ID(s): `SCREENREF-{tiberius_like,helixer_like}-s{0,1,2}`, `FLOOR-SCREEN-M1`, env `SETUP-CODINGRNA-TORCH`
- Path tested: same-budget reference baselines (random-init small-sample training) to set screen_anchor.
- Milestone: M1 (screen_anchor establishment). Track: baseline. Mode: submit-and-monitor. Profile: screen. Claim: cannot.

### Hypothesis being tested
A fair Track A bar must come from random-init reference architectures trained under our own small-sample budget, not from pretrained-model inference. The same-budget value should be well below the pretrained ceiling (0.9213).

### Architecture changes
- New harness `src/screen_anchor/{data,models,gff_io,train_screen_ref,floor_baseline}.py`: per-base 3-class labels from GFF, chromosome-level split, windowed one-hot dataset, two light random-init torch models (Tiberius-like CNN+biLSTM; Helixer-like dilated-conv, no recurrence = mechanism delta), ORF FLOOR, predict to CDS GFF, eval (--span-mode cds), aggregate.
- torch 2.5.1 installed into project env coding-rna.

### Data readiness
- 2 species (yeast+fly), chromosome-level split (check_data PASS: no seqid leakage). sample_fraction 0.3.

### Run status
- 6 GPU runs 8538949-8538954 COMPLETED on shared-gpu (8538912-17 FAILED@1s MKL set-u, fixed). FLOOR local.
- Output dirs: outputs/SCREENREF-*/ , outputs/FLOOR-SCREEN-M1/.

### Result summary
- Bracket (base-w CDS F1): FLOOR 0.3735 < screen_anchor 0.5579 (tiberius 0.5576 / helixer 0.5579 seed-means) < pretrained_ceiling 0.9213.
- screen_anchor set = 0.5579 in the goal contract. Confirms same-budget (0.56) far below pretrained ceiling (0.92).
- Semantic success: pass (6/6 runs completed, finite metrics, check_data + smoke gates passed).

### Tri-review consensus
- The contract decision behind this was tri-reviewed (M1-CONTRACT-REVIEW, 2/3 comparability-blocker). This execution iteration is the agreed build-out; per-run tri-review not separately required (non-claim screen references).

### Pivot decision
- Track A now unblocked: candidates must strictly exceed screen_anchor 0.5579. Pending: tighten completed_poor gate; re-derive on frozen anchor species; M2 sota_benchmark.

### Links
- result-log: `docs/06_results_log.md#result-m1-samebudget-screen-anchor`
- harness: `src/screen_anchor/`; protocol script under scripts/
- metrics: `outputs/SCREENREF-*/metrics/metrics.json`, `outputs/FLOOR-SCREEN-M1/metrics/metrics.json`

---

### ITER-B0-008 component seed runs (per-exp_id, for ledger reconciliation)
SCREENREF-tiberius_like-s0, SCREENREF-tiberius_like-s1, SCREENREF-tiberius_like-s2, SCREENREF-helixer_like-s0, SCREENREF-helixer_like-s1, SCREENREF-helixer_like-s2 — all COMPLETED; collectively analyzed in ITER-B0-008 / docs/06 Result: M1-SAMEBUDGET-SCREEN-ANCHOR.

---

## ITER-B0-009: TA-DECODER-M3 (Track A structured-decoder focused batch)
- Date 2026-06-11. Track A screen (M3), run-and-evaluate. Decoders on fixed tiberius_like backbone vs same-budget screen_anchor 0.5576.
- CONSTR (constrained-Viterbi) WON: seed-mean F1 0.5791 > gate 0.5676, ratio 2.74→1.12 (< 1.25 claim guardrail). primary_progress_gate MET; R5 satisfied (beat anchor + fixed coherence). But it is POST-PROCESSING, not learned structure.
- CRF dropped (correct, tests 5/5, but W=2048 forward+backward too slow); semi-CRF dropped (pure-python segment DP intractable). Both tractability, not correctness.
- tri-review 2/3 DEGRADED, 1-1 split (Claude continue/HOLD-Track-B vs Codex scale-to-Track-B); consensus = vectorize CRF/semi-CRF next + CONSTR is post-processing. Pivot DEFERRED the Track-B-now decision to user (autonomy exception: tie + resource).
- Components: SCREENREF-tiberius_like-constrained-s0, SCREENREF-tiberius_like-constrained-s1, SCREENREF-tiberius_like-constrained-s2, SCREENREF-tiberius_like-crf-s0, SCREENREF-tiberius_like-crf-s1, SCREENREF-tiberius_like-crf-s2.
- Links: docs/06#result-ta-decoder-m3, docs/07#tri-review-ta-decoder-m3, docs/08#pivot-decision-ta-decoder-m3.

---

## ITER-B0-010: TA-DECODER-VEC-M3 (vectorized LEARNED structured decoder)
- Date 2026-06-11. Track A screen (M3), run-and-evaluate. Decoders on FIXED tiberius_like backbone vs anchor 0.5576 + CONSTR 0.5791.
- WON: CRF-vec (vectorized linear-chain CRF) seed-mean F1 0.6186 (>gate 0.5676, >anchor, >CONSTR), ratio 0.88. Ladder: anchor 0.5576/2.74 < CONSTR 0.5791/1.12 < CRF-vec 0.6186/0.88. LEARNED structure beats post-processing beats per-base — core architecture bet validated.
- Vectorization: parallel-scan partition (O(log W)) + per-token NLL normalization + batched Viterbi predict; unit-tested vectorized==reference. semi-CRF DROPPED (vectorize segment DP deferred).
- CAVEAT: high seed variance (0.58-0.66; s2 loses to CONSTR). tri-review 2/3 DEGRADED (agy flaky), 2-0 consensus scale-to-track-b; Track B job#1 = seeds≥5-8 + CI + scale-data scalability test. LAUNCH pending user.
- Components: SCREENREF-tiberius_like-crf-s0, SCREENREF-tiberius_like-crf-s1, SCREENREF-tiberius_like-crf-s2.
- Links: docs/06#result-ta-decoder-vec-m3, docs/07#tri-review-ta-decoder-vec-m3, docs/08#pivot-decision-ta-decoder-vec-m3.

---


## ITER-FP-001 — FP-SEGMENTNT-PROBE-M1 (2026-06-11)
- Track: A screen (NON-CLAIM), foundation-probe #1 (first post-ruler-change architecture move).
- Architecture change: input = FROZEN SegmentNT(multi_species) 14 element logits (vs raw one-hot DNA), anchor-matched conv+biLSTM head -> clean INPUT-SIGNAL ablation. major_axis=training_signal/data_view; mechanism_delta=pretrained-foundation-features-as-input.
- Result: AXIS-2 F1 0.6888 (>> anchor 0.5576); AXIS-1 spec 0.8416 (< anchor 0.8710), macro 0.7543 (< 0.7978). not_yet, not Pareto-dominant.
- Pivot: iterate-probe / change-objective-or-loss (3/3). Parent of next round (FP-aware loss + raw-DNA fusion). See docs/08.
- Component exp_ids (3 seeds, all DONE/result-logged): FP-SEGMENTNT-PROBE-M1-convlstm-s0 FP-SEGMENTNT-PROBE-M1-convlstm-s1 FP-SEGMENTNT-PROBE-M1-convlstm-s2 (jobs 8548460-62).

## ITER-FP-002 -- TA-FOUNDATION-DECODER-M4 (2026-06-11)
- Track A screen (NON-CLAIM), MAIN architecture bet: foundation features -> structured decoder. 3 candidates x 5 seeds.
- FP-SEGNT-FPLOSS (loss_design): spec 0.9303 > anchor 0.8710, gbF1 0.6157 > anchor 0.5576, macro 0.8431 > gate -> PARETO-beats anchor -> PASS (Track-B candidate).
- FP-SEGNT-FUSION (data_view): spec 0.8615 < anchor -> no. FP-SEGNT-CRF (decoder): spec 0.8298 high-variance < anchor -> no (but best gene_count coherence 0.90).
- Parent: FP-SEGMENTNT-PROBE-M1 (ITER-FP-001). Pivot: see docs/08. Components: FP-SEGNT-{FPLOSS,FUSION,CRF}-s0..4 (jobs 8550151-66).
### Component exp_ids (ledger reconciliation)
FP-SEGNT-FPLOSS-s0 FP-SEGNT-FPLOSS-s1 FP-SEGNT-FPLOSS-s2 FP-SEGNT-FPLOSS-s3 FP-SEGNT-FPLOSS-s4 FP-SEGNT-FUSION-s0 FP-SEGNT-FUSION-s1 FP-SEGNT-FUSION-s2 FP-SEGNT-FUSION-s3 FP-SEGNT-FUSION-s4 FP-SEGNT-CRF-s0 FP-SEGNT-CRF-s1 FP-SEGNT-CRF-s2 FP-SEGNT-CRF-s3 FP-SEGNT-CRF-s4 — 15 seed-runs of ONE logical experiment TA-FOUNDATION-DECODER-M4, jobs 8550151-66, all COMPLETED.

## ITER-FP-003 -- TA-COHERENCE-FIX-M5 (2026-06-11)
- Track A screen (NON-CLAIM), M4 pivot follow-up: de-fragment FPLOSS + 5-seed anchor paired test.
- FP-FRAGFIX-CONSTR (FPLOSS + constrained post-proc): spec 0.9272 (paired +0.0836±0.037 vs 5-seed anchor 0.8436), gbF1 0.6581, macro 0.8555, gene_count 2.25->1.28. PARETO-beats anchor, paired-significant.
- 5-seed anchor 0.8436 (down from 3-seed 0.8710). Parent: TA-FOUNDATION-DECODER-M4 (ITER-FP-002). Pivot: docs/08. Components: FP-FRAGFIX-CONSTR-s0..4 + SCREENREF-tiberius_like-s3,s4.
- Component exp_ids: FP-FRAGFIX-CONSTR-s0 FP-FRAGFIX-CONSTR-s1 FP-FRAGFIX-CONSTR-s2 FP-FRAGFIX-CONSTR-s3 FP-FRAGFIX-CONSTR-s4 SCREENREF-tiberius_like-s3 SCREENREF-tiberius_like-s4.

## ITER-FP-004 -- TA-FRAGFIX-SWEEP-M6 (2026-06-11)
- Track A screen (NON-CLAIM), STEP-0 promote-gate: VAL-chosen constrained param sweep clears gene_count 1.28->0.939 (<=1.25).
- FP-FRAGFIX-CONSTR (rp, mfg=20/mcl=90): test spec 0.9262>anchor, gbF1 0.6376, macro 0.8389, gene_count 0.939 -> ALL 4 GATES PASS -> PROMOTE-READY.
- Parent: TA-COHERENCE-FIX-M5 (ITER-FP-003). Pivot: docs/08 (promote to Track B = user go-ahead). Component exp_ids: FP-FRAGFIX-CONSTR-rp-s0 FP-FRAGFIX-CONSTR-rp-s1 FP-FRAGFIX-CONSTR-rp-s2 FP-FRAGFIX-CONSTR-rp-s3 FP-FRAGFIX-CONSTR-rp-s4.

## ITER-FP-005 -- REANCHOR-HELDOUT-M7 (2026-06-12)
- Track A screen (NON-CLAIM), retrospective-derived re-anchor GATE before ③ Track-B. Held-out/UTR-rich cross-clade {Arabidopsis thaliana, Gallus gallus(subset)}.
- Held-out same-budget anchor (tiberius_like 3-seed): spec 0.8054 / macro 0.7804 / gbF1 0.7099. ANNEVO ceiling spec 0.9824.
- Candidate FP-FRAGFIX-CONSTR (5-seed, IDENTICAL config, VAL-selected mfg=20/mcl=60): spec 0.9604+-0.008 (all 5 > anchor, +0.155), macro 0.9621, gbF1 0.6664, gcount 0.9688 (<=1.25). HELD-OUT PARETO-PASS.
- Margin LARGER cross-clade (+0.155) than yeast+fly (+0.078); absolute HIGHER (0.9604 vs 0.9218). Retrospective concern REFUTED. Parent: TA-FRAGFIX-SWEEP-M6 (ITER-FP-004). Pivot: docs/08.
- Component exp_ids: FP-FRAGFIX-CONSTR-ho-s0 FP-FRAGFIX-CONSTR-ho-s1 FP-FRAGFIX-CONSTR-ho-s2 FP-FRAGFIX-CONSTR-ho-s3 FP-FRAGFIX-CONSTR-ho-s4 SCREENREF-tiberius_like-ho-s0 SCREENREF-tiberius_like-ho-s1 SCREENREF-tiberius_like-ho-s2 REANCHOR-CEILING-ANNEVO-M7 FP-SEGMENTNT-FEATCACHE-M7.

## ITER-FP-006 -- TB-GBF1-MULTICLASS-M8 (2026-06-12)
- Track B scale-up (NON-CLAIM), ③: multi-class structured output for gbF1 recovery on CLEAN held-out plants {arabidopsis,rice}.
- NEGATIVE (primary): mc-candidate gbF1 0.7189 NOT > 3c-candidate 0.7392 (multi-class did NOT recover gbF1; gcount 0.66 under-pred). M8 bet REFUTED.
- POSITIVE (clean): 3c-candidate Pareto-beats raw-DNA anchor on clean plants BOTH axes (spec 0.966>0.905, gbF1 0.739>0.696) — leakage-free (SegmentNT backbone excludes plants).
- Next axis: staged UNFREEZE/fine-tune SegmentNT or backbone-only self-train (frozen features cap gbF1). Parent: REANCHOR-HELDOUT-M7. Pivot: docs/08.
- Components: M8-MC-CAND-s0..4, M8-3C-CAND-s0/s2/s4 (s1/s3 transient-fail), SCREENREF-tiberius_like-m8clean-s0..2.

### ITER-FP-006 component exp_ids (ledger, verbatim): M8-MC-CAND-s0 M8-MC-CAND-s1 M8-MC-CAND-s2 M8-MC-CAND-s3 M8-MC-CAND-s4 M8-3C-CAND-s0 M8-3C-CAND-s1 M8-3C-CAND-s2 M8-3C-CAND-s3 M8-3C-CAND-s4 SCREENREF-tiberius_like-m8clean-s0 SCREENREF-tiberius_like-m8clean-s1 SCREENREF-tiberius_like-m8clean-s2 M8-MC-SMOKE

## ITER: TB-UNFREEZE-BACKBONE-M9 (CK3) — 2026-06-14

- **Track**: A screen (unfreeze-depth probe), NON-CLAIM. major_axis=training_signal, mechanism_delta=backbone top-N layer unfreeze (frozen vs L2 vs L4). FOCUSED ARCH BATCH on training_signal.
- **Component exp_ids (verbatim)**: M9-UNFREEZE-L0-s0, M9-UNFREEZE-L2-s0, M9-UNFREEZE-L4-s0. Jobs 8667188/8667189/8667190 COMPLETED ~9.5h, private-teodoro-gpu, single-species arabidopsis sample0.3 epochs4, generanno env.
- **Hypothesis**: M8 found frozen SegmentNT features cap gbF1; unfreezing NT-v2-500m top-N layers should lift gbF1 past frozen ceiling.
- **Result**: VALIDATED (dual-axis monotonic). gbF1 0.8284->0.8544->0.8759; spec 0.9656->0.9669->0.9754. L4 gbF1 nears ANNEVO 0.898. validate_goal=not_yet (constrained_gbF1=0 because L4 FPR 0.0246>0.02 screen threshold; eval line 237 hard gate).
- **Pivot**: continue-current-route -> deeper unfreeze (L6/L8/full) to push FPR<0.02 + multi-seed + cross-species (Track B, needs user gate >24h). See docs/08.
- **Prior-failure note**: 2 earlier batches (8575441-49, 8623290-92) TIMEOUT'd due to Write-tool non-persistence (docs/10 2026-06-13); fixed via Bash/ssh + sacct verification.

## ITER: TB-UNFREEZE-BACKBONE-M9-DEEP — 2026-06-14

- **Track**: A screen / Track-B preflight, NON-CLAIM. major_axis=training_signal, mechanism_delta=deeper NT-v2-500m top-layer unfreeze (L6/L8/L12) after L4 hit FPR 0.0246.
- **Execution mode**: submit-and-handoff; jobs 8751498/8751499/8751500 COMPLETED on RTX 3090, ~9.6-9.8h each.
- **Hypothesis**: deeper unfreeze should improve emissions enough to push FPR below 0.02 while preserving the L4 gbF1 lift.
- **Architecture change**: same 3-class intron-aware convLSTM head + FP-aware loss + constrained postproc as M9 CK3; only unfreeze depth changed from L4 to L6/L8/L12.
- **Result summary**: L6/L8/L12 all break FPR<=0.02. Best L12: intergenic_specificity 0.9810, FPR 0.0190, gbF1 0.9035, constrained_gbF1 0.9035, gene_count_ratio 0.792 on arabidopsis seed0. validate_goal=progress after fixing profile-aware guardrail regression.
- **Tri-review consensus**: 3/3 quorum; consensus = M9-L12 is primary next Track-B route, GENERanno LoRA+3-class is parallel challenger if budget allows.
- **Pivot decision**: scale-to-track-b. Primary next EXP=`M10-M9L12-CLEANPLANTS`; optional parallel EXP=`M10-GENERANNO-LORA-3C`. See docs/08.
- **Links**: docs/06 `Result: TB-UNFREEZE-BACKBONE-M9-DEEP`; sbatch `scripts/run_M9_arm.sbatch`; output dirs `outputs/M9-UNFREEZE-L{6,8,12}-s0`; logs `outputs/fp_segnt_logs/M9ARM_8751498-8751500.*`.

## ITER: M10-GENERANNO-LORA-3C-SMOKE — 2026-06-15

- **Track**: parallel challenger smoke, NON-CLAIM. major_axis=backbone_finetune, mechanism_delta=GENERanno 1.2b CDS-pretrained encoder + LoRA adapters + our intron-aware 3-class head.
- **Execution mode**: submit-and-handoff; job 8833070 completed on shared A100 40GB (`gpu022`), about 1h34m runtime after queue start.
- **Hypothesis**: replacing GENERanno's native binary CDS head with a 3-class intron-aware head should be technically feasible and is the correct path to repair native CDS-mask fragmentation.
- **Architecture change**: discarded released binary token classifier for prediction; froze base encoder, trained LoRA(q/k/v/o,r=8) plus 3-class FP-aware convLSTM head; smoke-limited to 8 train windows / 4 val windows / 1 arabidopsis test seqid.
- **Result summary**: runtime integration PASS, metric quality NOT READY. `intergenic_specificity=0.9491`, `intergenic_FPR=0.0509`, `gbF1=0.7525`, `constrained_gbF1=0.0`, gene_count_ratio `4.43`. Corrected validate=`progress` (smoke non-claim) after rerunning with status-file path.
- **Tri-review consensus**: combined M10 review 3/3 says park GENERanno LoRA for now; runtime-positive but metric-negative and not screen-ready.
- **Pivot decision**: park prepared GENERanno screen; only revisit with redesigned specificity-preserving schedule after M9-L12 calibration.
- **Links**: docs/06 `Result: M10-GENERANNO-LORA-3C-SMOKE`; sbatch `sbatch/M10-GENERANNO-LORA-3C-SMOKE.sbatch`; output dir `outputs/M10-GENERANNO-LORA-3C-SMOKE`; log `outputs/fp_segnt_logs/M10GENLORA_8833070.out`.

## ITER: M10-M9L12-CLEANPLANTS — 2026-06-15

- **Track**: screen / Track-B preflight, NON-CLAIM. major_axis=training_signal, mechanism_delta=NT-v2-500m top-12 unfreeze scaled from single arabidopsis seed to multi-seed clean plants `{arabidopsis,rice}`.
- **Component exp_ids (verbatim)**: `M10-M9L12-CLEANPLANTS-s0`, `M10-M9L12-CLEANPLANTS-s1`, `M10-M9L12-CLEANPLANTS-s2`. Slurm array `8833071_[0-2]`, private-teodoro-gpu `gpu034`, all COMPLETED.
- **Hypothesis**: M9-DEEP L12's arabidopsis seed0 win should generalize across clean plants and seeds, preserving high gbF1 while improving specificity enough to become the main route toward claim.
- **Architecture change**: same M9 L12 NT-v2 unfreeze + 3-class FP-aware convLSTM head + constrained postproc; scale axis is clean-plant pooled training and three seeds.
- **Result summary**: strong progress, not claim-ready. Seed mean `intergenic_specificity=0.9826`, `intergenic_FPR=0.0174`, `macro_specificity=0.9801`, `gbF1=0.8398`, gene_count_ratio `0.897`. validate=`progress`; full/scale FPR `<=0.01` remains unmet.
- **Tri-review consensus**: 3/3 `continue-current-route`; mainline remains M9-L12. Immediate blocker is FPR/constrained operating point, especially arabidopsis constrained zero due FPR>0.02.
- **Pivot decision**: targeted M9-L12 specificity calibration before full-scale promotion. Next primary EXP=`M11-L12-SPEC-CALIBRATION`; park GENERanno LoRA.
- **Links**: docs/06 `Result: M10-M9L12-CLEANPLANTS`; sbatch `sbatch/M10-M9L12-CLEANPLANTS.sbatch`; output dirs `outputs/M10-M9L12-CLEANPLANTS-s{0,1,2}`; logs `outputs/fp_segnt_logs/M10M9L12_8833071_{0,1,2}.out`.

## ITER: M11-L12-SPEC-CALIBRATION — 2026-06-16

- **Track**: screen / Track-B preflight, NON-CLAIM. major_axis=decode_calibration, mechanism_delta=validation-only intergenic-bias + constrained-param operating point selection on M9-L12 raw emissions.
- **Component exp_ids (verbatim)**: `M11-L12-SPEC-CALIBRATION-s0`, `M11-L12-SPEC-CALIBRATION-s1`, `M11-L12-SPEC-CALIBRATION-s2`. Slurm array `8934130_[0-2]`, private-teodoro-gpu `gpu034`, all COMPLETED.
- **Hypothesis**: M10's remaining FPR tail can be corrected by no-leakage validation-only decode calibration, preserving the strong M9-L12 gbF1 and gene-count coherence; stronger FP objective is only needed if calibration cannot reach `FPR<=0.01`.
- **Architecture change**: same M10 training (`fp_lambda=1.0`, L12 unfreeze), but trainer saves VAL/TEST raw emissions; calibrator sweeps intergenic bias, `min_cds_len`, and `max_fill_gap` on VAL only, then applies selected point once to TEST.
- **Result summary**: SUCCESS for stated blocker. Seed mean `intergenic_specificity=0.9913`, `intergenic_FPR=0.0087`, `macro_specificity=0.9909`, `gbF1=0.8178`, constrained_gbF1@0.01 `0.8178`, gene_count_ratio `1.003`. All seeds pass aggregate `FPR<=0.01` and `gene_count<=1.25`.
- **Tri-review consensus**: 2/3 `DEGRADED_REVIEW` (Codex + Antigravity) both `scale-to-track-b`; Claude failed the required marker heuristic and was not counted. Consensus: promote calibrated M9-L12 to full/scale/comparability preparation; no stronger FP objective now.
- **Pivot decision**: scale data/training via `M12-M9L12-FULLSCALE-CALIBRATED` preflight. Freeze VAL-only calibration protocol, close ANNEVO-compatible benchmark/pretraining-overlap blockers, save checkpoints, and report aggregate + per-species FPR sensitivity. Stronger FP objective remains fallback only if full/scale FPR fails.
- **Links**: docs/06 `Result: M11-L12-SPEC-CALIBRATION`; sbatch `sbatch/M11-L12-SPEC-CALIBRATION.sbatch`; output dirs `outputs/M11-L12-SPEC-CALIBRATION-s{0,1,2}`; logs `outputs/fp_segnt_logs/M11L12CAL_8934130_{0,1,2}.out`.

## ITER: M12-PUBLICATION-PREFLIGHT-TWOSEED — 2026-06-17

- **Track**: screen / publication-alignment preflight, NON-CLAIM. major_axis=paper_evidence, mechanism_delta=fixed-model cross-species evaluation + same-panel external baselines + bounded GENERanno challenger.
- **Component exp_ids (verbatim)**: `M12A-FIXEDMODEL-CROSSSPECIES-A2R-s0`, `M12A-FIXEDMODEL-CROSSSPECIES-A2R-s1`, `M12A-FIXEDMODEL-CROSSSPECIES-A2R-s2`, `M12B-SAMEPANEL-BASELINES-ANNEVO`, `M12B-SAMEPANEL-BASELINES-TIBERIUS`, `M12B-SAMEPANEL-BASELINES-HELIXER`, `M12C-GENERANNO-1P2B-CDS-SMOKE`, `M12C-GENERANNO-0P5B-BASE-SMOKE`.
- **Execution mode**: submit-and-monitor via smart-sbatch. M12A array `8974902_[0-2]`; M12B repaired array `8982048_[0-2]` after short `/tmp` TMPDIR repair; M12C array `8974903_[0-1]`.
- **Hypothesis**: M11's strong calibrated M9-L12 result must be tested as a paper-facing model: fixed train/calibrate on Arabidopsis should transfer to unseen rice; external gene callers should be compared on the same panel; GENERanno variants should test whether pretrained models generally solve the task.
- **Result summary (3-seed closure)**: M12A fixed Arabidopsis->rice is negative in all three seeds: mean `gbF1=0.6556`, `specificity=0.9689`, `FPR=0.0311`, constrained_gbF1@0.01 `0.0`, gene_count_ratio `1.755`. Seed2 completed after the initial two-seed user gate and confirmed the same direction.
- **External baseline summary**: Tiberius same panel is strongest on FPR (`specificity=0.9927`, `FPR=0.0073`, `gbF1=0.9252`) but under-calls genes (`gene_count_ratio=0.628`). ANNEVO has `gbF1=0.9269`, `FPR=0.0117`; Helixer has `gbF1=0.9220`, `FPR=0.0216`.
- **GENERanno summary**: 1.2B CDS-preview has signal but fails coherence (`gbF1=0.7527`, `FPR=0.0432`, gene_count_ratio `4.405`); 0.5B base collapses (`FPR~1.0`, gene_count_ratio `0.0002`).
- **Tri-review consensus**: formal M13-distance/generalization tri-review reached 3/3 quorum for `run-sanity-check-first`: a close Arabidopsis-relative plant scan is reasonable only as a bounded non-claim diagnostic after zero-GPU M12A failure-mode analysis and clean species/provenance selection.
- **Pivot decision**: pending formal pivot. Interim direction: stop M9-only micro-optimization as the mainline; first analyze M12A failure versus M11 pooled rice, then optionally run a frozen single-seed M13 distance scan.
- **Links**: docs/06 `Result: M12-PUBLICATION-PREFLIGHT-TWOSEED`; M12A outputs `outputs/M12A-FIXEDMODEL-CROSSSPECIES-A2R-s{0,1,2}`; M12B outputs `outputs/M12B-SAMEPANEL-BASELINES-{ANNEVO,TIBERIUS,HELIXER}`; M12C outputs `outputs/M12C-GENERANNO-{1P2B-CDS,0P5B-BASE}-SMOKE`.

## ITER: M13-DISTANCE-GENERALIZATION-SCAN-s0 — 2026-06-17

- **Track**: screen / distance-generalization diagnostic, NON-CLAIM. major_axis=data_view/generalization_distance, mechanism_delta=fixed Arabidopsis train+calibration tested on close Arabidopsis-relative plus far rice.
- **Component exp_ids (verbatim)**: `M13-DISTANCE-GENERALIZATION-SCAN-s0`.
- **Execution mode**: submit-and-monitor via smart-sbatch/pre-submit gate. Slurm job `9019532` COMPLETED on private-teodoro-gpu `gpu034` in `10:26:59`.
- **Prerequisite sanity**: `reports/M13_FAILURE_SANITY/report.md` shows M12A rice failure is emission/calibration shift plus fragmentation: selected-bias predecode intergenic false-genic rate `0.2320` vs M11 pooled `0.0161`; rice oracle grid 0/3 valid under `FPR<=0.01` and `gene_count<=1.25`.
- **Species freeze**: `Arabidopsis lyrata subsp. lyrata` RefSeq `GCF_000004255.2` (`v.1.0`, Annotation Release 101) prepared under `data/m1_screen/arabidopsis_lyrata`; top 8 large scaffolds, 194.2Mb, check_data PASS. Diagnostic-only caveat: assembly level is Scaffold.
- **Hypothesis**: if fixed Arabidopsis-trained M9-L12 transfers to A. lyrata but not rice, failure is distance-limited; if it fails A. lyrata too, abandon single-species fixed-model generalization as main claim route.
- **Architecture/config**: same M12A trainer/calibrator; train/val species `arabidopsis_thaliana`; test species `arabidopsis_lyrata`, `oryza_sativa`; seed `0`; no test-label decode/hyperparameter tuning.
- **Code review**: PASS_WITH_WARNINGS host-self screen review; machine gates at `outputs/M13-DISTANCE-GENERALIZATION-SCAN/code_review_gate.json` and `outputs/M13-DISTANCE-GENERALIZATION-SCAN-s0/code_review_gate.json`.
- **Result summary**: negative for fixed single-species generalization. VAL Arabidopsis selected point is valid (`FPR=0.0085`, gbF1 `0.8965`, gene_count_ratio `1.004`), but TEST aggregate fails (`FPR=0.0340`, gbF1 `0.7415`, constrained gbF1 `0.0`, gene_count_ratio `1.616`). A. lyrata close plant also fails (`FPR=0.0355`, gbF1 `0.8167`, gene_count_ratio `1.358`), so failure is not just rice distance.
- **Tri-review consensus**: completed in combined M13/M14/M16 review, 3/3 quorum. Consensus: comparability/SOTA-freeze blocker first; do not tune or promote current fixed-model route.
- **Pivot decision**: `Comparability audit first`. Next primary EXP=`M17-SAMEPANEL-GENERALIZATION-BASELINES`; parallel local audit=`M17-PRETRAINING-OVERLAP-AUDIT`.
- **Links**: docs/06 `Result: M13-DISTANCE-GENERALIZATION-SCAN-s0`; config `configs/M13-DISTANCE-GENERALIZATION-SCAN.yaml`; panel `configs/m13_distance_generalization_panel.yaml`; sbatch `sbatch/M13-DISTANCE-GENERALIZATION-SCAN.sbatch`; output `outputs/M13-DISTANCE-GENERALIZATION-SCAN-s0`; log `outputs/fp_segnt_logs/M13DIST_9019532.out`.

## ITER: M14/M15-PARALLEL-DIAGNOSTICS — 2026-06-17

- **Track**: screen / publication-alignment diagnostics, NON-CLAIM. Parallel portfolio opened after user clarified GPU budget is sufficient and waiting only on M13 is unnecessary.
- **Component exp_ids (verbatim)**: `M14-ANIMAL-DISTANCE-NEGCTRL-s0`, `M15-GENERANNO-1P2B-CDS-PANEL-SCREEN`, `M15-GENERANNO-0P5B-BASE-PANEL-SCREEN`.
- **Execution mode**: submit-and-monitor via smart-sbatch/pre-submit gate. M14 job `9022700` COMPLETED on private-teodoro-gpu after initial 0s guard false-positive job `9022458`; M15 was first submitted as shared-gpu `9022457_[0-1]`, then cancelled before start and rerouted to private-teodoro as `9023295_[0-1%1]` with array concurrency 1. Both M15 arms COMPLETED.
- **M14 hypothesis**: if fixed Arabidopsis-trained M9-L12 fails both gallus and drosophila badly, this supports distance-limited or multi-species/domain-adaptation framing; animal outputs are negative-control diagnostics and not clean final-claim evidence.
- **M15 hypothesis**: if GENERanno 1.2B CDS-preview improves under a bounded clean-plant panel screen while 0.5B base remains collapsed, official CDS pretraining has useful signal but still needs our intron-aware head/schedule; if both fail, pretrained backbones alone do not solve the paper task.
- **Architecture/config**: M14 reuses M13/M12A NT-v2 L12 trainer and M11 val-only calibrator; train/val species `arabidopsis_thaliana`, test species `gallus_gallus`, `drosophila_melanogaster`. M15 reuses existing GENERanno LoRA 3-class trainer with bounded train/val windows on clean plants `{arabidopsis,rice}`.
- **Code review**: PASS_WITH_WARNINGS host-self screen reviews. New configs/sbatch passed `bash -n`, YAML parse, `py_compile`, and `pre_submit_gate`.
- **M14 result summary**: negative-control failure overall. TEST aggregate `FPR=0.0406`, gbF1 `0.5448`, constrained gbF1 `0.0`, gene_count_ratio `2.254`. Drosophila has acceptable FPR (`0.0085`) but low gbF1/gene coherence; gallus is severe FP/fragmentation (`FPR=0.1156`, gene_count_ratio `4.987`).
- **M15 result summary**: 1.2B CDS-preview is better but not guardrail-valid (`gbF1=0.8510`, `FPR=0.0258`, gene_count_ratio `1.177`); 0.5B base recovers from M12C collapse with more bounded training but remains worse (`gbF1=0.7623`, `FPR=0.0562`, gene_count_ratio `1.326`).
- **Tri-review/pivot status**: M14 consumed by combined M13/M14/M16 review; M15 consumed as GENERanno challenger context. Pivot chooses `Comparability audit first` before further model-architecture GPU spend.
- **Links**: docs/06 `Result: M14-ANIMAL-DISTANCE-NEGCTRL-s0` and `Result: M15-GENERANNO-LORA-PANEL-SCREEN`; configs `configs/M14-ANIMAL-DISTANCE-NEGCTRL.yaml`, `configs/M15-GENERANNO-LORA-PANEL-SCREEN.yaml`; sbatch `sbatch/M14-ANIMAL-DISTANCE-NEGCTRL.sbatch`, `sbatch/M15-GENERANNO-LORA-PANEL-SCREEN.sbatch`; outputs under `outputs/M14-ANIMAL-DISTANCE-NEGCTRL-s0` and `outputs/M15-GENERANNO-*PANEL-SCREEN`.

## ITER: M16-MULTISPECIES-TRAIN-DISTANCE-SCAN-s0 — 2026-06-18

- **Track**: screen / generalization-mechanism diagnostic, NON-CLAIM. Opened because GPU budget is sufficient and the user's central concern is whether fixed-model generalization fails because training used too few species.
- **Component exp_ids (verbatim)**: `M16-MULTISPECIES-TRAIN-DISTANCE-SCAN-s0`.
- **Execution mode**: submit-and-monitor via smart-sbatch/pre-submit gate; submitted after M13/M14/M15 were already running.
- **Hypothesis**: if adding rice to Arabidopsis training/calibration improves A. lyrata and/or animal transfer versus M13/M14, species diversity becomes a credible mainline axis; if not, simply adding one more species is insufficient and the route needs clade adapters, richer architecture, or a broader curated panel.
- **Architecture/config**: same M9-L12 NT-v2 top-12 unfreeze + 3-class FP-aware convLSTM head + M11 validation-only decode calibration. Train/calibrate species `arabidopsis_thaliana`, `oryza_sativa`; test species `arabidopsis_lyrata`, `gallus_gallus`, `drosophila_melanogaster`; seed `0`.
- **Code review**: PASS_WITH_WARNINGS host-self screen review; static checks passed; machine gate at `outputs/M16-MULTISPECIES-TRAIN-DISTANCE-SCAN-s0/code_review_gate.json`.
- **Result summary**: completed/result-logged. VAL mixed plant selected point is valid (`FPR=0.0087`, gbF1 `0.8310`, gene_count_ratio `1.035`), but TEST aggregate still fails guardrails (`FPR=0.0197`, gbF1 `0.5615`, constrained gbF1 `0.5615`, gene_count_ratio `1.326`). A. lyrata still fails FPR (`0.0306`), gallus improves versus M14 but remains fragmented (`gene_count_ratio=2.368`), and drosophila has excellent FPR (`0.0012`) but poor gbF1 (`0.4466`).
- **Tri-review consensus**: completed in combined M13/M14/M16 review, 3/3 quorum. Consensus: fixed single-species route is not claimable, M16 partial rescue is insufficient, and same-panel baseline/comparability must be closed first.
- **Pivot decision**: `Comparability audit first`. Next primary EXP=`M17-SAMEPANEL-GENERALIZATION-BASELINES`; parallel local audit=`M17-PRETRAINING-OVERLAP-AUDIT`.
- **Links**: docs/06 `Result: M16-MULTISPECIES-TRAIN-DISTANCE-SCAN-s0`; config `configs/M16-MULTISPECIES-TRAIN-DISTANCE-SCAN.yaml`; sbatch `sbatch/M16-MULTISPECIES-TRAIN-DISTANCE-SCAN.sbatch`; output `outputs/M16-MULTISPECIES-TRAIN-DISTANCE-SCAN-s0`; log `outputs/fp_segnt_logs/M16MULTI_9065776.out`.

## ITER: M17-SAMEPANEL-GENERALIZATION-BASELINES — 2026-06-19

- **Track**: comparability / publication-alignment screen, NON-CLAIM. major_axis=baseline_comparability, mechanism_delta=released ANNEVO/Tiberius/Helixer on the same M13/M16 diagnostic panel and evaluator.
- **Component exp_ids (verbatim)**: `M17-SAMEPANEL-GENERALIZATION-BASELINES-ANNEVO`, `M17-SAMEPANEL-GENERALIZATION-BASELINES-TIBERIUS`, `M17-SAMEPANEL-GENERALIZATION-BASELINES-HELIXER`. Slurm array `9119473_[0-2%3]`, all COMPLETED.
- **Hypothesis**: before spending more GPU on our model, determine whether the diagnostic panel itself is broadly hard for released gene callers or whether our fixed/adapted model is specifically behind.
- **Result summary**: released callers remain strong but with distinct tradeoffs. ANNEVO: gbF1 `0.9115`, FPR `0.0240`, gene_count_ratio `0.840`; Tiberius: gbF1 `0.8791`, FPR `0.0173`, gene_count_ratio `0.556`; Helixer: gbF1 `0.8797`, FPR `0.0526`, gene_count_ratio `0.931`.
- **Interpretation**: M17 refutes the weak excuse that the M13/M16 panel is simply impossible. Our M16 broad fixed model is not competitive with released callers on gbF1. The next main question is how to approach the external-caller tradeoff frontier without losing calibrated specificity/gene-count coherence.
- **Next action**: keep M18 parallel diagnostics running; after M18, tri-review/pivot over M17+M18 rather than tuning M9 in isolation.
- **Links**: docs/06 `Result: M17-SAMEPANEL-GENERALIZATION-BASELINES`; config `configs/M17-SAMEPANEL-GENERALIZATION-BASELINES.yaml`; sbatch `sbatch/M17-SAMEPANEL-GENERALIZATION-BASELINES.sbatch`; outputs `outputs/M17-SAMEPANEL-GENERALIZATION-BASELINES-*`; logs `outputs/M17-SAMEPANEL-GENERALIZATION-BASELINES/logs/`.

## ITER: M18-MULTICLADE-TRAIN-DIAGNOSTIC-s0 — 2026-06-19

- **Track**: screen / generalization-mechanism diagnostic, NON-CLAIM. major_axis=species_breadth, mechanism_delta=NT-v2 L12 train/calibrate on a broader supervised panel.
- **Component exp_ids (verbatim)**: `M18-MULTICLADE-TRAIN-DIAGNOSTIC-s0`. Slurm job `9123661` COMPLETED.
- **Hypothesis**: if adding Drosophila to Arabidopsis+rice training/calibration materially improves A. lyrata, gallus, and yeast transfer, then fixed-model failure is partly a species-coverage problem; if only nearby/domain-related transfer improves, the paper route needs clade/domain adaptation or narrower claims.
- **Architecture/config**: NT-v2 top-12 unfreeze + 3-class FP-aware convLSTM head + validation-only decode calibration. Train/val species `arabidopsis_thaliana`, `oryza_sativa`, `drosophila_melanogaster`; test species `arabidopsis_lyrata`, `gallus_gallus`, `saccharomyces_cerevisiae`; seed `0`; train/val window caps `8192/4096`.
- **Code review**: PASS_WITH_WARNINGS host-self screen review; machine gate at `outputs/M18-MULTICLADE-TRAIN-DIAGNOSTIC-s0/code_review_gate.json`.
- **Result summary**: semantic-success PASS, validate=`progress` screen non-claim. TEST aggregate specificity `0.9556`, FPR `0.0444`, gbF1 `0.6170`, constrained gbF1 `0.0`, gene_count_ratio `1.572`. A. lyrata improves versus M13 (`FPR=0.0184`, gbF1 `0.7427`, gene_count_ratio `1.078`), but gallus fails severely (`FPR=0.1359`, gene_count_ratio `3.994`) and yeast under-calls genes (`0.591x` reference).
- **Interpretation**: broader supervised species coverage helps close/domain-related plant transfer but does not yield a broad fixed eukaryotic caller. This run is still far behind M17 released callers on broad-panel gbF1 and practical gene-count behavior.
- **Tri-review/pivot status**: deferred until the two parallel M18 GENERanno siblings finish, then M17+M18 should be reviewed together.
- **Links**: docs/06 `Result: M18-MULTICLADE-TRAIN-DIAGNOSTIC-s0`; config `configs/M18-MULTICLADE-TRAIN-DIAGNOSTIC.yaml`; sbatch `sbatch/M18-MULTICLADE-TRAIN-DIAGNOSTIC.sbatch`; output `outputs/M18-MULTICLADE-TRAIN-DIAGNOSTIC-s0`; log `outputs/fp_segnt_logs/M18MULTI_9123661.out`.

## ITER: M18-GENERANNO-0P5B-SPEC-OBJECTIVE-s0 — 2026-06-19

- **Track**: screen / GENERanno challenger control, NON-CLAIM. major_axis=pretraining_specialization, mechanism_delta=0.5B base with same stronger FP objective as 1.2B M18 sibling.
- **Component exp_ids (verbatim)**: `M18-GENERANNO-0P5B-SPEC-OBJECTIVE-s0`. Slurm job `9131867` COMPLETED.
- **Hypothesis**: if the 0.5B base model can be rescued by stronger FP objective and the same 3-class LoRA head, generic GENERanno pretraining may be sufficient; if not, official CDS-annotator adaptation is likely the relevant signal source.
- **Architecture/config**: `GenerTeam/GENERanno-eukaryote-0.5b-base` masked-LM backbone, k=1, LoRA r=8, 3-class convLSTM head, `fp_lambda=2.5`, `min_cds_len=90`, train/val/test on clean plants.
- **Result summary**: semantic-success PASS, validate=`progress` screen non-claim. Aggregate specificity `0.9033`, FPR `0.0967`, gbF1 `0.6561`, constrained gbF1 `0.0`, gene_count_ratio `1.617`. Arabidopsis is usable (`FPR=0.0182`, gbF1 `0.7813`), but rice fails badly (`FPR=0.1239`, gbF1 `0.4981`, gene_count_ratio `2.482`).
- **Interpretation**: stronger FP objective does not rescue 0.5B base; compared with M15 0.5B it is worse on aggregate. Treat as ablation evidence, not a scale candidate.
- **Tri-review/pivot status**: deferred until 1.2B M18 sibling finishes, then M17+M18 should be reviewed together.
- **Links**: docs/06 `Result: M18-GENERANNO-0P5B-SPEC-OBJECTIVE-s0`; config `configs/M18-GENERANNO-0P5B-SPEC-OBJECTIVE.yaml`; sbatch `sbatch/M18-GENERANNO-0P5B-SPEC-OBJECTIVE.sbatch`; output `outputs/M18-GENERANNO-0P5B-SPEC-OBJECTIVE-s0`; log `outputs/fp_segnt_logs/M18GENBASE_9131867.out`.

## ITER: M18-GENERANNO-1P2B-SPEC-OBJECTIVE-s0 — 2026-06-19

- **Track**: screen / GENERanno challenger, NON-CLAIM. major_axis=pretraining_specialization+objective, mechanism_delta=1.2B CDS-preview backbone with stronger FP objective and longer minimum CDS constraint.
- **Component exp_ids (verbatim)**: `M18-GENERANNO-1P2B-SPEC-OBJECTIVE-s0`. Slurm job `9122868` COMPLETED in `12:07:01`.
- **Hypothesis**: if M15 1.2B failed mainly because the schedule/objective did not suppress intergenic FP enough, stronger FP pressure should move the CDS-preview LoRA route into the FPR/gene-count guardrail region; if not, GENERanno remains only an ablation.
- **Architecture/config**: `GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview` token-classification backbone, k=6, LoRA r=8, 3-class convLSTM head, `fp_lambda=2.5`, `min_cds_len=90`, clean-plant panel `{arabidopsis,rice}`, seed `0`.
- **Result summary**: semantic-success PASS, validate=`progress` screen non-claim. Aggregate specificity `0.9929`, FPR `0.0071`, gbF1/constrained gbF1@0.01 `0.8494`, macro specificity `0.9943`, gene_count_ratio `0.864`. Arabidopsis is very strong (`FPR=0.0027`, gbF1 `0.9144`); rice is guardrail-valid but lower gbF1 (`FPR=0.0087`, gbF1 `0.7542`, gene_count_ratio `1.036`).
- **Interpretation**: this reverses the prior "GENERanno LoRA only fragments" result. Under the same stronger objective, 1.2B CDS-preview is far better than 0.5B base, so official CDS specialization is a real mechanism. It is still not a claim candidate because same-panel released callers have higher clean-plant gbF1 and pretraining overlap is unresolved.
- **Tri-review/pivot status**: pending combined M17+M18 review; likely candidate for next portfolio via raw-score validation-only calibration and two seeds.
- **Links**: docs/06 `Result: M18-GENERANNO-1P2B-SPEC-OBJECTIVE-s0`; config `configs/M18-GENERANNO-1P2B-SPEC-OBJECTIVE.yaml`; sbatch `sbatch/M18-GENERANNO-1P2B-SPEC-OBJECTIVE.sbatch`; output `outputs/M18-GENERANNO-1P2B-SPEC-OBJECTIVE-s0`; log `outputs/fp_segnt_logs/M18GENSPEC_9122868.out`.

## ITER: M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS — 2026-06-19

- **Track**: screen / Track-B-preflight, NON-CLAIM. major_axis=backbone+calibration, mechanism_delta=promote GENERanno 1.2B CDS-preview route with raw-score persistence and two seeds after M17+M18 pivot.
- **Component exp_ids (verbatim)**: `M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s0`, `M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s1`. Slurm array `9141356_[0-1%2]` COMPLETED on `gpu035`.
- **Hypothesis**: M18 1.2B's guardrail-valid single-seed result is stable enough to justify route promotion, and saved VAL/TEST raw logits will allow validation-only calibration and FPR sensitivity analysis without test leakage.
- **Architecture/config**: same successful M18 1.2B CDS-preview + LoRA r=8 + 3-class FP-aware head + `fp_lambda=2.5`, `min_cds_len=90`, clean plants `{arabidopsis,rice}`, now with `--save-raw-scores` and seeds `0/1`.
- **Code review**: PASS_WITH_WARNINGS host-self screen gates for top-level and seed arm exp_ids; static checks passed (`bash -n`, YAML parse, `py_compile`, `pre_submit_gate`).
- **Parallel cohort**: `M19-GENERANNO-PROVENANCE-AUDIT` completed locally. Public HF/GitHub sources do not expose a complete eukaryotic CDS-preview species/accession exclusion list, so GENERanno remains `overlap_unknown` for claim purposes.
- **Result summary**: completed/result-logged. Both seeds are semantically successful and remain aggregate FPR-valid/gene-count sane after validation-only calibration. Raw decode: s0 `gbF1=0.8390`, FPR `0.0088`, gene_count_ratio `0.967`; s1 `gbF1=0.8593`, FPR `0.0059`, gene_count_ratio `0.805`. Calibrated: s0 selected `b2p0_mcl60_mfg20` -> `gbF1=0.8421`, FPR `0.0083`, gene_count_ratio `1.083`; s1 selected `b0p0_mcl60_mfg20` -> `gbF1=0.8815`, FPR `0.0065`, gene_count_ratio `0.830`.
- **Comparability summary**: paper-facing table refreshed at `reports/M19-COMPARABILITY-EVIDENCE/comparison_tables.md`. M19 1.2B is stable and much stronger than the M18 0.5B base control, but clean-plant released callers still define the high-gbF1 frontier: Tiberius/ANNEVO/Helixer `0.922-0.927` gbF1. Tiberius is the closest practical comparator because it also passes aggregate FPR<=0.01, but it under-calls genes (`0.628x`).
- **Tri-review consensus**: M17+M18 combined review 2/3 degraded promoted GENERanno 1.2B into M19. M19 combined review then reached 3/3 quorum: M19 is stable and useful adaptation/comparability evidence, but current Arabidopsis/rice evidence cannot claim because GENERanno provenance is `overlap_unknown` and released callers still lead gbF1.
- **Pivot decision**: `Comparability audit first`: open `M20-CLAIM-CLEAN-PANEL-FREEZE` as the primary claim gate; allow `M20-STRUCTURED-DECODER-IMPL` and `M20-SOTA-ERROR-ANALYSIS` in parallel, but no claim-grade GENERanno GPU on the current panel.
- **Links**: docs/07 `Tri-Review: M17+M18 Combined Evidence` and `Tri-Review: M19-GENERANNO-COMBINED-DECISION`; docs/08 `Pivot Decision: M17+M18 Combined Evidence` and `Pivot Decision: M19-GENERANNO-COMBINED-DECISION`; config `configs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS.yaml`; sbatch `sbatch/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS.sbatch`; audit `refs/dossiers/m19_generanno_provenance_audit.md`; outputs `outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s{0,1}`.

## ITER: M23-NTV2-CLEAN-TRANSFER-s0 — 2026-07-01

- **Track**: screen / clean-provenance NT-v2 transfer-learning check, NON-CLAIM. major_axis=backbone provenance, mechanism_delta=return to public clean-provenance NT-v2 after M22 negative while forbidding M22 `gb_tversky`, trained CRF, and raw-score calibration.
- **Component exp_ids (verbatim)**: `M23-NTV2-CLEAN-TRANSFER-s0`. Slurm job `9854668` COMPLETED on private `gpu034` in `19:25:54`.
- **Hypothesis**: if clean-provenance NT-v2 direct transfer can recover a better gbF1/FPR tradeoff than the current GENERanno branch, it may become the claim-facing route despite GENERanno's stronger adaptation metrics being provenance-blocked.
- **Architecture/config**: `InstaDeepAI/nucleotide-transformer-v2-500m-multi-species`, top-12 unfreeze, 3-class convLSTM head, `fp_aware`, constrained postproc, clean plants `{arabidopsis_thaliana, oryza_sativa}`, seed `0`; no raw-score saving or calibration.
- **Code review**: PASS_WITH_WARNINGS host-self screen gate. Separate Codex review was attempted but blocked by bwrap namespace failure. Static checks, split checks, and `pre_submit_gate.py` passed.
- **Result summary**: semantic-success PASS, validate=`progress` screen non-claim. Aggregate gbF1 `0.8427`, FPR `0.01673`, specificity `0.98327`, macro specificity `0.98050`, gene_count_ratio `0.867`. Metrics exactly match historical `M10-M9L12-CLEANPLANTS-s0`, confirming this is a direct-transfer replication/claim-route checkpoint rather than a new performance improvement.
- **Interpretation**: clean-provenance NT-v2 remains cleaner for provenance but weaker than the best adapted GENERanno metric route (`M19 s1`: gbF1 `0.8815`, FPR `0.0065`). Direct NT-v2 does not pass hard FPR<=0.01 and should not be rerun as more seeds without a structural change.
- **Tri-review/pivot status**: pending if deciding next GPU route; recommended combined review over M22 negative + M23 clean NT-v2 result.
- **Links**: docs/06 `Result: M23-NTV2-CLEAN-TRANSFER-s0`; config `configs/M23-NTV2-CLEAN-TRANSFER.yaml`; sbatch `sbatch/M23-NTV2-CLEAN-TRANSFER-s0.sbatch`; output `outputs/M23-NTV2-CLEAN-TRANSFER-s0`; log `outputs/fp_segnt_logs/M23NTV2S0_9854668.out`.
