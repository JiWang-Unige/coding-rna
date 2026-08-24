# Code Review Log / 代码审前闸记录

> 由 `/code-review-gate` 维护。每次实现训练/评估/数据/配置/job 脚本后、提交真实训练前，必须记录一次审查。目标是把 label、metric、split、输出路径等会让结果作废的问题挡在运行前。

## Policy
- `full` / `scale` / claim-candidate 必须有 `PASS` 或 `PASS_WITH_WARNINGS`。
- `screen` 若改了共享训练/eval/data 代码，也必须有审查。
- `smoke` 可只做 host checklist，但若 smoke 之后直接提交 screen/full，仍需本 gate。
- `BLOCKED` 不得提交；用户强行豁免必须写 `WAIVED_BY_USER`。

## Review Entries

### Code Review Gate: M21-GENERANNO-1P2B-CRF-SCREEN 2026-06-21
- Reviewer mode: host-read-only.
- Scope: screen/non-claim M21 two-seed GENERanno 1.2B CRF decoder array; trainer path was smoke-proven in M20, this gate reviews the screen config/sbatch and current trainer/evaluator hashes.
- Verdict: PASS_WITH_WARNINGS
- Linked implementation: `src/foundation_probe/train_generanno_lora_3class.py`, `configs/M21-GENERANNO-1P2B-CRF-SCREEN.yaml`, `sbatch/M21-GENERANNO-1P2B-CRF-SCREEN.sbatch`
- Linked evaluator contract: docs/19_evaluator_contract.md
- Linked baseline reproduction: docs/20_baseline_reproduction.md

#### Blockers
- [x] No open blocker for screen/non-claim submission.

#### Warnings
- Review independence is `host_self`; acceptable for screen continuation but not sufficient for full/scale claim.
- GENERanno provenance remains `overlap_unknown`; all M21 outputs are adaptation/mechanism evidence only.
- CRF full-screen runtime is unmeasured and may exceed M19 due Viterbi decoding overhead; sbatch uses private partition and extended 71:50 walltime.

#### Confirmed OK
- M21 changes the structural axis only: same clean-plant panel, same M19 budget/objective, adds `--decoder crf --crf-aux-ce 1.0`.
- No smoke prediction/decode cap is present in the screen sbatch.
- `--decoder none` remains the trainer default, preserving M19 behavior for other runs.
- Evaluator uses `--span-mode cds`; aggregator emits validate_goal-required fields.
- Output roots are seed-isolated: `outputs/M21-GENERANNO-1P2B-CRF-SCREEN-s0` and `...-s1`.
- Static checks passed: `py_compile`, `bash -n`, YAML parse, and synthetic `LinearChainCRFVec` NLL/backprop/Viterbi test.
- `pre_submit_gate.py` passes for top-level and both seed exp_ids with non-empty reviewed file hashes.

#### Fix / waiver record
- No waiver. Machine gates written to `outputs/M21-GENERANNO-1P2B-CRF-SCREEN/code_review_gate.json` and `outputs/M21-GENERANNO-1P2B-CRF-SCREEN-s{0,1}/code_review_gate.json`.

### Code Review Gate: M21-GENERANNO-1P2B-CRF-SCREEN-s1-shared 2026-06-21
- Reviewer mode: host-read-only.
- Scope: shared-A100 wrapper reroute for `M21-GENERANNO-1P2B-CRF-SCREEN-s1`; training/evaluator/config files unchanged from the main M21 screen gate.
- Verdict: PASS_WITH_WARNINGS for screen/non-claim reroute.
- Linked implementation: `sbatch/M21-GENERANNO-1P2B-CRF-SCREEN-s1-shared.sbatch`, `configs/M21-GENERANNO-1P2B-CRF-SCREEN.yaml`, `src/foundation_probe/train_generanno_lora_3class.py`, `scripts/eval_gene_body_mask.py`, `scripts/aggregate_gene_body_metrics.py`, `scripts/pre_submit_gate.py`.

#### Blockers
- [x] No open blockers for reroute. Private seed1 task `9259965_1` was pending with estimated start `2026-06-24T12:19:00`; it was cancelled before start.
- [x] The wrapper changes only scheduler resources/job naming: shared-gpu, A100 40GB, 11:50 walltime, fixed `SEED=1`.
- [x] Static checks pass: `bash -n` for the wrapper sbatch and `pre_submit_gate.py` for both the wrapper exp_id and actual seed exp_id.

#### Warnings
- Review independence is `host_self`, acceptable only for screen/non-claim reroute.
- M21 remains non-claim because GENERanno provenance is `overlap_unknown`.

#### Confirmed OK
- Shared seed1 launched as Slurm job `9260587` on `gpu020`; seed0 continues as private job `9259965_0`.

#### Fix / waiver record
- No waiver. Machine gates written to `outputs/M21-GENERANNO-1P2B-CRF-SCREEN-s1/code_review_gate.json` and `outputs/M21-GENERANNO-1P2B-CRF-SCREEN-s1-shared/code_review_gate.json`.

### Code Review Gate: M21-GENERANNO-1P2B-CRF-DIAG-s1 2026-06-21
- Reviewer mode: host-read-only.
- Scope: runtime diagnostic only for M21 CRF throughput; not a model-quality run and not result-log eligible.
- Verdict: PASS_WITH_WARNINGS for smoke diagnostic.
- Linked implementation: `src/foundation_probe/train_generanno_lora_3class.py`, `sbatch/M21-GENERANNO-1P2B-CRF-DIAG-s1.sbatch`, `configs/M21-GENERANNO-1P2B-CRF-SCREEN.yaml`, `scripts/eval_gene_body_mask.py`, `scripts/aggregate_gene_body_metrics.py`, `scripts/pre_submit_gate.py`.

#### Blockers
- [x] No open blockers for diagnostic. Trainer patch adds optional `--pretokenize-windows` and `--train-progress-every`; default behavior remains unchanged for existing CLI calls.
- [x] Static checks pass: `python3 -m py_compile src/foundation_probe/train_generanno_lora_3class.py`, `bash -n sbatch/M21-GENERANNO-1P2B-CRF-DIAG-s1.sbatch`, and `pre_submit_gate.py`.

#### Warnings
- Diagnostic briefly landed on the same visible GPU as shared seed1; it was cancelled intentionally after 7 progress batches to avoid interfering with the main screen.
- No metrics or claim should be inferred from this run.

#### Confirmed OK
- Diagnostic stdout produced batch-level progress (`batch 1/8` through `7/8`) and estimated throughput at roughly `8.4-9.8 sec/batch` on A100.

#### Fix / waiver record
- No waiver. Machine gate written to `outputs/M21-GENERANNO-1P2B-CRF-DIAG-s1/code_review_gate.json`.

### Code Review Gate: M21-GENERANNO-1P2B-CRF-SCREEN-s1-opt 2026-06-21
- Reviewer mode: host-read-only.
- Scope: private rescue/hedge wrapper for `M21-GENERANNO-1P2B-CRF-SCREEN-s1` after shared-A100 walltime risk became clear.
- Verdict: PASS_WITH_WARNINGS for screen/non-claim rescue.
- Linked implementation: `sbatch/M21-GENERANNO-1P2B-CRF-SCREEN-s1-opt.sbatch`, `configs/M21-GENERANNO-1P2B-CRF-SCREEN.yaml`, `src/foundation_probe/train_generanno_lora_3class.py`, `scripts/eval_gene_body_mask.py`, `scripts/aggregate_gene_body_metrics.py`, `scripts/validate_goal.py`, `scripts/pre_submit_gate.py`.

#### Blockers
- [x] No open blockers for screen/non-claim submission.
- [x] Output root is unique: `outputs/M21-GENERANNO-1P2B-CRF-SCREEN-s1-opt`; it will not overwrite the original shared seed1 output.
- [x] Rescue is not a new model claim: same seed/objective/data/budget as M21 seed1, with only runtime instrumentation/throughput changes (`pretokenize`, progress logging, batched CRF decode).
- [x] Static checks pass: `python3 -m py_compile src/foundation_probe/train_generanno_lora_3class.py scripts/eval_gene_body_mask.py scripts/aggregate_gene_body_metrics.py scripts/validate_goal.py` and `pre_submit_gate.py`.

#### Warnings
- Review independence is `host_self`, acceptable only for screen/non-claim rescue.
- The original shared job `9260587` may still finish; if both seed1 variants complete, treat `s1-opt` as the runtime-rescue duplicate and compare outputs before choosing which enters the M21 result summary.

#### Confirmed OK
- Machine pre-submit gate passed for `M21-GENERANNO-1P2B-CRF-SCREEN-s1-opt`.
- Job submitted as Slurm `9298703` on `private-teodoro-gpu`.

#### Fix / waiver record
- No waiver. Machine gate written to `outputs/M21-GENERANNO-1P2B-CRF-SCREEN-s1-opt/code_review_gate.json`.

### Code Review Gate: M21-GENERANNO-1P2B-CRF-SCREEN-s1-fastval 2026-06-22
- Reviewer mode: host-read-only.
- Scope: shared-GPU fast-validation rescue for `M21-GENERANNO-1P2B-CRF-SCREEN-s1` after `s1-opt` showed healthy but slow CRF validation with no progress printing.
- Verdict: PASS_WITH_WARNINGS for screen/non-claim rescue.
- Linked implementation: `src/foundation_probe/train_generanno_lora_3class.py`, `sbatch/M21-GENERANNO-1P2B-CRF-SCREEN-s1-fastval.sbatch`, `configs/M21-GENERANNO-1P2B-CRF-SCREEN.yaml`, `scripts/eval_gene_body_mask.py`, `scripts/aggregate_gene_body_metrics.py`, `scripts/validate_goal.py`, `scripts/pre_submit_gate.py`.

#### Blockers
- [x] No open blockers for screen/non-claim submission.
- [x] Training semantics unchanged: seed, species, train/val window caps, optimizer, loss, CRF decoder, post-processing, and evaluator match M21 seed1.
- [x] Runtime delta is validation-only batching/progress: `--eval-batch 16 --eval-progress-every 10`; training DataLoader remains `--batch-size 1`.
- [x] Output root is unique: `outputs/M21-GENERANNO-1P2B-CRF-SCREEN-s1-fastval`; it will not overwrite original shared seed1 or `s1-opt`.
- [x] Static checks pass: `python3 -m py_compile src/foundation_probe/train_generanno_lora_3class.py`, `bash -n sbatch/M21-GENERANNO-1P2B-CRF-SCREEN-s1-fastval.sbatch`, and `pre_submit_gate.py`.

#### Warnings
- Review independence is `host_self`, acceptable only for screen/non-claim rescue.
- If both `s1-opt` and `s1-fastval` complete, treat them as duplicate runtime variants of the same seed and compare metrics before result-log; do not count them as independent seeds.

#### Confirmed OK
- Machine pre-submit gate passed for `M21-GENERANNO-1P2B-CRF-SCREEN-s1-fastval`.
- Private pending attempt `9343568` and shared over-time-limit attempt `9343632` were cancelled before start; active shared queue submission is job `9343635`.

#### Fix / waiver record
- No waiver. Machine gate written to `outputs/M21-GENERANNO-1P2B-CRF-SCREEN-s1-fastval/code_review_gate.json`.

### Code Review Gate: M20-STRUCTURED-DECODER-IMPL-SMOKE 2026-06-21
- Reviewer mode: separate Codex attempted first, then host-read-only fallback.
- Scope: smoke-only implementation of optional `--decoder crf` in `src/foundation_probe/train_generanno_lora_3class.py`, plus `configs/M20-STRUCTURED-DECODER-IMPL.yaml` and `sbatch/M20-STRUCTURED-DECODER-IMPL-SMOKE.sbatch`.
- Verdict: PASS_WITH_WARNINGS
- Linked implementation: `src/foundation_probe/train_generanno_lora_3class.py`, `configs/M20-STRUCTURED-DECODER-IMPL.yaml`, `sbatch/M20-STRUCTURED-DECODER-IMPL-SMOKE.sbatch`
- Linked evaluator contract: docs/19_evaluator_contract.md
- Linked baseline reproduction: docs/20_baseline_reproduction.md

#### Blockers
- [x] No open blocker for smoke. No full/screen/claim submission is authorized by this gate.

#### Warnings
- Separate Codex review was attempted with `codex exec --sandbox read-only`, but the reviewer process could not read any file because every shell command failed with `bwrap: Creating new namespace failed: No space left on device`; it correctly returned `BLOCKED` due unreadable inputs.
- Final machine gate therefore uses `host_self`, which is weaker than `separate_codex`. This is acceptable only because profile is `smoke`.
- The new CRF path is only smoke-gated. Any claim-grade GENERanno CRF run still needs a fresh review after the final config is frozen.

#### Confirmed OK
- Default behavior is preserved: `--decoder none` remains the parser default and keeps emission argmax plus optional constrained post-processing.
- CRF training path uses the existing 3-class labels, `mask = y != -100`, per-token CRF NLL, and auxiliary CE/FP-aware loss.
- Prediction path uses Viterbi only when `decoder is not None`; tail bases not covered by full windows remain intergenic, matching the existing full-window prediction contract.
- Smoke sbatch uses `--span-mode cds`, aggregates via `scripts/aggregate_gene_body_metrics.py`, writes `outputs/M20-STRUCTURED-DECODER-IMPL-SMOKE`, and calls `pre_submit_gate.py` with the matching exp_id before training.
- Static checks passed: `py_compile`, `bash -n`, YAML parse, trainer help, and a synthetic `LinearChainCRFVec` NLL/backprop/Viterbi test.
- `pre_submit_gate.py` consumes the machine gate with non-empty reviewed file hashes.

#### Fix / waiver record
- No waiver. Machine gate written to `outputs/M20-STRUCTURED-DECODER-IMPL-SMOKE/code_review_gate.json`.
- Post-review smoke runtime fix: first job `9246005` was manually cancelled after 10 minutes because smoke still attempted full-chromosome prediction (`3818` windows). Added default-off `--limit-predict-windows` to the trainer, changed the smoke exp_id/output to `M20-STRUCTURED-DECODER-IMPL-SMOKE2`, removed `--save-raw-scores`, and capped prediction to 8 windows. Static checks were repeated (`py_compile`, `bash -n`, YAML parse, trainer help). Refreshed machine gate: `outputs/M20-STRUCTURED-DECODER-IMPL-SMOKE2/code_review_gate.json`; the submit-facing Slurm file is `sbatch/M20-STRUCTURED-DECODER-IMPL-SMOKE2.sbatch` so hook-derived exp_id matches the gate.
- Second smoke runtime fix: job `9249307` was manually cancelled after 6 minutes because CRF decoding still iterated over the full seqid even though scoring was capped to 8 windows. The trainer now applies the same default-off `--limit-predict-windows` cap to CRF Viterbi decoding. Submit-facing exp_id is `M20-STRUCTURED-DECODER-IMPL-SMOKE3`; gate refreshed at `outputs/M20-STRUCTURED-DECODER-IMPL-SMOKE3/code_review_gate.json`.

### Code Review Gate: <exp_id> <date>
- Reviewer mode: code-plan-reviewer / host-read-only / cli-advisory
- Scope:
- Verdict: PASS / PASS_WITH_WARNINGS / BLOCKED / WAIVED_BY_USER
- Linked implementation:
- Linked evaluator contract: docs/19_evaluator_contract.md
- Linked baseline reproduction: docs/20_baseline_reproduction.md

#### Blockers
- [ ] 

#### Warnings
- 

#### Confirmed OK
- 

#### Fix / waiver record
- 

### Code Review Gate: M10 dual-track submit set 2026-06-14
- Reviewer mode: host-read-only
- Scope: `M10-M9L12-CLEANPLANTS-s{0,1,2}`, `M10-GENERANNO-LORA-3C-SMOKE`, prepared `M10-GENERANNO-LORA-3C-s0`
- Verdict: PASS_WITH_WARNINGS
- Linked implementation: `src/foundation_probe/train_generanno_lora_3class.py`, `src/foundation_probe/train_unfreeze_backbone.py`, M10 configs and sbatch files
- Linked evaluator contract: docs/19_evaluator_contract.md
- Linked baseline reproduction: docs/20_baseline_reproduction.md

#### Blockers
- [x] No open blockers. Fixed pre-submit partition issue before gate: replaced nonexistent `private-dpnc-gpu` with live Baobab partitions; LoRA smoke now explicitly requests A100 40GB.

#### Warnings
- `host_self` review independence is weaker than a separate external CLI review; acceptable for screen/smoke handoff here, not for full/scale claim.
- GENERanno LoRA 1.2B PEFT forward and memory behavior still require the bounded GPU smoke to prove runtime viability.
- All M10 outputs are screen/non-claim because published SOTA remains draft in `ACTIVE_GOAL.json`.

#### Confirmed OK
- Python files compile under `generanno`; M10 sbatch files pass `bash -n`.
- GENERanno tokenizer/config alignment verified: `6144 bp -> 1024 tokens`, `k=6`; M9 `2046 bp -> 341 tokens`, `k=6`.
- M10 clean-plant data split gate passed: train/val/test `seqid` groups disjoint for arabidopsis + rice.
- Mainline reuses already validated M9 trainer with pooled species support and evaluates arabidopsis/rice separately before base-weighted + macro aggregation.
- GENERanno LoRA trainer discards the released binary CDS head, freezes base encoder weights, trains LoRA adapters plus the 3-class FP-aware head, writes standard CDS GFF + eval_subsets.

#### Fix / waiver record
- No waiver. Machine gates written to `outputs/<exp_id>/code_review_gate.json` with reviewed file hashes for stale-code detection.
- Post-run fix 2026-06-15: patched all M10 sbatch scripts so `validate_goal.py --run-status` receives the STATUS file path and runs after `STATUS=COMPLETED` is written. This fixes result validation only; training/evaluation semantics unchanged. `bash -n` passes for `sbatch/M10-M9L12-CLEANPLANTS.sbatch`, `sbatch/M10-GENERANNO-LORA-3C-SMOKE.sbatch`, and `sbatch/M10-GENERANNO-LORA-3C.sbatch`.

### Code Review Gate: M11-L12-SPEC-CALIBRATION 2026-06-15
- Reviewer mode: external_cli (`claude` + separate `codex exec`), adversarial pre-submit pack
- Scope: `M11-L12-SPEC-CALIBRATION-s{0,1,2}` screen and `M11-L12-SPEC-CALIBRATION-SMOKE`
- Verdict: PASS_WITH_WARNINGS
- Linked implementation: `src/foundation_probe/train_unfreeze_backbone.py`, `scripts/experiments/M11-L12-SPEC-CALIBRATION/calibrate_decode.py`, `configs/M11-L12-SPEC-CALIBRATION.yaml`, `sbatch/M11-L12-SPEC-CALIBRATION*.sbatch`
- Linked evaluator contract: docs/19_evaluator_contract.md
- Linked baseline reproduction: docs/20_baseline_reproduction.md

#### Blockers
- [x] No open blockers. Reviewers found no label/split/metric-schema/test-set-calibration/output-overwrite issue that invalidates submission.

#### Warnings
- `validate_goal.py ... || true` is intentional for screen/smoke, but result-log must read `metrics/validate_goal.json.status`; Slurm exit or `STATUS=COMPLETED` alone is not semantic success.
- Calibration fallback may select the least-bad validation operating point if no candidate reaches `target_fpr=0.01`; if test FPR remains >0.01, M11 is a failed calibration attempt, not a passed FPR gate.
- Final metrics are written with `--profile screen`; for the M11 scientific question, report `intergenic_guardrail_pass_at_0.01` and `constrained_gene_body_F1_at_0.01`, not only screen `constrained_gene_body_F1`.
- The NT-v2 trainer predicts only full windows; short tail bases remain intergenic. This is inherited behavior and must be mentioned in result-log if relevant.

#### Confirmed OK
- Train/val/test split uses deterministic seqid groups; data gate passed with train=11, val=2, test=2 and no group overlap.
- Calibration dataflow is VAL-only: sweep raw `val_*` scores and val labels, select one operating point, then apply once to `test_*` scores.
- `metrics/metrics.json` is produced by `aggregate_gene_body_metrics.py` and is readable by `validate_goal.py`.
- Output dirs are exp-id isolated: `outputs/M11-L12-SPEC-CALIBRATION-s{0,1,2}` and `outputs/M11-L12-SPEC-CALIBRATION-SMOKE`.
- Static checks passed: `python3 -m py_compile`, `bash -n`, and a CPU toy calibration E2E.

#### Fix / waiver record
- No waiver. Machine gates written to `outputs/M11-L12-SPEC-CALIBRATION*/code_review_gate.json` with reviewed file hashes. `pre_submit_gate.py` must pass before submission.
- Post-review resource fix: smoke sbatch now requests `gpu:nvidia_geforce_rtx_3090:1` instead of generic `gpu:1`, preventing accidental placement on low-VRAM shared GPUs. Logic unchanged; `bash -n` and `pre_submit_gate.py --exp-id M11-L12-SPEC-CALIBRATION-SMOKE` pass.
- Post-smoke submit fix: main calibration grid widened from bias `0..1.5` / min CDS `{60,90}` to bias `0..4.0` / min CDS `{60,90,120}` after smoke showed the original grid cannot approach FPR `0.01`; strict selection constraints remain `gbF1>=0.70` and `gene_count<=1.25`. Added prediction progress logging only. `py_compile`, `bash -n`, and YAML parse pass before main submission.

### Code Review Gate: M12C-GENERANNO-FAIR-CHALLENGER-SMOKE 2026-06-17
- Reviewer mode: attempted separate_codex read-only, then host-read-only fallback.
- Scope: smoke-only array with `M12C-GENERANNO-1P2B-CDS-SMOKE` and `M12C-GENERANNO-0P5B-BASE-SMOKE`.
- Verdict: PASS_WITH_WARNINGS for smoke only.
- Linked implementation: `src/foundation_probe/train_generanno_lora_3class.py`, `configs/M12C-GENERANNO-FAIR-CHALLENGER-SMOKE.yaml`, `sbatch/M12C-GENERANNO-FAIR-CHALLENGER-SMOKE.sbatch`, `configs/m12_publication_panel.yaml`.
- Linked evaluator contract: docs/19_evaluator_contract.md.
- Linked baseline reproduction: docs/20_baseline_reproduction.md.

#### Blockers
- [x] separate_codex returned `BLOCKED` because its read-only sandbox failed with `bwrap: Creating new namespace failed: No space left on device`, so it could not read local files. This is an independence/tooling blocker, not a semantic bug in the changed files. Fallback host review is recorded and only permits smoke, not screen/full.
- [x] Masked-LM backbone contract was unproven. Fixed by adding a pre-training dummy forward probe that requires `backbone(...).last_hidden_state` and verifies `hidden_tokens * k == window`.
- [x] Token/base alignment could silently drift. Fixed by adding `_tokenize_window()` and using it in dataset and prediction paths; every full window now must satisfy `len(input_ids) * k == window`.
- [x] STATUS/validate ordering was too weak. Fixed sbatch so validate JSON is produced and atomically moved before writing `STATUS=COMPLETED`; validate exit code is recorded.
- [x] Output overwrite risk. Fixed sbatch to refuse non-empty arm output directories except pre-submit gate artifacts.

#### Warnings
- Review independence is downgraded to `host_self` after separate_codex sandbox ENOSPC. This is acceptable only for bounded smoke; screen/full requires a successful independent review.
- 0.5b-base uses k=1 and window 1024 for smoke to reduce attention length; this is not a full fairness setting yet.
- GENERanno overlap remains `unknown`; M12C results are mechanism/control evidence, not clean claim evidence.

#### Confirmed OK
- Static checks pass: `python3 -m py_compile src/foundation_probe/train_generanno_lora_3class.py`; `bash -n sbatch/M12C-GENERANNO-FAIR-CHALLENGER-SMOKE.sbatch`; YAML parse for `configs/M12C-GENERANNO-FAIR-CHALLENGER-SMOKE.yaml` and `configs/m12_publication_panel.yaml`.
- Output paths are arm-isolated and currently have no old run artifacts beyond code-review gate files.
- Evaluator path uses `--span-mode cds`, `--profile smoke`, `aggregate_gene_body_metrics.py`, and `validate_goal.py` with a STATUS file path.
- The smoke is explicitly non-claim and single-species; it cannot update M12 panel results or SOTA claims.

#### Fix / waiver record
- No waiver. Machine gates written to `outputs/M12C-GENERANNO-FAIR-CHALLENGER-SMOKE/code_review_gate.json`, `outputs/M12C-GENERANNO-1P2B-CDS-SMOKE/code_review_gate.json`, and `outputs/M12C-GENERANNO-0P5B-BASE-SMOKE/code_review_gate.json`.
- Refreshed M12C machine gates after `configs/m12_publication_panel.yaml` metadata was corrected for 0.5b-base (`k=1`, loader supported). `pre_submit_gate.py` passes again for both smoke arms.

### Code Review Gate: M12A-FIXEDMODEL-CROSSSPECIES 2026-06-17
- Reviewer mode: host-read-only
- Scope: `M12A-FIXEDMODEL-CROSSSPECIES-A2R-s{0,1,2}` screen preflight.
- Verdict: PASS_WITH_WARNINGS for screen preflight.
- Linked implementation: `src/foundation_probe/train_unfreeze_backbone.py`, `scripts/experiments/M11-L12-SPEC-CALIBRATION/calibrate_decode.py`, `configs/M12A-FIXEDMODEL-CROSSSPECIES.yaml`, `sbatch/M12A-FIXEDMODEL-CROSSSPECIES.sbatch`.
- Linked evaluator contract: docs/19_evaluator_contract.md.
- Linked baseline reproduction: docs/20_baseline_reproduction.md.

#### Blockers
- [x] Historical trainer pooled all provided species into train/val/test, which would leak rice into the fixed-model training pool. Fixed by adding optional train/val/test species allowlists; M10/M11 behavior is unchanged when the flags are omitted.
- [x] Historical calibrator assumed val/test species sets were identical, which would make M12A look for nonexistent `test_arabidopsis` or `val_oryza` raw scores. Fixed by adding optional `--val-species` and `--test-species` while keeping `--species` as the M11-compatible default.

#### Warnings
- Review independence is `host_self`, acceptable only for screen preflight. Full/scale or claim runs need independent review.
- M12A is still non-claim: `ACTIVE_GOAL.status=draft`, SOTA benchmark unresolved, and the panel is clean for NT-v2 but not yet a published-SOTA claim panel.

#### Confirmed OK
- Static checks pass: `python3 -m py_compile` for trainer/calibrator/M12B checker; `bash -n` for M12A and M12C sbatch; YAML parse for M12 configs.
- M12A dataflow is fixed-model: train species = Arabidopsis; validation/calibration species = Arabidopsis; test species = rice. Rice labels are not used for training, early stopping, or decode/FPR selection.
- `pre_submit_gate.py` passes for `M12A-FIXEDMODEL-CROSSSPECIES-A2R-s0/s1/s2`.

#### Fix / waiver record
- No waiver. Machine gates written to `outputs/M12A-FIXEDMODEL-CROSSSPECIES/code_review_gate.json` and `outputs/M12A-FIXEDMODEL-CROSSSPECIES-A2R-s{0,1,2}/code_review_gate.json`.

### Code Review Gate: M12B-SAMEPANEL-BASELINES-EXTERNAL 2026-06-17
- Reviewer mode: host-read-only
- Scope: `M12B-SAMEPANEL-BASELINES-{ANNEVO,TIBERIUS,HELIXER}` screen external baseline runner.
- Verdict: PASS_WITH_WARNINGS for screen preflight.
- Linked implementation: `configs/M12B-SAMEPANEL-BASELINES.yaml`, `sbatch/M12B-SAMEPANEL-BASELINES-EXTERNAL.sbatch`.
- Linked evaluator contract: docs/19_evaluator_contract.md.
- Linked baseline reproduction: docs/20_baseline_reproduction.md.

#### Blockers
- [x] No open blockers. The runner uses the same clean panel (`arabidopsis_thaliana`, `oryza_sativa`), `--span-mode cds`, and the active full-transcript intergenic ruler through `scripts/eval_gene_body_mask.py`.

#### Warnings
- Review independence is `host_self`, acceptable only for screen/preflight.
- Tiberius `angiosperms.yaml` may download/cache remote weights on first run; cache/version evidence must be recorded after the run before claim.
- These are external-tool same-panel screen results, not a published-SOTA claim.

#### Confirmed OK
- Arm-specific output dirs prevent cross-tool overwrite: `outputs/M12B-SAMEPANEL-BASELINES-ANNEVO`, `...-TIBERIUS`, `...-HELIXER`.
- ANNEVO uses pinned local `ANNEVO_Magnoliopsida.pt`; Helixer uses pinned local `land_plant_v0.3_m_0100.h5` with recorded sha256; Tiberius uses local SIF + `model_cfg/angiosperms.yaml`.
- Static checks pass: `bash -n sbatch/M12B-SAMEPANEL-BASELINES-EXTERNAL.sbatch`; YAML parse for `configs/M12B-SAMEPANEL-BASELINES.yaml`.
- `pre_submit_gate.py` passes for all three arm exp_ids.

#### Fix / waiver record
- No waiver. Machine gates written to `outputs/M12B-SAMEPANEL-BASELINES-{ANNEVO,TIBERIUS,HELIXER}/code_review_gate.json`.

### Code Review Gate: M12B-SAMEPANEL-BASELINES-EXTERNAL TMPDIR repair 2026-06-17
- Reviewer mode: attempted separate `codex exec --sandbox read-only`, then host-read-only fallback.
- Scope: runtime-only repair to `sbatch/M12B-SAMEPANEL-BASELINES-EXTERNAL.sbatch` after ANNEVO failed with Python multiprocessing `OSError: AF_UNIX path too long`.
- Verdict: PASS_WITH_WARNINGS for screen preflight rerun.
- Linked implementation: `sbatch/M12B-SAMEPANEL-BASELINES-EXTERNAL.sbatch`.

#### Blockers
- [x] ANNEVO multiprocessing used a project-depth `TMPDIR` under `outputs/.../tmp/jobtmp`, exceeding AF_UNIX socket path limits. Fixed by moving `TMP_DIR` to `/tmp/m12b_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}` and keeping `TMPDIR=${TMP_DIR}/jobtmp`.
- [x] Post-fix stale-code gate was refreshed for `M12B-SAMEPANEL-BASELINES-EXTERNAL` and all three arm exp_ids.

#### Warnings
- separate_codex review did not complete because the local read-only sandbox failed with `bwrap: Creating new namespace failed: No space left on device`; independence is therefore `host_self`.
- This repair does not change model/data/evaluator semantics. It only makes ANNEVO/Helixer temporary paths short enough for runtime IPC.

#### Confirmed OK
- Baseline commands, species panel, `--span-mode cds`, aggregate metrics, and `validate_goal.py` calls are unchanged.
- `pre_submit_gate.py` passes for `M12B-SAMEPANEL-BASELINES-EXTERNAL`, `...-ANNEVO`, `...-TIBERIUS`, and `...-HELIXER` after the hash refresh.

#### Fix / waiver record
- No waiver. Machine gates refreshed with sbatch hash `43b83ada203d7a9e3529f127d4c1dc9fdb95f759b1e134cc3dcbd23d08dae235`.

### Code Review Gate: M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS 2026-06-19
- Reviewer mode: host-read-only after M17+M18 2/3 degraded tri-review/pivot.
- Scope: `M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s{0,1}` screen / Track-B-preflight array.
- Verdict: PASS_WITH_WARNINGS for screen/non-claim preflight.
- Linked implementation: `configs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS.yaml`, `sbatch/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS.sbatch`, `src/foundation_probe/train_generanno_lora_3class.py`, existing evaluator/aggregate/validate scripts.
- Linked evaluator contract: docs/19_evaluator_contract.md.

#### Blockers
- [x] No open runtime/dataflow blocker for submitting the 2-seed raw-score run. Static checks pass: `bash -n` for the sbatch, YAML parse for the config, and `py_compile` for trainer/evaluator/gate scripts.
- [x] Output dirs were absent before gate creation; seed-specific outputs are isolated as `outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s0` and `...-s1`.
- [x] `pre_submit_gate.py` passes for both seed exp_ids with non-empty reviewed file hashes.

#### Warnings
- Review independence is `host_self`, acceptable only for bounded screen/non-claim preflight. Claim/full-scale work needs independent review.
- GENERanno 1.2B provenance/overlap remains a claim blocker.
- This run only saves raw scores and produces fixed constrained predictions; validation-only calibration must be performed after completion and must not use test labels for operating-point selection.

#### Confirmed OK
- The sbatch uses Slurm array `0-1%2`, one RTX3090 per seed, and preserves the successful M18 1.2B hyperparameters while adding `--save-raw-scores`.
- The trainer will write VAL/TEST raw logits under each run's `raw_scores/` directory for no-leak calibration.
- The run remains screen/non-claim and cannot update `ACTIVE_GOAL.sota_benchmark` or support SOTA language.

#### Fix / waiver record
- No waiver. Machine gates written to `outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS/code_review_gate.json` and `outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s{0,1}/code_review_gate.json`.

### Code Review Gate: M13-DISTANCE-GENERALIZATION-SCAN-s0 2026-06-17
- Reviewer mode: host-read-only.
- Scope: `M13-DISTANCE-GENERALIZATION-SCAN-s0` screen diagnostic.
- Verdict: PASS_WITH_WARNINGS for bounded non-claim distance scan.
- Linked implementation: `src/foundation_probe/train_unfreeze_backbone.py`, `scripts/experiments/M11-L12-SPEC-CALIBRATION/calibrate_decode.py`, `configs/M13-DISTANCE-GENERALIZATION-SCAN.yaml`, `configs/m13_distance_generalization_panel.yaml`, `sbatch/M13-DISTANCE-GENERALIZATION-SCAN.sbatch`, `scripts/experiments/M13-DISTANCE-GENERALIZATION-SCAN/prep_arabidopsis_lyrata.py`, `scripts/experiments/M13-DISTANCE-GENERALIZATION-SCAN/analyze_m12a_failure.py`.
- Linked evaluator contract: docs/19_evaluator_contract.md.

#### Blockers
- [x] No open blockers. M13 reuses the M12A train/val/test species allowlist path: train/calibrate = `arabidopsis_thaliana`; test = `arabidopsis_lyrata`, `oryza_sativa`.
- [x] Close-plant data is frozen and split-safe: `reports/M13_CLOSE_PLANT_FREEZE.json`; `reports/M13_CLOSE_PLANT_CHECK_DATA.json` status PASS with `--group-col id`.
- [x] Static checks pass: `python3 -m py_compile` for trainer/calibrator/M13 scripts; `bash -n sbatch/M13-DISTANCE-GENERALIZATION-SCAN.sbatch`; YAML parse for M13 configs and `configs/m1_data_manifest.yaml`.

#### Warnings
- Review independence is `host_self`, acceptable only for screen diagnostic. Claim/full-scale work needs independent review.
- `Arabidopsis lyrata` assembly level is Scaffold; use M13 as a distance diagnostic only, not final clean claim evidence.
- M13 is single seed by design. It can route the project but cannot support statistical/generalization claims.

#### Fix / waiver record
- No waiver. Machine gate written to `outputs/M13-DISTANCE-GENERALIZATION-SCAN-s0/code_review_gate.json`.

### Code Review Gate: M14/M15-PARALLEL-DIAGNOSTICS 2026-06-17
- Reviewer mode: host-read-only.
- Scope: `M14-ANIMAL-DISTANCE-NEGCTRL-s0` screen animal negative-control diagnostic and `M15-GENERANNO-LORA-PANEL-SCREEN` bounded GENERanno panel screen.
- Verdict: PASS_WITH_WARNINGS for screen/non-claim preflight.
- Linked implementation: `configs/M14-ANIMAL-DISTANCE-NEGCTRL.yaml`, `configs/M15-GENERANNO-LORA-PANEL-SCREEN.yaml`, `sbatch/M14-ANIMAL-DISTANCE-NEGCTRL.sbatch`, `sbatch/M15-GENERANNO-LORA-PANEL-SCREEN.sbatch`, existing trainers/evaluators.
- Linked evaluator contract: docs/19_evaluator_contract.md.

#### Blockers
- [x] No open blockers after the overwrite-guard repair. The first M14 submission `9022458` failed at 0s before training because the guard treated the newly-created empty `metrics/` directory as a non-empty output conflict.
- [x] Fixed guard in both M14 and M15 sbatch scripts to only treat existing non-gate files as overwrite conflicts (`find ... -type f ! -name code_review_gate.json ! -name STATUS`), then refreshed stale code-review hashes.

#### Warnings
- Review independence is `host_self`, acceptable only for screen diagnostics.
- M14 animals are negative-control diagnostics, not clean final-claim evidence.
- M15 is a bounded screen larger than smoke but still capped by train/val window limits; do not overinterpret negative or positive results as final GENERanno route closure.

#### Confirmed OK
- Static checks pass: `bash -n` for M14/M15 sbatch, YAML parse for M14/M15 configs, and `py_compile` for trainers/calibrator/evaluator scripts.
- `pre_submit_gate.py` passes for `M14-ANIMAL-DISTANCE-NEGCTRL-s0`, `M15-GENERANNO-1P2B-CDS-PANEL-SCREEN`, and `M15-GENERANNO-0P5B-BASE-PANEL-SCREEN`.
- M14 resubmitted as job `9022700` and entered RUNNING; M15 shared job `9022457_[0-1]` was cancelled before start and rerouted to private array `9023295_[0-1%1]`.

#### Fix / waiver record
- No waiver. Machine gates written to `outputs/M14-ANIMAL-DISTANCE-NEGCTRL-s0/code_review_gate.json`, `outputs/M15-GENERANNO-1P2B-CDS-PANEL-SCREEN/code_review_gate.json`, and `outputs/M15-GENERANNO-0P5B-BASE-PANEL-SCREEN/code_review_gate.json`.

### Code Review Gate: M16-MULTISPECIES-TRAIN-DISTANCE-SCAN-s0 2026-06-18
- Reviewer mode: host-read-only.
- Scope: `M16-MULTISPECIES-TRAIN-DISTANCE-SCAN-s0` screen/non-claim multi-species fixed-model diagnostic.
- Verdict: PASS_WITH_WARNINGS for screen preflight.
- Linked implementation: `configs/M16-MULTISPECIES-TRAIN-DISTANCE-SCAN.yaml`, `sbatch/M16-MULTISPECIES-TRAIN-DISTANCE-SCAN.sbatch`, existing M9-L12 trainer/calibrator/evaluator.
- Linked evaluator contract: docs/19_evaluator_contract.md.

#### Blockers
- [x] No open blockers. M16 changes the species policy only: train/calibrate species are `arabidopsis_thaliana` + `oryza_sativa`; test species are `arabidopsis_lyrata`, `gallus_gallus`, and `drosophila_melanogaster`.
- [x] Calibration remains validation-only: selected decode parameters come from Arabidopsis+rice VAL raw scores, then apply once to test species.
- [x] Output path is unique: `outputs/M16-MULTISPECIES-TRAIN-DISTANCE-SCAN-s0`; no overwrite conflict with M13/M14/M15.
- [x] Static checks pass: `bash -n` for sbatch, YAML parse for config, and `py_compile` for trainer/calibrator/evaluator/goal scripts.

#### Warnings
- Review independence is `host_self`, acceptable only for screen diagnostics.
- M16 is single seed and includes animal negative controls, so it can route the project but cannot support final claim language.
- Drosophila has many small seqids; runtime may be dominated by prediction/evaluation rather than training.

#### Confirmed OK
- `pre_submit_gate.py` consumes the machine gate with non-empty reviewed file hashes.
- `validate_goal.py` reads `metrics/metrics.json` under the same screen profile used by M13/M14.

#### Fix / waiver record
- No waiver. Machine gates written to `outputs/M16-MULTISPECIES-TRAIN-DISTANCE-SCAN-s0/code_review_gate.json` and `outputs/M16-MULTISPECIES-TRAIN-DISTANCE-SCAN/code_review_gate.json`.

### Code Review Gate: M17-SAMEPANEL-GENERALIZATION-BASELINES 2026-06-18
- Reviewer mode: host-read-only.
- Scope: `M17-SAMEPANEL-GENERALIZATION-BASELINES-{ANNEVO,TIBERIUS,HELIXER}` screen/non-claim baseline comparability audit.
- Verdict: PASS_WITH_WARNINGS for screen baseline inference.
- Linked implementation: `configs/M17-SAMEPANEL-GENERALIZATION-BASELINES.yaml`, `sbatch/M17-SAMEPANEL-GENERALIZATION-BASELINES.sbatch`, existing released baseline runners/evaluator scripts.
- Linked evaluator contract: docs/19_evaluator_contract.md.

#### Blockers
- [x] No open blockers for submission. ANNEVO has local `Magnoliopsida`, `Aves`, and `Insecta` weights; Tiberius has local model configs for `angiosperms`, `vertebrates`, and `insecta`; Helixer has local `land_plant`, `invertebrate`, and newly pinned `vertebrate` weights.
- [x] M17 is explicitly NON-CLAIM: released-weight baseline inference under the same project evaluator, not same-budget retraining.
- [x] Static checks pass: YAML parse, `bash -n` for sbatch, and `py_compile` for evaluator/gate scripts.

#### Warnings
- Review independence is `host_self`, acceptable only for screen baseline audit.
- A. lyrata remains scaffold-level; animals remain diagnostic until overlap audit is complete.
- ANNEVO `Magnoliopsida` is reused for rice for continuity with M12B; record this caveat in result-log.

#### Confirmed OK
- `pre_submit_gate.py` passes for all three arms: `M17-SAMEPANEL-GENERALIZATION-BASELINES-ANNEVO`, `...-TIBERIUS`, and `...-HELIXER`.
- Helixer vertebrate weight pinned at `refs/weights/helixer-2025/vertebrate/vertebrate_v0.3_a_0400.h5`, sha256 `4cfa6290d4162db51370eb66d02e6ac7fd448fb9c5535e7e4d508788975d7778`.

#### Fix / waiver record
- No waiver. Machine gates written under `outputs/M17-SAMEPANEL-GENERALIZATION-BASELINES-{ANNEVO,TIBERIUS,HELIXER}/code_review_gate.json`.

### Code Review Gate: M18-PARALLEL-DIAGNOSTICS 2026-06-19
- Reviewer mode: attempted `separate_codex` read-only review first; environment blocked it before file reads (`bwrap: Creating new namespace failed: No space left on device`). Fallback reviewer mode: host-read-only.
- Scope: `M18-MULTICLADE-TRAIN-DIAGNOSTIC-s0` and `M18-GENERANNO-1P2B-SPEC-OBJECTIVE-s0`, both screen/non-claim.
- Verdict: PASS_WITH_WARNINGS for screen preflight.
- Linked implementation: `configs/M18-MULTICLADE-TRAIN-DIAGNOSTIC.yaml`, `sbatch/M18-MULTICLADE-TRAIN-DIAGNOSTIC.sbatch`, `configs/M18-GENERANNO-1P2B-SPEC-OBJECTIVE.yaml`, `sbatch/M18-GENERANNO-1P2B-SPEC-OBJECTIVE.sbatch`, existing NT-v2/GENERanno trainers and M11 calibrator/evaluator.
- Linked evaluator contract: docs/19_evaluator_contract.md.

#### Blockers
- [x] No open blockers for screen submission.
- [x] `M18-MULTICLADE` keeps train/calibrate species (`arabidopsis_thaliana`, `oryza_sativa`, `drosophila_melanogaster`) disjoint from test species (`arabidopsis_lyrata`, `gallus_gallus`, `saccharomyces_cerevisiae`) and uses saved raw scores with M11 validation-only decode selection on VAL species only.
- [x] `M18-GENERANNO` does not claim VAL-only calibration; it is an objective/postproc pressure test (`fp_lambda=2.5`, `min_cds_len=90`) using existing constrained predictions.
- [x] Output paths are unique and empty before submission: `outputs/M18-MULTICLADE-TRAIN-DIAGNOSTIC-s0` and `outputs/M18-GENERANNO-1P2B-SPEC-OBJECTIVE-s0`.
- [x] Static checks pass: YAML parse, `bash -n` for both sbatch scripts, and `py_compile` for the trainer/calibrator/evaluator/gate scripts.

#### Warnings
- Review independence is `host_self` because the separate Codex read-only sandbox was unavailable; acceptable only for screen/non-claim diagnostics.
- `M18-MULTICLADE` includes overlap-caveat species and is not final clean held-out claim evidence.
- `M18-GENERANNO` cannot answer whether decode calibration would rescue GENERanno; it only tests stronger FP objective and conservative post-processing.

#### Confirmed OK
- Both sbatch scripts call `pre_submit_gate.py` with exp_ids matching their output roots.
- Both experiments write `metrics/metrics.json` and run `validate_goal.py --profile screen`.
- Cross-tool/screen evaluator remains `scripts/eval_gene_body_mask.py --span-mode cds` plus `scripts/aggregate_gene_body_metrics.py`, matching docs/19 for screen comparability.
- Scheduler fix after submission attempts: `M18-MULTICLADE` walltime reduced from `71:50:00` to `47:50:00`, then route changed from `shared-gpu` to `private-teodoro-gpu` because shared-gpu kept returning `PartitionTimeLimit`; gate hashes refreshed.
- Runtime fix before final resubmit: first live start showed `train_windows=32037` and `val_windows=28274` because Drosophila has many small seqids. Added optional `--limit-train-windows/--limit-val-windows` to `train_unfreeze_backbone.py` and set M18MULTI to `8192/4096`; gate hashes refreshed.

#### Fix / waiver record
- No waiver. Machine gates written to both parent and component paths: `outputs/M18-MULTICLADE-TRAIN-DIAGNOSTIC/code_review_gate.json`, `outputs/M18-MULTICLADE-TRAIN-DIAGNOSTIC-s0/code_review_gate.json`, `outputs/M18-GENERANNO-1P2B-SPEC-OBJECTIVE/code_review_gate.json`, and `outputs/M18-GENERANNO-1P2B-SPEC-OBJECTIVE-s0/code_review_gate.json`.

### Code Review Gate: M18-GENERANNO-0P5B-SPEC-OBJECTIVE 2026-06-19
- Reviewer mode: host-read-only.
- Scope: `M18-GENERANNO-0P5B-SPEC-OBJECTIVE-s0`, screen/non-claim sibling of the running 1.2B GENERanno specificity-objective run.
- Verdict: PASS_WITH_WARNINGS for screen preflight.
- Linked implementation: `configs/M18-GENERANNO-0P5B-SPEC-OBJECTIVE.yaml`, `sbatch/M18-GENERANNO-0P5B-SPEC-OBJECTIVE.sbatch`, existing `src/foundation_probe/train_generanno_lora_3class.py` and evaluator scripts.

#### Blockers
- [x] No open blockers for screen submission.
- [x] Output path is unique and empty except for the machine gate: `outputs/M18-GENERANNO-0P5B-SPEC-OBJECTIVE-s0`.
- [x] Static checks pass: YAML parse, `bash -n` for the sbatch script, and `py_compile` for the GENERanno trainer.

#### Warnings
- Review independence is `host_self`, acceptable only for screen/non-claim diagnostics.
- This is an objective-control sibling, not a strict same-context-length comparison: 0.5B base uses the M15-proven 1024 bp window because it is 1 bp/token.
- GENERanno trainer still does not emit raw scores for VAL-only decode calibration.

#### Confirmed OK
- `pre_submit_gate.py` passes for both parent and component gates.
- The run uses the same clean-plant species set, `fp_lambda=2.5`, `min_cds_len=90`, and CDS-span evaluator as the 1.2B M18 sibling.

#### Fix / waiver record
- No waiver. Machine gates written to `outputs/M18-GENERANNO-0P5B-SPEC-OBJECTIVE/code_review_gate.json` and `outputs/M18-GENERANNO-0P5B-SPEC-OBJECTIVE-s0/code_review_gate.json`.

### Code Review Gate: M22-GENERANNO-1P2B-NONCRF-GBTVERSKY 2026-06-23
- Reviewer mode: separate Codex read-only pre-submit pack. Initial review returned `BLOCKED`; after fixing `--decoder none` and adding M22-specific promotion gates, final status is `PASS_WITH_WARNINGS`.
- Scope: `M22-GENERANNO-1P2B-NONCRF-GBTVERSKY-s{0,1}`, screen/non-claim non-CRF objective experiment.
- Linked implementation: `src/foundation_probe/train_generanno_lora_3class.py`, `configs/M22-GENERANNO-1P2B-NONCRF-GBTVERSKY.yaml`, `sbatch/M22-GENERANNO-1P2B-NONCRF-GBTVERSKY.sbatch`, `scripts/experiments/M22-GENERANNO-1P2B-NONCRF-GBTVERSKY/aggregate_promotion_gate.py`.

#### Blockers
- [x] No open blockers for screen preflight.
- [x] Real blocker fixed: sbatch now explicitly passes `--decoder none`, so M22 cannot silently fall back into the abandoned CRF route.
- [x] M22-specific hard gate added: each seed writes `metrics/m22_promotion_gate.json` requiring `FPR<=0.01`, `gbF1 > M19 s1`, and `gene_count<=1.25`.
- [x] Seed aggregation added: `aggregate_promotion_gate.py` writes `any_seed_promote` / `continue_route` after both seed outputs exist.

#### Warnings
- GENERanno provenance remains `overlap_unknown`; screen results are adaptation/challenger evidence only.
- `validate_goal.py --profile screen` still uses the project-wide screen threshold, so M22 route decisions must consume `m22_promotion_gate.json` and the seed-level `promotion_summary.json`.
- The final independent review remained `PASS_WITH_WARNINGS` because it reviewed the aggregation strategy before the last hardening patch; host self-review confirmed the patch recomputes required checks directly and handles malformed gate JSON as non-promoting.

#### Confirmed OK
- Static checks pass: `py_compile` for trainer + aggregation script, `bash -n` for sbatch, YAML parse.
- Data leakage checks pass for Arabidopsis and rice split JSONL with `id` as id/group: no train/val/test overlap.
- Evaluator path remains `scripts/eval_gene_body_mask.py --span-mode cds` plus `scripts/aggregate_gene_body_metrics.py`, matching `docs/19_evaluator_contract.md`.

#### Fix / waiver record
- No waiver. Machine gates written to `outputs/M22-GENERANNO-1P2B-NONCRF-GBTVERSKY-s0/code_review_gate.json` and `outputs/M22-GENERANNO-1P2B-NONCRF-GBTVERSKY-s1/code_review_gate.json`.

### Code Review Gate: M23-NTV2-CLEAN-TRANSFER-s0 2026-06-30
- Reviewer mode: attempted separate Codex read-only review first; environment blocked file reads with `bwrap: Creating new namespace failed: No space left on device`. Fallback reviewer mode: host-read-only.
- Scope: `M23-NTV2-CLEAN-TRANSFER-s0`, single-seed clean-provenance NT-v2 transfer-learning screen after M22 negative.
- Verdict: PASS_WITH_WARNINGS for screen/non-claim preflight.
- Linked implementation: `configs/M23-NTV2-CLEAN-TRANSFER.yaml`, `sbatch/M23-NTV2-CLEAN-TRANSFER-s0.sbatch`, existing `src/foundation_probe/train_unfreeze_backbone.py`, CDS-span evaluator, aggregator, and `validate_goal.py`.
- Linked evaluator contract: docs/19_evaluator_contract.md.

#### Blockers
- [x] No open blockers for screen submission.
- [x] Forbidden M22 continuations are not used. The sbatch command does not pass `--save-raw-scores`, does not call `calibrate_decode.py`, does not use a CRF decoder, and the NT-v2 trainer only supports `loss in {ce, fp_aware}`.
- [x] Output path is unique: `outputs/M23-NTV2-CLEAN-TRANSFER-s0`; it does not overwrite M10/M11/M19/M22 artifacts.
- [x] Static checks pass: YAML parse, `bash -n` for the sbatch script, and `py_compile` for trainer/evaluator/aggregator/goal scripts.
- [x] Data split checks pass for Arabidopsis and rice split JSONL with `id` as both split id and seqid group: no train/val/test overlap.

#### Warnings
- Review independence is `host_self` because the separate Codex read-only sandbox failed before file reads; acceptable only because this is screen/non-claim.
- This run intentionally duplicates the M10 direct NT-v2 L12 recipe under a claim-route/provenance decision label; interpret against M10/M11/M19/M22 rather than as a new architecture.

#### Confirmed OK
- The run uses `InstaDeepAI/nucleotide-transformer-v2-500m-multi-species` via `train_unfreeze_backbone.py`, `--unfreeze-layers 12`, `--loss fp_aware`, constrained post-processing, and clean plant species `{arabidopsis_thaliana, oryza_sativa}`.
- The evaluator remains `scripts/eval_gene_body_mask.py --span-mode cds` followed by `scripts/aggregate_gene_body_metrics.py`, producing `metrics/metrics.json` for `validate_goal.py --profile screen`.
- Runtime budget is compatible with historical M10 single-seed elapsed times around 20h on private 3090; sbatch requests private 1 GPU, 80G, and 35:50h.

#### Fix / waiver record
- No waiver. Machine gate written to `outputs/M23-NTV2-CLEAN-TRANSFER-s0/code_review_gate.json`.
