# TA-DECODER-M3 · Protocol (Track A · structured-decoder focused batch)

## Permissions
Full built-in tools. May extend `src/screen_anchor/*` + add `configs/`, `scripts/run_screen_ref.sbatch` analog. May spawn read-only `code-plan-reviewer` before submit and `experiment-implementer` per-candidate (each writes ONLY its own exp_id config/notes). CLAUDE.md / docs/03 / docs/09 = draft patch text only, no direct Edit. tri-review = CLI, not subagent.

## Final goal (M3 milestone)
Test the project's PRIMARY architecture bet — a STRUCTURED DECODER on top of the existing per-base backbone — against the same-budget `screen_anchor = 0.5576` (gene_body_F1_unconstrained, CDS span). Motivation: the per-base reference models fragment severely (helixer ratio 51-153, tiberius 1.8-4.1); a structured decoder should improve BOTH base-level F1 AND gene-level coherence (lower gene_count_ratio).
- **primary_progress_gate**: a candidate's seed-mean `gene_body_F1_unconstrained` (CDS) > `screen_anchor + 0.01` = **0.5676**.
- **Coherence (R5)**: report seed-wise `predicted_gene_count_ratio_vs_reference`; a candidate that beats 0.5676 AND brings ratio toward the coherent envelope (≪ the per-base baseline; ideally → claim guardrail 1.25) is the promotion winner.
- Screen NEVER claims SOTA. status stays draft. Compare ONLY vs screen_anchor (same-budget), NEVER vs pretrained_ceiling 0.9213.

## Track + resource (frozen screen protocol — IDENTICAL to the anchor refs)
2 species (yeast+fly), chromosome-level split (check_data PASS already on this split), window 2048, sample_fraction 0.3, 8 epochs, patience 3, **3 seeds {0,1,2}**, class-weighted CE (sqrt-inv), Adam lr 1e-3, metric = gene_body_F1_unconstrained under `--span-mode cds`. Backbone = `tiberius_like` (the coherent reference backbone) held FIXED; only the decoder/head varies. Resource: shared-gpu RTX3090 (private if free), ~15-40 min/run; reuse `scripts/run_screen_ref.sbatch` extended with a `--decoder` arg. conda env `coding-rna` (torch 2.5.1); sbatch must `set +u` around `conda activate` (MKL) + short node-local TMPDIR.

## Orthogonality declaration (Track A batch — SOFT_WARN focused arch batch on `decoder`)
| exp_id | major_axis | mechanism_delta | why structural | not merely hyperparam? |
|---|---|---|---|---|
| TA-CRF | decoder | linear-chain CRF: learned KxK transition matrix + Viterbi decode over {intergenic,CDS,intron} | new transition param tensor + forward (CRF NLL) + Viterbi decode replace per-base softmax | yes — new params + loss + decode |
| TA-SEMICRF | decoder | semi-CRF / semi-Markov segment decoder: segment-level scoring + duration model, segment DP | segment-level DP + duration potentials; models segments not per-base | yes — new layer + DP inference |
| TA-CONSTR | decoder | constrained-Viterbi biological decode: CDS phase/length-%3 + start/stop + fragment-merge post-processing on per-base logits | structured decode constraints + merge step; reduces fragmentation deterministically | yes — new decode algorithm (no per-base softmax argmax) |
Verdict: **SOFT_WARN — focused arch batch on `decoder`**. Allowed (CLAUDE.md §8): distinct mechanism_delta, all attack the structured-decoding hypothesis. NOT a hyperparam sweep.

## Execution mode
run-and-evaluate (Track A screen, ≤12h). If implementation+training trends >12h (e.g. semi-CRF debugging), switch the REMAINING training to submit-and-handoff and result-process in a follow-up goal. First /implement + sanity-smoke each decoder BEFORE the full 3-seed screen.

## Implementation plan (/implement, with sanity smoke + bounded debug each)
1. Extend `src/screen_anchor/models.py`: add a CRF head (e.g. linear-chain CRF layer: transition matrix, forward-backward NLL, Viterbi) wrapping the tiberius_like backbone; add a semi-CRF segment decoder; keep per-base softmax as the baseline path.
2. Extend `src/screen_anchor/train_screen_ref.py`: `--decoder {softmax,crf,semicrf,constrained}`; CRF/semi-CRF change the loss (structured NLL) + the predict path (Viterbi / segment-DP) → per-base labels → existing `labels_to_cds_gff`. `constrained` reuses softmax training but a constrained-Viterbi predict + fragment-merge.
3. Sanity smoke each decoder (tiny: sample 0.02, ≤300 windows, 1-2 epochs) on srun: assert no crash, predicted_genes>0, CDS F1>0, AND gene_count_ratio LOWER than the per-base baseline (the decoder's whole point). Bounded auto-debug (repair_advisor) ≤3 tries; semi-CRF DP is the highest-risk — if it cannot be made correct/bounded in the smoke, drop TA-SEMICRF from this batch (note it), keep CRF+CONSTR.
4. check_data: reuse the existing chromosome-split manifest (already PASS); re-run if the split changes (it must NOT).

## Pre-submit gate (HARD)
Each candidate must pass its sanity smoke (runs end-to-end, predicted_genes>0, gene_count_ratio improved vs per-base baseline) before its 3-seed screen submits. If a candidate's smoke fails after bounded debug, that candidate's final pivot is `fix_eval`/drop; do NOT submit its full screen.

## Screen + eval (per candidate)
Run the frozen protocol (3 seeds) per surviving candidate via the extended sbatch; eval each species `--span-mode cds`; aggregate (base-weighted + macro + per-species); compute seed-mean base-w gene_body_F1_unconstrained + seed-wise gene_count_ratio. Compare each candidate's seed-mean vs 0.5676 (progress gate) and report coherence.

## Promotion discipline (ACTIVE_GOAL.track_a_promotion / R5)
- screen_pass = seed-mean base-F1 > 0.5676. Beating the anchor is a DIRECTION signal, NOT auto Track B.
- Report seed-wise gene_count_ratio + macro/per-species. A high-F1 but still-fragmented candidate (ratio ≫ 1.25) is NOT promoted to Track B unless its plan fixes coherence. The structured decoders are EXPECTED to reduce ratio — that improvement is the key promotion evidence, not just base-F1.

## Pivot menu (track-A-screen)
promote-to-track-B (a decoder beats 0.5676 AND improves coherence) / focused-rebatch-decoder (try another decoder mechanism) / change-axis (decoder didn't help → backbone or foundation-probe path) / fix-eval / return-to-literature. Anti-tuning: if best candidate gap to anchor ≥0.05 and decoders don't help, change axis (NOT tune lr/batch).

## Decision autonomy
Per ACTIVE_GOAL/0.8: multi-option decisions (which decoder to promote, drop semi-CRF or not, rebatch vs change-axis) → list options + tri-review + cost-adjusted ROI auto-pick + write docs/08, no pause. Exceptions (pause): destructive ops, route abandon (docs/09), >24h net new compute, tied 3-way reviewer split.

## Skill invocation chain
/implement (decoders + smoke) → /smart-sbatch Phase 1 (reuse extended sbatch) → submit screen → /result-log → validate_goal → /tri-review → /pivot → /exp-log. /retrospective advisory (triggered: ≥5 iters, no prior retrospective) — advisory only, does not block.

## Constraints (full)
- Backbone fixed (tiberius_like); only decoder varies — keeps the batch a clean decoder-axis test.
- Same frozen protocol as the anchor refs (no protocol drift) — else not comparable to screen_anchor.
- screen never claims; compare only vs screen_anchor 0.5576/+0.01; never vs pretrained_ceiling.
- Report seed-wise gene_count_ratio (R5); coherence improvement is required promotion evidence.
- semi-CRF is highest-impl-risk; bounded-debug then drop if needed, don't blow the budget.
- ANNEVO env / SegmentNT / GENERanno NOT needed this batch (foundation-probe path is a later batch).
