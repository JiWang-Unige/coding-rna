# TA-DECODER-VEC-M3 · Protocol (Track A · vectorized LEARNED structured decoders)

## Permissions
Full built-in tools. Extend `src/screen_anchor/decoders.py` + train_screen_ref + the screen launcher. May spawn read-only `code-plan-reviewer` pre-submit. CLAUDE.md / docs/03 / docs/09 = draft patch only. tri-review = CLI (Reviewer C/agy now fixed: concise-prompt directive in reviewer_c_antigravity.sh).

## Final goal (the CORE learned-structure bet, now made tractable)
TA-DECODER-M3 validated structured decoding via CONSTR (post-processing) but the LEARNED decoders (CRF/semi-CRF) were dropped for W=2048 tractability — the core hypothesis is UNTESTED. This iteration makes them tractable and tests them fairly vs the same-budget `screen_anchor=0.5576` (base-weighted gene_body_F1_unconstrained, CDS span).
- primary_progress_gate: candidate seed-mean gene_body_F1_unconstrained (CDS) > 0.5676 (anchor + 0.01).
- Required (Codex/Claude both): report PER-SEED PAIRED delta vs softmax same-seed + seed-wise gene_count_ratio + CI/spread (seed variance was high: 0.5319/0.5779/0.6275).
- Compare learned decoders to BOTH the softmax anchor AND CONSTR (0.5791, ratio 1.12) — does LEARNED structure beat cheap post-processing? That is the scientific question.
- Screen NEVER claims SOTA; never compare vs pretrained_ceiling 0.9213. status stays draft.

## Track + resource (frozen screen protocol — IDENTICAL to anchor refs)
yeast+fly, chrom-split, window 2048, sample_fraction 0.3, 8 epochs, patience 3, 3 seeds, class-weighted CE (aux for CRF/semicrf), metric gene_body_F1_unconstrained --span-mode cds. Backbone FIXED = tiberius_like. shared-gpu RTX3090. env coding-rna (sbatch `set +u` around conda activate; short node-local TMPDIR).

## Vectorization plan (the implementation crux — /implement)
Why TA-DECODER-M3's CRF/semi-CRF were too slow: pure-python loops — CRF forward/Viterbi over W=2048 sequential steps × thousands of windows; semi-CRF triple-loop segment DP O(W·L·C²) per sample. Fixes to TRY (bounded; drop a method if still intractable = documented negative result):
- CRF: keep the time-loop sequential (CRF is inherently so) but minimize python/kernel overhead — operate in LOG space with batched (B,C,C) tensor ops only (already small C=3); the killer was per-epoch Viterbi on ALL val windows — already removed (eval=argmax). For TRAINING forward (NLL over 2048 steps × ~1000 batches): vectorize by processing the whole batch at once (done) and reduce overhead; if still too slow, CHUNK the sequence into fixed blocks (e.g. 256) with block-boundary continuity, or use an associative log-space scan (parallel-prefix) for the forward recursion. Final Viterbi predict runs once per test set (acceptable if minutes).
- semi-CRF: bound max_seg_len small (e.g. 16-32) AND vectorize the segment-emission via prefix sums (already) + replace the python (t,e,c) triple loop with batched tensor ops over the (t,e) grid; if not vectorizable within budget, drop semi-CRF and keep vectorized CRF.
- Acceptance for "tractable": a full 3-seed screen run completes in << 6h (target ≤1-2h/run) so it fits the same-budget screen and run-and-evaluate.

## Orthogonality declaration (Track A batch — SOFT_WARN focused arch batch on `decoder`)
| exp_id | major_axis | mechanism_delta | why structural |
|---|---|---|---|
| TA-CRFVEC | decoder | vectorized linear-chain CRF (learned KxK transitions + Viterbi) | learned transition tensor + structured NLL + Viterbi decode |
| TA-SEMICRFVEC | decoder | vectorized bounded semi-Markov CRF (segment scoring + duration) | segment-level DP + duration potentials |
Verdict: SOFT_WARN focused arch batch on `decoder`. Distinct mechanisms (per-base-transition vs segment-level); both attack the LEARNED-structure hypothesis. Drop semi-CRF if not bounded-fixable (keep CRF).

## Execution mode
run-and-evaluate (Track A screen, target ≤12h once vectorized). If a decoder still trends >12h after vectorization attempts, drop it (documented) or hand off; do NOT break the frozen protocol to force speed.

## Pre-submit gate (HARD)
Each candidate: (1) decoder unit-tested (extend tests/test_decoders.py: vectorized == reference on tiny inputs — partition≥gold, Viterbi correct); (2) sanity smoke end-to-end on a fast subset, predicted_genes>0, CDS F1>0; (3) TRACTABILITY check: a tiny-but-representative timing must extrapolate to << 6h/run. Fail after bounded debug (≤3) → drop that candidate.

## Required chain
1. /implement: vectorize CRF (+ semi-CRF) in src/screen_anchor/decoders.py; extend tests + train wiring.
2. Unit test (vectorized==reference) + sanity smoke + tractability timing; list survivors.
3. /smart-sbatch Phase 1 on the screen launcher; submit frozen 3-seed screen per survivor.
4. Eval --span-mode cds → aggregate → seed-mean + PER-SEED PAIRED delta vs softmax + gene_count_ratio + spread/CI.
5. /result-log → validate_goal → /tri-review (3/3 now, Reviewer C fixed) → /pivot → /exp-log.

## Pivot menu (track-A-screen)
promote-learned-decoder-to-Track-B (a learned decoder beats 0.5676 AND beats/ties CONSTR with better coherence) / keep-CONSTR-as-the-decoder (learned decoders don't beat cheap post-processing → important negative result) / change-axis (decoder axis exhausted → backbone / foundation-probe) / fix-eval. Anti-tuning: gap ≥0.05 → change axis, not lr/batch.

## Decision autonomy
Per ACTIVE_GOAL/0.8: multi-option decisions → list + tri-review (3 CLI, Reviewer C via agy with concise prompt) + ROI auto-pick + docs/08, no pause. Exceptions (pause): destructive ops, route abandon, >24h net-new compute / new long sub-iteration, ≥2 reviewers oppose / tie-no-leader.

## Parallel strand (separate, lower priority — not this goal's chain)
CONSTR → Track B scale-up (more data/seeds/CI) can run independently; tracked in docs/05. This goal focuses on the learned-decoder core bet.

## Skill chain
/implement → /smart-sbatch P1 → submit → /result-log → validate_goal → /tri-review → /pivot → /exp-log. /retrospective advisory (triggered ≥5 iters; advisory only).

## Constraints (full)
- Backbone FIXED (tiberius_like); only decoder varies; same frozen protocol as anchor refs (no drift).
- screen never claims; compare vs screen_anchor 0.5576/+0.01 AND vs CONSTR 0.5791 (the bar to beat to justify learned structure over post-processing); never vs pretrained_ceiling.
- Report per-seed PAIRED delta + gene_count_ratio + spread (high seed variance).
- "Learned decoders don't beat CONSTR" is a VALID, valuable negative result — record it, don't force a win.
- Drop a decoder that stays intractable after bounded vectorization (documented); keep the tractable ones.
