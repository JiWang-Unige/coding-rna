# Active Goal

## Council 2026-06-17 — 辩题: M12 是否应从 M9-only full/scale 改为发表证据收束 portfolio
- 立场分配: Proponent=Claude, Opponent=Codex, Referee=Antigravity · Quorum: 3/3 in Round 1 and 3/3 in Round 2.
- 辩题: 下一阶段是否应把单线 `M12-M9L12-FULLSCALE-CALIBRATED` 改造成发表证据收束组合：fixed-model cross-species evaluation、same-panel Tiberius/Helixer/ANNEVO baseline comparison、budget-capped fair GENERanno challenger；并停止 M9-only 性能微修补。
- 已夯实(共识+依据):
  - M11 已经证明 M9-L12+VAL-only calibration 在 arabidopsis/rice 同物种池上能达到 strong internal screen result：mean specificity `0.9913`, FPR `0.0087`, gbF1/constrained_gbF1@0.01 `0.8178`, gene_count `1.003`; 继续只修 M9 的 FP objective 或 decode 不是当前最高发表收益。
  - 当前最大投稿缺口不是再提高 internal metric，而是三类证据: 固定模型未见物种泛化、同面板 Tiberius/Helixer/ANNEVO 直观比较、以及证明“不是任意 pretrained model 加足训练都可以”的 GENERanno 对照。
  - `ACTIVE_GOAL.status=draft`、published `sota_benchmark` 未冻结、NT-v2/GENERanno pretraining-overlap audit 未完成，是任何 claim 前的硬 blocker。
  - GENERanno 不应平权成为 full/scale 主线；它应作为一次预算封顶、预注册判停条件的 challenger/ablation。若公平 screen 后 FPR/gene_count 仍失控，即可归档为 “pretraining alone insufficient / GENERanno CDS prior lacks intron coherence” 的机制证据。
- 仍争议 → 验证:
  - 哪些物种可以作为 clean same-panel held-out set? → 先做 NT-v2 与 GENERanno pretraining-overlap audit，再冻结 species panel。
  - M9 fixed model 在未见物种上是否还能守住 FPR/gene-count? → `M12A-FIXEDMODEL-CROSSSPECIES` zero-shot/fixed-calibration preflight。
  - Baseline runners在 same-panel 上能否可比运行? → `M12B-SAMEPANEL-BASELINES` inference dry-run first, then decide full panel.
  - GENERanno 是否只是 smoke 不公平? → `M12C-GENERANNO-FAIR-CHALLENGER` with bounded budget and explicit stop criteria.
- 否决/降级假设 + 理由:
  - 否决 “继续 M9-only 性能微修补作为下一步主线”：M11 后再优化同一物种池对发表主张边际收益低。
  - 降级 “GENERanno 与 M9 平权并行 full/scale”：M10 smoke 证据弱且工程成本高；先做 bounded ablation/challenger。
  - 暂缓 “stronger FP objective”：只有 fixed-model held-out/full panel 再次 FPR 失败时才触发。
- 用户裁决: 用户已明确指出当前推进偏离发表目标，并要求并行收束；本 council 将其转化为 M12 方向建议。正式 GPU 提交前仍需用户确认具体 panel/预算。
- 最强反方论证(保留警惕): portfolio 若不先完成 overlap audit 和 benchmark contract，会同时扩大变量并产出不可 claim 的数字；因此 M12 必须先做零/低 GPU 的 overlap audit + same-panel dry-run，再投入 full/scale。
- Raw outputs:
  - `/tmp/council_publication_alignment_20260617/round1_A_proponent.md`
  - `/tmp/council_publication_alignment_20260617/round1_B_opponent.md`
  - `/tmp/council_publication_alignment_20260617/round1_C_referee.md`
  - `/tmp/council_publication_alignment_20260617/round2_A_proponent.md`
  - `/tmp/council_publication_alignment_20260617/round2_B_opponent.md`
  - `/tmp/council_publication_alignment_20260617/round2_C_referee.md`

## last_result_summary
- exp_id: `M23-NTV2-CLEAN-TRANSFER-s0`.
- date: 2026-07-01 CEST.
- track: screen / NON-CLAIM clean-provenance NT-v2 transfer-learning check.
- primary_metric: intergenic_specificity `0.98327`; gene_body_F1_unconstrained `0.84266`; intergenic_FPR `0.01673`; gene_count_ratio `0.867`.
- **M23 FINDING**: clean-provenance NT-v2 direct transfer is semantically valid but exactly reproduces historical `M10-M9L12-CLEANPLANTS-s0`. It passes the screen FPR<=0.02 profile but fails the hard FPR<=0.01 claim guardrail (`constrained_gbF1@0.01=0.0`).
- **REFERENCE COMPARISON**: M23 is weaker than M19 GENERanno calibrated s1 on both decisive axes: gbF1 `0.8427 < 0.8815`, FPR `0.0167 > 0.0065`. Its advantage is clean public provenance, not current metric frontier.
- semantic_success: pass; Slurm job `9854668` completed `0:0` in `19:25:54`; metrics finite/parseable; loss decreased `0.7453 -> 0.4823`; no OOM/traceback.
- tri_review_status: pending.
- pivot_status: pending.
- recommended_next: run combined `$tri-review`/`$pivot` over M22 negative + M23 clean NT-v2. Do not rerun more direct M10/M23-style NT-v2 seeds; a clean-provenance continuation needs a structural change or a claim-panel/comparability purpose.

## prior_result_summary_M17_M18
- exp_id: `M17-SAMEPANEL-GENERALIZATION-BASELINES` (screen / NON-CLAIM released baseline comparability; Slurm array `9119473_[0-2%3]` COMPLETED) plus completed M18 diagnostics (`9123661`, `9131867`, and `9122868`).
- date: 2026-06-19 CEST.
- **M17 COMPARABILITY FINDING**: the A. lyrata/rice/gallus/drosophila diagnostic panel is not uniformly impossible for released gene callers. ANNEVO reaches gbF1 `0.9115` but FPR `0.0240`; Tiberius reaches gbF1 `0.8791`, specificity `0.9827`, FPR `0.0173`, but under-calls genes (`0.556x` reference); Helixer reaches gbF1 `0.8797` but FPR `0.0526`, with drosophila good and A. lyrata/gallus poor.
- **M18 MAINLINE FINDING**: broader NT-v2 train/calibration on Arabidopsis+rice+drosophila helps close-plant A. lyrata (`FPR=0.0184`, gbF1 `0.7427`, gene_count_ratio `1.078`) but aggregate held-out A. lyrata+gallus+yeast remains non-claim (`FPR=0.0444`, gbF1 `0.6170`, gene_count_ratio `1.572`) because gallus fails severely (`FPR=0.1359`, gene_count_ratio `3.994`) and yeast under-calls genes.
- **M18 GALLUS ORACLE DIAGNOSTIC**: test-label oracle calibration cannot rescue gallus. The best `FPR<=0.01` + gene-count-valid point has gbF1 `0.0060`; the best gbF1 with sane gene count has FPR `0.7347`. This points to emission/coherence failure, not a transferable decode-threshold bug.
- **M18 GENERANNO 0.5B FINDING**: stronger FP objective does not rescue `GenerTeam/GENERanno-eukaryote-0.5b-base`: aggregate FPR `0.0967`, gbF1 `0.6561`, gene_count_ratio `1.617`; rice fails badly (`FPR=0.1239`, gene_count_ratio `2.482`). 0.5B is ablation evidence only.
- **M18 GENERANNO 1.2B FINDING**: stronger FP objective makes `GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview` guardrail-valid on clean plants: aggregate FPR `0.0071`, specificity `0.9929`, macro specificity `0.9943`, gbF1/constrained gbF1@0.01 `0.8494`, gene_count_ratio `0.864`. This is much better than the M18 0.5B base sibling and promotes 1.2B CDS-preview to a serious challenger, though not a claim candidate.
- **M17+M18 TRI-REVIEW/PIVOT**: 2/3 degraded review converges operationally on switching primary next work from broad fixed NT-v2 to GENERanno 1.2B CDS-preview preflight. Decision=`change_backbone` to GENERanno 1.2B primary challenger; no generic M9 tuning and no 0.5B scale-up.
- **M19 RUN STATUS**: superseded by current `last_result_summary`; the array `9141356_[0-1%2]` later completed and is result-logged.
- **M19 PROVENANCE AUDIT**: public HF/GitHub sources do not expose a complete species/accession exclusion list for GENERanno eukaryotic CDS-preview. Arabidopsis/rice remain `overlap_unknown`; M18/M19 are mechanism/challenger evidence, not clean no-overlap claim evidence.
- **Implication for our current route**: M16 (`gbF1=0.5615`, FPR `0.0197`, gene_count_ratio `1.326`) is not competitive as a broad fixed released-caller replacement. The paper route cannot be "current M9-L12 fixed model beats Tiberius/Helixer/ANNEVO"; it must either improve toward the external-caller tradeoff frontier or narrow the claim.
- Prior generalization diagnostics remain valid: fixed Arabidopsis-only M9-L12 fails A. lyrata/rice/animal transfer; adding rice partially rescues coherence but not enough. GENERanno 1.2B CDS-preview is stronger than 0.5B base but still not guardrail-valid under M15.
- semantic_success: pass for M17 and all M18 arms. Metrics parse and are finite; no OOM/Traceback. Helixer internal logs contain NaN in empty-class diagnostic tables, but our evaluator metrics are finite.
- current_running: M19 GENERanno 1.2B raw-score 2-seed array `9141356_[0-1%2]`.
- next_decision: monitor M19 to completion, result-log both seeds, then run validation-only calibration from saved raw scores. If stable/calibratable, next pivot should consider structured decoder/segment-aware head or a claim-clean species panel; if unstable, treat M18 as single-seed positive and diagnose ceiling/provenance.
- --- prior results kept below for trend ---

## prior_result_summary
- exp_id: `M12-PUBLICATION-PREFLIGHT-TWOSEED` (publication-alignment preflight, screen / NON-CLAIM; now closed with M12A seed2 included)
- date: 2026-06-17 UTC
- **FIXED-MODEL CROSS-SPECIES NEGATIVE / EXTERNAL BASELINES STRONG**: M12A trains and calibrates on Arabidopsis only, then tests rice as unseen species. M12A three-seed mean on rice: intergenic_specificity **0.9689**, FPR **0.0311**, gene_body_F1 **0.6556**, constrained_gbF1@0.01 **0.0**, gene_count_ratio **1.755**. Same-panel external baselines are much stronger: Tiberius `gbF1=0.9252/spec=0.9927/FPR=0.0073/gcount=0.628`; ANNEVO `gbF1=0.9269/spec=0.9883/FPR=0.0117/gcount=0.726`; Helixer `gbF1=0.9220/spec=0.9784/FPR=0.0216/gcount=0.820`.
- GENERanno fair challenger smoke: 1.2B CDS-preview has signal but severe fragmentation (`gbF1=0.7527`, `FPR=0.0432`, `gene_count_ratio=4.405`); 0.5B base collapses (`FPR~1.0`, `gene_count_ratio~0.0002`).
- --- prior results kept below for trend ---
- exp_id: `M11-L12-SPEC-CALIBRATION` (primary M11 mainline, screen / Track-B preflight, NON-CLAIM)
- date: 2026-06-16 UTC
- **FPR BLOCKER CLEARED / CLAIM STILL NON-CLAIM SCREEN**: M9-L12 raw VAL/TEST emissions were saved and decoded via validation-only operating-point calibration. Corrected validate_goal=`progress` for all seeds.
- Seed-mean metrics: intergenic_specificity **0.9913**, FPR **0.0087**, macro_specificity **0.9909**, gene_body_F1 **0.8178**, constrained_gbF1@0.01 **0.8178**, gene_count_ratio **1.003**. All three seeds pass aggregate `FPR<=0.01` and `gene_count<=1.25`.
- Selected no-leakage decode points: s0 `b2p5_mcl60_mfg20`, s1 `b3p0_mcl60_mfg20`, s2 `b1p5_mcl60_mfg0`. This fixes M10's FPR tail (`0.0174 -> 0.0087`) without launching a stronger FP objective.
- Claim blockers remaining: screen/non-claim profile, `ACTIVE_GOAL.status=draft`, and published SOTA benchmark still not frozen under the same full-transcript intergenic ruler. Per-species caveat: arabidopsis seed2 FPR `0.0111`, so full/scale must report per-species sensitivity.
- TRI-REVIEW/PIVOT (2/3 DEGRADED quorum): originally `scale-to-track-b` with primary next `M12-M9L12-FULLSCALE-CALIBRATED`; this has been **reframed by 2026-06-17 retrospective+council+tri-review** into a publication-alignment portfolio before GPU submission: `M12-PREREQ-AUDIT`, `M12A-FIXEDMODEL-CROSSSPECIES`, `M12B-SAMEPANEL-BASELINES`, and bounded `M12C-GENERANNO-FAIR-CHALLENGER`. Stronger FP objective remains fallback only if fixed-model/full-panel FPR reopens the blocker.
- --- prior results kept below for trend ---
- exp_id: `M10-M9L12-CLEANPLANTS` (primary M10 mainline, screen / Track-B preflight, NON-CLAIM)
- date: 2026-06-15 UTC
- **STRONG MAINLINE PROGRESS / CLAIM BLOCKED BY FPR TAIL**: NT-v2-500m top-12 unfreeze + 3-class FP-aware convLSTM head ran across clean plants `{arabidopsis,rice}` and seeds `0/1/2`. Corrected validate_goal=`progress` for all seeds after rerunning with `--run-status outputs/.../STATUS`.
- Seed-mean metrics: intergenic_specificity **0.9826**, FPR **0.0174**, macro_specificity **0.9801**, gene_body_F1 **0.8398**, gene_count_ratio **0.897**. This is far above the same-budget screen anchor and confirms M9-L12 as the empirical mainline.
- Claim blockers: full/scale hard FPR `<=0.01` is still unmet (`0.0142-0.0212` aggregate; arabidopsis `0.022-0.0285`), and `ACTIVE_GOAL.status=draft` because the published SOTA benchmark is not frozen.
- Mechanism conclusion: deeper NT-v2 unfreeze solves the frozen-feature gbF1 ceiling and keeps gene counts coherent, but needs explicit specificity calibration / FP suppression before claim-scale promotion.
- TRI-REVIEW/PIVOT (3/3 quorum): continue-current-route. **Primary next = `M11-L12-SPEC-CALIBRATION`** (M9-L12 validation-only decode/FPR calibration + optional stronger FP objective). **GENERanno LoRA parked** until redesigned because smoke damaged specificity and did not learn intron/coherence.
- --- prior results kept below for trend ---
- exp_id: `M10-GENERANNO-LORA-3C-SMOKE` (parallel challenger smoke, NON-CLAIM)
- date: 2026-06-15 UTC
- **RUNTIME POSITIVE / METRIC NEGATIVE**: GENERanno 1.2b encoder + LoRA(q/k/v/o,r=8) + our 3-class FP-aware head ran successfully on A100 40GB and produced finite predictions/metrics. Corrected validate_goal=`progress` after rerunning with `--run-status outputs/.../STATUS`.
- Metrics on arabidopsis test smoke: intergenic_specificity **0.9491**, FPR **0.0509**, gene_body_F1 **0.7525**, constrained_gbF1 **0.0**, gene_count_ratio **4.43**. This beats the screen specificity anchor but fails the practical smoke/screen FPR expectation and remains heavily fragmented.
- Mechanism conclusion: the LoRA+3class route is technically viable, but the 8-window smoke did not learn intron/gene-body-nc (`val class2 F1=0`) and damaged GENERanno's native specificity. Do **not** submit the prepared GENERanno screen until combined M10 tri-review/pivot decides whether to redesign the schedule or park the route.
- Engineering note: current M10 sbatch validate calls pass literal `COMPLETED`; `validate_goal.py --run-status` expects a STATUS file path. Manually rerun validate for smoke; do the same for mainline after completion and fix sbatch before new submissions.
- --- prior results kept below for trend ---
- exp_id: `TB-UNFREEZE-BACKBONE-M9-DEEP` (Track A screen / Track-B preflight, NON-CLAIM)
- date: 2026-06-14 UTC
- **M9-DEEP STRONG POSITIVE**: deeper NT-v2-500m unfreeze broke the L4 FPR barrier. On arabidopsis seed0, all deeper arms L6/L8/L12 reached FPR<=0.02; best L12 = intergenic_specificity **0.9810**, FPR **0.0190**, gene_body_F1 **0.9035**, constrained_gbF1 **0.9035**, gene_count_ratio **0.792**.
- **Mechanism conclusion**: unfreezing more backbone layers improves emissions enough to lift BOTH axes: L4 gbF1 0.8759/spec 0.9754/FPR 0.0246 -> L12 gbF1 0.9035/spec 0.9810/FPR 0.0190. This validates M8's diagnosis that frozen features capped gbF1 and supports NT-v2 unfreeze as a main route.
- **Gate status**: validate_goal=`progress` after fixing a profile-aware guardrail regression in `scripts/validate_goal.py`; screen profile cannot claim SOTA. Full/scale claim remains blocked by draft SOTA benchmark, single species/seed, and full/scale FPR<=0.01.
- TRI-REVIEW/PIVOT (3/3 quorum): scale-to-track-b. **Primary next = `M10-M9L12-CLEANPLANTS`** (M9-L12 multi-seed + clean plants {arabidopsis,rice}); **parallel challenger = `M10-GENERANNO-LORA-3C`** if GPU/shared-code isolation permits, because native GENERanno has exceptional specificity but no intron concept/coherence.
- --- prior results kept below for trend ---
- exp_id: `TB-GBF1-MULTICLASS-M8` (③ Track-B scale-up, NON-CLAIM)
- date: 2026-06-12 UTC
- **M8 PRIMARY BET FAILED (key negative result)**: multi-class structured output did NOT recover gbF1 on CLEAN held-out plants {arabidopsis,rice}: mc-candidate gbF1 0.7189 NOT > 3c-candidate 0.7392 (−0.020, gcount 0.66 under-pred). gbF1->ceiling gap (~0.16) NOT closed by richer decoder labels -> structural, frozen features likely cap gbF1.
- **CLEAN POSITIVE side-finding**: 3c-candidate (frozen SegmentNT+FP-aware+constrained) PARETO-beats raw-DNA anchor on clean plants BOTH axes (spec 0.966 vs 0.905 +0.062; gbF1 0.739 vs 0.696 +0.043), leakage-free (SegmentNT backbone excludes plants) — honest replacement for M7's chicken-contaminated +0.155 headline.
- NEXT (pending tri-review/pivot): multi-class NOT scaled; next gbF1 axis = staged UNFREEZE/fine-tune SegmentNT OR backbone-only self-trained head (route-level, >24h, USER go-ahead).
- --- prior results kept below for trend ---
- **CRITICAL 2026-06-12 (M8-CK4 SegmentNT audit, surfaced early)**: SegmentNT(multi_species) segmentation head was FINE-TUNED on {human,mouse,chicken,fly,zebrafish,worm} -> our **chicken + fly** evals are LABEL-LEAKAGE contaminated; **only arabidopsis (plant) is truly-clean held-out** (candidate still wins there: spec 0.954 vs 0.892, +0.06 = the honest signal). M7 headline +0.155 is inflated by chicken leakage. PRE-CLAIM GATE: claim only on segmentation-clean species (not the 6). M8 redirected: clean held-out = {arabidopsis, rice}; chicken = contaminated robustness stratum (labeled). See docs/10.
- exp_id: `REANCHOR-HELDOUT-M7` (Track A screen, NON-CLAIM, retrospective-derived re-anchor GATE before ③ Track-B)
- date: 2026-06-12 UTC
- **HELD-OUT RE-ANCHOR GATE PASSED (Pareto-ADMISSIBLE) — SPECIFICITY axis reinforced cross-clade; gbF1 axis LOSES to anchor (tri-review correction).** On held-out/UTR-rich {Arabidopsis(plant), Gallus(vertebrate)} (UTR 42%/62% of exon vs yeast+fly ~0): candidate FP-FRAGFIX-CONSTR (5-seed, IDENTICAL config) intergenic_specificity **0.9604+-0.008** (all 5 > held-out anchor 0.8054, +0.155) > macro 0.9621 (>0.7804) > gbF1 0.6664 (>floor) ; gene_count 0.9688 (<=1.25, mild 3% under-pred, VAL-selected no leakage). ANNEVO ceiling 0.9824.
- KEY: margin over same-budget anchor LARGER cross-clade (+0.155) than yeast+fly (+0.078); absolute spec HIGHER (0.9604 vs 0.9218), near ceiling. Retrospective worry (numbers on low-UTR outliers) POSITIVELY REFUTED.
- CAVEAT (Track-B): mild held-out gene_count under-prediction (per-clade band calib); chicken subset gene-dense (macrochromosome intergenic untested); candidate gbF1<anchor (trades CDS-F1 for spec; multi-class output planned).
- NEXT: ③ Track-B promotion now de-risked on held-out — USER GO-AHEAD required (>24h compute). Optional /revise-goal to record held-out anchor 0.8054 / ceiling 0.9824 alongside yeast+fly.
- --- prior results kept below for trend ---
- exp_id: `TA-FRAGFIX-SWEEP-M6` (Track A screen, NON-CLAIM, STEP-0 promote-gate)
- date: 2026-06-11 UTC
- **STEP-0 GATE CLEARED -> FP-FRAGFIX-CONSTR is PROMOTE-READY.** VAL-chosen (no test leakage) constrained params mfg=20/mcl=90 -> TEST 5-seed: intergenic_specificity **0.9262** (>anchor 0.8710/0.8436, all 5 seeds), gene_body_F1 0.6376 (>floor/anchor), macro 0.8389 (>gate), **gene_count_ratio 1.28 -> 0.939 (<=1.25 guardrail)** -> ALL 4 gates PASS.
- CAVEAT: gene_count seed variance 0.55-1.35 (mcl=90 aggressive, some under-predict); milder mfg=20/mcl=30 (gcount~1.0) is a Track-B option.
- PIVOT (pending tri-review): PROMOTE FP-FRAGFIX-CONSTR to Track B = ③ (USER GO-AHEAD): scale data/epochs/seeds + Tiberius multi-class (CDS/intron/intergenic/phase/splice) + staged SegmentNT unfreeze.
- ladder (new ruler): FLOOR 0.8805 / anchor 0.8436(5s) / FP-FRAGFIX-CONSTR 0.926 / ceiling 0.9917. status draft.
- --- prior results kept below for trend ---
- exp_id: `TA-COHERENCE-FIX-M5` (Track A screen, NON-CLAIM, M4 pivot follow-up: de-fragment FPLOSS + 5-seed anchor)
- date: 2026-06-11 UTC
- **FP-FRAGFIX-CONSTR (FPLOSS + deterministic constrained post-proc) = paired-significant Pareto winner**: intergenic_specificity 0.9272 ± 0.036 (paired +0.0836 ± 0.037 vs 5-seed anchor 0.8436, ALL 5 seeds positive), gene_body_F1 0.6581 > anchor 0.5768, macro 0.8555 > gate, gene_count_ratio 2.25 -> **1.28** (fragmentation 95% fixed; 0.03 above full/scale guardrail 1.25 -> trivially tunable). KEEPS M4 FPLOSS spec (0.927) + IMPROVES F1 (0.658). Strong Track-B promotion candidate.
- 5-seed anchor 0.8436 < old 3-seed 0.8710 (new seeds weaker) -> screen_anchor candidate for /revise-goal update; doesn't change conclusion (CONSTR beats both).
- PIVOT (pending tri-review): promote-ready (gene_count 1.28 trivially closeable). ③ = Track-B scale + richer multi-class (phase/splice) + maybe unfreeze SegmentNT -> USER GO-AHEAD.
- screen_anchor=0.8710(3-seed, 5-seed=0.8436); ceiling 0.9917; status draft.
- --- prior results kept below for trend ---
- exp_id: `TA-FOUNDATION-DECODER-M4` (Track A screen, NON-CLAIM, MAIN architecture bet: foundation -> structured decoder)
- date: 2026-06-11 UTC
- 3 candidates x 5 seeds (8 epochs) on frozen SegmentNT features (reused FEATCACHE), NEW ruler, vs anchor 0.8710.
- **WINNER FP-SEGNT-FPLOSS** (FP-aware specificity-targeted loss): intergenic_specificity **0.9303 +/- 0.036 > anchor 0.8710** (ALL 5 seeds > anchor mean) AND gene_body_F1 **0.6157 > anchor 0.5576** AND macro 0.8431 > gate 0.7978 -> **PARETO-beats the same-budget anchor on the dual co-primary**. First candidate to strictly exceed the anchor on the new ruler. MAIN bet validated at screen.
- FUSION 0.8615 (just below anchor, no); CRF 0.8298 (high var ±0.119, no; but best gene_count coherence 0.90).
- PIVOT (pending tri-review): FPLOSS = Track-B promotion candidate (scale-up = new long sub-iteration -> user go-ahead). Obvious synthesis next: FP-aware loss + CRF decoder (specificity + coherence). 2-epoch smoke was misleading -> 5-seed mandate validated.
- screen_anchor=0.8710 (PROVISIONAL); ceiling 0.9917; status draft.
- --- prior results kept below for trend ---
- exp_id: `FP-SEGMENTNT-PROBE-M1` (foundation-probe #1, Track A screen, NON-CLAIM)
- date: 2026-06-11 UTC
- WHAT: frozen SegmentNT(multi_species) 14 element features -> anchor-matched conv+biLSTM head (clean input-signal ablation), same-budget, NEW dual co-primary ruler.
- RESULT (3 seeds, bw): AXIS-2 gene_body_F1 **0.6888 >> anchor 0.5576 (+0.13, PASS)**; AXIS-1 intergenic_specificity **0.8416 < anchor 0.8710**, macro 0.7543 < gate 0.7978 (FAIL). Per-species: fly spec ~0.85 GOOD, yeast(fungus) ~0.65 POOR (over-predict 1.8-2.1x). Not Pareto-dominant -> not_yet. High spec seed variance (s1=0.897 > anchor).
- FINDING: frozen human/vertebrate foundation features improve gene DETECTION but not cross-clade intergenic specificity (weak transfer to divergent fungus). Same ↑recall/↓specificity trade-off as structured decoders -> intergenic spillover is the central obstacle.
- PIVOT (3/3 tri-review consensus): ITERATE-PROBE via change-objective-or-loss. Next (new goal): (A) FP-aware/specificity-targeted loss on frozen features; (B) raw-DNA ⊕ SegmentNT fusion; ≥5 seeds+CI+paired test. DEFER semi-CRF (until FP controlled), unfreeze (Track B), GENERanno (parallel). Pre-claim: verify test clade not in SegmentNT pretraining.
- screen_anchor=0.8710 (PROVISIONAL); status draft. CK1-CK6 all complete.
- --- prior results kept below for trend ---
- exp_id: `REVISE-INTERGENIC-PRIMARY-M1` (`/revise-goal` — FOUNDATIONAL eval-ruler change + DUAL co-primary)
- date: 2026-06-11 UTC
- WHAT: intergenic redefined = genome − FULL-transcript span (incl UTR), decoupled from gene-body-F1 span_mode. NEW primary = `intergenic_specificity` (=1−FPR). Re-scored existing same-budget predictions (no GPU/retrain). tri-review 3/3 approve-with-modifications + user human gate.
- **RANKING FLIPS** (new ruler, 3-seed base-weighted spec / gene_body_F1): FLOOR(ORF) 0.8805/0.3735 (BLOCKED by F1 floor) · **tiberius_like ANCHOR 0.8710/0.5576** · CONSTR 0.8369/0.5791 · helixer_like 0.7954/0.5579(frag) · **CRF-vec 0.7138/0.6186 (was old-ruler WINNER → now WORST)**. Structured decoders raise recall by spilling into intergenic → worse specificity (highest FPR 0.2862).
- CONTRACT: DUAL CO-PRIMARY (Pareto). AXIS-1 headline=intergenic_specificity (screen_anchor 0.8710 PROVISIONAL, macro 0.8278 gate); AXIS-2=gene-level F1 (SOTA-comparable claim). Promotable iff specificity>anchor AND gene_body_F1>=floor(0.5276 screen/0.5576 promo) AND macro_spec>=0.7978. Both M1 gaming modes blocked.
- **CRF-vec Track-B promotion INVALIDATED/PAUSED** (never launched → caught before GPU spent, per user's "白跑" warning). CRF-vec → ablation, not scaled.
- NEXT: foundation-probe (SegmentNT/GENERanno features → cut intergenic FP w/o losing recall) → semi-CRF + FP-aware objective. Re-derive anchor on UTR-rich species before Track-B HARD gating; recompute pretrained_ceiling under new ruler.
- screen_anchor=0.8710 (intergenic_specificity, PROVISIONAL); status draft (full-claim needs sota_benchmark, M2). eval tests 2 + validate_goal tests 6 pass.
- --- prior results kept below for trend ---
- TA-DECODER-M3: CONSTR (post-processing) won batch 1 (0.5791, ratio 1.12) but learned decoders were untested → this iteration tested them. (2026-06-11)
- exp_id: `M1-SAMEBUDGET-SCREEN-ANCHOR` — established screen_anchor 0.5576 (FLOOR 0.3735 < anchor 0.5576 < pretrained_ceiling 0.9213); completed_poor gate tightened; helixer fragmentation finding. (2026-06-10)
- exp_id: `BASE-ANNEVO-SAC-DMEL-SMOKE-M1` (completes the reproduce-baselines trio)
- date: 2026-06-10 UTC
- track: baseline / M1 ANNEVO two-species smoke — THIRD full gene-caller
- execution: dedicated `annevo` env (Py3.10/torch2.1; setuptools<81 fix); ran on shared-gpu gpu021 (private full → rerouted), job 8537422 COMPLETED after bounded debug (set-u activate wrap, AF_UNIX short TMPDIR). eval --span-mode cds.
- gate: `not_yet`/`completed_poor` — completed_poor verified a 3rd time in production
- ANNEVO CDS: base-w F1 0.9197 / macro 0.9429; S.cer F1 0.9735 FPR 0.0072 (lowest FPR + highest precision of the 3); D.mel F1 0.9122
- **TRIO COMPLETE** (CDS base-weighted): Tiberius 0.8608 / Helixer 0.9213 / ANNEVO 0.9197 → screen_anchor = max = **Helixer 0.9213** (ANNEVO ties but does not raise → confirms prior provisional). NOTE: this screen bar ≠ published-SOTA story (ANNEVO strength is broad-clade locus/exon gffcompare; yeast/fly are gene-dense pilots).
- POST tri-review+pivot (M1-CONTRACT-REVIEW, 2/3 DEGRADED, both comparability-blocker High; user confirmed): the 0.9213 trio max is a PRETRAINED-INFERENCE CEILING, NOT a same-budget anchor → moved to non-gating `pretrained_ceiling`; `screen_anchor` reset to PENDING (same-budget definition).
- semantic_success: pass
- recommended_next: **M1-SAMEBUDGET-SCREEN-ANCHOR sub-stage** (user scope=2 family refs): freeze one unified small-sample screen protocol → implement+train random-init Tiberius-like + Helixer-like (ANNEVO-light if tractable) + cheap FLOOR → screen_anchor = seed-avg max → /revise-goal → status active. Tighten completed_poor exemption (semantic_success flag is a constant; require unconstrained≥~0.05 + count ratio). Track A PAUSED until anchor exists. This is a substantial /implement phase (2 light architectures) — awaiting user go-ahead to start.

## 当前研究方向
构建一个用于 eukaryotic protein-coding RNA gene body 的 ab initio 跨物种逐碱基多分类基因组注释模型，输入 DNA sequence，输出每个碱基的注释标签，并严格超越已发表深度学习 SOTA。

## 任务边界
- 输入: 原始 genomic DNA sequence，不依赖 RNA-seq、蛋白同源比对、已有转录本证据或物种专属外部注释特征。
- 输出: 逐碱基多分类标签，至少区分 protein-coding gene body 与 intergenic；deep research 需要核实 SOTA 使用的标签体系，例如 exon / intron / intergenic、CDS / UTR / intron / intergenic、strand-aware labels 等。
- 算这个任务: ab initio protein-coding gene annotation；跨物种泛化；gene body 与 intergenic 边界识别；在序列层面预测基因区域和相关子区域。
- 不算这个任务: RNA folding / RNA 3D structure prediction；RNA-seq 或 protein homology evidence-based annotation；只做 transcript-level 分类；只做 promoter / enhancer / TF binding / ncRNA 注释；论文写作本身。
- 应用场景: 学术 benchmark 和可发表模型开发，最终目标是达到 Nature Methods 级别的模型贡献；SOTA claim 必须建立在可比数据、split、metric、preprocessing 和 test-time inference 协议上。

## 候选数据集
- 暂不预设物种和数据集；优先对齐已验证 SOTA 研究使用的数据、物种、split 和评估协议。
- deep research 需重点查明 Tiberius、Helixer 2025、ANNEVO、SegmentNT、GENERanno-eukaryote-1.2b-cds-annotator-preview 使用的数据来源、物种集合、训练/测试划分、跨物种泛化设置和标签定义。
- 若已有 SOTA 数据不可直接复用，后续 benchmark-roadmap 需定义公开可复现的 Ensembl / GENCODE / RefSeq 等多物种 benchmark，并按物种/染色体/同源簇去冗余切分。

## 评估指标
- Primary 候选: 在控制 intergenic false positives 的前提下，优先关注 gene-body F1 或 coding-gene-body F1。
- Secondary 候选: gene-level F1；segment-level F1；exon/intron/intergenic macro-F1；precision/recall by class；boundary accuracy；cross-species held-out performance；inference cost。
- 指标待 deep research 核实: SOTA 是否使用 base-level、segment-level、gene-level、transcript-level 或 CDS-specific 指标；gene-level F1 的 exact matching / overlap threshold / strand consistency 规则；intergenic FP 如何计算。

## 既有 SOTA(用户感知,待 /sota-inventory 验证)
- Tiberius
- Helixer 2025 version
- ANNEVO
- SegmentNT
- GENERanno-eukaryote-1.2b-cds-annotator-preview
- 传统 gene prediction 工具如 AUGUSTUS / BRAKER / GeneMark 等不作为主要指标锚点，除非 deep research 发现其仍在某些 benchmark 上有强可比性；可从传统方法中提取机制启发，例如 GT/AG splice donor/acceptor motifs 可用于边界过滤或后处理。

## 差异化假设
- 用户觉得现有 SOTA 的薄弱点: 即使深度学习方法已超越传统 gene prediction，现有模型可能仍在跨物种泛化、gene body/intergenic 边界、长序列上下文、gene-level consistency、低 FP intergenic 控制和可复现 benchmark 对齐上存在缺口。
- 我们打算如何不同: 采用组合型风险策略，一条稳健复现/改进 SOTA-like baseline，配合两条以上结构性高风险路线并行筛选；优先考虑 backbone、decoder/head、objective/loss、boundary-aware post-processing 或 data view 的真实架构创新，而不是只调学习率、batch size、dropout。
- 机制启发: 传统 gene prediction 的显式生物序列约束可作为神经模型边界校正线索，例如 GT/AG splice signals、start/stop codon consistency、ORF/gene-body consistency，但必须通过 ablation 和可比评估证明有效。

## 资源约束
- GPU 量级: 大资源配额；允许多 GPU、多天训练、Track A 并行筛查、Track B full/scale 验证，以及必要的预训练或大模型微调。
- 时间预算: 可以承担较长研究迭代，但每次真实训练仍需遵守 Slurm、smoke/screen 先行、job_watch 对账和 validate_goal 确定性闸门。
- 风险偏好: 组合型。最终目标是 Nature Methods 级别；没有架构或机制层面的显著创新很难达成目标。因此需要稳健 baseline + 中高风险结构性创新路线并行推进。

## Handoff to /research-synthesize

- 用户的核心 motivation: 做出一个跨物种泛化能力强、逐碱基标注 protein-coding gene body/intergenic 区域、严格超越深度学习 SOTA 的 ab initio gene annotation 模型。
- 隐含假设(用户没明说但贯穿对话的): 当前深度学习 SOTA 已经明显强于传统 gene prediction，主要竞争对象应是 Tiberius、Helixer、ANNEVO、SegmentNT、GENERanno 等现代模型；传统方法更适合作为机制启发，而不是主要指标锚点。
- 期望 deep research 回答的问题:
  1. 当前 eukaryotic ab initio protein-coding gene annotation 的真实深度学习 SOTA 是哪些，尤其是 Tiberius、Helixer 2025、ANNEVO、SegmentNT、GENERanno-eukaryote-1.2b-cds-annotator-preview 的 paper/code/weights/datasets/metrics。
  2. 这些 SOTA 的标签体系、数据来源、物种集合、split、跨物种泛化协议和 primary/secondary metrics 如何定义，哪些可以作为严格可比 benchmark。
  3. 哪些架构或机制缺口最可能带来 Nature Methods 级别创新，包括 long-context backbone、structured decoder、boundary-aware objective、gene-level consistency、cross-species training 和 biological constraint post-processing。
- 用户**不**关心的方向(防止 synthesize 跑偏): RNA folding / RNA 3D；RNA-seq/protein homology evidence-based annotation；传统工具排行榜式罗列；只做小型窗口分类；只靠调参或 ensemble 的路线；不可复现或私有数据 benchmark。

## direction_clarified_2026-06-09

> Generated by `$grill` after `$sota-inventory`; input for `$benchmark-roadmap`, not a final benchmark contract.

### Conceptual boundaries

- Roadmap/screen stage target: exceed the best baseline reproduced under our unified small-budget protocol (`screen_anchor`). Screen evidence is for architecture selection only and can never support a SOTA claim.
- Published claim stage target: full/scale runs must strictly exceed a frozen published SOTA anchor. The anchor is not yet final; `$benchmark-roadmap` should treat ANNEVO as provisional primary anchor and Helixer/Tiberius as complementary baselines until reproduction/metric alignment freezes `sota_benchmark`.
- Tiberius 2024 is not a broad-eukaryote SOTA anchor. It is a mammal-trained, mammal-tested structured-decoder baseline and a strong mechanism source for HMM/phase/border labels/objective design.
- Helixer is a broad-lineage open baseline; ANNEVO is a provisional cross-clade published SOTA anchor candidate. ANNEVO paper values and current GitHub release values must remain separate.
- Use the term `structured neural gene callers` for Tiberius/Helixer/ANNEVO. Do not frame them as merely "traditional models"; the relevant distinction is structured gene grammar vs foundation-model segmentation.

### Core architecture bet

- Primary track: `foundation probe -> structured decoder`.
- Mechanism delta: test whether SegmentNT/GENERanno-style pretrained base/CDS/splice signals can feed a segment-level structured decoder, with semi-CRF / segment-level decoding as the main candidate, to improve interval/gene-level consistency over Tiberius-like or Helixer-like baselines.
- Semi-CRF is the main decoder candidate for roadmap drafting; HMM/CRF/constrained Viterbi should be explicit ablations or fallback directions rather than an equally weighted first choice.

### Metric and story

- Primary paper story: interval/gene-level F1 improvement. Nucleotide-level F1 may be secondary, but cannot be sacrificed without guardrails.
- Screen metric proposal: `constrained_gene_body_F1` with profile-aware intergenic FPR guardrail. During M1 evaluator freeze, smoke/screen use `intergenic_FPR <= 0.02`; full/scale claim candidates keep `intergenic_FPR <= 0.01`. Always report sensitivity at `0.005`, `0.01`, and `0.02`.
- Additional guardrails for roadmap: nucleotide gene-body F1 should not drop by more than roughly 2-3 points, and predicted gene count should not inflate far beyond reference (initial warning threshold: >1.25x reference gene count).
- BUSCO is support-only, never the primary metric, because it can hide false-positive genes/exons.

### Data and benchmark sequencing

- First reproduce/align Tiberius-like mammal protocol to learn exact split, preprocessing, transcript collapsing, metric implementation and screen mechanics.
- Then add cross-clade screen to target held-out clade generalization within eukaryotes/animals. Broad eukaryote cross-domain generalization is a stronger follow-up/Phase 8 target, not the first claim definition.
- Binary `intergenic/gene-body` is allowed only as smoke/sanity and fast signal probing. It is not sufficient as a main publication task for full gene annotation, because it lacks exon/intron/phase/splice-border/strand gene structure.
- Formal screen should use structure-aware labels: strand-aware CDS/intron/intergenic plus phase/splice-border or equivalent information close to Tiberius-style 15-class supervision.

### Foundation-model gate and stop rules

- SegmentNT/GENERanno raw or frozen-feature probes are valid early tests but not claimable full gene-caller baselines.
- If foundation raw/frozen probe has primary screen gap `>= 0.05` versus Tiberius-like screen_anchor and shows no near-baseline signal on key labels (CDS, splice donor/acceptor, exon boundary), downgrade the foundation route and do not proceed to high-cost large-scale fine-tuning.
- If `foundation + semi-CRF` fails to exceed Tiberius-like screen_anchor and gap remains `>= 0.05`, the anti-tuning gate applies: do not continue with hyperparameter tuning; pivot architecture axis.
- If gene-level gains are caused by gene-count inflation or intergenic FPR guardrail violation, treat as invalid progress even if nominal gene-body F1 improves.

### Downgraded / secondary routes

- RMT / long-context backbone: secondary track; activate if error stratification shows long genes/long introns are dominant failure modes.
- MoE / clade routing: secondary track; ANNEVO already occupies this mechanism space, so we need a sharper novelty claim before making it primary.
- Raw foundation model direct claim: disallowed unless converted into a comparable full gene-caller with structured decoding and benchmark alignment.

### Open verification items

- Freeze published SOTA anchor and exact metric value after benchmark-roadmap/reproduce-baselines. Candidate anchors: ANNEVO paper/supplement, Helixer broad-lineage metrics, Tiberius mammal mechanism baseline.
- Verify exact species lists, transcript collapsing, preprocessing and gffcompare settings before any full/scale claim.
- Quantify whether foundation probes retain useful label-specific signal even if global primary metric lags.

## last_result_summary (2026-06-14)
TB-UNFREEZE-BACKBONE-M9 CK3 (screen, non-claim): unfreeze NT-v2-500m backbone VALIDATED — gbF1 0.828(frozen)->0.854(L2)->0.876(L4) AND spec 0.966->0.967->0.975, dual-axis monotonic, L4 nears ANNEVO ceiling 0.898. validate_goal=not_yet: constrained_gbF1=0 because L4 FPR 0.0246 just exceeds screen 0.02 threshold. Next lever (Track B, user gate): deeper unfreeze L6/L8/full to push FPR<0.02 + multi-seed + cross-species CI.
