# Decisions Log (read before each new iteration)

每次 /goal-prompt 生成新迭代前,Claude 必须先读完整个本文件,确认新方向与任何 abandoned route 都没有 unexplained overlap。

如果新方向落在某个 cousin 列表里,必须在 /goal command 的「差异化说明」段明确写"这次为什么不同",或考虑放弃。

**注意:单次实验失败不进本文件**,进 docs/06_results_log.md 就够。只有 /tri-review + /pivot 决定 abandon 整条 route 才进这里。

每个 entry 用 ## DEC-<NNN>: <route name> 开头。模板见 /decisions-log SKILL.md。

---

## DEC-001: GENERanno 1.2B trained CRF decoder route

- Date: 2026-06-22 CEST.
- Decision source: `$tri-review`/`$pivot` for `M21-GENERANNO-1P2B-CRF-SCREEN`.
- Route status: abandoned / do not scale / do not tune as current component.

### Path tried
`GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview` with our 3-class FP-aware LoRA adaptation plus a trained CRF decoder and Viterbi prediction (`--decoder crf`, auxiliary CE), evaluated on the clean-plant Arabidopsis/rice screen panel.

### Evidence why failed
- M21 CRF seed0: gbF1 `0.8544`, FPR `0.0273`, gene_count_ratio `0.956`.
- M21 CRF seed1 rescue: gbF1 `0.8744`, FPR `0.0192`, gene_count_ratio `0.690`.
- Best M19 non-CRF seed on the same panel remains better on both core axes: gbF1 `0.8815`, FPR `0.0065`, gene_count_ratio `0.830`.
- Released same-panel callers remain at higher gbF1: Tiberius `0.9252`, ANNEVO `0.9269`, Helixer `0.9220`.
- Tri-review quorum 3/3: A/B `replace-component`, C `abandon-route`; all reject CRF scale-up or CRF tuning.

### What we now believe
The trained CRF decoder does not solve GENERanno's recall/gene-recovery gap on this setup. It removes the main advantage of the M19 route, namely hard-FPR-valid specificity, while failing to improve gbF1 over the non-CRF baseline. The problem is not a small tuning miss; it is a component-level mismatch for the current emissions/objective.

### Cousins to avoid unless re-entry criteria are met
- CRF transition regularization / temperature / LR sweeps on the same M21 formulation.
- Full/scale or claim runs of `GENERanno 1.2B + trained CRF decoder` on the current panel.
- HMM/CRF-style trained transition decoders that do not first demonstrate a mechanism for preserving M19-like FPR.

### Re-entry criteria
Only revisit a CRF-like route if at least one condition is true:
- A distinct emission model first reaches M19-like hard FPR (`<=0.01`) and a clear gbF1 ceiling where structure is the only remaining blocker.
- The decoder is fundamentally different from this trained per-base CRF, e.g. segment-level/semi-CRF with explicit FP guardrails and a validation-only FPR constraint.
- A local smoke/screen proves, before full training, that the structured decoder preserves FPR within `+0.002` of a matched non-CRF baseline while improving gbF1.
