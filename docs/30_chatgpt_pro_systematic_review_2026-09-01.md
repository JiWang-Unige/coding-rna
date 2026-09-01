# ChatGPT Pro systematic repository review — 2026-09-01

Status: completed two-round external AI review and Codex synthesis. This is a decision record, not an experimental result and not authorization to run Stage 1 or any later branch.

## Access verification and review boundary

The review was conducted in the already logged-in Ji Wang project ChatGPT Pro conversation. Pro reported `ACCESS_VERIFIED` for the private repository `JiWang-Unige/coding-rna` at commit `fe38a675f2d4a0e9349046c226a5679b850a82b0` (`Archive M25R result and research roadmap`). Its opening verification reproduced facts that were not included in the access-verification prompt:

- exact `docs/28_current_research_state_2026-09-01.md` title: `Current research state — 2026-09-01`;
- M25R Slurm job `12116383`, `COMPLETED`, exit `0:0`, experiment status `STOP_M25_BRANCH`;
- best-ranked non-admissible tuple: interval F1 `0.1204141`, chain F1 `0.3249883`, FPR `0.0124683`, predicted-gene-count ratio `0.3252713`.

This is consistent with direct access to the authorized private snapshot and was the access check used for this review; no independent OAuth/connector audit log was available to Codex. The snapshot deliberately does not contain the three large checkpoints or the full `5,625`-row grid, so Pro could audit code, configs, summaries and logic but could not independently rerun checkpoint inference from GitHub. It did not read or search Setaria annotation. Setaria remains entirely unmeasured.

The principal repository evidence reviewed was:

- `docs/28_current_research_state_2026-09-01.md`;
- `docs/29_end_to_end_technical_roadmap_2026-09-01.md`;
- `reports/M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0/` compact result;
- M24, M19 and M20 reports;
- M25/M25R execution and repair plans;
- `src/foundation_probe/train_generanno_structural_heads.py`;
- `scripts/eval_m25_structure.py` and focused tests.

## Joint current verdict

| Question | Joint verdict | Evidence boundary |
|---|---|---|
| Is there a positive machine-learning result? | Yes, but only at coarse adaptation level. M19 shows that the released GENERanno 1.2B CDS-preview package can be adapted to a low-FP coarse coding/gene-body regime on Arabidopsis/rice development data. | M19 is not strand-aware exact CDS or a complete gene model. |
| Did M24 change the scientific target? | Yes. It showed that coarse gbF1 does not imply exact structure: M19 exact interval/coordinate-pseudo-chain results are far below released callers on the same held-out ranges. | M24 is the direct reason to require exact interval, chain/gene, strand, phase and independent validity. |
| Is the initial M25 result usable? | No. M25 is implementation-invalid because empty structural masks caused non-finite boundary loss. | M25 checkpoints cannot support scientific inference. |
| Is M25R a valid result? | Yes. It is a valid development no-go for the frozen combined representation, supervision and decoder contract. | It does not isolate GENERanno backbone quality. |
| Is the project ready for a Nature Communications performance claim? | No. The project remains in `structural feasibility diagnosis`. | No development-passing complete annotator, Setaria result, paired blind baseline table, multi-species transfer or mechanism proof exists. |

The central interpretation is therefore narrow but useful: CDS-specialized foundation-model signal can support coarse low-FP adaptation, but the current M25R route has not converted that signal into sufficiently many exact complete gene models.

## Why M25R is a decisive no-go

All `5,625` validation tuples pass intergenic FPR `<=0.020`, while none passes gene-count ratio `0.80–1.20`; the full observed count range is only `0.0953–0.3305`. This rules out excess intergenic overcalling as the immediate gate failure. The system is obtaining low FPR by emitting far too few complete genes.

The implementation explains why the boundary-threshold grid could not rescue recovery. `validation_candidates()` first runs the complete decoder with boundary thresholds set to zero. Region blocks, motif availability, phase, ORF and internal-stop rules have already discarded models before `_filter_models()` applies the four searched boundary thresholds. That filter can only remove existing models; it cannot create missing ones.

Training loss falls from `0.8205` to `0.6669` to `0.6177`, while exact recovery and gene count worsen across epochs. This supports objective-versus-complete-gene misalignment and increasing conservatism. It does not identify whether the cause is the backbone, LoRA, structural heads, targets/loss, windowing or decoder.

The development metric `structurally_valid_complete_fraction = 1.0 if predictions exist else 0.0` is tautological. It does not mean that all emitted transcripts are valid. This flaw does not overturn the no-go because every tuple already fails gene count, but future diagnostics must audit each emitted transcript independently.

## Revised immediate decision

`GO` only for planning and, after separate execution approval, implementing:

`M25R-DEV-REDECODE-ERROR-DECOMPOSITION`

`NO-GO` for new training, threshold or decoder search, SegmentNT promotion, extra seeds, Setaria inference or annotation access, baseline execution, and external provenance contact.

The exact Stage 1 contract is now in `docs/29_end_to_end_technical_roadmap_2026-09-01.md`. Its essential constraints are:

1. replay exactly one original highest-ranked tuple per epoch, with no new ranking or search;
2. reproduce all four frozen aggregates within `1e-5`;
3. use the exact `6,450` reference-chain universe and the complete M25R training sample set actually observed as `1,536` windows, without constructing a new subset;
4. calculate fixed transition, motif-candidate and truth-assisted exact-chain upper bounds;
5. trace the real decoder order and reconcile every reference, candidate lineage and emitted model;
6. compare train and validation raw region, boundary and phase outputs;
7. independently audit every emitted transcript;
8. write the report and stop.

Stage 1 can identify a learned-emission/candidate, structural-head, development-generalization, decoder or mixed failure. It cannot by itself call a `backbone representation bottleneck` because the representation, LoRA, heads, targets, loss and training budget remain coupled.

## Hard requirements versus operational heuristics

The second review corrected the first review's overstatement of numerical gates:

- `<=1e-5` aggregate reproduction, complete accounting, `100%` validity-audit coverage, unchanged model/decoder and zero Setaria reads are hard diagnostic-integrity requirements.
- A single stage explaining `>=50%` of missing chains is only an operational rule for using the word `dominant`.
- Truth-assisted headroom `>=0.15` is only an operational information-value rule for proposing one experiment.
- Development interval `>=0.60` and chain `>=0.50` were withdrawn; they were not in the frozen M25R contract and are not justified by the non-paired M24 baseline ranges.
- Strongest-baseline `+2 pp` is a publication-planning target, not a Nature Communications rule.

The next branch's absolute FPR/count/validity guardrails and exact structural competitiveness gate must be frozen only after Stage 1 and paired A/rice development baseline scoring, and before any branch experiment is submitted.

## Post-Stage-1 decision tree

Stage 1 always terminates before this table can authorize work.

| Observed pattern | Accurate interpretation | One discriminating next experiment to consider |
|---|---|---|
| Most chains disappear before usable region/CDS-run/transition candidates | learned-emission/candidate-generation bottleneck | One released SegmentNT native-long-context, fixed-context, development-only emission probe; no fine-tuning and no full caller. |
| Transition/motif reachability is high, but boundary/phase heads are weak on train and validation | structural-head/supervision-fit bottleneck | Freeze backbone, LoRA and region head; one boundary/phase-head-only fit. |
| Raw heads fit training data but fail validation | development data/distribution bottleneck | One species-balanced fit with all other contracts fixed. |
| Raw candidate and truth-assisted ceilings are high, but final emission is low | decoder/grammar bottleneck | Existing scores, one implicated grammar rule, read-only re-decode. |
| No dominant stage or no useful single-stage headroom | mixed bottleneck | Stop incremental M25-family patching and redesign the task representation. |

A true backbone-versus-supervision test would require a later matched two-arm experiment in which GENERanno and SegmentNT representations share the same chromosomes, training examples, targets, head capacity, optimization budget, decoder and evaluator. Even that would compare representation packages, not isolate one pretraining ingredient.

## Baseline and blind-test timing

The agreed order is:

`Stage 1 -> review -> A/rice same-range baseline re-score -> freeze one branch gate -> one branch experiment`

The development baseline re-score uses the exact M25R validation chromosomes, reference release, longest-CDS primary-transcript policy and evaluator. It defines structural competitiveness but must not select a checkpoint or decoder.

Only after a later development candidate passes its frozen gate should the project:

1. freeze candidate config/checkpoint/decoder/ablation and baseline versions/commands;
2. generate candidate, ablation, Helixer, Tiberius and ANNEVO Setaria predictions from FASTA while annotation remains closed;
3. record the exact commit, configs and output paths;
4. open Setaria annotation once and score all predictions together;
5. end the branch on blind failure without target-label tuning.

No additional checksum or prediction-hash layer is added for Stage 1 because it would not change a decision. The existing Git commit, resolved config, checkpoint/grid-row identity, exact paths and read-only execution record are sufficient for this diagnostic.

## SegmentNT boundary

SegmentNT is conditional, not the next mainline by default. If Stage 1 indicates learned-emission/candidate loss, the sole proposed probe is the released SegmentNT model at one fixed native long-context setting on the same A/rice development chromosomes. It reports exon/intron and donor/acceptor emissions plus ordered splice-event reachability. It does not report phase, exact start/stop, full coding-gene F1 or structural GFF3 validity.

The existing independent 6-kb cache only supports the statement that this extraction has moderate exon signal and weak donor/acceptor signal. It neither validates nor rejects a native-long-context SegmentNT route. The exact released artifact and context must be verified before any future probe is frozen.

## Provenance and permissible claim

GENERanno generic pretraining and CDS-specific post-training require separate species/accession audits. Until an exact manifest is obtained, status remains `overlap_unknown`. Preparing an audit checklist or a contact draft is internal work; sending a message to GenerTeam requires separate authorization.

If the manifest remains unavailable, the acceptable claim is:

> We adapted a publicly released CDS-specialized GENERanno checkpoint for direct structural gene annotation and evaluated downstream generalization to species excluded from our task-specific training and calibration. The species and accession composition of the upstream generic pretraining and CDS-specific post-training corpora is not fully disclosed; therefore, these experiments do not establish pretraining-held-out, no-overlap or zero-shot generalization.

Claims that must then be removed include `clean pretraining-held-out`, `zero-shot`, `no-overlap unseen species`, `Setaria unseen to the backbone`, universal gene grammar without prior exposure, and attribution of the entire 1.2B-versus-0.5B difference to CDS post-training.

## Route to a Nature Communications decision

Necessary evidence:

- a development-passing exact structural model under a pre-registered gate;
- frozen candidate, decoder, ablation and one-time untouched Setaria blind evaluation;
- Helixer, Tiberius and ANNEVO on identical ranges and evaluator;
- final-architecture three-seed robustness and a multi-species held-out panel;
- a matched mechanism ablation;
- exact interval, chain/gene, strand, phase, FPR, count and independent validity;
- long-gene, tile-edge, non-canonical and repeat-rich strata;
- independent RNA, protein-homology and synteny support;
- runtime, resource and reproducibility evidence;
- provenance closure or an explicitly narrowed adaptation claim.

Enhancing evidence includes a larger species panel, low-homology subsets, blinded disagreement loci, reference-correction analysis and orthogroup-level utility. Low-value work to delete includes further coarse gbF1 threshold tuning, GENERanno 0.5B performance development, generic CRF/Tversky sweeps, early extra seeds, premature MoE/context expansion, Setaria access, organelles, unplaced scaffolds, alternative isoforms and BUSCO as a substitute for exact structure.

The final paper does not need an arbitrary `+2 pp` accuracy margin, but it does need a cross-species paired advantage with uncertainty support and no unacceptable FPR/count cost, or an equally strong pre-registered Pareto advantage with independent biological relevance. A single-species metric win is insufficient.

## Codex synthesis and recommended decisions

Codex agrees with Pro's second-round corrections and recommends the following defaults:

1. accept the revised `docs/29` as the current planning contract;
2. preserve Setaria as a one-time blind reserve and prohibit every Stage 1 Setaria read;
3. if Stage 1 implementation is later approved, limit it to one experiment-specific diagnostic script, one focused test and one report; do not modify the trainer or create a framework;
4. use the complete frozen M25R training sample set for read-only raw-head inference, not a new subset;
5. prepare provenance questions internally, but do not contact GenerTeam without separate approval.

No new experiment, baseline job, SegmentNT job, Setaria operation or external provenance message was performed during this review.
