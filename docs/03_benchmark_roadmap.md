# Benchmark + Roadmap

> Finalized by `$benchmark-roadmap` on 2026-06-09 after user confirmation.
> This document defines the benchmark contract, dataset readiness plan, and technical roadmap.
> This stage does not download data; data preparation happens inside a concrete `$reproduce-baselines` or `$goal-prompt` iteration.

## 1. Target task

- Input: raw eukaryotic genomic DNA sequence, with no RNA-seq, protein homology, transcript evidence, or species-specific external annotation features as model input.
- Output: structure-aware per-base labels for ab initio protein-coding gene annotation. Formal screen/full runs must support strand-aware CDS/intron/intergenic labels plus phase/splice-border/start-stop or an equivalent representation that can be decoded into valid gene models.
- Evaluation setting: held-out species or held-out clade evaluation with chromosome/species/homology-aware leakage checks. Random sequence split is disallowed.
- Explicitly out of scope: RNA folding / RNA 3D, evidence-based RNA-seq/protein-homology pipelines, transcript-only classification, promoter/enhancer/TF-binding/ncRNA-only tasks, binary-only gene-body/intergenic as a publication claim.

Binary `gene-body/intergenic` is allowed only as a smoke or fast sanity signal. It is not sufficient for the main paper task.

## 2. Metrics

### Primary metric

- Name: `constrained_gene_body_F1`
- Implementation: convert reference and prediction GFF/GTF to protein-coding gene-body base masks under the frozen transcript-collapsing rule, compute gene-body precision/recall/F1, then accept the score only if the intergenic false-positive guardrail passes.
- Direction: higher_is_better
- Profile-aware guardrail after `BASE-TIBERIUS-MINISMOKE` tri-review/pivot: smoke/screen use `intergenic_FPR <= 0.02` while the M1 evaluator is being frozen; full/scale claim candidates keep the stricter `intergenic_FPR <= 0.01`.
- Sensitivity report: always report the same primary metric under `intergenic_FPR <= 0.005`, `<= 0.01`, and `<= 0.02`.
- Claim restriction: screen/smoke values can never claim SOTA; full/scale values can claim only after strict exceedance of `sota_benchmark` and comparability review.

The exact implementation should be created during M1 and should cross-check:
- ANNEVO repo's `gffcompare -r <ref> <pred> --no-exon-merge --strict-match` notes.
- Tiberius config tables using `gffcompare v0.12.10`.
- Helixer repo scripts, especially annotation-level gffcompare and `scripts/accs_genic_intergenic.py`.

### Secondary metrics

- Gene/locus-level F1 after strict transcript/locus matching.
- Exon F1 and CDS F1.
- Intron-chain F1.
- Coding phase F1.
- Splice donor/acceptor precision, recall, F1.
- Boundary error distribution.
- Predicted gene count ratio versus reference; warning threshold `>1.25x`.
- Nucleotide gene-body F1 drop versus anchor; warning threshold `>0.03`.
- BUSCO completeness as support-only, never primary, because it can hide false-positive genes/exons.

## 3. Three-layer gates

| Gate | Threshold | Trigger |
|---|---:|---|
| `primary_progress_gate` | screen profile: `constrained_gene_body_F1 > screen_anchor + 0.01`, smoke/screen `intergenic_FPR <= 0.02` during M1 evaluator freeze, nucleotide gene-body F1 drop <= 0.03 when anchor exists, predicted gene count ratio <= 1.25 | Progress signal only; not a SOTA claim |
| `sota_claim_gate` | full/scale profile: strict `constrained_gene_body_F1 > sota_benchmark`; equal is failure | Claim candidate only if comparability contract is all PASS and human gate passes |
| `review_decision_gate` | any primary pass without strict SOTA exceedance, any gap < 0.05, any guardrail near/failing threshold, or any metric/preprocessing mismatch | Force `$tri-review` + `$pivot` |

Anti-tuning rule: if the gap to the active anchor is `>= 0.05`, hyperparameter tuning is disallowed; pivot to a different architecture axis.

## 4. SOTA reference table

| Model | Dataset | Split | Metric | Value | Source | Comparable? | Notes |
|---|---|---|---|---:|---|---|---|
| ANNEVO paper | RefSeq/Ensembl-derived broad eukaryotic benchmark; 566 species overall; 12 model-species table | Test species intentionally absent from train/val; exact supplement species table still needs expansion | nucleotide mean F1 | 0.92 | `refs/dossiers/annevo-2026.md`; Nature Methods 2026 | Partial | Provisional published SOTA anchor, but metric is not exactly our primary |
| ANNEVO current repo | 12 model species across mammal/insect/fungi/plant/bird/fish groups | Repo notes say model species absent from ANNEVO train/val | exon recall/precision; locus recall/precision; BUSCO | locus R/P 76.3/74.3; exon R/P 91.4/90.2 | ANNEVO repo @ `37bdd9a`; `docs/performance.md` | Likely after pinning | Best candidate for deriving our frozen full-eval `sota_benchmark` |
| Tiberius 2024 | 37 mammalian RefSeq species | Held-out human/cow/beluga; close clades excluded from train | gene F1; exon F1 | avg gene 55.1; exon 89.7 | `refs/dossiers/tiberius-2024.md`; local PDF | Mammal-only | First baseline to reproduce and main structured-decoder mechanism source |
| Tiberius current configs | Multi-clade configs in repo | Per-clade held-out species listed in configs | gffcompare gene/exon F1 | vertebrates avg gene 55.3; exon 85.8 | Tiberius repo @ `8c49fd0` | Partial | Useful for broad-lineage screen anchor, but not 2024 paper claim |
| Helixer 2025 | RefSeq/Phytozome lineage datasets | Species-level lineage splits | phase F1; exon/intron/intron-chain/transcript F1 | vertebrate phase F1 0.8829; feature table varies by lineage | `refs/dossiers/helixer-2025.md` | Partial | Broad-lineage open baseline; original metrics must not be mixed with ANNEVO comparisons |
| SegmentNT 2025 | Human GENCODE/ENCODE-style element segmentation | chr20/21 test, chr22 validation, other chromosomes train | MCC/auPRC/Jaccard/F1/SOV | no comparable full gene-caller value | `refs/dossiers/segmentnt-2025.md` | No for claim | Probe or feature source only |
| GENERanno eukaryote CDS preview | Public example fly parquet and HF model card; full eukaryote benchmark unclear | Unknown | model card lists F1 but no public table | no comparable value | `refs/dossiers/generanno-2025.md` | No for claim | Probe only until schema/split/benchmark are public |

## 5. Comparability contract

| Model | Dataset version | Split scheme | Metric impl | Preprocessing | Weights version | Test-time inference | Verdict |
|---|---|---|---|---|---|---|---|
| ANNEVO paper | Need supplement expansion | Need exact train/val/test species | Paper nucleotide mean F1; repo gffcompare values differ | Excludes non-coding transcripts/UTRs; invalid start/stop and intron lengths filtered | Paper frozen weights unclear; repo `saved_model/*.pt` available | Viterbi decoder | PENDING; likely published anchor after value mapping |
| ANNEVO repo/current | Repo-pinned 12 species table | Test species absent from train/val per repo | gffcompare strict locus/exon; must add our gene-body mask metric | Same as repo notes | repo @ `37bdd9a` | Repo inference + Viterbi | PENDING; preferred `sota_benchmark` reproduction target |
| Tiberius 2024 | RefSeq mammal data; exact accessions need supplement | Human/cow/beluga held out; close clades excluded | gffcompare v0.12.10 gene/exon F1 plus our mask metric | Softmasking with RepeatModeler2/RepeatMasker/TRF | paper snapshot vs current weights must be separated | Forward/reverse strand, Viterbi, tile boundary second pass | PARTIAL VERIFIED by mini-smoke; full paper split still pending |
| Helixer 2025 | RefSeq/Phytozome/Zenodo model versions pinned | Species-level lineage splits | phase/feature F1; repo genic/intergenic script | Lineage-specific preprocessing; UTR/CDS/intron classes | Zenodo 17850139 v2 unless frozen otherwise | HelixerPost HMM | PENDING; broad-lineage baseline |
| SegmentNT | HF model card; no full gene benchmark | Human chromosome split | element segmentation metrics only | 6-mer path does not handle `N` cleanly | HF `InstaDeepAI/segment_nt` | segmentation head only | NON-COMPARABLE for claim |
| GENERanno preview | Dataset card incomplete | Unknown | CDS/token F1 unclear | Unknown | HF preview model | token classification CLI | NON-COMPARABLE for claim |

## 6. Two-layer anchors

| Anchor | How to obtain | Value status | Use |
|---|---|---|---|
| `screen_anchor` | User-confirmed option A: run `max(Tiberius-like, Helixer-like, ANNEVO-light/available)` under one unified small-budget screen protocol with the same sample fraction, split, epochs, patience, seed policy, preprocessing, and metric implementation | Pending M1 | Track A screen compares only against this value; never claim |
| `sota_benchmark` | User-confirmed option A: reproduce ANNEVO under the frozen full evaluation and use the ANNEVO-compatible `constrained_gene_body_F1` as the published-anchor value; paper value remains provenance, not directly substituted if metric differs | Pending M1/M2 baseline reproduction | Full/scale strict SOTA claim only |

Do not put ANNEVO paper nucleotide mean F1 `0.92` directly into `sota_benchmark` unless the M1 metric audit proves it is the same quantity as `constrained_gene_body_F1`.

## Baseline Reproduction Report: Tiberius mini-smoke 2026-06-10

- Experiment ID: `BASE-TIBERIUS-MINISMOKE`
- Status: semantic success for mini-smoke; not a screen_anchor and not a SOTA claim.
- Reported metric: Tiberius repo integration test requires CDS F1 >= 0.75 and transcript-chain F1 >= 0.28 on bundled `Panthera_pardus` test data.
- Reproduced: CDS exact F1 = 0.8594; transcript-chain exact F1 = 0.3124; both integration thresholds pass.
- Our provisional metric on the same mini data: unconstrained gene-body F1 = 0.9196, intergenic FPR = 0.0187, so `constrained_gene_body_F1 = 0.0` under the roadmap guardrail `intergenic_FPR <= 0.01`.
- Metric implementation VERIFIED for mini-smoke: `refs/repos/tiberius-2024/tests/integration_tests/test_prediction.py` parses CDS features from GTF, compares exact `(chrom,start,end,strand)` CDS intervals, and compares each transcript as a frozenset of exact CDS intervals.
- Dataset rawness VERIFIED for mini-smoke: bundled repo data at `refs/repos/tiberius-2024/test_data/Panthera_pardus/inp.tar.gz`; extracted genome has 42,253,745 bp and reference GTF has CDS/intron/start/stop features but no `gene`/`transcript` features.
- Split scheme VERIFIED for mini-smoke: not a claim split; this is Tiberius bundled test data used for workflow validation only.
- Preprocessing / test-time: Tiberius 2.0.5 launcher, official Singularity image `tiberius_2.0.5.sif`, `mammalia_softmasking_v2`, `seq_len=259992`, `batch_size=8`, dual-strand current CLI inference.
- Output paths: `outputs/BASE-TIBERIUS-MINISMOKE/tiberius_prediction.gtf`, `outputs/BASE-TIBERIUS-MINISMOKE/metrics/metrics.json`, `outputs/BASE-TIBERIUS-MINISMOKE/logs/srun_tiberius_predict_4090_bound.log`.
- Report-vs-reproduce gap: no gap for integration threshold; paper-level human/cow/beluga metrics remain unreproduced.
- Comparability implication: this validates the container and exact-CDS/transcript-chain mini metric, but `screen_anchor` still requires M1 unified Tiberius-like/Helixer-like/ANNEVO-light screen runs.

## M1 Evaluator Freeze: gene-body/FPR screen contract 2026-06-10

- Contract: `docs/11_evaluator_contract.md`
- Implementation: `scripts/eval_gene_body_mask.py`
- Validator support: `scripts/validate_goal.py` supports `threshold_by_profile` and profile-scoped guardrails.
- Regression tests: `tests/test_eval_gene_body_mask.py`, `tests/test_validate_goal_profiles.py`
- Screen baseline manifest: `configs/m1_screen_baselines.yaml`
- Data manifest: `configs/m1_data_manifest.yaml`
- Test status: `pytest -q tests/test_eval_gene_body_mask.py tests/test_validate_goal_profiles.py` passes.
- Gene-body mask rule: symmetric transcript-span masks from `CDS/exon/intron/start_codon/stop_codon`, grouped by `transcript_id` with `gene_id` fallback, merged per sequence before base overlap.
- FPR rule: `intergenic_FPR = predicted_intergenic_false_positive_bases / reference_intergenic_bases`.
- Profile thresholds: smoke/screen `<=0.02`; full/scale `<=0.01`; sensitivity always reports `0.005/0.01/0.02`.
- Count guardrail rule: `predicted_gene_count_ratio_vs_reference = unique predicted gene_id / unique reference gene_id`; transcript multiplicity is separate as `predicted_transcript_count_ratio_vs_reference`.
- Tiberius mini-smoke evalfix: under smoke threshold, `constrained_gene_body_F1=0.9196`; under full threshold, the same artifacts remain `not_yet` because `intergenic_FPR=0.0187 > 0.01`.
- Next non-training step: fill the data manifest with frozen source URLs/accessions/checksums, then verify Helixer and ANNEVO runner environments before screen submissions.

## 7. Dataset readiness plan

| Dataset | Purpose | Required by | Timing | URL | Size | Split source | Hash needed? | Notes |
|---|---|---|---|---|---:|---|---|---|
| Tiberius mammal RefSeq protocol | First metric/split/preprocessing reproduction; screen mechanics | M1 | now | Tiberius paper/repo; RefSeq accessions from supplement/configs | TBD | Paper held-out species + close-clade exclusions | yes | Learn gffcompare, transcript collapsing, softmasking, tile inference |
| ANNEVO 12 model-species benchmark | Published-anchor reproduction and full-eval candidate | M1-M2 | now | ANNEVO paper/repo/supplement | TBD | ANNEVO supplement/repo model species tables | yes | Freeze paper vs repo/current release separately |
| Helixer lineage benchmark subset | Broad-lineage baseline and screen anchor component | M1-M2 | now/on-demand | Helixer paper/repo/Zenodo | TBD | `docs/model_overview.md`; Zenodo version | yes | Use for cross-lineage behavior, not as primary published anchor unless ANNEVO fails comparability |
| SegmentNT chr20/21/22 or selected windows | Foundation signal probe | Path 1 | on-demand | HF `InstaDeepAI/segment_nt` | model/data TBD | HF chromosome split | yes | Probe logits/embeddings, not claim |
| GENERanno fly/example CDS data | Foundation CDS behavior probe | Path 1 | on-demand | HF model/card and repo examples | TBD | Unknown public split | yes | Probe only; stop if weights/schema download fails |
| Wider eukaryote held-out clades | Phase 8 generalization | M5 | later | RefSeq/Ensembl/Phytozome as selected | TBD | Frozen clade split + homology dedup | yes | Only after primary full/scale success |

## 8. Technical Roadmap

### 8.1 SOTA weaknesses (mechanism-level)

| SOTA model | Weakness | Mechanism | Evidence | Exploitable? |
|---|---|---|---|---|
| ANNEVO | Paper/repo metric drift and MoE novelty space already occupied | Clade/evolution expert routing plus Viterbi; published metric does not directly equal our FP-controlled gene-body metric | `refs/dossiers/annevo-2026.md` | Yes; attack evaluation alignment and decoder/objective, not generic MoE |
| Tiberius | Mammal-centric 2024 claim and limited broad-eukaryote scope | Strong 15-class + CCE-F1 + differentiable HMM, but trained/evaluated mainly on mammals in 2024 paper | `refs/dossiers/tiberius-2024.md` | Yes; generalize structured decoder beyond mammal protocol |
| Helixer | Two-stage objective mismatch | Base-wise neural predictor is not trained end-to-end through final HMM decoder | `refs/dossiers/helixer-2025.md` | Yes; end-to-end or segment-level decoder can improve consistency |
| SegmentNT | Not a full gene caller | Foundation segmentation logits are not constrained into valid gene models | `refs/dossiers/segmentnt-2025.md` | Yes as representation source |
| GENERanno preview | Public benchmark/schema incomplete | CDS token classifier does not expose full eukaryotic gene-model protocol | `refs/dossiers/generanno-2025.md` | Yes as behavior probe only |

### 8.2 Differentiated paths

#### Path 1: Foundation Features + Semi-CRF Segment Decoder

- Hypothesis: SegmentNT/GENERanno-style pretrained signals contain useful CDS/splice/gene-body evidence, but need a segment-level structured decoder to become valid gene models.
- Architecture change: frozen or lightly adapted foundation embeddings/logits plus a trainable semi-CRF segment decoder over gene-body/intergenic/CDS/intron/splice-border states.
- Why this attacks SOTA weakness: converts non-claimable foundation segmentation/CDS signal into a gene grammar-compatible predictor, targeting ANNEVO/Tiberius/Helixer's structured consistency while testing the foundation-model route cheaply.
- Track A screen design: sample_fraction initially 5-10%, screen epochs 3-5, patience 2, seed 1; expected walltime <=12h per candidate on private RTX3090 or shared high-memory GPU if needed.
- Track B scale-up rule: promote if screen exceeds `screen_anchor + 0.01` and all guardrails pass; expand to 25-50% data, 2-3 seeds, longer patience.
- Required data: on-demand for SegmentNT/GENERanno probes; M1 screen data is now.
- Expected gain: better interval/gene consistency if foundation signals reduce local ambiguity at splice/CDS boundaries.
- Risk: foundation logits may be weak or non-comparable; if gap `>=0.05` and key label probes fail, downgrade and do not large-scale fine-tune.
- Failure detection: no CDS/splice/boundary signal, intergenic FPR violation, gene-count inflation, or no improvement over Tiberius-like screen anchor.

#### Path 2: Tiberius-Like Structured Baseline with CRF/Semi-CRF Head Replacement

- Hypothesis: keeping Tiberius-style 15-class supervision and metric-aligned losses while replacing/augmenting HMM with CRF or semi-CRF transitions can improve gene-level consistency and boundary precision.
- Architecture change: decoder/head replacement; compare differentiable HMM, linear-chain CRF, semi-CRF, and constrained Viterbi over the same backbone and label schema.
- Why this attacks SOTA weakness: directly tests whether the structured decoder, not just neural backbone, is the main bottleneck.
- Track A screen design: same M1 small-budget split; 3-5 epochs, patience 2, seed 1; no tuning-only variants.
- Track B scale-up rule: promote if `screen_anchor + 0.01` and guardrails pass; otherwise pivot architecture axis if gap `>=0.05`.
- Required data: now.
- Expected gain: improved boundary legality, fewer fragmented/fused genes, better locus/gene-level F1.
- Risk: HMM constraints may already be near-optimal; CRF/semi-CRF may add compute without gain.
- Failure detection: gene-body F1 flat, gene-level F1 flat/down, invalid transitions, or runtime/memory dominated by decoder.

#### Path 3: Long-Context / RMT Backbone + Structured Decoder

- Hypothesis: long genes, long introns, and tile-boundary artifacts cause fragmentation/fusion errors that local CNN/LSTM windows cannot solve robustly.
- Architecture change: RMT/hierarchical memory/state-space or other long-context backbone feeding the same structured decoder family.
- Why this attacks SOTA weakness: Tiberius and Helixer use tiled inference; long-context memory may reduce boundary stitching and long-intron false calls.
- Track A screen design: include long-gene and long-intron strata in the screen set; 3-5 epochs, patience 2, seed 1; expected walltime <=12h if memory allows.
- Track B scale-up rule: promote only if aggregate screen passes and long-gene strata improve without intergenic FPR inflation.
- Required data: now for stratification; extra long-gene OOD sets later.
- Expected gain: fewer gene fragmentation/fusion failures and stronger gene-level F1 in long-locus strata.
- Risk: high compute and unclear benefit if errors are mostly decoder/objective-driven rather than context-limited.
- Failure detection: no stratum-specific gain, memory blow-up, or worse short-gene calibration.

#### Path 4: FP-Aware Objective + Biological Constraint Decoder

- Hypothesis: intergenic false positives and biologically illegal structures are primary blockers for a genome-usable predictor.
- Architecture change: FP-aware loss/hard-negative mining plus ablated biological legality constraints: GT/AG splice motifs, start/stop codon consistency, phase legality, minimum exon/intron lengths, no in-frame stop codons in CDS.
- Why this attacks SOTA weakness: Tiberius ablations show CCE-F1 and border/phase classes matter; BUSCO/local F1 can hide FP-heavy gene inflation.
- Track A screen design: same M1 screen set; compare objective/constraint module against the same backbone/decoder.
- Track B scale-up rule: promote if it improves or preserves primary F1 while reducing intergenic FPR and gene-count inflation.
- Required data: now.
- Expected gain: lower FP under whole-genome inference and cleaner gene models.
- Risk: can become an ad hoc post-filter that sacrifices recall or overfits motif assumptions.
- Failure detection: primary gain caused only by deleting predictions, recall collapse, invalid comparison to baseline, or no ablation-supported benefit.

### 8.3 Track A / Track B strategy

#### Track A: small-sample parallel architecture screening

- Default sample_fraction: 5-10% of the frozen M1 screen data, stratified by species/clade and gene/intergenic content.
- Default epochs: 3-5.
- Default patience: 2.
- Default seeds: 1 for screen; add seeds only near promotion.
- Candidate set: user selected 4 directions: Path 1, Path 2, Path 3, Path 4.
- Execution rule: project hard limit remains `max_parallel_directions=3`, so the 4 directions run as two waves or require an explicit later config revision.
- Claim policy: never claim SOTA from Track A / screen.

#### Track B: scale-up promising candidates

- Promotion criteria from Track A: user-confirmed rule A, `constrained_gene_body_F1 > screen_anchor + 0.01` and all guardrails pass.
- Data expansion rule: expand promoted candidate to 25-50% data first, then full data if the gain holds.
- Epoch / patience expansion: increase patience and walltime only after M2 promotion; avoid hyperparameter search when gap `>=0.05`.
- Seed expansion: 2-3 seeds for Track B; more only near claim.
- Failure interpretation: if screen gap `>=0.05`, tune is disallowed and the next pivot must change architecture axis.

#### Parallel progression rule

Track A continues exploring orthogonal architecture changes while Track B scales candidates that passed promotion criteria. No more than three directions run concurrently unless `ACTIVE_GOAL.json` and `cluster_config.yaml` are explicitly revised.

### 8.4 Milestones

| ID | Milestone | Threshold | Track | Expected runs | Completion evidence |
|---|---|---|---|---|---|
| M1 | Freeze eval implementation and establish `screen_anchor` | unified screen metric implemented; max of Tiberius-like/Helixer-like/ANNEVO-light recorded as `screen_anchor`; guardrail metrics parseable | baseline/screen | 3 baseline screen runs | result-log + metrics JSON + `ACTIVE_GOAL.json` screen_anchor update |
| M2 | Reproduce ANNEVO-compatible full benchmark value | ANNEVO full/frozen eval produces `constrained_gene_body_F1`; value recorded as `sota_benchmark` if comparability passes | baseline/full | 1-2 full runs | result-log + comparability contract all PASS or documented gaps |
| M3 | One architecture path hits primary_progress_gate | screen score > `screen_anchor + 0.01` and guardrails pass | Track A | 4 screen directions, executed <=3 concurrently | result-log + tri-review + pivot |
| M4 | Best M3 candidate holds on larger data | gain persists on 25-50% data and 2-3 seeds without FP/gene-count inflation | Track B | 1-3 full-ish runs | full result-log + pivot |
| M5 | First strict exceed on full/scale | full/scale `constrained_gene_body_F1 > sota_benchmark`, equal fails | Track B | 1-3 full/scale runs | comparability all PASS + human gate before claim |
| M6 | Phase 8 comprehensive superiority | primary, secondary, OOD, robustness, cost, ablation all reviewed | generalization | varies | `$generalization` output |

### 8.5 Risk dial

| Path | Risk | Reasoning | User priority |
|---|---|---|---|
| Path 1 Foundation + semi-CRF | High | New combination and possible foundation-signal mismatch, but largest novelty upside | Primary |
| Path 2 Structured head replacement | Medium | Mechanism-close to Tiberius; strongest controlled architecture test | Co-primary/control |
| Path 3 Long-context/RMT | High | Could matter greatly for long genes but compute-heavy and failure-mode dependent | Included in first Track A set |
| Path 4 FP-aware objective/constraints | Medium | Directly attacks user's FP concern; needs strict ablation to avoid heuristic overfitting | Co-primary |

### 8.6 Resource budget

| Stage | Resource profile | Compute estimate | Wall-clock estimate | Notes |
|---|---|---|---|---|
| M1 baseline screen anchor | screen | 3 reference runs, 1 GPU each | <=12h each if small screen succeeds | Must use `$smart-sbatch`; no claim |
| M2 ANNEVO-compatible full baseline | full | 1-2 full inference/eval or light reproduction runs | likely submit-and-handoff if >12h | Sets `sota_benchmark` if comparable |
| M3 Track A portfolio | screen | 4 architecture directions, but <=3 concurrent | one or two waves, <=12h per screen target | User selected 4 directions; hard parallel cap remains 3 |
| M4 Track B scale-up | full | 1-3 promoted candidates, 2-3 seeds | multi-day possible | Use private-teodoro-gpu when possible |
| M5 final claim candidate | scale | multi-seed full/scale evidence | multi-day | Human gate before claim |
| M6 generalization | scale/generalization | varies by OOD matrix | varies | After strict primary exceed only |

## 9. User-confirmed technical choices

- `screen_anchor`: option A, max of Tiberius-like, Helixer-like, and ANNEVO-light/available under one unified small-budget protocol.
- `sota_benchmark`: option A, ANNEVO method reproduced under our frozen full evaluation as `constrained_gene_body_F1`.
- First Track A candidate set: option C, all 4 directions included; execute in waves because project hard limit is 3 parallel directions.
- Track B promotion rule: option A, screen exceeds `screen_anchor + 0.01` and guardrails pass.
- Main route: Path 1 as primary, Path 2 as reproduction/structured control, Path 4 as low-FP constraint route; Path 3 included early because user selected 4-direction Track A.

## 10. Open uncertainties

- [ ] Expand ANNEVO supplementary species tables and freeze the exact species/split/provenance for M1/M2.
- [ ] Decide transcript-collapsing rule before metric implementation: longest CDS per gene is the default candidate because ANNEVO/Tiberius both use longest/strict transcript handling in evaluation notes.
- [ ] Decide whether softmasking is part of model input, preprocessing only, or excluded from our architecture; Tiberius uses a masked-repeat track.
- [ ] Build and test a single metrics JSON schema consumed by `scripts/validate_goal.py`.
- [ ] Keep GENERanno and SegmentNT as non-claim probes unless they are wrapped in a comparable structured gene caller.

## 11. TODO for docs/05_todo.md

- [ ] Run `$reproduce-baselines` for M1 before writing our own model code.
- [ ] Implement/freeze the `constrained_gene_body_F1` evaluator and intergenic FPR guardrail.
- [ ] Establish `screen_anchor` and write the value to `ACTIVE_GOAL.json`.
- [ ] Reproduce ANNEVO-compatible full eval and write `sota_benchmark` to `ACTIVE_GOAL.json`.
- [ ] Only after M1/M2, use `$goal-prompt` for a concrete Track A portfolio iteration.
