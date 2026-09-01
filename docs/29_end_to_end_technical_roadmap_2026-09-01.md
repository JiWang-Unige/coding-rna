# End-to-end technical roadmap — Pro-reviewed 2026-09-01

Status: reviewed decision contract. This document does not authorize an experiment or a Slurm submission. The sole scientifically admissible next analysis is Stage 1; its implementation and execution still require an explicit follow-up decision.

## Target outcome

Build a raw-FASTA, primary-nuclear-chromosome gene annotator that emits complete protein-coding gene models with exact CDS coordinates, strand and phase; demonstrates a clear mechanistic gain from the proposed model; transfers without target-label tuning; and is compared with Helixer, Tiberius and ANNEVO on identical genome ranges and metrics.

The present evidence supports a method-development story, not yet a Nature Communications result. The paper must eventually establish all four layers:

1. exact and valid gene structures;
2. fair advantage or a clearly useful trade-off against released callers;
3. genuine transfer with frozen inference decisions and audited provenance;
4. biological and technical robustness beyond one species and one favorable metric.

## Stage 0 — preserve the current evidence boundary

Current state: complete.

- Keep M25R checkpoints, raw validation grid and runtime outputs on Baobab only.
- Version only code, frozen configs, compact summaries and small logs.
- Keep Setaria annotation inaccessible while the model/checkpoint/decoder/predictions are not frozen.
- Keep the analysis restricted to primary nuclear chromosomes.
- Do not submit more M25-like training, extra seeds, SegmentNT training or baseline reruns.

## Stage 1 — locate the M25R attrition point without retraining

Experiment/analysis ID: `M25R-DEV-REDECODE-ERROR-DECOMPOSITION`.

Use only the three existing M25R checkpoints, the original Arabidopsis/rice development primary nuclear chromosomes, the original validation grid and the complete set of training examples actually used by M25R. That training set is deterministically defined by the frozen config, seed, split, allowlist and sampling order; its observed size is `1,536`, but Stage 1 must not create a new sampled subset. No Setaria file, including FASTA, is needed or allowed.

For each epoch, replay exactly one tuple: the highest-ranked row for that epoch under the original validation-grid rank and tie-break. Do not re-rank epochs, search thresholds or create a new selected tuple. The global epoch-1 row must reproduce interval F1 `0.1204141`, chain F1 `0.3249883`, FPR `0.0124683` and count ratio `0.3252713`.

Define two reference views without changing the core denominator:

- `R_all`: the exact `6,450` pooled development reference chains used by the frozen M25R evaluator;
- `R_canonical`: the tagged subset representable by the frozen canonical start/stop and GT-AG decoder.

All accounting uses `R_all`. Non-canonical chains remain visible as `unsupported_by_frozen_canonical_decoder`; they must not disappear from the denominator.

The diagnostic must:

1. reproduce interval F1, chain F1, FPR and gene-count ratio for all three replayed tuples within absolute tolerance `1e-5`;
2. report train-versus-validation region metrics, start/stop/donor/acceptor AUCPR, exact/`±1`/`±3`/`±6` event recall, CDS-base phase accuracy and phase accuracy at decoded CDS starts;
3. compute three fixed post-hoc upper bounds without altering any parameter: transition reachability within the decoder's frozen `±6 bp` radius, motif-candidate reachability under the frozen motif dictionary, and truth-assisted exact-chain recovery when an exact truth coordinate already exists in the frozen candidate set;
4. trace attrition in the production order: region state path, non-intergenic block, ordered CDS runs, legal terminal transitions, start/stop motif candidates, ordered donor/acceptor candidates, learned boundary choice, phase check, complete-ORF/internal-stop check, frozen boundary-threshold filter and emitted GFF3;
5. use deterministic same-chromosome/same-strand one-to-one matching so every reference chain enters exactly one earliest-failure category and every candidate lineage is reconciled with the final GFF3 count;
6. report boundary offsets, strand and phase errors, split/fusion, missing and false-positive models, single-exon versus multi-exon strata, CDS-count strata, genes spanning more than `6,144 bp`, and events within `6 bp` of a production tile edge;
7. independently check every emitted transcript for parent linkage, coordinates, strand, CDS order, phase continuity, start/stop completeness, in-frame ORF, internal stops and the declared splice grammar;
8. stop for scientific review.

Development truth is permitted only for aggregate reproduction, diagnostics and upper bounds. It must not select a checkpoint, tuple, threshold, radius, motif dictionary, grammar rule or decoder. Oracle results are not model performance.

Stage success means diagnostic integrity, not model success. It requires exact reference and prediction accounting, `100%` independent-validity audit coverage, unchanged weights/thresholds/decoder, and zero Setaria reads. If no transcript is emitted, report `n_emitted=0`, `n_checked=0`, `complete_empty_audit=true` and model validity fraction `not_applicable`; count and recall then carry the model failure. Any aggregate mismatch, unassigned or double-counted reference, unexplained candidate count or production-trace mismatch blocks interpretation.

Stage 1 can localize a learned-emission/candidate, structural-head, data-generalization or frozen-decoder failure. It cannot by itself prove that the GENERanno backbone representation is intrinsically insufficient, because backbone, LoRA, heads, targets, loss and training budget remain coupled.

### Threshold taxonomy

| Rule | Status | Decision use |
|---|---|---|
| M25R validation FPR `<=0.020`, count `0.80–1.20` and the original validity gate | frozen historical M25R contract | explains why M25R stopped; not automatically inherited by a redesigned model |
| Aggregate reproduction error `<=1e-5`; all `6,450` references and all candidates/predictions reconciled; validity audit coverage `100%`; no Setaria reads; no model or decoder change | hard Stage 1 integrity requirements | violation makes the diagnostic uninterpretable |
| One stage explains `>=50%` of missing chains | operational attribution heuristic | permits the word `dominant`; below it, report a mixed failure |
| Removing one stage has truth-assisted headroom `>=0.15` absolute | operational information-value heuristic | may justify proposing one discriminating experiment; it does not authorize it |
| Development interval `>=0.60` and chain `>=0.50` | withdrawn first-round suggestion | not frozen; the next experiment's structural gate is set only after Stage 1 and same-range development baseline scoring |
| Strongest baseline `+2` percentage points | publication-planning heuristic | useful target, not a journal requirement; a defensible Pareto advantage may substitute |

## Stage 2 — choose exactly one bottleneck branch

Stage 1 ends before any branch is authorized. After independent review, re-score Helixer, Tiberius and ANNEVO on the exact same Arabidopsis/rice development chromosomes and evaluator. Use this paired comparison to freeze the exact-structure gate for one later branch experiment; never use baseline results to select a checkpoint, decoder or threshold.

| Observed evidence | Interpretation | Minimal next branch | Stop condition |
|---|---|---|---|
| Most chains disappear before usable region, CDS-run or transition candidates exist | learned-emission/candidate-generation bottleneck, not proven backbone failure | Propose one released SegmentNT native-long-context, fixed-context, development-only emission probe. Do not train or emit full GFF3. | Stop SegmentNT promotion if donor/acceptor and ordered splice-chain reachability do not improve consistently in both development species. |
| Transition/motif reachability is high, but boundary/phase heads are weak on both train and validation | structural-head/supervision-fit bottleneck | Freeze backbone, LoRA and region head; propose one boundary/phase-head-only experiment. | Stop if the heads still cannot fit training targets or raw gains do not become exact structures. |
| Train raw heads are strong but validation is weak | development data/distribution bottleneck | Propose one species-balanced fit with every other contract fixed. | Stop if only training metrics improve or one development species improves by collapsing the other. |
| Raw candidates, phase and truth-assisted chain ceiling are high, but final emission is low | decoder/grammar bottleneck | Re-decode existing scores after changing only the single rule responsible for the largest attrition. | Stop if multiple rules are needed or the one change cannot realize its registered headroom without FPR/validity damage. |
| No stage is dominant or no single stage has useful headroom | mixed bottleneck | Do not patch M25 incrementally; redesign the task representation and review again. | No automatic experiment. |

The SegmentNT probe, if later approved, uses the same development chromosomes, transcript policy, event coordinates, one-to-one matching and upper-bound evaluator as GENERanno. It fixes one released native long-context setting, does not sweep context, does not fine-tune and reports only exon/intron and donor/acceptor emissions plus splice-event/chain reachability. It does not claim phase, start/stop, full CDS chains or complete gene models.

The sequence is fixed: `Stage 1 -> review -> A/rice same-range baseline re-score -> freeze one branch gate -> consider one branch experiment`. Nothing after Stage 1 is automatic.

## Stage 3 — pre-register one development experiment

Before submission, freeze:

- biological question and single causal change;
- training and validation species/chromosomes;
- one fit/seed for discovery;
- model, tokenizer, window/context policy and supervision;
- checkpoint selection and decoder search space;
- exact interval, chain, coding-gene, strand, phase, validity, FPR and count metrics;
- unchanged-input or component ablation that tests the claimed mechanism;
- numeric success and stop conditions;
- all files needed to reproduce the decision.

Before a later branch is submitted, separately freeze (a) absolute usability guardrails for FPR, gene count and independent validity and (b) exact interval/chain/gene competitiveness targets informed by the paired development baselines. M25R's FPR `<=0.020`, count `0.80–1.20` and validity `>=0.99` are reasonable starting proposals, not automatically inherited law. A passing FPR with low complete-gene recovery is never success.

## Stage 4 — freeze the candidate before opening the blind test

Only after a development candidate passes its registered gate:

1. freeze the candidate resolved config, model checkpoint, exact decoder parameters, chromosome allowlist and ablation;
2. freeze Helixer, Tiberius and ANNEVO versions, weights and commands;
3. generate candidate full/ablation and all baseline Setaria predictions from FASTA without reading Setaria annotation;
4. record the Git commit, resolved configs, exact output paths and the absence of target-label choices;
5. release the Setaria annotation once and score every frozen prediction with one evaluator.

For historical clarity, the frozen M25/M25R Setaria discovery target was:

- strand-aware exact CDS interval F1 `>=0.80`;
- exact CDS-chain/transcript F1 `>=0.55`;
- exact coding-gene F1 `>=0.50`;
- matched-gene strand accuracy `>=0.98`;
- exact-matched CDS phase accuracy `>=0.90`;
- structurally valid complete-transcript fraction `>=0.99`;
- intergenic FPR `<=0.020`;
- predicted-gene-count ratio `0.80–1.20`;
- full minus ablation interval and chain F1 each `>=+0.10`;
- full minus ablation FPR `<=+0.005`.

Failure ends that frozen branch. It does not authorize tuning on Setaria.

A redesigned future branch must register its own blind gate before training, without consulting Setaria annotation. It may retain the historical gate, but any change must be justified from development evidence and same-range baselines rather than target labels.

## Stage 5 — score released callers under the identical blind contract

Use the Helixer, Tiberius and ANNEVO predictions generated before annotation release on the same frozen Setaria primary chromosomes. Score them with the same reference transcript policy and evaluator. Report:

- exact strand-aware CDS interval precision/recall/F1;
- exact CDS-chain/transcript and coding-gene precision/recall/F1;
- strand and phase accuracy where represented;
- independent structural validity;
- intergenic FPR and predicted-gene-count ratio;
- runtime, memory and hardware.

Do not compare full-genome baseline numbers with subset candidate numbers. Do not tune released callers on Setaria annotation. Tool-specific unsupported metrics remain `not_applicable`, not zero.

## Stage 6 — establish reproducibility and transfer

If the one-fit candidate passes the blind gate:

- run two additional training seeds, for three total;
- plan a publication panel targeting at least six held-out nuclear-genome species in total, with accessions and scope frozen before their annotations are used for model selection;
- keep one global decoder unless clade specialization is itself the registered method;
- report per-species and macro metrics, chromosome-level uncertainty, gene-count behavior and failure strata;
- distinguish same-clade transfer from broader eukaryotic transfer.

If a global model fails but evidence supports clade specialization, an adapter/MoE route is publishable only if gating is determined from sequence or declared taxonomy rather than target labels, each expert has sufficient training support, and the comparison includes a parameter/data-matched global model.

## Stage 7 — prove mechanism rather than checkpoint luck

Required ablations depend on the selected branch, but should minimally separate:

- CDS-specialized GENERanno 1.2B from an appropriate generic-pretraining control without claiming that the current 1.2B-versus-0.5B comparison isolates post-training alone;
- structural heads/supervision from unchanged coarse-region input;
- learned structural emissions from decoder-only motif repair;
- full decoder from the implicated decoder component;
- context length if long-context extraction is claimed;
- optional SegmentNT backbone from the same structural decoder if SegmentNT becomes mainline.

The GENERanno generic-pretraining and CDS-post-training species/accession audit must be completed or the paper must explicitly use adaptation, not clean zero-shot language.

## Stage 8 — robustness and biological utility

After the core method passes, evaluate failure modes that are reachable in the supported nuclear-genome use case:

- single-exon versus multi-exon genes;
- gene length and CDS count, including genes longer than one model window;
- canonical versus supported non-canonical splice sites;
- low- versus high-gene-density regions;
- repeat-rich and low-complexity nuclear regions if they occur in the frozen chromosomes;
- gene split, merge, missing-model and spurious-model rates;
- agreement with independent RNA-seq splice junctions and protein-homology evidence on at least one held-out genome.

Organelles and unplaced scaffolds remain outside this route unless a later paper scope explicitly adds them.

## Stage 9 — Nature Communications decision package

A credible submission package should contain:

1. a clearly defined direct gene-annotation method and mechanistic diagram;
2. a frozen multi-species benchmark with exact structural metrics and released callers on identical ranges;
3. one untouched blind-test result plus broader held-out transfer;
4. three-seed robustness for trained candidates;
5. causal ablations showing why the method works;
6. provenance and leakage audit;
7. biological evidence that exact predictions recover real genes/splice junctions;
8. practical runtime/resource comparison and reproducible code/config release;
9. transparent negative results and failure boundaries.

The manuscript is **no-go** if the best model only improves coarse gene-body F1, if exact structures remain far below released callers, if gains require target-label tuning, or if the provenance prevents the intended transfer claim. No fixed `+2 pp` margin is a journal rule: the central result must instead show a cross-species paired advantage with uncertainty support and no unacceptable FPR/count cost, or an equally strong pre-registered Pareto value such as lower false positives, better long-gene recovery, runtime or independently supported reference correction. A narrower adaptation/methods paper may still be viable, but it would be a different claim and possibly a different journal target.

## Review disposition

ChatGPT Pro verified access to private commit `fe38a675f2d4a0e9349046c226a5679b850a82b0` and reviewed the frozen evidence in two rounds. The first round reconstructed M9-M25R and returned a Stage 1 go/no-go; the second round challenged unsupported thresholds, label use, upper-bound definitions, baseline timing and implementation scope. The resulting joint decision is recorded in `docs/30_chatgpt_pro_systematic_review_2026-09-01.md`.

The review does not authorize Stage 1 execution, a baseline run, SegmentNT, training or external provenance contact. It freezes the questions and decision rules so a later implementation can remain minimal and auditable.
