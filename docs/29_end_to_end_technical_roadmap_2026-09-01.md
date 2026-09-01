# End-to-end technical roadmap — draft for ChatGPT Pro review

Status: planning only. No experiment is authorized by this document.

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

Use the three existing M25R checkpoints and Arabidopsis/rice development data only. For each checkpoint and selected representative decoder tuples:

1. reproduce the frozen interval F1, chain F1, FPR and gene-count aggregates within absolute tolerance `1e-5`;
2. record the number of reference coding genes surviving each stage: usable region emission, candidate region, oriented boundary candidates, legal chain, complete start/stop model, unique emitted model;
3. assign each of the `6,450` pooled reference chains exactly once to its earliest failure stage;
4. maintain a corresponding ledger for all emitted predictions and reconcile every stage count with the final GFF3 count;
5. measure raw score recall around true start, stop, donor and acceptor sites before decoder filtering;
6. separate strand mistakes, phase mistakes, boundary-offset errors, chain split/merge errors, missing models and false-positive models;
7. validate every emitted transcript independently for strand-consistent order, legal phase continuity, start/stop completeness, non-overlapping CDS ordering and expected splice grammar;
8. stop for scientific review.

Stage success means the diagnosis is internally reconciled, not that the model passes. Any aggregate mismatch, unassigned reference, double assignment or unexplained prediction-count difference blocks interpretation and must be fixed before proceeding.

## Stage 2 — choose exactly one bottleneck branch

The Stage 1 evidence determines the next branch. Do not blend branches in one experiment.

| Observed evidence | Interpretation | Minimal next branch | Stop condition |
|---|---|---|---|
| Low raw recall at true CDS/boundary sites before decoding | representation/emission bottleneck | Test one stronger structural representation: either longer-context/directly adapted SegmentNT or a revised GENERanno structural emission head, selected by Pro after reviewing the score upper bounds. | Stop if raw-event recall and exact-interval upper bound do not materially exceed M25R. |
| Adequate raw event recall, but large loss at structural-head scoring | structural-supervision bottleneck | Change only the supervision/target formulation implicated by the diagnostic; retain backbone and decoder. | Stop if boundary/phase validation improves without recovering complete-gene count. |
| Adequate emissions and head scores, but large loss in snapping/transitions/completeness/uniqueness | decoder bottleneck | Re-decode existing scores or make one minimal grammar/chain reconstruction correction; no backbone retraining first. | Stop if gene-count recovery requires FPR above the frozen limit or exact chains do not improve. |
| Different failure profiles across Arabidopsis and rice after all internal stages are sound | data/domain bottleneck | Freeze a broader development panel or a clade-aware adaptation design; do not use Setaria labels. | Stop universal-model framing if gains require target-specific labels or thresholds. |
| Multiple comparable losses with no dominant stage | mixed bottleneck | Choose the single stage with the largest recoverable exact-chain upper bound and test it first. | Stop if the intervention does not move the predicted bottleneck. |

SegmentNT is therefore conditional, not automatically the next mainline. The present 6-kb cache shows moderate exon/intron signal but weak rare-boundary signal. It becomes the next backbone only if Stage 1 shows that GENERanno emissions are the limiting factor and a longer-context SegmentNT extraction has a plausible structural upper bound. If Stage 1 shows decoder attrition despite adequate emissions, changing the backbone would not answer the causal question.

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

At minimum, development admission must retain M25R's non-negotiable constraints: intergenic FPR `<=0.020`, predicted-gene-count ratio `0.80–1.20`, and independently measured structurally valid complete transcripts `>=0.99`. Pro must approve the exact structural improvement threshold before execution. A passing FPR with low complete-gene recovery is not success.

## Stage 4 — freeze the candidate before opening the blind test

Only after a development candidate passes its registered gate:

1. freeze the resolved config, model checkpoint, exact decoder parameters and chromosome allowlist;
2. generate both full and unchanged-input/component-ablation Setaria predictions without reading Setaria annotation;
3. freeze the complete prediction GFF3 files;
4. record that no target-label choice occurred;
5. release the Setaria annotation once for evaluation.

The registered discovery target remains:

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

## Stage 5 — run released callers under the identical blind contract

Run Helixer, Tiberius and ANNEVO on the same frozen Setaria primary chromosomes and score them with the same reference transcript policy and evaluator. Report:

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
- freeze at least two additional held-out nuclear-genome species chosen before viewing their annotations for model selection;
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

The manuscript is **no-go** if the best model only improves coarse gene-body F1, if exact structures remain far below released callers, if gains require target-label tuning, or if the provenance prevents the intended transfer claim. A narrower adaptation/methods paper may still be viable, but it would be a different claim and possibly a different journal target.

## Immediate decision requested from ChatGPT Pro

Pro should review the repository evidence and return:

- a go/no-go on Stage 1 as the sole immediate action;
- whether the proposed error ledger can distinguish representation, supervision and decoder failures;
- the one next branch for each possible Stage 1 outcome;
- exact success and stop criteria for the next development experiment;
- whether the Setaria blind target and baseline protocol are sufficient;
- what minimum extra evidence is required for a Nature Communications-level story;
- which stages can be removed or simplified without weakening the claim.
