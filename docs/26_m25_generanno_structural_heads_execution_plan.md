# M25 GENERanno structural-heads execution plan

Status: `approved_by_chatgpt_pro`
Experiment: `M25-GENERANNO-1P2B-STRUCTURAL-HEADS-s0`
Frozen branch: `GENERANNO_STRUCTURAL_HEADS`
Product: raw genome FASTA to one primary protein-coding gene model; no candidate GFF at inference

ChatGPT Pro reviewed the corrected M24 result and returned `AUDIT_VERDICT: VALID_FOR_BRANCH_DECISION`. It selected this branch over SegmentNT structural emission and `STOP_NO_GPU`. After two required-change rounds, Pro returned `VERDICT: APPROVED` for this implementation document.

## 1. Evidence and mechanism under test

M19 has useful coarse CDS/gene-body emission but not usable exact structure. Across Arabidopsis and rice, exact CDS-run F1 is `0.0531–0.1498` and pseudo-chain F1 is `0.0082–0.0123`. Exact boundary recall is low, while ±6 bp boundary recall is `0.7614–0.8617`. The single hypothesis tested by M25 is that nucleotide identity plus explicit boundary and phase supervision can convert these near-boundary emissions into exact coding structures.

The current 6-kb SegmentNT cache has moderate exon signal but weak donor/acceptor signal: primary-view exon AUCPR is `0.5866–0.6569`, while donor/acceptor AUCPR is `0.0314–0.0443`. It is not the first GPU branch. This does not reject SegmentNT under longer context or direct checkpoint adaptation.

The direct-route blockers are frozen in this order:

1. single-nucleotide coordinate decoding;
2. missing strand, phase, start, stop, donor and acceptor semantics;
3. CDS-chain grouping;
4. coarse region emission.

## 2. Data contract and Setaria embargo

- Train: existing Arabidopsis and rice train seqids.
- Validation, checkpoint selection and decode calibration: existing Arabidopsis and rice validation seqids, pooled.
- Test: `Setaria viridis`, assembly `GCF_005286985.2` (`Setaria_viridis_v4.0`), annotation `RS_2025_03`.
- Seed: `0`; one fit only.

Scope amendment approved after the first GPU smoke exposed organellar trans-splicing and circular-origin records that the frozen decoder cannot represent. Per the user decision, M25 trains, validates, infers and evaluates only explicit primary nuclear chromosome seqids. Chloroplast, mitochondrion and unplaced scaffolds are excluded as whole FASTA records; individual GFF lines are never silently skipped. The original full-sequence list still determines train/validation splits before the chromosome allowlist is applied. ChatGPT Pro returned `AMENDMENT_VERDICT: APPROVED`, `REQUIRED_CHANGES: NONE`, `FORMAL_FIT_ALLOWED: YES`.

- Arabidopsis: `NC_003070.9`, `NC_003071.7`, `NC_003074.8`, `NC_003075.7`, `NC_003076.8`.
- Rice: `NC_089041.1`, `NC_089042.1`, `NC_089043.1`, `NC_089044.1`, `NC_089046.1`.
- Setaria: `NC_048263.2` through `NC_048271.2`.

The exact allowlist is written and hashed in the embargo-release marker before the Setaria annotation is acquired. Post-release reference parsing is restricted to the identical seqids.

Only the Setaria assembly FASTA may be acquired before inference. The Setaria annotation must not be downloaded or parsed until all of the following are immutable and recorded: resolved training config, checkpoint hash, validation-selected decode parameters, success gate, Setaria FASTA hash and both full/ablation prediction GFF3 hashes. The annotation is then acquired once for evaluation. Setaria labels, coding fraction, gene count and external-caller results cannot change the checkpoint, thresholds, decoder or gate.

The embargo is enforced by the execution path:

- before training and inference, sbatch fails if any Setaria GFF, GTF, GBFF or reference-annotation path/file exists in the project data or is passed to a process;
- the trainer receives only the Setaria FASTA path;
- after config, checkpoint, validation decode and both prediction hashes are written, sbatch creates `SETARIA_EMBARGO_RELEASED.json` containing those hashes;
- annotation download/parsing is forbidden before this marker;
- `eval_m25_structure.py` fails unless the marker exists and all recorded hashes still match.

Arabidopsis and rice are development species after M24. They cannot be described as M25 unseen tests.

## 3. Frozen model and training change

Inherit unchanged from M19:

- `GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview`;
- official 6-mer tokenizer and 6,144-bp windows;
- LoRA `r=8`, `alpha=16`, targets `q_proj,k_proj,v_proj,o_proj`, dropout `0.05`;
- head LR `8e-4`, LoRA LR `2e-5`, batch size `1` and exactly three completed epochs with no early stopping;
- train-window cap `1536`, validation-window cap `768`, bf16, seed `0`;
- existing 3-class region head and FP-aware region loss with `fp_lambda=2.5`.

Make one structural change after the per-base upsampled 6-mer representation:

- concatenate A/C/G/T/N nucleotide identity so the head can distinguish the six positions represented by one token;
- add four boundary logits: start, stop, donor and acceptor;
- add four phase logits: none, phase 0, phase 1 and phase 2;
- infer forward sequence and reverse complement separately, then map the reverse-complement predictions back to genomic coordinates; no placeholder strand field remains.

Training is explicitly strand-conditioned. For every sampled window, the forward example receives region, boundary and phase targets only from `+` primary transcripts. Its reverse-complement example receives transformed targets only from `-` primary transcripts. Opposite-strand structures are absent from that orientation's region and structural targets; no undifferentiated two-strand gene-body mask is reused.

Use longest-CDS primary transcripts for structural labels. Mask partial or ambiguous transcripts from boundary/phase loss. Keep the existing region loss and add, without a sweep:

`loss = M19_region_FP + 2.0 * boundary_focal(gamma=2) + 0.5 * phase_CE`.

Boundary focal loss is computed independently for start, stop, donor and acceptor. Each type uses fixed train-set class balancing: positive weight `N_negative / (N_positive + N_negative)` and negative weight `N_positive / (N_positive + N_negative)`, computed once before training; the four type losses are averaged. Phase CE is evaluated only at CDS-labelled bases with fixed inverse-square-root balancing among phases 0, 1 and 2; the `none` class does not dominate the reduction. The coefficients above do not change.

The only unchanged-input ablation uses the same trained checkpoint, Setaria FASTA, forward/reverse-complement inference, region logits, region transitions and canonical motif grammar. It ignores learned boundary/phase scores, bypasses their learned threshold checks, selects the nearest valid motif with the deterministic distance/tie rule below and derives phase only from decoded CDS lengths. It is not a second fit and cannot fail merely because constant logits fall below learned thresholds.

## 4. Minimal decoder

Select one global parameter set on pooled Arabidopsis/rice validation only. Training always completes exactly three epochs and saves the end of epochs 1, 2 and 3; early stopping cannot remove a candidate. Candidate thresholds are evaluated in the listed config order: region genic probability `[0.40, 0.50, 0.60]`; independently for each start, stop, donor and acceptor sigmoid probability `[0.10, 0.20, 0.30, 0.40, 0.50]`. Boundary snapping radius is fixed at ±6 bp.

Choose the checkpoint/threshold tuple that maximizes pooled validation exact CDS-chain F1 subject to validation intergenic FPR `<=0.020`, gene-count ratio `0.80–1.20` and structurally valid complete transcript fraction `>=0.99`. Tie-break by higher exact CDS-interval F1, then lower FPR, then the fixed epoch/threshold enumeration order. If no tuple satisfies all constraints, write `STOP_M25_BRANCH` before any Setaria inference.

For each strand independently:

- softmax the three region logits; set state `I` when `1 - P(intergenic)` is below the selected genic threshold, otherwise set `C` when `P(CDS) >= P(gene-body non-CDS)` and `G` otherwise;
- use these exact transition candidates in orientation coordinates: `I→C` permits start, `G→C` permits start or acceptor, `C→G` permits stop or donor, and `C→I` permits stop; `I↔G` alone emits no coding boundary;
- a single-exon gene is the fixed `I→C ... C→I` case; UTR-bearing boundaries use `I→G→C` for start and `C→G→I` for stop, while internal `C↔G` transitions use donor/acceptor under the phase and alternation constraints;
- choose the highest-scoring motif-valid boundary within ±6 bp of a region transition; equal scores tie by smallest absolute offset and then smallest genomic coordinate;
- require donor/acceptor alternation and continuous CDS phase;
- require a valid start and an in-frame stop with no internal in-frame stop;
- emit only uniquely decoded complete canonical primary coding models;
- do not tune or repair rules after Setaria inference.

The first run supports canonical GT–AG splicing. Noncanonical reference genes remain in the overall evaluation and are also reported as a fixed stratum. No HMM, CRF, MoE, long-context change, candidate GFF, UTR model or alternative-isoform model is added.

## 5. Minimum implementation surface

Implement only the M25-specific training, decode, evaluation and submission path. Reuse current FASTA/GFF parsing, M19 backbone/LoRA loading, M24 structural metrics and Baobab environment where their contracts match. Do not introduce a general gene-annotation framework.

Planned files:

- `configs/M25-GENERANNO-1P2B-STRUCTURAL-HEADS-s0.yaml`: frozen model, data, loss, decode, embargo and gate values;
- `src/foundation_probe/train_generanno_structural_heads.py`: M19-derived single experiment trainer/inference path with nucleotide, boundary and phase outputs;
- `scripts/eval_m25_structure.py`: thin use of the M24 parser/metrics plus validity, signed-boundary-offset and modulo-6 diagnostics;
- `tests/test_generanno_structural_heads.py`: primary-label, reverse-complement coordinate, phase-continuity, decode and unchanged-input-ablation fixtures;
- `sbatch/M25-GENERANNO-1P2B-STRUCTURAL-HEADS-s0.sbatch`: one dynamic-routed GPU job that enforces and releases the Setaria annotation embargo.

Do not modify the historical M19 trainer or its outputs. Do not add retries, fallback parsers, exception swallowing, metric-triggered resubmission or generic abstractions.

## 6. Required outputs

Under `outputs/M25-GENERANNO-1P2B-STRUCTURAL-HEADS-s0/` write only:

- `config_resolved.yaml`, the frozen threshold grid/selection rule and data/checkpoint/prediction hashes;
- `train_summary.json` and validation-selected decode parameters;
- the selected checkpoint;
- Setaria raw structural scores required for the fixed ablation;
- full and ablation Setaria GFF3 predictions and `SETARIA_EMBARGO_RELEASED.json`;
- one structure-metrics JSON and Slurm logs;
- terminal `COMPLETED`, `STOP_M25_BRANCH` or `PASSED_DISCOVERY_GATE` status.

GFF3 output is limited to `gene`, `mRNA`, `CDS`, `start_codon` and `stop_codon`. M25 makes no UTR, alternative-isoform, noncoding-RNA, partial-gene or reference-correction claim.

## 7. Frozen success and stop gates

All success conditions must hold on the once-opened Setaria evaluation:

| Metric | Success threshold |
|---|---:|
| strand-aware exact CDS interval F1 | `>= 0.80` |
| exact CDS-chain / coding-transcript F1 | `>= 0.55` |
| exact coding-gene F1 | `>= 0.50` |
| matched-gene strand accuracy | `>= 0.98` |
| exact-matched CDS phase accuracy | `>= 0.90` |
| structurally valid complete transcript fraction | `>= 0.99` |
| intergenic FPR | `<= 0.020` |
| predicted gene-count ratio | `0.80–1.20` |
| full minus ablation exact interval F1 | `>= +0.10` |
| full minus ablation exact chain F1 | `>= +0.10` |
| full minus ablation FPR | `<= +0.005` |

The interval and chain floors sit just below the M24 same-scope external minima (`0.8117` and `0.5850`). The FPR floor is near the external median and excludes the high-FP range; the gene-count range contains the external `0.893–1.114` range with bounded transfer margin.

Any failed success condition ends this branch after reporting. In addition, mark `STOP_M25_BRANCH` immediately if annotation embargo is violated, placeholder strand/phase remains, structurally valid fraction is below `0.95`, exact interval F1 is below `0.60`, exact chain F1 is below `0.30`, exact gene F1 is below `0.25`, FPR exceeds `0.030`, gene-count ratio falls outside `0.50–1.50`, or either exact-structure ablation gain is below `0.10`.

Do not change loss weights, thresholds or decoder rules; do not submit another fit. Only if every success condition passes may seeds 1 and 2 be proposed in a new reviewed document.

## 8. Uncertainties that remain explicit

- GENERanno training-species/accession provenance remains `overlap_unknown`; a successful Setaria result is not clean pretraining-held-out evidence.
- M24 cannot distinguish whether ±6 errors originate from 6-mer upsampling, 3-class supervision, constrained postprocessing or motif ambiguity.
- M25 is a coding-structure MVP, not full transcript annotation.
- Canonical-only decoding may suppress real noncanonical genes; they remain in total evaluation and a fixed stratum.
- A 6,144-bp window may fail on long genes or window edges. Report fixed strata; do not enlarge context in this experiment.
- Setaria external-baseline performance is unknown and cannot change the frozen gate.
- `RS_2025_03` is an evaluator reference, not absolute biological truth.
- Peak memory and runtime of forward plus reverse-complement prediction are not yet measured. Use one bounded smoke on development-only windows only if required to verify the new code path, then submit the one real fit through live Baobab routing.

Tests must freeze 0-based internal and 1-based GFF3 coordinate conventions for start, stop, donor, acceptor and CDS phase on both strands. They must also prove that reverse-complement mapping returns the original genomic coordinates without duplicate strand calls, and that the fixed ablation emits motif-derived models without consulting learned boundary/phase thresholds.

## 9. Execution record

- GPU smoke `12078092` stopped before model execution because strict parsing exposed 10 organellar trans-splicing/circular-origin records; no record was skipped.
- GPU smoke `12078382` completed but was rejected because its one-window class statistics produced a non-finite loss.
- GPU smoke `12078563` passed on a representative primary-chromosome development window: exit `0:0`, 61 seconds, finite loss `3.440891981124878`, and no Setaria inference.
- Formal job `12094731` is frozen for `private-teodoro-gpu`, one 24-GB RTX 3090, 8 CPUs, 96 GB RAM and a six-day limit. It is held until `2026-08-28T10:05:00`, after the scheduled Baobab maintenance window.

Runtime remains uncertain. M19 required about 17 hours for two epochs plus a smaller inference scope; M25 adds a third epoch, bidirectional full-chromosome validation after every epoch and bidirectional inference over nine Setaria chromosomes. Six days is an operational upper bound, not a measured runtime claim.
