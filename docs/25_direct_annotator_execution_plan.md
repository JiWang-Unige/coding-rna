# Direct annotator execution plan

Status: `approved_after_required_changes`
Decision owner: Ji Wang
Frozen product choice: **A — raw genome FASTA to direct gene annotation**

ChatGPT Pro verdict: `APPROVED_WITH_REQUIRED_CHANGES` on 2026-08-24. The five required corrections below are incorporated; no optional expansion was added.

No research code or Slurm job may be changed or submitted from this plan until ChatGPT Pro has reviewed it. Review approval applies first to the 48-hour diagnostic. The first GPU experiment remains gated by that diagnostic result.

## 1. Product and claim contract

The paper method must infer a primary protein-coding gene annotation from raw genome sequence. Existing annotations and external caller GFF files may be used as references or baselines during training/evaluation, but they are not inference inputs.

The current strongest internal result, M19, is a downstream adaptation of `GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview`: the official token-classification head is discarded, LoRA is applied to the released backbone, and the project trains its own 3-class ConvLSTM head with an FP-aware objective and validation-only calibration. M19 is not a complete gene annotator.

The official SegmentNT checkpoint and the generic NT-v2 backbone are separate routes:

- SegmentNT supplies 14 supervised genomic-element probabilities, including exon, intron, splice donor and splice acceptor. Its plant evaluations are cleaner than its fly/chicken evaluations because plants were absent from both disclosed NT-v2 pretraining and SegmentNT supervised species.
- M9–M12 trained the underlying NT-v2-500M backbone directly. Those runs are not evaluations of the unchanged official SegmentNT segmentation checkpoint.
- GENERanno 1.2B is the strongest span-level development backbone, but its pretraining and CDS-post-training species/accession overlap is unknown. `overlap_unknown` must not be rewritten as confirmed overlap or confirmed exclusion.

## 2. Uncertainties that remain live

1. M19 predictions contain grouped CDS runs but force `strand=+` and `phase=0`. It is unknown how much exact CDS-boundary and chain signal remains before those unsupported fields are considered.
2. M8 multiclass predictions were collapsed to the existing 3-class representation before GFF and raw-prediction persistence. They cannot be used as saved strand/phase predictions.
3. The cached SegmentNT 14-element features do retain exon/intron/donor/acceptor probabilities, but their plant test-set boundary accuracy has not been measured.
4. Existing external baseline predictions cover full genomes, while M19 evaluation uses held-out seqids. Their published-looking aggregate rows are not yet an exact same-scope structural comparison.
5. No threshold for promoting a direct structural backbone is justified before the first structural diagnostics exist. The diagnostic will report effect sizes; it will not manufacture a pass threshold from the current span metric.

## 3. M24: 48-hour CPU diagnostic

Experiment ID: `M24-DIRECT-STRUCTURE-DIAGNOSTIC`

### Question

Before changing the model task, what structural information is already present in M19 CDS runs and the official SegmentNT feature cache, and how far is it from the exact structures emitted by released callers on the identical held-out seqids?

### Existing inputs

Reference/scope source:

- `outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s1/eval_subsets/<species>/genome.fa`
- `outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s1/eval_subsets/<species>/reference.gff3`

Candidate predictions:

- M19 seeds 0 and 1 under `outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s*/predictions/`
- clean-plant SegmentNT-derived M8 3-class predictions from the three completed artifacts: seeds 0, 2 and 4 under `outputs/M8-3C-CAND-s{0,2,4}/predictions/`
- `outputs/FP-SEGMENTNT-FEATCACHE/segment_nt_multi_species/{arabidopsis_thaliana,oryza_sativa}.npz`

Released baselines:

- `outputs/M12B-SAMEPANEL-BASELINES-ANNEVO/predictions/`
- `outputs/M12B-SAMEPANEL-BASELINES-HELIXER/predictions/`
- `outputs/M12B-SAMEPANEL-BASELINES-TIBERIUS/predictions/`

All rows are filtered in memory to the seqids present in the M19 evaluation FASTA. No filtered GFF copies or new data manifests are created.

### Frozen normalization and scope

- Evaluate protein-coding genes only.
- Choose one primary transcript per reference/predicted gene by longest total CDS length; break equal-CDS-length ties by deterministic transcript ID ordering.
- Exclude partial genes from exact complete-structure metrics and report their counts separately.
- Derive introns from adjacent CDS/exon features only when gene–transcript–feature parent mapping is unambiguous.
- For every method row, assert identical FASTA seqids and lengths, reference primary-transcript count, reference CDS-interval count and genome bases. Also report retained prediction records and require zero out-of-bounds coordinates.

### Minimal implementation after review approval

Add only:

- `scripts/eval_structure_diagnostic.py`
- `tests/test_eval_structure_diagnostic.py`
- `sbatch/M24-DIRECT-STRUCTURE-DIAGNOSTIC.sbatch`

The evaluator will reuse the repository's existing GFF/GTF attribute parsing where it is correct. It will fail with the original exception on malformed required input rather than silently skipping records.

The test file will cover exact CDS interval/chain matching, supported versus unsupported structure fields, and one minimal fixture for each real baseline syntax: Tiberius GTF, ANNEVO GFF and Helixer GFF3. The fixtures verify gene–transcript–feature grouping and the frozen primary-transcript rule. It will not create a general GFF framework.

### Required outputs

Write runtime artifacts under `outputs/M24-DIRECT-STRUCTURE-DIAGNOSTIC/` and one tracked result report at `reports/M24-DIRECT-STRUCTURE-DIAGNOSTIC/report.md`.

For M19 and M8 3-class outputs, report:

- strand-agnostic `exact_CDS_run` interval precision, recall and F1;
- genomic left/right CDS-run boundary precision and recall at exact, ±1 bp, ±3 bp and ±6 bp; these are coordinate boundaries, not biological CDS start/stop codons;
- strand-agnostic coordinate pseudo-chain precision, recall and F1 using the prediction's existing `gene_id` grouping;
- split/merge overlap counts between predicted groups and reference primary coding transcripts;
- current same-scope CDS-span gbF1, intergenic FPR and gene-count ratio.

For M19/M8, strand-aware exact interval/chain, phase accuracy, exact intron-chain, exact transcript and exact gene metrics are `not_applicable` because strand, phase and complete transcript semantics are absent. The fixed `strand=+`, `phase=0` writer may be shown only as a separately labelled placeholder-writer lower-bound diagnostic and never enters model ranking.

For ANNEVO, Helixer and Tiberius, report on the same seqids and the same reference primary-transcript policy, subject to a small explicit metric-capability table:

- exact CDS exon, intron and CDS-chain metrics when the artifact supplies or unambiguously supports them;
- distinguish coding-only exact CDS-transcript/gene from full exon/UTR transcript/gene metrics;
- strand and phase accuracy only when represented by the artifact;
- the existing CDS-span/FPR/count metrics for continuity with earlier reports.

Unsupported baseline metrics are `not_applicable`, not zero.

For cached SegmentNT features:

- record the actual cache path, feature order, seqids/splits and extraction `tile_bp` before scoring;
- compute test AUCPR directly from continuous probabilities;
- select each element's F1 threshold on that species' existing validation seqid only and apply it once to the test seqid;
- separately report exon, intron, splice-donor and splice-acceptor AUCPR/F1 against (a) the frozen primary-transcript reference and (b) an all-isoform multilabel-union reference;
- report the positive-class prevalence beside AUCPR so rare-boundary performance is interpretable.

The existing cache was extracted in independent 6,000 bp tiles. A weak result supports only the statement that the current 6-kb cached extraction has weak signal; it does not reject the official SegmentNT checkpoint or a longer-context extraction. This is a diagnostic of an already released checkpoint, not a trained candidate and not a paper claim.

### Execution

Initiate the submission from local Codex over SSH as explicitly requested by the user. On the cluster, obey the execution-time `cluster_config.yaml` submission mode and live scheduler state. The project-specific Baobab skill must select the CPU partition and walltime at submission time; the plan does not hard-code a partition. The job requests no GPU.

The job must terminate non-zero on a parser or metric failure. No catch-all exception handler, fallback parser, checksum layer or retry wrapper is permitted.

### M24 decision output

M24 must answer, without changing the product definition:

1. Does M19 already localize exact CDS boundaries often enough to justify adding structural heads to the GENERanno route?
2. Do official SegmentNT donor/acceptor and exon/intron outputs carry usable clean-plant structural signal?
3. Which apparent differences versus ANNEVO/Helixer/Tiberius remain after identical seqid filtering?
4. Is the direct route blocked primarily by emissions, grouping/decoding, strand/phase representation, or all three?

M24 itself has no invented numeric promotion gate. After its result table exists, one short Pro review will choose the first GPU branch and freeze a metric gate against the observed baseline distribution.

## 4. First GPU experiment: conditional branch, one fit

Experiment family: `M25-DIRECT-STRUCTURE-EMISSION-s0`

Before viewing the M24 result, freeze the new held-out M25 test species as `Setaria viridis`, RefSeq `Setaria_viridis_v4.0`, assembly `GCF_005286985.2`, annotation release `RS_2025_03`. It is absent from the current repository's train/calibration/test records. M24's Arabidopsis/rice labels may guide branch selection, so those species become development-only after M24. Setaria must not enter training, validation, threshold calibration or branch selection. Its disclosed clean-provenance status applies to the NT-v2/SegmentNT route; GENERanno remains `overlap_unknown` even on this species.

Only one branch is submitted first, with fixed seed 0:

- If SegmentNT has useful donor/acceptor and exon/intron test signal, use the official SegmentNT checkpoint/features as the first clean-provenance structural-emission route.
- If SegmentNT boundary signal is weak but M19 has strong exact CDS-boundary localization, minimally extend the existing GENERanno trainer with strand/phase/boundary auxiliary outputs.
- If neither condition holds, do not submit a GPU job. Report that the direct structural route lacks an emission basis and return to scientific design; do not switch to candidate-GFF refinement without a new user decision.

Whichever branch is selected must:

- use raw FASTA only at inference;
- train/calibrate only on predefined non-Setaria species and evaluate once on frozen Setaria without target-label calibration;
- preserve exact nucleotide coordinates through writing;
- produce a strand-aware, phase-aware diagnostic GFF rather than collapsing to 3 classes;
- use one fit/seed for discovery;
- compare against its own unchanged-input ablation and the same-scope released baselines;
- declare success and stop criteria before submission.

Longer context, additional training species and clade-specific adapters/MoE are deferred. They are not tested until a single direct structural-emission branch demonstrates usable exact boundaries and valid GFF output. Candidate conditioning remains out of scope under choice A.

## 5. Fourteen-day sequence after Pro approval

- Days 1–2: implement the single diagnostic evaluator/test/sbatch; run M24 on existing artifacts.
- Day 3: write the M24 result report and expose all unsupported metrics and provenance limitations.
- Day 4: Pro reviews the result table and chooses exactly one M25 backbone branch with a frozen metric gate.
- Days 5–8: make the smallest branch-specific model/writer change and run one bounded smoke only if it exercises a code path not exercised by the real run.
- Days 9–12: submit one M25 seed/fit using live Baobab routing.
- Days 13–14: evaluate exact structures and decide whether the winner merits two additional seeds.

In parallel with M24, perform only read-only work that does not change the experiment: inspect public GENERanno provenance material and prepare a concise species/accession question for GenerTeam. Do not contact the authors without separate user authorization.

## 6. Seed and stopping policy

- New architecture discovery: one fixed seed.
- Failed gate: stop; no extra seeds.
- Passed gate: add two seeds later, giving three total for a claim-facing estimate.
- Do not run five seeds unless observed variance or reviewer requirements justify it.
- Passing FPR alone is insufficient for M25; structural metrics and valid strand/phase output are required.

Metric failure is a research result. It must not trigger automatic threshold tuning, loss replacement or resubmission.
