# M28 B2: implemented methods and interpretation boundaries

Working Methods supplement, 2026-10-01. This document describes the implementation fixed for the B2 first fit and paired decoding experiment; it is not a results or novelty claim. The completed fit, full paired DEV evaluation and reference-level diagnosis are recorded in [b2_execution.md](../reports/M28-METHOD-RESTART/b2_execution.md). All 8,690 paired direction-windows and the 7,728-reference denominator were validated; the predeclared early-versus-late continuation criterion was not met. This is a bounded result for this version, not a general impossibility claim. The pre-fit [B2 contract](../reports/M28-METHOD-RESTART/b2_first_fit_contract.md) remains the authority for the data, budget, stopping rules and continuation criterion. This supplement does not amend it or replace the historical [R4 Methods](M28_implemented_methods_R4.md).

The B2 runtime was fixed before the first optimizer update at commit [63730f07d4c19a0850e9e1813a7c89747bd34403](https://github.com/JiWang-Unige/coding-rna/tree/63730f07d4c19a0850e9e1813a7c89747bd34403). Source pointers below refer to that implementation, not to a proposed alternative.

## Study question and scope

B2 asks whether applying a learned, history-dependent compatibility score before fixed search caps improves complete coding-chain recovery compared with applying the same learned score only after the original additive search. Both decoding arms use the same final new checkpoint, frozen base model and candidate-budget settings. Thus this is a single-seed exploratory comparison of scoring position, not two independently trained replicates.

The output target is a primary protein-coding CDS exon chain. Full UTR annotation, recovery of every alternative isoform, unrestricted genome-wide applicability, independent zero-shot generalization and superiority to current external tools are not established by this experiment. A successful B2 continuation criterion would justify further investigation on its stated DEV scope, not a publication-level claim.

## Frozen representation and data

The local representation, segmentation, endpoint and link parameters are inherited from the B1 checkpoint `outputs/M28-PILOT-R5/B1/step_004608.pt`. This entire base is kept in evaluation mode with gradients disabled. Its existing cached GENERanno features use revision `b0483c23b6b63787b61a6d3a204a9b517d6ba345`; neither the language model nor its cache is retrained or replaced. The frozen M28 core supplies a 128-dimensional position feature, denoted h_i, within each valid oriented window.

TRAIN and DEV retain their existing Arabidopsis and rice allowlists, 24,576-bp windows, 12,288-bp stride, terminal anchoring and orientation rules. Training metadata are filtered to TRAIN before reference-label parsing, and R4 supervision is joined by window identity rather than taken from old cached labels. DEV features and references do not enter the fit or checkpoint selection. Test and Setaria remain sealed.

The frozen backbone is the released `GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview`, an annotation-specialized release, not a verified exposure-free foundation checkpoint. In the [model card at the exact revision](https://huggingface.co/GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview/raw/b0483c23b6b63787b61a6d3a204a9b517d6ba345/README.md) inspected on 2026-10-01, no species or assembly-accession training manifest is supplied. The [author repository](https://github.com/GenerTeam/GENERanno) links a [CDS dataset card](https://huggingface.co/datasets/GenerTeam/cds-annotation/raw/main/README.md) that describes prokaryotic evaluation datasets; that card does not resolve training exposure for this eukaryotic checkpoint. This bounded source review leaves upstream genome and annotation exposure unknown. It establishes neither absence of overlap nor actual leakage.

The saved development census records the following allowed chromosome identities:

| Species | TRAIN | DEV |
| --- | --- | --- |
| Arabidopsis | NC_003070.9, NC_003071.7, NC_003076.8 | NC_003074.8 |
| Rice | NC_089043.1, NC_089044.1, NC_089046.1 | NC_089041.1 |

These sets are chromosome-disjoint within each species, not species-disjoint. The [split implementation](../src/screen_anchor/data.py) assigns sorted sequence identities by index modulo five; it does not group homologs. Its historical “leakage-safe” wording therefore applies only to whole-sequence coordinate separation, not to homology or upstream annotation exposure. Both [Arabidopsis](../reports/M28-METHOD-RESTART/arabidopsis_census.json) and [rice](../reports/M28-METHOD-RESTART/rice_census.json) censuses explicitly retain the homology-isolation limitation. No new sequence-level overlap or homology analysis was performed for this source review.

An early-versus-late result can still answer the frozen within-DEV scoring-position question because the upstream backbone and final head are shared. It cannot by itself establish independent homology-isolated or upstream-unseen generalization. Such a future claim requires a separately authorized evaluation protocol and adequate model-exposure evidence; it is not obtained by renaming the existing DEV panel or automatically opening Setaria.

The original 1,536 paired TRAIN draws, order, sampling probabilities and importance weights are repeated three times, with batch size one, for 4,608 new-head updates. Each epoch contains 768 draws from each species and 1,491 distinct window identities in total. Replacement sampling and repeated epochs are not independent observations. The frozen B1 base already received 4,608 updates before B2; the additional 4,608 head updates must not be reported as the total training exposure of the complete method or as matched exposure against C0/B1.

## Learned prefix value

The only trainable module is PrefixValue128, with 184,369 parameters. It has an exon encoder, a GRUCell over the exon history, a three-stage event embedding and an absolute scalar readout. There is no residual addition of the old cumulative proposal score.

Coordinates in this module are half-open intervals in the oriented window. For exon j with interval [a_j, b_j), its input concatenates the mean of h_i over the exon, h_(a_j), h_(b_j-1), and three scalars:

- log(1 + exon length);
- log(1 + gap from the preceding exon end), or zero for the first exon;
- the total coding length of the preceding exons modulo three, divided by two.

This 387-dimensional input is projected to 128 dimensions with a linear layer and GELU. A 128-to-128 GRUCell updates the history state u_j from that exon encoding and u_(j-1). The empty history is the zero vector.

The scalar value at an event is computed from the concatenation of the history state (128), a local event feature (128), a learned entry/donor/final embedding (16), and four geometry scalars (4). The 276-dimensional concatenation passes through linear-128, GELU and linear-1 layers. The geometry is log(1 + exon count), log(1 + chain span), log(1 + total coding length), and log(1 + distance from the last exon end to the entry position). Span is zero for the empty chain; the last distance is nonzero only for a nonempty entry history.

Entry uses h_(position), donor uses h_(position-1), and final readout uses h_(last exon end-1). An empty history is scored only at entry. Its learned entry value is distinct from the fixed zero final null score.

The search grammar separately tracks the actual unfinished codon string across introns and uses it for ORF validation, including junction-spanning stops. That string is not explicitly concatenated into the new network. The network instead sees the inherited sequence features and the preceding-exon coding-phase scalar. A description that the new head explicitly encodes the full grammar carry would therefore overstate the implementation.

The final value is a logit learned with balanced candidate supervision. It is not calibrated to the natural probability that a generated chain is correct, and it is not an admissible future-score bound. The implementation does not provide an A* or global optimality guarantee.

Source: [prefix_value.py](../src/m28/prefix_value.py).

## Search and scoring position

The inherited search uses its frozen endpoint logits, CDS probabilities, link head, geometry and ORF rules. Starts are ATG; terminal stops are TAA, TAG or TGA. Supported introns use GT/GC donors and AG acceptors, with a physical minimum length of 5 bp.

The fixed budget retains eight sites of each endpoint type per 1,024-bp block, eight ORF-valid exon edges per start/carry/end-kind, eight links per donor, two beam histories per event and actual carry, and at most 128 exons per chain. The final cap is 8 × ceil(valid window length / 1,024) chains over the whole window. It is a global final cap, not a separate spatial quota in each block.

A history consists of the old additive score, its exon tuple and its actual unfinished codon string. Duplicate chain/carry identities retain the larger old score. Histories compete within an event/carry group, including a fresh empty entry history when applicable. The late arm uses the original additive ordering at these cuts. The early arm orders by the new absolute prefix value, with a deterministic lexicographic identity tie-break, before the event beam and final cap.

Both arms preserve the endpoint, conditional-edge, link, grammar and budget rules. They need not visit the same histories, perform identical amounts of work or return the same final candidate pool: changing pruning order is intended to change reachability. In particular, the exon-edge limit is conditional on a start and an actual carry; it does not imply a single globally identical set of evaluated edges in both arms.

Discrete candidate generation is outside autograd. The recorded `proposal_scores` remain the old additive sums in both arms, for provenance and replay. They must not be interpreted as the new selection gains.

Source: [prefix_search.py](../src/m28/prefix_search.py), [prefix_inference.py](../src/m28/prefix_inference.py).

## Training supervision

The first fit is off-policy: all free candidates and first-cut observations come from the frozen original additive search, not from the evolving early policy. FirstCutObserver only observes actual search groups; it does not change their order, budget or candidates.

The complete-chain loss uses the frozen free pool augmented, for this loss only, with missing complete positive reference chains. Injected chains keep a separate origin label and never become free inference candidates. Known positives and known reference-inconsistent negatives enter a class-balanced binary cross-entropy loss. When both classes exist, each contributes one half; when only one known class exists, it has unit weight. Unknown labels (-1) are ignored, and a pool with no known examples contributes differentiable zero. An empty background window does not receive an invented positive chain.

The prefix term uses the first actual beam or final cut reached by an eligible positive reference path under the control search. Entry supervision matches the existing exon prefix and the next reference exon start; donor supervision matches a proper nonempty prefix; final supervision matches the complete chain. If an endpoint, exon edge, link or ORF restriction blocks the path earlier, no fictitious reachable prefix is injected.

For a beam cut, comparisons use the two actual survivors. For a final cut, comparisons use the top two survivors and the survivor at the cutoff, with duplicate identities removed. Only known-negative competitors form ranking pairs. Exact alternative-isoform prefixes and unknown or non-callable spans retain the inherited -1 protection; they are not forced into the negative class.

For each usable cut, the pair loss is the mean of softplus(v_negative - v_positive) over its known-negative competitors. The window prefix loss is the mean over cuts with at least one such competitor. Cuts without usable negatives do not dilute that mean; if none are usable, this term is differentiable zero.

The per-window objective is:

`L_window = importance_weight × (L_complete + L_prefix)`

The two component coefficients are both one. Each component first follows its stated within-window normalization; the original importance weight is not divided away by batch normalization. Logged component averages are unweighted, whereas the logged total is importance-weighted. Consequently, adding the reported component means does not generally reproduce the reported total.

Frozen-base features are obtained with no_grad rather than inference_mode, so they can be used by the new head's backward pass. PrefixWindow detaches the base features and uses cumulative sums for exon means. Its hidden-state cache is local to one window and one optimizer step, with separate no-grad search and gradient-bearing loss objects. New-head hidden states, logits and graphs are not cached across updates.

The fixed free pool and first-cut observations may be reused by window identity during this fit because their inputs, base, budget and labels do not change. New-head values and gradients are recomputed for every draw. This reuse is not evidence of on-policy training.

Source: [prefix_training.py](../src/m28/prefix_training.py), [fit_b2.py](../scripts/experiments/M28-METHOD-RESTART/fit_b2.py).

## Optimization and checkpoint

After constructing the frozen base, the new head is initialized with torch.manual_seed(0). Only its parameters are passed to AdamW: learning rate 3e-4, weight decay 0.01, betas (0.9, 0.999), epsilon 1e-8, FP32 and global gradient-norm clipping at 1.0. The learning rate is fixed, with no seed or optimizer sweep.

The official checkpoint is `step_004608.pt`, after the final planned update. Epoch checkpoints are retained as records, not selected by DEV performance. The completed fit and observed costs are recorded separately in [b2_fit_result.json](../reports/M28-METHOD-RESTART/b2_fit_result.json) and [b2_execution.md](../reports/M28-METHOD-RESTART/b2_execution.md). A falling TRAIN loss alone does not establish improved annotation or passage of the B2 criterion.

## Paired inference and controls

The paired driver evaluates all 8,690 frozen DEV direction-windows across both species and orientations. Both arms load the same final new-head checkpoint and frozen base. The late arm calls the literal original additive generator and applies the new head only to complete chains; the early arm uses the new prefix ranking before the fixed cuts. Inference does not read reference annotations or inject chains.

Within each window, the two candidate sets are merged into a sorted, identity-deduplicated union. A fresh PrefixWindow reads out each complete union chain once, and that value is mapped back into each arm. Common chains therefore receive exactly the same final gain in both arms. The early search's no-grad cache is not reused for this readout. This controls final batch/cache numerical differences while allowing candidate formation to differ.

For every window, the late arm must replay the saved R5 B1 candidate identities and order, retained links, counts, budget and reference-free flag exactly. Its old proposal scores must be finite and differ by at most 1e-5 in maximum absolute value. A mismatch invalidates the paired control; affected windows are not silently excluded and the tolerance is not relaxed.

The comparison is accepted as complete only after the entire frozen grid, both arms and all replay/common-value controls finish. Partial-window results are not used to rank variants or change the checkpoint, K, seed or gain threshold.

Source: [infer_b2.py](../scripts/experiments/M28-METHOD-RESTART/infer_b2.py).

## Assembly and evaluation

The assembly and evaluator reuse the existing frozen ruler. Chain ownership is assigned by span midpoint to the nearest window center, with smaller window start breaking ties. Duplicate chains use the maximum gain within the same arm. The unchanged weighted interval selector chooses nonoverlapping chains with positive gain on each strand; the final null score and gain threshold remain zero.

The primary reference contains 7,728 chains: 5,437 Arabidopsis and 2,291 rice. It is reconstructed against the frozen identity ledger, preserving its supported-target definitions and the retained long, overlapping or otherwise unsupported records rather than switching to an easier denominator. The historical 6,450-chain table remains a separately labelled secondary view.

Reported measurements include exact-chain precision, recall and F1 per species, arithmetic macro summaries and pooled-count micro summaries. Exact coverage is tracked through free generation, ownership and selection to expose where reference chains are lost. The final selected TP count must agree with the corresponding stage accounting.

CDS-base metrics use the final emitted CDS union, not the frozen local head's raw predictions. Background span coverage, CDS coverage and chains per Mb use the complement of reference gene spans on either strand. These are reference-relative burdens, not independently adjudicated biological false positives.

The final GFF3 is parsed back and compared with emitted chain identities and records. CDS-structure support is reported as a diagnostic, not used as a post hoc output filter. Reverse-oriented coordinates are converted back by the existing assembly rules.

C0/B1 and ANNEVO, Helixer and Tiberius are evaluated from the existing fixed caches with the same ruler; they are not new native runs or a claim to the latest available versions. Their known provenance and remaining version/exposure uncertainties are recorded in [baseline_cache_provenance.md](../reports/M28-METHOD-RESTART/baseline_cache_provenance.md). In particular, a cached baseline row must not automatically be described as independently rice-unseen zero-shot evidence.

Source: [evaluate_b2.py](../scripts/experiments/M28-METHOD-RESTART/evaluate_b2.py).

## Continuation criterion, resources and claim limits

The frozen criterion compares early against the same-head late arm. Both species must show strictly more exact free chains, strictly more exact owner chains and improved final exact-chain F1, while neither species increases background span coverage. A partial or single-species improvement does not pass. All results, including failures, are retained; no additional threshold, K, seed or training adjustment is used to rescue this version.

The formal first-fit, paired-inference and CPU-evaluation allocations are capped at 11 allocated GPU-hours, 4 independent allocated CPU-hours and 2 GiB of combined new artifacts. GPU-stage I/O and CPU waiting remain part of GPU allocation cost. Preparation and R6, R4/R5 and other historical costs remain separate records, not erased by a new round.

The paired driver shares base forward passes, union readout, I/O and control loading across arms. Its instrumented stage timings are not two independent standalone early/late wall-clock runs, and do not by themselves establish an inference speedup. Actual allocated resources, program timings and inherited training exposure must be distinguished in any later efficiency claim.

Timeout, OOM, nonfinite values, unexpected base gradients, input/label mismatch, failed replay, reference use in free inference or artifact-budget violations invoke the existing stop policy. This Methods supplement does not authorize retries, extension, data release or new experiments.

Learned prefix ranking and search-error training are established ideas. B2's frozen-control first-cut softplus/BCE objective does not implement the on-policy update/reset procedure of LaSO or beam-search optimization, and no corresponding guarantee is claimed. The task-specific contribution would need to be demonstrated by defensible complete-chain gains and broader evidence, not inferred from naming a new head. See [related work and open evidence requirements](M28_learned_search_related_work_20261001.md).
