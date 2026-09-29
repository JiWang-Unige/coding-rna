# M28 implemented Methods — R4 working draft

2026-09-30. This is an implementation-grounded Methods draft, not a results section
or a claim of novelty. The implementation is fixed at
[2bb6d3a](https://github.com/JiWang-Unige/coding-rna/tree/2bb6d3a249d81f0af3a97b5454fdd71072affb40).
The [R4 pre-fit contract](../reports/M28-METHOD-RESTART/r4_pilot_contract.md)
controls the running pilot. At drafting, complete paired fit/DEV results are still
pending; the procedures below describe the configured analysis, not its success.

## Task, data and reference policy

The model takes genomic DNA and fixed model weights as inference inputs. The
current output target is a primary protein-coding CDS chain, not a complete
UTR annotation, all transcript isoforms or a biological proof of gene function.
Both orientations of each allowed chromosome are processed independently.
Internal intervals are zero-based and half-open.

The pilot uses the existing Arabidopsis thaliana and Oryza sativa TRAIN/DEV
allowlists, reconstructed and checked against the saved chromosome split files by
[census.py](../scripts/experiments/M28-METHOD-RESTART/census.py). It does not use the
legacy M25R training sampler or its fitted LoRA/head. The DEV chromosomes are
NC_003074.8 (23,459,830 bp) and NC_089041.1 (29,936,421 bp). Test chromosomes and
Setaria remain sealed. Chromosome separation is not asserted to provide homology
isolation or an unseen-species test.

Within each gene, the complete-chain reference is the longest CDS-assessable
isoform, with transcript-ID tie-breaking, selected before applying M28 training
support restrictions. The frozen primary DEV reference contains 5,437 and 2,291
unique chains, respectively. Unsupported splice structures, overlapping primary
genes and genes longer than the model window are retained in this evaluation
denominator. Other assessable isoforms are recorded separately, not substituted
for the selected primary when scoring. Genes without an assessable complete
isoform can provide certain local supervision from the longest available
isoform but cannot become complete-chain positives. Reference construction and
training-label alignment are implemented in
[common_ruler.py](../scripts/experiments/M28-METHOD-RESTART/common_ruler.py) and
[align_labels.py](../scripts/experiments/M28-METHOD-RESTART/align_labels.py).

## Shared sequence representation

Windows span 24,576 bp with a 12,288-bp stride and a terminal anchored window.
Short sequences are oriented before right-padding. Each window is divided into
four 6,144-bp blocks and tokenized without special tokens into nonoverlapping
six-mers. Whole-window and concatenated block token IDs must agree.

The frozen GENERanno 1.2B CDS-annotator backbone uses revision
b0483c23b6b63787b61a6d3a204a9b517d6ba345, with its original released weights,
evaluation mode, BF16 and SDPA. The 4,096 by 2,048 last-hidden-state matrix per
window is cached in BF16 and supplied identically to both arms. No backbone
optimizer updates are made. Backbone attention remains block-local; cross-block
context is supplied by the trainable head, not by a claim of full-window
foundation-model attention.

The head projects each token vector to 128 dimensions and applies a bidirectional
GRU with 64 units per direction. Its outputs are repeated across the six
corresponding bases. A separate five-channel A/C/G/T/N branch applies two
32-channel one-dimensional convolutions (kernel widths 9 and 5), each followed
by GELU. The concatenated 160-dimensional representation is mapped to 128
dimensions by a linear layer, GELU and layer normalization. Padded bases are
masked. A shared linear head produces 15 local-state logits: intergenic,
phase-specific CDS/intron/donor/acceptor states, start and stop.
See [core.py](../src/m28/core.py) and
[cache_features.py](../scripts/experiments/M28-METHOD-RESTART/cache_features.py).

## B1 proposals and complete-chain scores

B1 adds four endpoint logits and a donor–acceptor link network. Motif-compatible
sites comprise ATG starts, TAA/TAG/TGA stops and GT/GC–AG splice junctions.
The generator keeps eight sites per endpoint type per 1,024-bp bin, eight
ORF-valid exon edges per start/carry/end-kind combination, and eight links per
donor. Links use the two boundary feature vectors and log intron length.
Intron length is at least five bases.

The exon graph carries the actual unfinished codon string across splice
junctions, rather than phase alone, to reject internal stops created across
exons. Approximate proposal search retains two beam states per event/carry and
at most 128 exons per chain. Exon proposal scores combine an endpoint logit with
mean CDS log-odds; paths also accumulate start, acceptor and link scores. The
final proposal limit is a window-wide cap of 8 times ceil(valid length / 1,024)
chains, not a spatial quota in each bin. These are actual pruned
free proposals, not an oracle enumeration of all feasible chains.

For each exon, B1 concatenates the mean base representation, its first and last
base representations, log exon length, log preceding intron length (zero for the
first exon), and cumulative coding phase divided by two. A linear/GELU map gives
an exon vector z_j. The complete-chain gain is

g(c) = sum_j a(z_j) + r(GRU(z_1, ..., z_m)) - n(mean_i h_i),

where a and n are linear maps and r is a two-layer tanh network. The null term
is shared by candidates in the same window. Multiple nonoverlapping genes can
coexist; candidates are not normalized by a single window-wide softmax.
A positive gain is the fixed calling operating point, not a calibrated claim
that natural-population correctness probability exceeds 0.5.

Discrete site selection, pruning and proposal generation are outside autograd.
Loss gradients reach the neural scoring modules and shared features for the
selected/injected structures. The nonadditive network reranks a bounded
additive proposal pool; it is not a globally optimal nonadditive decoder.
See [candidates.py](../src/m28/candidates.py),
[core.py](../src/m28/core.py) and [training.py](../src/m28/training.py).

## Paired training and C0 control

Both arms use the same 1,536 saved TRAIN draws (768 per species, 1,491 unique
windows), the same draw order repeated three times, and seed 0. Within each
chromosome/strand stratum, sampling mixes 75% uniform windows with 25%
gene-first sampling; it chooses an original-manifest eligible gene uniformly
then one of its supporting grid windows uniformly. Strata without eligible
genes use uniform sampling. The recorded marginal probability q(w) is retained
after label-policy alignment. Each window loss is multiplied directly by
1 / (N q(w)), where N is the number of natural windows in that stratum.
The batch-one weight is not divided away by normalization to a batch weight
sum. The three passes are repeated exposure, not independent data replicates.

C0 uses local fine-state likelihood, with grouped region likelihood where only
the broader local class is resolved. B1 adds endpoint, link and complete-chain
binary losses, all with coefficient one. Each binary output column averages its
known positive and negative class means with equal mass when both are present;
one observed class has unit mass and no known examples contribute zero.
Alternative-isoform conflicts, partial/unsupported structures and noncallable
bases are masked according to [labels.py](../src/m28/labels.py).

For each B1 training window, free proposals are generated and saved using only
DNA and neural scores before reference-derived targets or injections are built.
Complete positive chains and known positive
links can then be injected for training, with origin labels retained. Injected
positives must not be counted as inference candidate recovery. Nonexact
supervised candidates are reference-inconsistent examples, not experimentally
proven nongenes.

The configured fit is 4,608 batch-one updates per arm, FP32 trainable heads,
AdamW (learning rate 0.0003, weight decay 0.01, betas 0.9/0.999, epsilon 1e-8),
and global gradient clipping at one. C0 and B1 have 366,191 and 565,495 trainable
parameters, respectively, with matched shared/local-head initialization.
The only main-evaluation checkpoint is the last update; no DEV epoch selection
or seed/threshold/proposal-budget sweep is part of this pilot.

C0 decodes its local probabilities through the separately obtained, unmodified
ANNEVO HMM source at commit 37bdd9aa62ddf24fa55941fb827061f7ed49ce53. The bridge
explicitly permutes phase columns, preserves native probability clipping and
zero transition penalties, and uses internal minimum intron length one
(physical length five bases). A decoded gene's gain is the native path objective
difference obtained by replacing that gene's entire span with intergenic states,
including incoming/outgoing transitions. This is neither CDS-length-normalized
nor a full-chromosome HMM run. C0 does not call the ANNEVO neural model.
See [c0_decode.py](../src/m28/c0_decode.py) and
[fit_pilot.py](../scripts/experiments/M28-METHOD-RESTART/fit_pilot.py).

## Chromosome output and evaluation

Both arms infer all 8,690 oriented DEV windows without reference input or
reference-chain injection. Coordinates are mapped back to each chromosome and
strand. The owner of a candidate is the window whose center is nearest its
CDS-span midpoint, breaking ties toward the lower window start. Only
owner-generated candidates compete; chains seen only in nonowner windows are
counted as lost. Identical chains are deduplicated with the highest same-arm
score. Weighted interval scheduling selects a maximum-total-positive-gain,
same-strand nonoverlapping set; zero gain retains null and ties retain the
previous dynamic-programming solution. Opposite strands are assembled separately.

GFF3 export converts the interval start to one-based coordinates and assigns CDS
phase in transcription order. Exported chains are parsed back and checked for
coordinate/chain preservation. Structural diagnostics are reported without
post-hoc removal of predictions.

Exact-chain precision, recall and F1 use unique strand-aware CDS coordinate
chains. Species-wise metrics, arithmetic macro averages and count-pooled micro
metrics are distinct outputs. CDS-base scores use an unstranded union and cannot
replace chain accuracy. A separate historical table retains the 6,450-chain
parent reference. Baseline ANNEVO, Helixer, Tiberius and M25R predictions are
unchanged historical caches rescored on the common ruler, not new runs of the
latest tool releases.

Background is the complement of all annotated gene-feature spans on either
strand. Predicted-span/CDS coverage and wholly-background chains per Mb describe
reference-relative burden; they do not establish biological false-positive
rates. Actual primary-chain recovery is decomposed across the free proposal,
owner and final selection stages. B1's additive-only readout uses the same
fitted model and free candidates; it is a diagnostic, not a separately trained
causal ablation. B1 versus C0 is a whole-system comparison, not an isolated
effect of the chain GRU.
See [assembly.py](../src/m28/assembly.py),
[export.py](../src/m28/export.py),
[infer_pilot.py](../scripts/experiments/M28-METHOD-RESTART/infer_pilot.py) and
[evaluate_pilot.py](../scripts/experiments/M28-METHOD-RESTART/evaluate_pilot.py).

## Execution and claim limits

The pilot is capped at eight allocated GPU hours and 200 GiB of new artifacts
(180 GiB for features). Allocated time includes CPU/I/O waits inside GPU jobs.
Head/backbone/proposal/decoding timings and scheduler allocation are reported
separately. Incomplete arms cannot support a formal paired ranking.

These two species have already informed development. Backbone pretraining/CDS
training exposure is unresolved; the pilot cannot support a clean zero-shot or
independent-generalization claim. Mechanistic novelty, superiority to current
tools and biological support for reference-discordant calls remain empirical
questions. This draft supplies reproducible implementation text; it does not
replace the missing paired results, independent evaluation or publication.
