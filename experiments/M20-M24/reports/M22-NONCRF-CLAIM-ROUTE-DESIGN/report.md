# M22-NONCRF-CLAIM-ROUTE-DESIGN

Date: 2026-06-23

## Question

After M21 refuted `GENERanno 1.2B + trained CRF decoder`, should M22 spend GPU on more calibration/post-processing, on a non-CRF emission/objective change, or on clean-provenance backbone transfer?

## Inputs

- M21 pivot: abandon trained GENERanno-CRF decoder route; do not tune CRF/HMM-style transition decoders unless the route is materially different and preserves M19 FPR.
- M19 best non-CRF adaptation point: seed1 `gbF1=0.8815`, `intergenic_FPR=0.0065`, `gene_count_ratio=0.830`.
- M20 same-panel clean-plant comparators: Tiberius `gbF1=0.9252`, `FPR=0.0073`, `gene_count=0.628`; ANNEVO `gbF1=0.9269`, `FPR=0.0117`; Helixer `gbF1=0.9220`, `FPR=0.0216`.
- M21 provenance scout: NT-v2 remains the cleanest public-provenance plant backbone; GENERanno 1.2B remains `overlap_unknown` and can only support adaptation/challenger evidence unless provenance clears.

## Local diagnostic

I replayed the saved M19 validation calibration grid on the clean-plant test split as a diagnostic oracle. This is not a valid production selector; it is a ceiling check for whether decode/calibration still hides enough gbF1 under the hard FPR/gene-count constraints to justify another calibration iteration.

Hard-valid definition:

- `test_intergenic_FPR <= 0.01`
- `test_gene_count_ratio <= 1.25`

The replay evaluated 48 candidate decode settings per seed, 96 total. Metrics were written under each M19 output directory:

- `outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s0/calibration/test/`
- `outputs/M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS-s1/calibration/test/`

## Results

| Seed | Selected tag | Selected test gbF1 | Selected FPR | Selected gene_count | Hard-valid candidates | Best hard-valid tag | Best hard-valid gbF1 | Best hard-valid FPR | Best hard-valid gene_count | Best FPR tag | Best FPR | Best-FPR gbF1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| s0 | `b2p0_mcl60_mfg20` | 0.8421 | 0.0083 | 1.083 | 48 | `b0p0_mcl60_mfg20` | 0.8603 | 0.0094 | 1.011 | `b2p0_mcl120_mfg0` | 0.0070 | 0.7714 |
| s1 | `b0p0_mcl60_mfg20` | 0.8815 | 0.0065 | 0.830 | 48 | `b0p0_mcl60_mfg20` | 0.8815 | 0.0065 | 0.830 | `b2p0_mcl120_mfg0` | 0.0046 | 0.7950 |

## Interpretation

The calibration/post-processing frontier is effectively exhausted for the useful M19 seed. Seed1's validation-selected setting is also the test-oracle hard-valid best setting. Seed0 has a small hidden gain (`0.8421 -> 0.8603`) from a different threshold, but it remains well below seed1 and still far below the released fixed-caller gbF1 frontier.

The low-FPR behavior is real and valuable, but the remaining gap is not a decode-grid problem. Pushing FPR even lower with longer minimum CDS length collapses gbF1 (`0.7950` on seed1 best-FPR), confirming that simple filtering trades away recall rather than recovering genes.

## Decision

M22 should not launch another calibration/post-processing GPU or a CRF tuning run.

Primary next GPU direction:

1. `M22-GENERANNO-1P2B-NONCRF-EMISSION-OBJECTIVE`: keep the non-CRF head/argmax decode that preserved low FPR in M19, but change the training objective or emissions to target recall/gene recovery without relaxing `FPR<=0.01`.

Parallel claim route:

2. `M22-NTV2-CLEAN-BACKBONE-TRANSFER`: translate the strongest FP-aware recipe to NT-v2 or another clean-provenance backbone. This is more claimable but may have weaker emissions; keep it as a parallel/claim-boundary line rather than replacing the immediate non-CRF mechanism test.

Do not run:

- trained CRF retuning
- HMM/transition decoder tuning
- decode-only calibration sweeps as the main M22 GPU action

## Gate For Implementation

Before submitting new GPU jobs:

- implement a distinct non-CRF objective/emission change, not a decoder retune;
- run code-review gate against label mapping, FPR metric, raw-score output, and calibration contract;
- treat GENERanno results as adaptation/challenger evidence unless provenance is cleared;
- keep `M19 s1` as the local guardrail: the new run must preserve `FPR<=0.01` and aim to improve gbF1 beyond `0.8815` without gene-count explosion.
