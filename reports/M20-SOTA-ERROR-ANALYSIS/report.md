# M20-SOTA-ERROR-ANALYSIS

Same clean-plant panel, same CDS-span evaluator. GENERanno rows are our adapted models; ANNEVO/Tiberius/Helixer rows are released fixed-model baselines.

## Aggregate Metrics

| Model | Kind | ordinary gbF1 | historical constrained gbF1 | Precision | Recall | Spec | FPR | Gene count ratio | FPR<=0.01 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| ANNEVO-Magnoliopsida | released_fixed_model | 0.9269 | 0.9269 | 0.9563 | 0.8993 | 0.9883 | 0.0117 | 0.7263 | False |
| Tiberius-angiosperm | released_fixed_model | 0.9252 | 0.9252 | 0.9667 | 0.8871 | 0.9927 | 0.0073 | 0.6280 | True |
| GENERanno-1.2B-LoRA-s1 | adapted_pretrained | 0.8815 | 0.8815 | 0.9611 | 0.8141 | 0.9935 | 0.0065 | 0.8299 | True |
| GENERanno-1.2B-LoRA-s0 | adapted_pretrained | 0.8421 | 0.8421 | 0.9534 | 0.7541 | 0.9917 | 0.0083 | 1.0827 | True |
| Helixer-land_plant | released_fixed_model | 0.9220 | 0.0000 | 0.9194 | 0.9246 | 0.9784 | 0.0216 | 0.8204 | False |

2026-09-26 reporting correction: the old `gbF1` column contained `constrained_gene_body_F1`, not ordinary accuracy. Both fields above are copied from the preserved `summary.json`; no metrics were rerun. The historical constrained field is not recomputed from the separate FPR<=0.01 flag (ANNEVO demonstrates that they are different). Helixer's zero must not be described as zero prediction accuracy.

Scope warning: the separate historical `interval_overlap` diagnostic did not restrict full-genome prediction intervals to the reference FASTA seqids. Those rows/CSV are quarantined from scientific comparison pending same-scope rescoring; this does not by itself invalidate the independently loaded aggregate metrics above. M24 and M25R also use different chromosome scopes and denominators and must not be merged into one performance table.

## Interpretation

- Tiberius is the strongest released fixed-model comparator under the hard FPR guardrail, but it under-calls gene count relative to reference.
- ANNEVO has the best gbF1 among released fixed baselines on this panel, but aggregate FPR exceeds the `0.01` claim guardrail.
- Helixer strongly over-calls intergenic bases on this panel under the current evaluator, which makes it useful as a practical-specificity contrast.
- GENERanno LoRA is stable across two seeds and keeps FPR under `0.01`, but its remaining weakness is recall/gene recovery rather than specificity. The structured-decoder line should target this exact error mode.

## Artifacts

- `summary.json`
- `aggregate_metrics.csv`
- `per_species_metrics.csv`
- `interval_overlap.csv`
