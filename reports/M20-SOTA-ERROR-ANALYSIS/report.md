# M20-SOTA-ERROR-ANALYSIS

Same clean-plant panel, same CDS-span evaluator. GENERanno rows are our adapted models; ANNEVO/Tiberius/Helixer rows are released fixed-model baselines.

## Aggregate Metrics

| Model | Kind | gbF1 | Precision | Recall | Spec | FPR | Gene count ratio | FPR<=0.01 |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| ANNEVO-Magnoliopsida | released_fixed_model | 0.9269 | 0.9563 | 0.8993 | 0.9883 | 0.0117 | 0.7263 | False |
| Tiberius-angiosperm | released_fixed_model | 0.9252 | 0.9667 | 0.8871 | 0.9927 | 0.0073 | 0.6280 | True |
| GENERanno-1.2B-LoRA-s1 | adapted_pretrained | 0.8815 | 0.9611 | 0.8141 | 0.9935 | 0.0065 | 0.8299 | True |
| GENERanno-1.2B-LoRA-s0 | adapted_pretrained | 0.8421 | 0.9534 | 0.7541 | 0.9917 | 0.0083 | 1.0827 | True |
| Helixer-land_plant | released_fixed_model | 0.0000 | 0.9194 | 0.9246 | 0.9784 | 0.0216 | 0.8204 | False |

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
