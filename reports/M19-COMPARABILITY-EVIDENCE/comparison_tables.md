# M19 Comparability Evidence

Scope: same-evaluator, same clean-plant panel comparison for paper-facing utility tables. Rows marked pending will be filled automatically after M19 seed metrics exist.

## Clean Plant Aggregate
| Model | Type | Status | gbF1 | gbF1@0.005 | gbF1@0.01 | gbF1@0.02 | Spec | FPR | Macro spec | Gene count ratio | Utility |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GENERanno-1.2B-s0 | our adapted/fine-tuned pretrained-CDS backbone | available | 0.8494 | 0.0000 | 0.8494 | 0.8494 | 0.9929 | 0.0071 | 0.9943 | 0.864 | FPR<=0.01; gene count sane |
| GENERanno-0.5B-base-s0 | our adapted/fine-tuned generic pretrained backbone | available | 0.6561 | 0.0000 | 0.0000 | 0.0000 | 0.9033 | 0.0967 | 0.9290 | 1.617 | FPR>0.02; over-calls genes |
| M19-GENERanno-1.2B-s0 | our adapted/fine-tuned pretrained-CDS backbone | available | 0.8421 | 0.0000 | 0.8421 | 0.8421 | 0.9917 | 0.0083 | 0.9936 | 1.083 | FPR<=0.01; gene count sane |
| M19-GENERanno-1.2B-s1 | our adapted/fine-tuned pretrained-CDS backbone | available | 0.8815 | 0.0000 | 0.8815 | 0.8815 | 0.9935 | 0.0065 | 0.9947 | 0.830 | FPR<=0.01; gene count sane |
| Tiberius | released fixed model | available | 0.9252 | 0.0000 | 0.9252 | 0.9252 | 0.9927 | 0.0073 | 0.9936 | 0.628 | FPR<=0.01; under-calls genes |
| ANNEVO | released fixed model | available | 0.9269 | 0.0000 | 0.0000 | 0.9269 | 0.9883 | 0.0117 | 0.9903 | 0.726 | FPR 0.01-0.02; under-calls genes |
| Helixer | released fixed model | available | 0.9220 | 0.0000 | 0.0000 | 0.0000 | 0.9784 | 0.0216 | 0.9793 | 0.820 | FPR>0.02; gene count sane |

## Clean Plant Per Species
| Model | Species | gbF1 | gbF1@0.01 | Spec | FPR | Gene count ratio | Utility |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GENERanno-1.2B-s0 | arabidopsis_thaliana | 0.9144 | 0.9144 | 0.9973 | 0.0027 | 0.779 | FPR<=0.01; gene count sane |
| GENERanno-1.2B-s0 | oryza_sativa | 0.7542 | 0.7542 | 0.9913 | 0.0087 | 1.036 | FPR<=0.01; gene count sane |
| GENERanno-0.5B-base-s0 | arabidopsis_thaliana | 0.7813 | 0.0000 | 0.9818 | 0.0182 | 1.187 | FPR 0.01-0.02; gene count sane |
| GENERanno-0.5B-base-s0 | oryza_sativa | 0.4981 | 0.0000 | 0.8761 | 0.1239 | 2.482 | FPR>0.02; over-calls genes |
| M19-GENERanno-1.2B-s0 | arabidopsis_thaliana | 0.9222 | 0.9222 | 0.9974 | 0.0026 | 0.930 | FPR<=0.01; gene count sane |
| M19-GENERanno-1.2B-s0 | oryza_sativa | 0.7226 | 0.0000 | 0.9897 | 0.0103 | 1.389 | FPR 0.01-0.02; over-calls genes |
| M19-GENERanno-1.2B-s1 | arabidopsis_thaliana | 0.9358 | 0.9358 | 0.9971 | 0.0029 | 0.761 | FPR<=0.01; gene count sane |
| M19-GENERanno-1.2B-s1 | oryza_sativa | 0.8038 | 0.8038 | 0.9923 | 0.0077 | 0.968 | FPR<=0.01; gene count sane |
| Tiberius | arabidopsis_thaliana | 0.9553 | 0.9553 | 0.9959 | 0.0041 | 0.634 | FPR<=0.01; under-calls genes |
| Tiberius | oryza_sativa | 0.8687 | 0.8687 | 0.9913 | 0.0087 | 0.613 | FPR<=0.01; under-calls genes |
| ANNEVO | arabidopsis_thaliana | 0.9473 | 0.9473 | 0.9956 | 0.0044 | 0.736 | FPR<=0.01; under-calls genes |
| ANNEVO | oryza_sativa | 0.8902 | 0.0000 | 0.9850 | 0.0150 | 0.702 | FPR 0.01-0.02; under-calls genes |
| Helixer | arabidopsis_thaliana | 0.9506 | 0.0000 | 0.9816 | 0.0184 | 0.811 | FPR 0.01-0.02; gene count sane |
| Helixer | oryza_sativa | 0.8708 | 0.0000 | 0.9769 | 0.0231 | 0.845 | FPR>0.02; gene count sane |

## Broad Diagnostic Aggregate (M17)
| Model | Type | gbF1 | gbF1@0.01 | Spec | FPR | Macro spec | Gene count ratio | Utility |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ANNEVO-M17 | released fixed/clade-matched models | 0.9115 | 0.0000 | 0.9760 | 0.0240 | 0.9721 | 0.840 | FPR>0.02; gene count sane |
| Tiberius-M17 | released fixed model | 0.8791 | 0.0000 | 0.9827 | 0.0173 | 0.9786 | 0.556 | FPR 0.01-0.02; under-calls genes |
| Helixer-M17 | released fixed/clade-matched models | 0.8797 | 0.0000 | 0.9474 | 0.0526 | 0.9453 | 0.931 | FPR>0.02; gene count sane |

## Interpretation

- M18 GENERanno 1.2B is not a random baseline: it is FPR-valid and gene-count sane on clean plants, unlike the 0.5B base result.
- Clean-plant released callers still define the high-gbF1 frontier: ANNEVO/Tiberius/Helixer are around 0.922-0.927 gbF1; Tiberius is the closest practical comparator because it also passes FPR<=0.01, but it under-calls genes.
- Current GENERanno evidence should be written as pretrained-CDS backbone adaptation/comparability, not clean no-overlap held-out SOTA, until provenance clears or a cleaner species panel is selected.

## Artifacts

- `clean_plant_aggregate.csv`
- `clean_plant_per_species.csv`
- `broad_panel_aggregate.csv`
- `broad_panel_per_species.csv`
- `comparison_tables.json`
