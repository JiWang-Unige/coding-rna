# M13 failure sanity: M12A Arabidopsis->rice

## Verdict

- M12A failure is not just final constrained-decoder post-processing noise. The selected Arabidopsis-calibrated operating points transfer poorly to rice before final GFF scoring: pre-decode false genic rate on true rice intergenic bases is much higher than M11 pooled training.
- The dominant pattern is cross-species emission/calibration shift plus fragmentation. M12A keeps reasonable true-CDS genic sensitivity, but it marks too many rice intergenic bases as genic and produces too many gene-body runs/genes after decode.
- Test-oracle diagnostic grid found valid `FPR<=0.01` and `gene_count<=1.25` rice operating points in 0/3 M12A seeds. With this grid, the failure is not merely an Arabidopsis-to-rice calibration-transfer bug; rice emissions and fragmentation remain incompatible with the hard guardrails even under diagnostic test-label oracle selection.
- Therefore a close Arabidopsis-relative scan is justified as a distance diagnostic: if a near plant also fails, stop the single-species fixed-model generalization route; if it succeeds while rice fails, reframe as distance-limited transfer.

## Aggregate comparison

| family | test rice FPR | test rice gbF1 | test rice gene_count_ratio | predecode intergenic false-genic rate | predecode true-CDS genic rate |
|---|---:|---:|---:|---:|---:|
| M11 pooled Arabidopsis+rice | 0.0082±0.0013 | 0.6965±0.0142 | 1.156±0.114 | 0.0161±0.0042 | 0.6752±0.0097 |
| M12A fixed Arabidopsis->rice | 0.0311±0.0049 | 0.6556±0.0125 | 1.755±0.207 | 0.2320±0.0604 | 0.7619±0.0481 |

## Per-seed selected calibration transfer

| family | seed | selected | VAL FPR | VAL gbF1 | VAL gcount | rice TEST FPR | rice TEST gbF1 | rice TEST gcount |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| M11 | 0 | b2p5_mcl60_mfg20 | 0.0087 | 0.8310 | 1.035 | 0.0088 | 0.6954 | 1.224 |
| M11 | 1 | b3p0_mcl60_mfg20 | 0.0095 | 0.8423 | 1.061 | 0.0092 | 0.7112 | 1.219 |
| M11 | 2 | b1p5_mcl60_mfg0 | 0.0095 | 0.8214 | 0.925 | 0.0067 | 0.6828 | 1.025 |
| M12A | 0 | b3p0_mcl60_mfg20 | 0.0085 | 0.8965 | 1.004 | 0.0328 | 0.6521 | 1.990 |
| M12A | 1 | b2p0_mcl60_mfg0 | 0.0098 | 0.9088 | 0.853 | 0.0349 | 0.6694 | 1.675 |
| M12A | 2 | b3p0_mcl60_mfg20 | 0.0096 | 0.8998 | 0.952 | 0.0256 | 0.6452 | 1.599 |

## M12A rice test-oracle diagnostic

| seed | valid FPR<=0.01 & gcount<=1.25? | best valid / best sane point | FPR | gbF1 | gcount |
|---:|---|---|---:|---:|---:|
| 0 | False | b4p0_mcl120_mfg0 | 0.0127 | 0.4959 | 1.596 |
| 1 | False | b4p0_mcl120_mfg0 | 0.0190 | 0.5562 | 1.545 |
| 2 | False | b0p0_mcl120_mfg20 | 0.0261 | 0.6297 | 1.239 |

## Next action

1. Freeze one close Brassicaceae/near-dicot species with high-quality genome+GFF provenance.
2. Run M13 only as a bounded single-seed distance scan: train/calibrate on Arabidopsis, test close plant and rice; no test-label calibration.
3. Treat fly/chicken only as diagnostic/negative controls unless overlap-clean status is resolved.

Machine-readable summary: `/srv/beegfs/scratch/shares/ds4dh/common/coding-rna/reports/M13_FAILURE_SANITY/summary.json`
Per-seed table: `/srv/beegfs/scratch/shares/ds4dh/common/coding-rna/reports/M13_FAILURE_SANITY/per_seed.tsv`
