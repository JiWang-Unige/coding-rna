# M25R terminal summary

- Experiment: `M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0`
- Slurm job: `12116383`
- Terminal state: `COMPLETED`, exit `0:0`, elapsed `1-13:03:39`
- Experiment status: `STOP_M25_BRANCH`

## Scope and embargo

M25R used Arabidopsis/rice development splits and predefined primary nuclear chromosome allowlists. Chloroplast, mitochondrion and unplaced scaffolds were excluded. The frozen Setaria input was FASTA-only. No Setaria annotation was read, no Setaria prediction was generated, and no embargo release marker exists.

## Training validity

- Three epoch checkpoints exist on Baobab and all checkpoint tensors are finite.
- Training loss decreased from `0.8205138235` to `0.6669163239` to `0.6176625575`.
- The run used `1,536` training windows, `768` validation-loss windows, seed `0`, one fit and three epochs.
- The only stderr content is a PyTorch `use_reentrant` warning; no traceback, OOM or non-finite-loss error occurred.
- This run contains the minimal empty-structural-mask repair, so it supersedes the implementation-invalid M25 run for scientific interpretation.

## Frozen validation grid

The complete grid contains `5,625` finite tuples: three checkpoints, three genic-region thresholds and `5^4` start/stop/donor/acceptor thresholds.

| Quantity | Result |
|---|---:|
| admissible tuples | `0 / 5,625` |
| pass intergenic FPR `<=0.020` | `5,625 / 5,625` |
| pass gene-count ratio `0.80–1.20` | `0 / 5,625` |
| recorded structural-valid fraction pass | `5,625 / 5,625` |
| exact CDS interval F1 range | `0.0350929–0.1205736` |
| exact CDS-chain F1 range | `0.1273345–0.3249883` |
| intergenic FPR range | `0.0030302–0.0126640` |
| predicted-gene-count ratio range | `0.0953488–0.3305426` |

The structural-valid fraction is recorded as `1.0` for every tuple, but the current function makes this effectively tautological for non-empty output. It cannot be treated as independent structural-validity evidence.

## Best-ranked non-admissible tuple

- epoch `1`
- region threshold `0.4`
- start `0.5`, stop `0.5`, donor `0.1`, acceptor `0.1`
- exact CDS interval F1 `0.1204141`
- exact CDS-chain F1 `0.3249883`
- intergenic FPR `0.0124683`
- predicted-gene-count ratio `0.3252713` (`2,098 / 6,450` pooled development reference chains)

The best exact-interval tuple reaches `0.1205736` interval F1, `0.3240849` chain F1, FPR `0.0126630`, and count ratio `0.3299225`; it is also non-admissible.

Best-ranked complete-gene recovery declines with training epoch:

| Epoch | Exact interval F1 | Exact chain F1 | FPR | Gene-count ratio |
|---:|---:|---:|---:|---:|
| 1 | `0.1204141` | `0.3249883` | `0.0124683` | `0.3252713` |
| 2 | `0.0783022` | `0.2377457` | `0.0050774` | `0.1829457` |
| 3 | `0.0382618` | `0.1402276` | `0.0034152` | `0.1034109` |

## Outputs intentionally absent

No admissible tuple means there is no selected checkpoint, frozen decoder parameter file, Setaria raw-score directory, full GFF3, unchanged-input ablation GFF3, embargo release marker or Setaria evaluation JSON. Consequently strand-aware exact interval F1, exact chain/transcript F1, coding-gene F1, strand accuracy, phase accuracy, independent structural validity and full-minus-ablation metrics are **not measured**, not zero.

## Scientific verdict

M25R is a valid development no-go for the frozen combined system. Its dominant observed failure is severe under-recovery of complete genes, not intergenic false positives. Current artifacts cannot determine whether missing genes arise primarily from backbone representation, structural-head supervision, strand/phase scores, decoder thresholds, transition grammar or complete/unique-model filtering.

The only approved next analysis is a read-only development re-decode/error decomposition using these checkpoints. It must reproduce aggregates, account for every reference and prediction at each decoder stage, replace the tautological validity statistic with independent per-transcript validation, avoid Setaria annotation, and stop for review.

## Archived compact artifacts

- `config_resolved.yaml`
- `primary_chromosome_allowlist.json`
- `train_summary.json`
- `validation_grid_summary.json`
- `env.txt`
- `M25RGENSTRUCTs0_12116383.out`
- `M25RGENSTRUCTs0_12116383.err`
- `STATUS`
- `JOBID`

Epoch checkpoints and the full `validation_grid_diagnostics.json` remain on Baobab and are intentionally excluded from Git.
