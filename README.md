# coding-rna

Cross-species ab initio protein-coding gene annotation experiments on the UNIGE Baobab cluster.

This private repository is a lean research snapshot. It contains runnable source code, experiment configurations, Slurm submission scripts, tests, compact metrics, and result summaries. Raw genomes, reference annotations, model weights, caches, full runtime outputs, logs, and generated prediction dumps are intentionally excluded.

## Contents

- `src/` — model, data, decoding, and evaluation code.
- `scripts/` — experiment, aggregation, validation, and analysis scripts.
- `configs/` — experiment and benchmark configurations.
- `sbatch/` — Baobab Slurm submission scripts.
- `tests/` — focused regression tests.
- `reports/` — compact JSON/CSV/Markdown/HTML result artifacts.
- `docs/06_results_log.md` — chronological experiment results.
- `docs/10_findings.md` — consolidated research findings.
- `docs/12_publication_strategy.md` — publication-oriented assessment.
- `docs/14_validation_matrix.md` — validation coverage and remaining gaps.
- `docs/15_evidence_register.md` — evidence and provenance register.
- `docs/experiments/` — selected per-experiment summaries.

## HPC workspace

The full working directory is available after SSH login at:

```text
login1.baobab.hpc.unige.ch
/home/users/j/jwang/coding-rna
```

Large data and runtime artifacts remain only on the cluster and are ignored by Git.
