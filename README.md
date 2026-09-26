# coding-rna

Cross-species ab initio protein-coding gene annotation experiments on the UNIGE Baobab cluster.

This public repository is a lean research snapshot. It contains source code, experiment configurations, Slurm submission scripts, tests, compact metrics, and result summaries. Raw genomes, reference annotations, model weights, caches, full runtime outputs, logs, and generated prediction dumps are intentionally excluded.

## Current evidence and research status — 2026-09-26

The current fixed M25R route is **NO-GO**, not a publication-ready general-purpose annotator. The broader research objective remains open. Engineering validity, reference-match accuracy, biological validity, and independent generalization are separate judgments.

Start with the [current research roadmap](docs/research_roadmap_20260926.md). Read the evidence by scientific question, not experiment numbering:

The new [M27 WT replay result](reports/M27-ALLELIC-INTEGRITY/wt_result.md) verifies exact whole-chr22 native/cache replay for ANNEVO (530 CDS transcripts) and Tiberius (575), after [reference-only preparation](reports/M27-ALLELIC-INTEGRITY/preparation_result.md) fixed 21 paired sites. Unique WT-exact targets are 12 and 10, with 9 shared. **No mutant inference has run**: this is a measurement prerequisite, not an allelic effect, biological validation, or publication claim.

| Question | Evidence | What it establishes |
|---|---|---|
| How does the fixed decoder compare with cached callers on the same scope? | [Same-scope comparison](reports/M26-SAME-SCOPE-MECHANISM/result.md) | B improves over A but trails the three cached baselines in exact-chain F1. |
| Where are remaining reference chains lost? | [Candidate support partition](reports/M26-SAME-SCOPE-MECHANISM/support_result.md) | Reference-assisted support bounds for this fixed graph, not attainable model accuracy. |
| Do reference-source differences survive aligned coding eligibility? | [Semantic source pairing](reports/M26-SAME-SCOPE-MECHANISM/semantic_reference_pair_result.md) | The two Arabidopsis sources share 7,970 eligible CDS isoform chains; the local source-dispute hypothesis is closed. |
| Does fixed-budget global boundary ranking retain more complete chains? | [Endpoint-budget experiment](reports/M26-SAME-SCOPE-MECHANISM/endpoint_budget_result.md) | 685 gains and 910 losses, net −225, with opposite species directions; this fixed proposal is closed. |

The historical 6,450-chain mechanism cohort and the later CDS-assessable source comparison have different denominators and must not be merged. Baseline caches are historical, not a claim to have rerun every latest release under perfectly reconstructed provenance. GENERanno training exposure remains unresolved; adaptation results are not clean zero-shot evidence. Setaria remains sealed.

Each linked result identifies its contract, scripts, compact metrics, and execution evidence. The full runtime inputs are not included here: cloning this snapshot alone is insufficient to reproduce genome-scale jobs. Historical publication strategies and framework-era records are retained as history; the dated current roadmap takes precedence. ChatGPT Pro discussions are AI methodological consultation, not independent experiments or peer review.

## Research workflow (2026-09-08)

Research is driven by the user's current conversation and the local AGENTS.md/Skills, not by Auto Research. No cluster-side ACTIVE_GOAL value, reviewer quorum, external CLI review, approval sentinel, or framework ledger is required. Goals are discussed with the user; an absent machine goal is not a blocker. Complete authorized implementation, verification, and result handling without asking again for routine steps. External review (including ChatGPT Pro) is optional and only used when the user requests it.

Historical docs, goals, review records and experiment protocols remain evidence. Their generic Auto Research instructions (mandatory tri-review/pursue/pivot/context-pack/goal signatures) are retired, not current operating rules. Preserve concrete scientific controls: data isolation, held-out embargoes, frozen metrics/thresholds and original results. Changes to those controls or materially different compute/research scope still require user direction; removing this framework does not authorize training, test-set access or publication.

Use cluster_config.yaml only for compute configuration, verifying live scheduling facts when submitting. Research conda environments, shared user tools, source code, configs, data and outputs are retained. Legacy launchers which required the removed framework were archived unchanged outside this project; current M25/M25R launchers remain in place. Consult docs/auto_research_retirement_2026-09-08.md for the archive and exact scope. Do not reinstall or regenerate the retired framework as a startup step.

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
