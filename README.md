# coding-rna

Cross-species ab initio protein-coding gene annotation experiments on the UNIGE Baobab cluster.

This public repository is a lean research snapshot. It contains source code, experiment configurations, Slurm submission scripts, tests, compact metrics, and result summaries. Raw genomes, reference annotations, model weights, caches, full runtime outputs, logs, and generated prediction dumps are intentionally excluded.

## Current evidence and research status — 2026-09-29

The current fixed M25R route is **NO-GO**, not a publication-ready general-purpose annotator. The broader research objective remains open. Engineering validity, reference-match accuracy, biological validity, and independent generalization are separate judgments.

Start with the [current research roadmap](docs/research_roadmap_20260926.md). Read the evidence by scientific question, not experiment numbering:

The [M27 paired-native result](reports/M27-ALLELIC-INTEGRITY/paired_result.md) is complete: 44/44 whole-chr22 single-SNV runs, 3.7222 allocated GPUh. Each caller has two confirmed PTC-specific strict bypass responses and two unassignable pairs; fixed-denominator identification bounds are [0.1667, 0.3333] for ANNEVO (n=12) and [0.2, 0.4] for Tiberius (n=10). These are not confidence intervals. The shared nine loci have positive lower bounds but no shared confirmed bypass locus. The [frozen contract](reports/M27-ALLELIC-INTEGRITY/paired_native_contract.md), [raw compact results](reports/M27-ALLELIC-INTEGRITY/paired_result.json), and [execution record](reports/M27-ALLELIC-INTEGRITY/paired_native_execution.md) distinguish native response from unproven component attribution and biological validity. The preregistered next candidates are ANNEVO/POTEH and Tiberius/CRYBA4; no publication claim is established.

The [M27 novelty boundary](reports/M27-ALLELIC-INTEGRITY/novelty_scope.md) distinguishes this diagnostic from prior ACE/SGRF premature-stop experiments, early Helixer mutagenesis, and published pangenome annotation-consistency studies. A positive response in a modern caller would not by itself establish a new mechanism or a publication-level contribution.

The [native confirmation](reports/M27-ALLELIC-INTEGRITY/confirmation_result.md) reproduced all four preregistered target and whole-chromosome results, including ANNEVO/POTEH exon skipping and the Tiberius/CRYBA4 15-nt GC–AG intron. Tiberius filter snapshots and request metadata also match. These are computational repeats, not new independent observations.

The [single-edge test](reports/M27-ALLELIC-INTEGRITY/single_edge_result.md) is complete and negative for its frozen restoration endpoint: ANNEVO/POTEH same-PTC cache replay exactly reproduced all 530 chains, but releasing the sole registered CDS1_TA→CDS2 edge changed neither the bypass chain nor its score; the edge was not selected. Both arms used the same native min-intron pass and no new neural calls. This fixed pilot is closed without rescue combinations. The result does not establish neural-only causation, global ORF-constraint irrelevance, selective risk diagnosis, or biological truth. Total M27 usage remains 4.3133/8 GPUh and 34.30/50 GiB; the publication objective is open.

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
