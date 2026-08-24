# Baseline Reproduction Ledger / SOTA 复现中央账本

> v4.1 central ledger. Migrated on 2026-06-14 from legacy `docs/12_baseline_reproduction.md`; the legacy file is preserved for compatibility, but this file is now the preferred read-first ledger for SOTA/baseline reproduction state.

## 0. Status
- Required before own model iteration: yes, unless a waiver is explicitly recorded in `docs/03_benchmark_roadmap.md` or `ACTIVE_GOAL.json`.
- Current gate: partial_migrated
- Last migrated: 2026-06-14

## 1. Migration Notes
- Source of truth before v4.1: `docs/12_baseline_reproduction.md`.
- v4.1 source of truth from now on: `docs/20_baseline_reproduction.md`.
- The migrated record below preserves the M1 baseline reproduction details and warnings, including that pretrained-inference ceilings are not same-budget screen anchors.
- Future `/reproduce-baselines` or `/sota-randomized` work should append structured entries above or below the migrated record and register evidence in `docs/15_evidence_register.md`.

## 2. Migrated Legacy Record

# Baseline Reproduction Record (M1)

> Authoritative recipe to reproduce the three full gene-caller baselines (Tiberius, Helixer,
> ANNEVO) used in M1. Written 2026-06-10. Goal: any future session / colleague can re-run
> these exactly. Pairs with `docs/11_evaluator_contract.md` (the evaluator) and the verified
> ledger `docs/02_sota_model_inventory.md`.

## ⚠️ READ FIRST — what these runs ARE and ARE NOT

All three baselines here are **PRETRAINED models doing inference/annotation** on our pilot
genomes. They were trained by their authors on large multi-species data.

- These results = a **pretrained-inference reference ceiling** (what a fully-trained published
  model achieves on our pilot species under our fair CDS evaluator). Useful as: runner/env/
  evaluator validation, an upper reference, and a comparability anchor for metric definitions.
- These results are **NOT** a same-budget `screen_anchor`. The `screen_anchor` (per CLAUDE.md
  §2 / ACTIVE_GOAL) must be the **`-like` / `-light` reference architectures RANDOMLY
  INITIALIZED and trained under the SAME small-sample screen protocol** (sample fraction,
  epochs, patience, seeds, preprocessing, metric) that our own from-scratch candidate
  architectures use. Comparing a from-scratch small-sample model against a fully-pretrained
  caller is the exact unfair "small-sample-vs-large-sample-SOTA" comparison the two-tier
  anchor system exists to prevent.
- **TODO (true screen_anchor):** implement small-sample training of Tiberius-like / Helixer-like
  / ANNEVO-light reference architectures (random init), train under the unified screen
  protocol, take the max. The current `screen_anchor=0.9213` is a placeholder = pretrained
  ceiling and must be replaced (or explicitly redefined) before it gates Track A promotion.

## Pilot species (runner/evaluator validation; gene-dense outliers, NOT final anchor species)

| Species | Clade | RefSeq accession | genome | reference | local path |
|---|---|---|---|---|---|
| Saccharomyces cerevisiae | Fungi | GCF_000146045.2_R64 | genome.fa | reference.gff3 | `data/m1_screen/saccharomyces_cerevisiae/` |
| Drosophila melanogaster | Insecta | GCF_000001215.4_Release_6_plus_ISO1_MT | genome.fa | reference.gff3 | `data/m1_screen/drosophila_melanogaster/` |

Downloaded via `scripts/download_refseq_accessions.py`; integrity in `data/m1_screen/check_data_report.json`.

## Evaluator (cross-tool fair)

Always `--span-mode cds` for cross-tool comparison (CDS-only gene-body span; see docs/11).
```
python scripts/eval_gene_body_mask.py --reference-gtf <ref.gff3> --prediction-gtf <pred> \
  --genome-fasta <genome.fa> --output-json <out.json> --experiment-id <id> --profile screen --span-mode cds
python scripts/aggregate_gene_body_metrics.py --metrics <sp1.json> <sp2.json> \
  --output-json <agg.json> --experiment-id <id> --profile screen      # base-weighted + macro + per-species
```

## Tiberius  (exp: BASE-TIBERIUS-PILOT-M1)

- Container: `refs/repos/tiberius-2024/singularity/tiberius_2.0.5.sif` (built under srun; see docs/10).
- Env: `coding-rna` (launcher). Bind: `APPTAINER_BINDPATH=/srv:/srv,/home:/home`.
- Models: `fungi` (S.cer), `insecta` (D.mel) — current multi-clade release configs.
- Output: GTF (CDS/exon/intron/start/stop/transcript; **no UTR**).
- Run: `scripts/run_BASE-TIBERIUS-PILOT-M1.sbatch` (private-teodoro-gpu RTX3090).
- VERIFIED metric (cds span): S.cer F1 0.9888 / FPR 0.0186; D.mel F1 0.8413 / FPR 0.0227. Base-weighted F1 0.8608 / macro 0.9150.

## Helixer  (exp: BASE-HELIXER-SAC-DMEL-SMOKE-M1)

- Container: `refs/containers/helixer-docker_latest.sif` (Apptainer). Weights: `refs/weights/helixer-2025/{fungi/fungi_v0.3_a_0100.h5, invertebrate/invertebrate_v0.3_m_0100.h5}`.
- Env: `coding-rna` (launcher).
- Run command (per species): `apptainer run --nv $SIF Helixer.py --lineage <fungi|invertebrate> --model-filepath <weight.h5> --subsequence-length <21384 fungi|213840 animal> --fasta-path <genome> --species <Name> --gff-output-path <pred.gff3> --batch-size 16`. NOTE `--model-filepath` (offline-safe) REQUIRES `--subsequence-length`.
- Output: GFF3 (gene/mRNA/exon/CDS/**five_prime_UTR/three_prime_UTR**) — UTR present → MUST eval with `--span-mode cds`.
- Run: `scripts/run_BASE-HELIXER-SAC-DMEL-SMOKE-M1.sbatch`. (transcript-span eval gives a misleading S.cer FPR 0.654; cds span gives 0.033.)
- VERIFIED metric (cds span): S.cer F1 0.9869 / FPR 0.0333; D.mel F1 0.9118 / FPR 0.0224. Base-weighted F1 0.9213 / macro 0.9494.

## ANNEVO  (exp: BASE-ANNEVO-SAC-DMEL-SMOKE-M1)

- Env: dedicated `annevo` conda env (sanctioned exception; NOT coding-rna). Build: `scripts/setup_annevo_env.sh` (= `mamba env create -f refs/repos/annevo-2026/ANNEVO.yml -n annevo`), then `mamba install -n annevo -c conda-forge "setuptools<81"` (torchmetrics 0.8.2 needs pkg_resources). Py3.10 / torch 2.1.0 / cu12.1.
- Models: `refs/repos/annevo-2026/saved_model/ANNEVO_Fungi.pt` (S.cer), `ANNEVO_Insecta.pt` (D.mel).
- Run command (from repo dir): `python annotation.py -g <genome> -m saved_model/ANNEVO_<Clade>.pt -l <Fungi|Insecta> -o <pred.gff> --batch_size 32 -t <cpus> --show_log`.
- **Gotchas**: (1) wrap `conda activate annevo` in `set +u/-u` (MKL activate.d unbound var). (2) set a SHORT node-local `TMPDIR=/tmp/annevo_$JOBID` (decoding multiprocessing pymp sockets hit AF_UNIX 108-char limit on long beegfs paths). (3) run from the ANNEVO repo dir.
- Output: GFF (gene/mRNA/exon/CDS/UTR) — eval with `--span-mode cds`.
- Run: `scripts/run_BASE-ANNEVO-SAC-DMEL-SMOKE-M1.sbatch` (routed to shared-gpu when private full; ANNEVO is light: 3.8GB, ~18min two species).
- VERIFIED metric (cds span): S.cer F1 0.9735 / FPR 0.0072 (lowest FPR, prec 0.9971); D.mel F1 0.9122 / FPR 0.0352. Base-weighted F1 0.9197 / macro 0.9429.

## Three-tool summary (pretrained inference, CDS span, pilot species)

| Tool | base-w F1 | macro F1 | base-w FPR | note |
|---|---:|---:|---:|---|
| Tiberius | 0.8608 | 0.9150 | 0.0225 | CDS-only output |
| Helixer | 0.9213 | 0.9494 | 0.0228 | max on this metric |
| ANNEVO | 0.9197 | 0.9429 | 0.0341 | lowest S.cer FPR; published-SOTA strength is broad-clade locus/exon, not this |

`pretrained-inference ceiling (max) = Helixer 0.9213`. This is NOT the same-budget screen_anchor (see READ FIRST).

## Reproduction status
- Runner + env + evaluator: VERIFIED reproducible for all three.
- Metric implementation + dataset rawness + span definition: VERIFIED (docs/11).
- Same-budget screen_anchor (random-init small-sample training of `-like` refs): NOT YET DONE — required before Track A promotion gating.

