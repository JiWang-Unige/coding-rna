# Dossier: GENERanno: A Genomic Foundation Model for Metagenomic Annotation / eukaryote 1.2B CDS annotator preview

- slug: `generanno-2025` · type: sota · added: 2026-06-09
- Paper/preprint: https://doi.org/10.1101/2025.06.04.656517
- Direct PDF: blocked by Cloudflare during archive.
- Repo: `refs/repos/generanno-2025/` cloned @ `0e3cb65`
- GitHub: https://github.com/GenerTeam/GENERanno
- HF model: https://huggingface.co/GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview
- Why relevant: user-requested early probe for eukaryotic CDS annotation behavior.

## Dataset source

- Repo lists `GenerTeam/cds-annotation`, `GenerTeam/gener-tasks`, and `GenerTeam/prokaryotic-gener-tasks`.
- Eukaryotic CDS example uses `hf://datasets/GenerTeam/cds-annotation/examples/fly_GCF_000001215.4.parquet`.
- Source audit found `GenerTeam/annotation_data_eukaryote` visible but without a complete public dataset card.

## Metric implementation

- HF model card metadata lists `f1`.
- No public eukaryote CDS/gene benchmark table, exact label schema, or split protocol was found.
- Public statement says more technical details are coming soon.

## Split scheme

- Unknown for eukaryote CDS preview.
- Do not use as claimable SOTA until split and schema are published or independently reconstructed under our own benchmark.

## Weights / license

- Repo license: MIT.
- HF model card license: MIT.
- Model card: `GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview`, pipeline tag `token-classification`.

## Architecture / training

- Repo overview describes GENERanno as a genomic foundation model with context length 8k bp and 500M base model parameters trained on 386B bp eukaryotic DNA.
- Eukaryote 1.2B CDS annotator preview is listed as trained on 1T data.
- CLI: `python src/tasks/downstream/cds_annotation.py --organism eukaryote`; BF16 recommended when supported; multi-GPU optional.

## Reproducibility notes

- Reproducibility: `unknown` for benchmark reproduction; `moderate` for a simple inference probe if weights download works.
- Worth reproducing: `partial`.

## Relevance to our project

- Probe raw CDS annotation behavior early, especially false positives in intergenic regions.
- Not a published SOTA anchor and not sufficient for gene-level claim without our own decoder/evaluation.
