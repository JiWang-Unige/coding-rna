# Dossier: Annotating the genome at single-nucleotide resolution with DNA foundation models / SegmentNT

- slug: `segmentnt-2025` · type: sota · added: 2026-06-09
- Paper: https://www.nature.com/articles/s41592-025-02881-2
- PDF: `refs/pdfs/segmentnt-2025.pdf`
- Repo: `refs/repos/segmentnt-2025/` cloned @ `2dc37b8`
- GitHub: https://github.com/instadeepai/nucleotide-transformer
- HF: https://huggingface.co/InstaDeepAI/segment_nt
- Why relevant: foundation-model segmentation reference and early probe, not a full-gene claim anchor.

## Dataset source

- HF model card: trained on all human chromosomes except chromosomes 20 and 21 for test and chromosome 22 for validation.
- Uses randomly sampled training sequences and fixed 30 kb sliding-window validation/test sequences.
- Labels derive from GENCODE / ENCODE-style gene and regulatory annotations; no separate public dataset card was found.

## Metric implementation

- Paper/reports use MCC, auPRC, Jaccard, F1, and SOV for element segmentation.
- Labels are multi-label genomic elements, not mutually exclusive gene grammar states.
- Full gene-model/locus comparison is not directly comparable to Tiberius/Helixer/ANNEVO without extra decoding.

## Split scheme

- Human chromosome split: chr20/21 test, chr22 validation, remaining chromosomes train.
- Homology exclusion is reported by source audit, but potential distal regulatory homology leakage is still noted in paper discussion.

## Weights / license

- HF model card license: CC BY-NC-SA 4.0.
- Repo license: CC BY-NC-SA 4.0.
- Paper license: CC BY-NC-ND 4.0.
- Model: `InstaDeepAI/segment_nt`; repo docs also support `segment_nt_multi_species`.

## Architecture / training

- Nucleotide Transformer v2 500M encoder with language-model head removed.
- 1D U-Net segmentation head with two downsampling and two upsampling convolutional blocks.
- HF card reports 53M segmentation-head parameters, 562M total.
- Trained on 3 kb, 10 kb, 20 kb, then 30 kb sequences; effective batch size 256; DGXH100 8 GPUs for 3 days; 23B tokens.
- 6-mer tokenization path does not handle `N` cleanly for SegmentNT docs.

## Reproducibility notes

- Reproducibility: `moderate` for inference/probing.
- Worth reproducing: `partial`; use as base-wise representation/logit source, not as full-gene SOTA.

## Relevance to our project

- Early probe for gene-body, intron, splice donor/acceptor, exon and UTR signals.
- If useful, attach a structured decoder under our benchmark rather than claim raw SegmentNT output as a gene caller.
