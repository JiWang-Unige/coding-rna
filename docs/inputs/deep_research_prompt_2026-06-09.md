# Deep Research Prompt: Cross-species ab initio protein-coding gene annotation (FRESH)

把下面这段复制到 ChatGPT (Deep Research mode) / Gemini (Deep Research) / Claude (Deep Research)。建议同一份 prompt 跑三个平台，以便后续 `/research-synthesize` 处理冲突。

---

I am working on the following research problem:

Build a cross-species ab initio deep learning model for eukaryotic protein-coding gene annotation. The input is raw genomic DNA sequence. The output is a per-base multi-class annotation, at minimum distinguishing protein-coding gene body from intergenic sequence, and ideally using the label system used by current SOTA models (for example exon / intron / intergenic, CDS / UTR / intron / intergenic, or strand-aware labels if that is the SOTA protocol).

Task scope:
- Input: raw genomic DNA sequence only.
- Output: per-base multi-class labels for protein-coding gene annotation, including gene-body / intergenic distinction and any SOTA-relevant sublabels.
- Application setting: academic benchmark and model development, with the eventual goal of a Nature Methods-level model that strictly exceeds published deep learning SOTA under a comparable benchmark.
- Explicitly NOT in scope: RNA folding or RNA 3D structure prediction; RNA-seq/protein-homology evidence-based annotation; transcript-level classification only; promoter/enhancer/TF-binding/ncRNA annotation; private non-reproducible benchmarks.

Please write a Chinese-language research-status report covering:

1. Method families. Identify 4-6 method families used for eukaryotic ab initio protein-coding gene annotation. For each family give the core idea, representative works (paper title + year + venue/preprint + first author), strengths, weaknesses, and typical performance level. Include both modern deep learning methods and traditional gene prediction methods only where they remain informative.

2. Current SOTA candidates. Focus especially on Tiberius, Helixer 2025 version, ANNEVO, SegmentNT, and GENERanno-eukaryote-1.2b-cds-annotator-preview. For each candidate, provide paper URL, official GitHub URL, pretrained weights URL (HuggingFace / Zenodo / other if available), reported metric values, dataset/split, species used for training and testing, and whether the result is cross-species or same-species.

3. Datasets and benchmarks. List the main datasets and benchmark protocols used in this area. For each: source (Ensembl / GENCODE / RefSeq / custom), release/version, species, genome build, label definitions, train/val/test split scheme, whether held-out species evaluation is used, download URLs, and known issues such as annotation quality, homology leakage, chromosome-level leakage, or inconsistent transcript collapsing.

4. Evaluation metrics. What metrics are standard for this problem? Prioritize metrics that control intergenic false positives while measuring gene-body F1. Also report gene-level F1, segment-level F1, exon/intron/intergenic macro-F1, base-level precision/recall/F1, boundary accuracy, and any SOTA-specific scoring rules. For gene-level F1, specify exact matching rules: overlap threshold, strand consistency, transcript/gene collapsing, CDS-only versus full gene body, and handling of partial genes.

5. Known limitations of current SOTA. Specifically identify where current deep learning SOTA fails: cross-species generalization, long-context modeling, intergenic false positives, boundary precision, gene-level consistency, rare gene structures, fragmented predictions, low-quality assemblies, annotation noise, inference cost, or dependence on external pretraining.

6. Open research opportunities. Identify gaps explicitly called out in papers or visible from method comparisons. Focus on opportunities that could plausibly support a Nature Methods-level contribution: architecture-level innovation, long-context backbone, structured decoder/head, objective/loss design, gene-level consistency, cross-species representation learning, biological-constraint-aware decoding, or foundation-model adaptation.

7. Biological constraints and post-processing. Compare whether current neural models use or ignore classical gene prediction signals such as GT/AG splice donor/acceptor motifs, start/stop codon consistency, ORF constraints, phase consistency, and strand-aware decoding. Discuss whether these constraints are integrated into the model, objective, decoder, or post-processing, and where they could be useful for boundary filtering.

8. Reproducibility status. Which SOTA models have working code and weights? Which only have paper/code/weights partially available? Which are hard to reproduce due to missing data, private weights, ambiguous preprocessing, or undocumented metric scripts? For each model, estimate setup difficulty and whether it is worth reproducing first.

Constraints on your answer:
- Cite specific papers with full URLs. Do not say "studies have shown" without a citation.
- If you are uncertain about a number, metric, dataset, or split, say so. Do not fabricate.
- If different sources report conflicting numbers, list both and flag the conflict.
- Clearly distinguish primary-source facts from your inference.
- Use plain markdown. Do not use HTML.

Length: aim for 3000-6000 words.

---

## 怎么用这份 prompt

1. **三家固化**：把上面三横线之间的内容整段复制到以下三个 deep research 平台：
   - Claude (Deep Research)
   - Gemini (Deep Research)
   - ChatGPT (Deep Research mode)
2. 报告回来后直接粘贴进已经预生成好的占位文件：
   - `docs/inputs/deep_research_claude_20260609.md`
   - `docs/inputs/deep_research_gemini_20260609.md`
   - `docs/inputs/deep_research_chatgpt_20260609.md`
3. 三份报告全部粘贴完毕后调用 `$research-synthesize`。
