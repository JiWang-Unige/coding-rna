# Dossier: Tiberius: end-to-end deep learning with an HMM for gene prediction

- slug: `tiberius-2024` · type: sota · added: 2026-06-09
- Paper: https://academic.oup.com/bioinformatics/article/40/12/btae685/7903281
- PDF: `refs/pdfs/tiberius-2024.pdf` (added/read 2026-06-09; sha256 `72478a3bf04c71b8b9a6476188035487c72d3c6df8d9abd6d4c2eee1c93f6d31`)
- Repo: `refs/repos/tiberius-2024/` cloned @ `8c49fd0`
- GitHub: https://github.com/Gaius-Augustus/Tiberius
- Why relevant: strongest open differentiable-HMM structured-decoder baseline.

## Dataset source

- 2024 paper: 37 mammalian species; RefSeq annotations from NCBI.
- Genomes were soft-masked with RepeatModeler2, RepeatMasker, and Tandem Repeats Finder.
- Species were included only if annotation BUSCO completeness was >90%, to reduce apparent false negatives in reference annotations.
- The dataset covers model species such as human and mouse, but the test split excludes close taxonomic groups around the test species.
- Current repo configs extend to Mammalia, Vertebrata, Insecta, Angiosperms, Fungi, Diatoms, Chlorophyta and others.

## Metric implementation

- Paper emphasizes gene-level F1.
- Paper evaluates predicted gene structures against reference annotation at two levels: exon and gene. For each level it computes TP, FP, FN, recall = TP/(TP+FN), precision = TP/(TP+FP), and F1 as harmonic mean.
- BUSCO completeness is computed as an auxiliary annotation metric, but the paper explicitly warns it is not suitable for gene-structure accuracy and can hide false-positive genes/exons.
- Repo config benchmark tables use `gffcompare v0.12.10` and report exon Sn/Pr/F1 and gene Sn/Pr/F1.
- Training labels retain only one transcript isoform per gene; docs warn training works only on genes with one transcript isoform per gene.
- For fair gene-level comparison against methods that can emit isoforms, the paper uses only the transcript variant with the longest coding sequence for each predicted/reference gene set.

## Split scheme

- 2024 paper test species: human, cow, beluga whale.
- Validation species include Panthera pardus and Rattus norvegicus.
- Close clades Hominidae, Ruminantia, and Cetacea were excluded from training for held-out tests.
- The paper states the minimal evolutionary distance from a training species is 43 MYA to human and 64 MYA to cow/beluga.
- De novo human mode uses ClaMSA signals from a 64-species Zoonomia alignment; only human de novo results are reported because inference alignment was human-referenced.
- Current per-clade configs list test species and warn when species overlap with ANNEVO/Helixer training.

## Architecture details from paper

- Input: one-hot nucleotide sequence over A/C/G/T/N stacked with a masked-repeat track.
- Model: CNN + LSTM backbone plus differentiable, GPU-vectorized HMM layer; approximately 8M trainable parameters.
- Output: during training, probabilities for 15 gene-structure classes; during inference, full label/gene-structure sequences.
- Strand handling: one strand at a time; reverse strand is predicted by feeding the reverse complement.
- HMM: 15 hidden states: 1 intergenic, 3 intron, and 11 exon-related states. Exon states include reading-frame states plus start/stop codon and donor/acceptor splice-site border states.
- Biological constraints: reading-frame consistency, start/stop codons, canonical splice patterns, and no in-frame stop codons within exons; some intron-spanning in-frame stop cases are later filtered.
- Inference: parallel Viterbi on both strands, GTF output. Genome is processed in roughly 500 kb tiles, with a second pass over ~1 Mb windows around tile boundaries that appear genic.
- De novo mode: adds four ClaMSA-derived codon-start logits per position to the input.

## Training setup from paper

- Training sequence length: seamless tiles of 9,999 bp.
- Inference sequence length: paper figure reports 500,004 bp; current mammal config default is 400,050 bp.
- Hardware/time: 15 days on 4x NVIDIA A100 80 GB.
- Two-phase training: 6 days pre-HMM training, then 9 days end-to-end fine-tuning with HMM.
- Optimizer: Adam, learning rate 1e-4.
- Batch size per GPU: 250 during pre-HMM training, 128 during fine-tuning.
- Validation: every 1,000 training steps on selected leopard/rat genome segments; final model chosen by combined validation exon/gene F1.
- Loss: CCE-F1 loss, combining categorical cross entropy with an F1-loss over exon classes; lambda = 2. If a sequence has no exon labels, an exon false-positive-rate term is used instead.

## Ablations and design signals

- End-to-end HMM integration improves test F1 by +2.6 percentage points at gene level and +1.1 at exon level over pre-HMM.
- Increasing model size from ~2M to ~8M parameters improves gene/exon F1 by +4.9/+1.2 percentage points.
- Removing softmasking reduces average exon/gene F1 by only 0.2/0.9 percentage points on mammalian tests, but larger effects appear for more distant species.
- Replacing CCE-F1 with CCE reduces gene/exon F1 by 8.8/2.1 percentage points.
- Mapping the 15-class output to a 5-class label system reduces gene/exon F1 by 11.8/5.4 points.
- Removing both CCE-F1 and separate exon-border classes reduces gene/exon F1 by 21.1/10.3 points.

## Weights / license

- License: MIT.
- Example current mammalia config: `model_cfg/mammalia_softmasking_v2.yaml`, weights URL `https://bioinf.uni-greifswald.de/bioinf/tiberius/models/tiberius_weights_v2.tar.gz`, `tiberius_version: 1.1.5`, date 2025-05-12, default sequence length 400050.
- Current configs expose public `weights_url` entries for multiple clades.

## Reported values

- 2024 paper abstract: human gene-level F1 62% for Tiberius vs 21% for the next best ab initio method.
- 2024 paper average over human/cow/beluga ab initio comparison: Tiberius exon F1 89.7%, gene F1 55.1%; Helixer 72.9% / 19.3%; AUGUSTUS 67.3% / 12.4%.
- 2024 paper BUSCO average over the three test species: Tiberius 96.0%, Helixer 92.1%, AUGUSTUS 74.2%.
- 2024 paper comparison to extrinsic-evidence methods: Tiberius gene/exon F1 55.1% / 89.7%, BRAKER3 53.7% / 83.2%, GALBA 41.8% / 86.2%, BRAKER2 14.9% / 48.7%.
- 2024 paper human de novo mode: gene F1 65.5%, exon F1 92.6%; paper states two out of three human genes have exactly correct exon-intron structure in de novo mode.
- 2024 paper runtime table: Tiberius 1:39 h, Helixer 8:54 h, AUGUSTUS 2:25 h, BRAKER3 48:53 h.
- 2024 paper runtime: human de novo mode took 2:05 h, only 9 min longer than ab initio mode.
- Current `vertebrates.yaml`: average exon F1 85.8 and gene F1 55.3 across six non-mammalian vertebrate test species.

## Reproducibility notes

- Reproducibility: `moderate`.
- Code, config, weights, and containers are available, but paper snapshot and current repo version differ.
- PDF text has been read from local archive; no longer blocked for the 2024 paper.
- Preferred first baseline reproduction because license and setup are comparatively clean and metric tables are explicit.
- Need supplementary tables for exact species accession list, RNA-seq library list for BRAKER3, and full per-species metric tables.

## Baseline reproduction verified 2026-06-10

- Mini-smoke experiment: `BASE-TIBERIUS-MINISMOKE`.
- Execution: `srun` on baobab shared-gpu RTX4090 (`job_id=8527962`) using the official Tiberius Singularity image and repo-bundled `Panthera_pardus` data.
- Environment: project conda env `coding-rna`; Tiberius launcher installed editable from `refs/repos/tiberius-2024`; official SIF cached at `refs/repos/tiberius-2024/singularity/tiberius_2.0.5.sif` (14 GB).
- Command core: `python tiberius.py --singularity --genome outputs/BASE-TIBERIUS-MINISMOKE/data/inp/genome.fa --model_cfg mammalia_softmasking_v2 --out outputs/BASE-TIBERIUS-MINISMOKE/tiberius_prediction.gtf --seq_len 259992 --batch_size 8`.
- Critical runtime detail: the container must bind `/srv` explicitly via `APPTAINER_BINDPATH=/srv:/srv,/home:/home` or it cannot see project files.
- Metric implementation verified for mini-smoke: `tests/integration_tests/test_prediction.py` computes exact CDS F1 from `(chrom,start,end,strand)` CDS intervals and exact transcript-chain F1 from frozensets of CDS intervals.
- Reproduced values: CDS exact F1 `0.8594`; transcript-chain exact F1 `0.3124`; both pass the repo integration thresholds (`>=0.75` and `>=0.28`).
- Project metric probe: unconstrained gene-body F1 `0.9196`, intergenic FPR `0.0187`, constrained gene-body F1 `0.0` under `intergenic_FPR <= 0.01`.
- Limitation: this does not reproduce the 2024 paper's human/cow/beluga benchmark and must not be used as `screen_anchor` or `sota_benchmark`.

## Limitations / benchmark cautions

- 2024 paper model is mammal-trained; the paper says non-vertebrate annotation is not recommended without retraining.
- It predicts one label sequence per position and does not naturally represent alternative splicing.
- It cannot predict gene structures with spliced start codons.
- It does not prevent intron-spanning in-frame stop codons inside the HMM itself; such transcripts are removed by a simple post-filter.
- BUSCO is explicitly inadequate as a primary metric because false-positive genes/exons may not affect BUSCO.

## Relevance to our project

- Primary mechanism baseline for trainable structured decoding.
- Inspect differentiable HMM state design, long input length, transcript collapsing, softmasking, and gffcompare pipeline.
- Most actionable design lesson: gene prediction needs both grammar-aware decoding and metric-aligned losses/border labels; the large CCE-F1/15-class ablation gap is more relevant than ordinary hyperparameter tuning.
