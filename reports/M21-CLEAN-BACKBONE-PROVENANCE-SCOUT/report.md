# M21-CLEAN-BACKBONE-PROVENANCE-SCOUT

Date: 2026-06-21
Status: local provenance/claim-boundary report

## Question

M21 is running a non-claim GENERanno 1.2B CRF screen. In parallel, we need to decide which backbone can support a paper-grade claim on the clean plant panel, and which backbones should remain only adaptation/challenger evidence.

## Source checks

| Backbone/source | Public provenance signal | Claim implication |
|---|---|---|
| `InstaDeepAI/nucleotide-transformer-v2-500m-multi-species` | HF model card states the model was pretrained on 850 NCBI genomes and explicitly says plants and viruses were not included. HF dataset card for `InstaDeepAI/multi_species_genomes` repeats that plant and virus genomes were not taken into account. | Strongest current clean-provenance fit for `{Arabidopsis thaliana, Oryza sativa}` as held-out plant evidence, assuming no downstream supervised plant tuning is introduced. |
| `InstaDeepAI/segment_nt` | HF model card states the released SegmentNT model was trained on human chromosomes, with chr20/21 test and chr22 validation. It uses the NT-v2 encoder plus a segmentation head. | For plant claim, supervised SegmentNT head is not plant-trained. However, SegmentNT itself predicts human genomic element classes and is less aligned with our target than the NT-v2 unfreeze route. |
| `GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview` | HF card gives only a short CDS-annotator description and says more technical details are coming soon. GitHub lists this checkpoint as 1.2B / 1T / eukaryote, but does not expose a complete species/accession exclusion manifest. | `overlap_unknown`; current Arabidopsis/rice results cannot be claimed as clean held-out. Use as pretrained-CDS adaptation/challenger evidence only. |
| `GenerTeam/GENERanno-eukaryote-0.5b-base` | HF card/GitHub describe 386B bp eukaryotic pretraining, no complete species/accession exclusion manifest. Our runs show it is materially weaker than 1.2B CDS-preview. | Not claim-clean on current panel and not empirically competitive. Keep as ablation/negative control; do not allocate main GPU beyond targeted controls. |
| External released callers: Tiberius / ANNEVO / Helixer | Baseline provenance varies by lineage/weights and must be documented per panel. They are comparator tools, not candidate backbones for our adapted model. | Use same evaluator tables and mark released fixed model vs adapted/fine-tuned model explicitly. |

## Decision

The cleanest near-term paper backbone for plant held-out provenance is still NT-v2/SegmentNT-family, because public NT-v2 documentation explicitly excludes plants from the multi-species pretraining dataset. But the fixed-model NT-v2 route has already failed the key publication question: Arabidopsis-only or small-panel calibration does not generalize robustly across even close plants, rice, and animals.

GENERanno 1.2B CDS-preview is currently the strongest adaptation/challenger backbone empirically: M18/M19 show stable hard-FPR-valid clean-plant adaptation, and M21 is testing whether a real CRF head can recover gene-body F1/coherence. It should not be promoted to clean held-out claim on Arabidopsis/rice unless the GenerTeam provenance blocker is resolved or a new demonstrably clean panel is selected.

## Recommended next action

1. Let M21 finish as a non-claim mechanism screen.
2. If M21 CRF improves gbF1/gene count without breaking FPR, position GENERanno as "CDS-pretrained backbone adapted to low-label gene annotation", not as clean held-out pretraining generalization.
3. Keep the clean-claim route anchored on NT-v2 plant-exclusion provenance, but stop M9-only tuning unless a new data/architecture plan directly attacks fixed-model transfer.
4. For paper tables, separate:
   - released fixed callers: Tiberius / ANNEVO / Helixer;
   - clean-provenance adapted backbone: NT-v2 on plants;
   - high-signal overlap-unknown challenger: GENERanno 1.2B CDS-preview.

## Sources

- NT-v2 HF model card: https://huggingface.co/InstaDeepAI/nucleotide-transformer-v2-500m-multi-species
- NT-v2 multi-species dataset card: https://huggingface.co/datasets/InstaDeepAI/multi_species_genomes
- SegmentNT HF model card: https://huggingface.co/InstaDeepAI/segment_nt
- GENERanno 1.2B CDS-preview HF model card: https://huggingface.co/GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview
- GENERanno 0.5B base HF model card: https://huggingface.co/GenerTeam/GENERanno-eukaryote-0.5b-base
- GENERanno GitHub: https://github.com/GenerTeam/GENERanno
