# M20-CLAIM-CLEAN-PANEL-FREEZE

- Date: 2026-06-21
- Status: `BLOCKED_FOR_GENERANNO_CLEAN_HELDOUT_CLAIM`
- Scope: GENERanno 1.2B CDS-preview / 0.5B base as candidate backbones for clean held-out gene-annotation claims.

## Verdict

The current Arabidopsis/rice GENERanno panel is **not claim-clean**. Public documentation still does not expose a complete species/accession manifest for eukaryote pretraining plus CDS-preview annotation tuning, so absence from training cannot be certified.

Do not spend claim-grade GENERanno GPU on the current Arabidopsis/rice clean-plant panel until a public or author-provided species/accession manifest excludes the claim species. Use M18/M19/M20 as adaptation/comparability evidence, or move the clean held-out claim to a backbone/training protocol with certifiable provenance.

## Panel Freeze

| Species | Data path | Path exists | Project role | GENERanno overlap status | Claim use |
|---|---:|---:|---|---|---|
| `arabidopsis_thaliana` | `data/m1_screen/arabidopsis_thaliana` | True | current clean-plant train/val/test screen species | `overlap_unknown` | adaptation/comparability only |
| `oryza_sativa` | `data/m1_screen/oryza_sativa` | True | current clean-plant train/val/test screen species | `overlap_unknown` | adaptation/comparability only |
| `arabidopsis_lyrata` | `data/m13_distance_screen/arabidopsis_lyrata` | False | close-plant diagnostic candidate | `overlap_unknown` | diagnostic only unless manifest exclusion is obtained |
| `gallus_gallus` | `data/m1_screen/gallus_gallus` | True | animal negative-control diagnostic | `overlap_unknown` | diagnostic only |
| `drosophila_melanogaster` | `data/m1_screen/drosophila_melanogaster` | True | animal negative-control diagnostic | `overlap_unknown` | diagnostic only |

## Claim Consequence

- `allowed`: present GENERanno results as pretrained-backbone adaptation, same-panel comparability, and mechanism evidence.
- `blocked`: present Arabidopsis/rice M18/M19/M20 as clean no-overlap held-out SOTA evidence.
- `unlock condition`: a complete public/author-provided training and fine-tuning manifest excludes the exact claim species/accessions, or the claim is moved to a backbone/protocol with controlled provenance.

## Sources Used

- GENERanno 1.2B CDS annotator preview model card: https://huggingface.co/GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview — Confirms model identity/task; no complete eukaryotic species/accession exclusion list.
- GENERanno 0.5B eukaryote base model card: https://huggingface.co/GenerTeam/GENERanno-eukaryote-0.5b-base — Confirms broad eukaryotic pretraining; does not certify plant exclusion.
- GENERanno GitHub README: https://github.com/GenerTeam/GENERanno — Confirms release lineage and eukaryotic annotation focus; no full training manifest.
- GENERanno eukaryote pretraining dataset card: https://huggingface.co/datasets/GenerTeam/pretrain_data_eukaryote — RefSeq-derived broad eukaryotic corpus family; public card is not a species exclusion manifest.
- M19 local provenance audit: refs/dossiers/m19_generanno_provenance_audit.md — Prior project audit already marked Arabidopsis/rice overlap as unknown.
