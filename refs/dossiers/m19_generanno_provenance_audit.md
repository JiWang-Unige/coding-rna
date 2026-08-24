# M19-GENERANNO-PROVENANCE-AUDIT

- Date: 2026-06-19
- Scope: `GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview` and `GenerTeam/GENERanno-eukaryote-0.5b-base` provenance/overlap status for interpreting M18/M19 clean-plant results.
- Status: preliminary public-source audit; claim blocker remains open.

## Verdict

`GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview` should remain `overlap_unknown` for Arabidopsis/rice claim purposes. Public model cards and the GitHub README confirm the eukaryotic base/pretraining scale and the existence of the 1.2B CDS annotator preview, but do not provide a complete species/accession exclusion list for the CDS-preview model or its annotation fine-tuning data.

This does not invalidate M18/M19 as mechanism/challenger evidence. It means the current clean-plant numbers cannot be used as clean held-out SOTA claim evidence unless one of the following happens:

- GenerTeam releases a full species/accession training/fine-tuning manifest that excludes our claim species/splits.
- We choose a new claim panel with species/accessions provably absent from the released GENERanno training/fine-tuning data.
- We frame GENERanno results as pretrained-model adaptation/comparative evidence, not clean no-overlap generalization.

## Public facts checked

| Source | Public fact | Implication |
|---|---|---|
| HF model card `GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview` | Token-classification CDS annotator; says more technical details are coming soon; links GitHub for end-to-end genome annotation. | Confirms model identity and task, but not species/accession provenance. |
| HF model card `GenerTeam/GENERanno-eukaryote-0.5b-base` | 0.5B base, 8k context, trained on 386B bp of eukaryotic DNA; developmental; technical details pending. | Confirms broad eukaryotic pretraining, not no-overlap for plants. |
| GitHub `GenerTeam/GENERanno` README | News says eukaryotic annotation expert model released 2026-02-10; overview table lists `GENERanno-eukaryote-1.2b-cds-annotator-preview` with `1T` eukaryote data. | Confirms 1.2B CDS-preview uses larger eukaryotic data scale, but no species list. |
| HF dataset `GenerTeam/pretrain_data_eukaryote` | Gene-centric RefSeq-derived corpus; schema includes taxonomy and species_type; broad categories include `<pln>` plant; raw files are partitioned by broad clade categories, not a published exclusion list in the model card. | Plant overlap cannot be excluded from public docs alone. |
| HF dataset `GenerTeam/cds-annotation` | Public README describes prokaryotic CDS annotation summaries; repository includes prokaryotic summaries and example fly files. | Does not provide a claim-grade eukaryotic CDS-preview fine-tuning manifest. |
| Existing M12/M17 dossiers | M12 already marked GENERanno as `overlap_unknown`; M17 kept GENERanno as mechanism/challenger evidence. | M19 does not change claim status; it reinforces the blocker. |

## Species-level status for current panel

| Species | Current role | Public GENERanno overlap status | Claim use |
|---|---|---|---|
| `arabidopsis_thaliana` | M18/M19 train/val/test subset in clean-plant screen | unknown; plant category exists in released pretraining corpus family, no accession exclusion list | mechanism/challenger only |
| `oryza_sativa` | M18/M19 train/val/test subset in clean-plant screen | unknown; plant category exists in released pretraining corpus family, no accession exclusion list | mechanism/challenger only |
| `arabidopsis_lyrata` | M13/M17 diagnostic close plant | unknown | diagnostic only |
| `gallus_gallus`, `drosophila_melanogaster` | M17/M18 diagnostic animals | unknown for GENERanno; known/likely overlap caveats for some external baselines | diagnostic only |

## Consequence for M19

M19 can proceed as a screen / Track-B-preflight because the scientific question is now: "Does the 1.2B CDS-specialized backbone plus our 3-class FP-aware route produce stable, calibratable, high-specificity gene annotation behavior?" It cannot by itself close the paper claim. If M19 is strong, the next gate is to either find a claim-clean species panel or make the paper story explicitly about pretrained adaptation/practical calibration rather than no-overlap cross-species SOTA.

## Sources

- Hugging Face model card: `https://huggingface.co/GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview`
- Hugging Face model card: `https://huggingface.co/GenerTeam/GENERanno-eukaryote-0.5b-base`
- GitHub README: `https://github.com/GenerTeam/GENERanno`
- Hugging Face dataset card: `https://huggingface.co/datasets/GenerTeam/pretrain_data_eukaryote`
- Hugging Face dataset card: `https://huggingface.co/datasets/GenerTeam/cds-annotation`
- Local prior dossier: `refs/dossiers/m12_prereq_audit.md`
- Local prior dossier: `refs/dossiers/m17_pretraining_overlap_audit.md`
