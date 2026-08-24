# M12-PREREQ-AUDIT · publication-alignment preflight

- Date: 2026-06-17
- Scope: freeze the first publication-facing panel contract before M12A/M12B/M12C GPU work.
- Status: v0.1 frozen for implementation planning; not a SOTA claim.

## Decision summary

Use a two-layer M12 panel:

1. **Claim-clean lead panel for M9/NT-v2**: plants `{Arabidopsis thaliana, Oryza sativa}` because the official NT-v2 500M multi-species model card states that plants were not included in the 850-genome pretraining collection. This makes the existing clean plant panel the safest immediate fixed-model/generalization panel for the M9 lead route.
2. **Mechanistic challenger panel for GENERanno**: same plant panel, but marked `overlap_unknown`, because GenerTeam's public eukaryotic pretraining / annotation datasets are gene-centric RefSeq-derived and do not expose a complete species/accession exclusion list. GENERanno results can support mechanism/ablation and practical comparison, but should not be used as clean held-out claim evidence until the overlap is resolved.

Do **not** submit M9-only full/scale before this panel contract is implemented in M12A/B/C.

## Model / data overlap audit

| Model | Public training-data fact | M12 panel implication | Claim status |
|---|---|---|---|
| `InstaDeepAI/nucleotide-transformer-v2-500m-multi-species` | HF model card: pretrained on 850 NCBI genomes; plants and viruses are not included. | `{arabidopsis,rice}` is admissible as a clean plant panel for NT-v2/M9. | claim-clean for plant panel, pending same-panel baselines and published SOTA freeze |
| `GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview` | HF model card: token-classification CDS annotator; repo/dossier says split/schema are not fully public; eukaryote preview is developmental. | Same plant panel is useful for fair challenger, but overlap cannot be excluded from public data alone. | mechanism/control only |
| `GenerTeam/GENERanno-eukaryote-0.5b-base` | HF model card: 500M fill-mask/base model, 8k bp context, trained on 386B bp eukaryotic DNA; developmental, technical details pending. | Include as M12C-base to test whether an unadjusted GENERanno base model can learn our 3-class intron-aware task. | mechanism/control only |
| `Tiberius` | Current configs cover angiosperms; published mammal SOTA is not the plant claim anchor. | Use as same-panel external tool if angiosperm weights run on the plant panel. | comparator, not published SOTA anchor by itself |
| `Helixer` | Plant model exists in current model set but local cache currently only has fungi/invertebrate weights. | Need plant weight download/pin before M12B. | comparator pending weight pin |
| `ANNEVO` | Repo includes `ANNEVO_Magnoliopsida.pt`; repo notes say listed model species were intentionally excluded from ANNEVO train/validation, but full supplementary table still must be frozen before claim. | Best near-term published-anchor comparator on `{arabidopsis,rice}`. | primary comparator; claim anchor still pending M2 freeze |

## Frozen M12 v0.1 species panel

| Species | Accession | Clade | Role | Existing path | Clean for NT-v2? | Clean for GENERanno? |
|---|---|---|---|---|---|---|
| Arabidopsis thaliana | GCF_000001735.4 | Magnoliopsida | train/calibrate candidate; same-panel baseline test subset | `data/m1_screen/arabidopsis_thaliana` | yes | unknown |
| Oryza sativa Japonica Group | GCF_034140825.1 | Magnoliopsida | unseen-species fixed-model test; same-panel baseline test subset | `data/m1_screen/oryza_sativa` | yes | unknown |

Rationale:
- Both species already have genome/reference checksums in `configs/m1_data_manifest.yaml` and download reports.
- Both are plant species, which avoids known NT-v2 850-genome multi-species pretraining overlap.
- Both were already used in M8-M11 clean plant evidence, so data loading and evaluation are de-risked.

Limitations:
- This is a plant-panel first pass, not a broad-eukaryote final claim.
- A later M12+/generalization panel should add a non-plant clean species only after overlap is audited for the chosen backbone(s).

## M12A fixed-model protocol

Primary pilot:
- Train and calibrate M9-L12 on `Arabidopsis thaliana` train/val only.
- Freeze checkpoint and validation-selected decode/calibration parameters.
- Evaluate once on `Oryza sativa` test seqids without rice labels in training or calibration.

Optional symmetric check:
- Train/calibrate on rice and evaluate Arabidopsis if the primary pilot is ambiguous or reviewer pressure requires symmetry.

Completion gate:
- Report aggregate and per-species intergenic specificity/FPR, macro specificity, gbF1, constrained gbF1, and gene_count ratio.
- This remains non-claim until published SOTA benchmark and same-panel baseline contract are frozen.

## M12B same-panel baseline protocol

Run/evaluate on the same frozen test subsets:
- M9-L12 fixed model from M12A.
- ANNEVO `Magnoliopsida` model.
- Tiberius angiosperm model/config if the public weights resolve.
- Helixer land_plant model after downloading/pinning the current weight record.

All outputs must be evaluated with:
- CDS-span gene-body F1 for cross-tool common comparability.
- Full-transcript-span complement for intergenic specificity where reference supports it.
- FPR sensitivity at `0.005/0.01/0.02`.
- Gene-count ratio and runtime/resource notes.

## M12C GENERanno fair challenger protocol

Two backbones:

1. `GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview`
   - Loader: `AutoModelForTokenClassification`; use `full_model.model`.
   - Existing M10 trainer mostly applies.
   - Purpose: tests whether the official CDS-adjusted model can become coherent with our 3-class head.

2. `GenerTeam/GENERanno-eukaryote-0.5b-base`
   - Loader: `AutoModelForMaskedLM`; use `full_model.model`.
   - Needs a small trainer change because existing code currently hardcodes `AutoModelForTokenClassification`.
   - Purpose: tests whether an unadjusted GENERanno base model can learn the task, separating base pretraining from official CDS adjustment.

Shared design:
- Same plant panel, same 3-class head, same FP-aware loss, same constrained post-processing, same evaluator.
- Start with bounded smoke/screen; do not equal-priority full/scale this route.

Stop criteria:
- Stop M12C if FPR remains `>0.02` or gene_count ratio remains `>1.25` after a fair bounded screen.
- Stop M12C if intron/gene-body-nc class remains effectively unlearned on validation.
- If 0.5b-base performs close to or better than 1.2b CDS preview, the mechanism story shifts toward "GENERanno base representations transfer well"; if only M9 remains strong, the story supports NT-v2/intron-aware calibration specificity.

## Immediate implementation blockers

- Add a base-loader branch to `src/foundation_probe/train_generanno_lora_3class.py` or create a narrowly scoped M12C trainer.
- Pin/download Helixer plant weights before M12B.
- Confirm Tiberius angiosperm config/weights run on the plant panel.
- Expand ANNEVO supplementary/model-species overlap audit before any full published-SOTA claim.
- Reconcile stale `FP-SEGMENTNT-FEATCACHE-M7` tracker row; do not rely on it for M12.

## Sources

- HF `InstaDeepAI/nucleotide-transformer-v2-500m-multi-species` model card, read 2026-06-17.
- HF `GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview` model card, read 2026-06-17.
- HF `GenerTeam/GENERanno-eukaryote-0.5b-base` model card and config, read 2026-06-17.
- HF `GenerTeam/pretrain_data_eukaryote` and `GenerTeam/cds-annotation` dataset cards, read 2026-06-17.
- Local dossiers: `refs/dossiers/{generanno-2025,segmentnt-2025,annevo-2026,helixer-2025,tiberius-2024}.md`.
- Local manifest: `configs/m1_data_manifest.yaml`.
