# M17-PRETRAINING-OVERLAP-AUDIT

- Date: 2026-06-18
- Scope: M13/M14/M16 diagnostic species `{arabidopsis_lyrata, oryza_sativa, gallus_gallus, drosophila_melanogaster}` and baseline/backbone caveats for interpreting `M17-SAMEPANEL-GENERALIZATION-BASELINES`.
- Status: local dossier, non-claim.

## Summary

M17 should be interpreted as a same-evaluator diagnostic comparison, not a clean held-out SOTA claim. The plant panel remains the cleanest for NT-v2/M9 based on the existing M12 audit. The animal panel is useful as a negative/control stress test, but several released external baselines explicitly include gallus/drosophila as benchmark/test species or note Helixer overlap. Therefore animal results can answer "how released tools behave on the same evaluator" but cannot by themselves establish a clean no-overlap claim.

## Species-level caveats

| Species | Role in M13/M16 | Current caveat | Claim status |
|---|---|---|---|
| `arabidopsis_lyrata` | close Arabidopsis-relative diagnostic | scaffold-level assembly in local manifest; not part of M12 clean plant claim panel | diagnostic only |
| `oryza_sativa` | far plant / M12A fixed-model test | M12 audit treats plants as clean for NT-v2; ANNEVO Magnoliopsida reused for rice for M12B continuity | best current plant comparator, still non-claim until benchmark frozen |
| `gallus_gallus` | animal negative control | Tiberius `vertebrates.yaml` lists Gallus gallus as a test species and notes Gallus was used during Helixer training/validation; animal overlap status is not clean | diagnostic only |
| `drosophila_melanogaster` | animal negative control | Tiberius `insecta.yaml` lists Drosophila melanogaster as a test species and notes Drosophila was used during Helixer training/validation; NT-v2 animal overlap not fully audited | diagnostic only |

## Model / baseline caveats

| Model / tool | Local fact checked | Implication for M17 |
|---|---|---|
| NT-v2 500M multi-species | Existing M12 audit records the HF model-card statement that plants were not included in the 850-genome collection. | Plant results are the cleanest NT-v2/M9 evidence; animal diagnostics need more overlap scrutiny. |
| GENERanno 1.2B / 0.5B | Existing M12 audit records incomplete public species/accession disclosure for eukaryote pretraining/annotation data. | GENERanno remains mechanism/challenger evidence, not clean held-out claim evidence. |
| ANNEVO | Local weights cover `Magnoliopsida`, `Aves`, `Insecta`; M12 audit notes model-species exclusion claims still need full supplementary freeze before claim. | M17 can compare released ANNEVO behavior on the panel, but SOTA benchmark freeze remains pending. |
| Tiberius | `model_cfg/vertebrates.yaml` lists Gallus as a test species; `model_cfg/insecta.yaml` lists Drosophila as a test species. | Tiberius animal results are released-tool diagnostic references, not clean unseen-species claims. |
| Helixer | Local `model_list.csv` provides land_plant/invertebrate/vertebrate weights; Tiberius configs state Gallus and Drosophila were used during Helixer training/validation. | Helixer animal results are expected to be potentially advantaged and should be marked overlap-contaminated. |

## M17 interpretation rule

- If released baselines outperform M9 on A. lyrata/gallus/drosophila, conclude the current M9 fixed-model route is practically weak on this diagnostic panel, but do not claim a clean SOTA gap from animal rows alone.
- If released baselines also degrade sharply on A. lyrata/gallus/drosophila, conclude the diagnostic panel is intrinsically hard or annotation/evaluator mismatched; follow with a cleaner held-out species panel before architecture scale-up.
- Use `oryza_sativa` and future clean plant species as the primary near-term claim-safe route until broader no-overlap species are frozen.

## Sources

- `refs/dossiers/m12_prereq_audit.md`
- `refs/repos/tiberius-2024/model_cfg/vertebrates.yaml`
- `refs/repos/tiberius-2024/model_cfg/insecta.yaml`
- `refs/weights/helixer-2025/model_list.csv`
- `configs/m1_data_manifest.yaml`
