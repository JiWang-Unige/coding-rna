#!/usr/bin/env python3
"""Freeze the M20 claim-clean panel/provenance gate for GENERanno.

This is a local evidence synthesis script: it does not infer hidden training
membership from performance.  It records whether public provenance is strong
enough to support a clean held-out claim.
"""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT_MD = ROOT / "refs/dossiers/M20-CLAIM-CLEAN-PANEL-FREEZE.md"
OUT_JSON = ROOT / "reports/M20-CLAIM-CLEAN-PANEL-FREEZE/summary.json"

SOURCES = [
    {
        "name": "GENERanno 1.2B CDS annotator preview model card",
        "url": "https://huggingface.co/GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview",
        "local_use": "Confirms model identity/task; no complete eukaryotic species/accession exclusion list.",
    },
    {
        "name": "GENERanno 0.5B eukaryote base model card",
        "url": "https://huggingface.co/GenerTeam/GENERanno-eukaryote-0.5b-base",
        "local_use": "Confirms broad eukaryotic pretraining; does not certify plant exclusion.",
    },
    {
        "name": "GENERanno GitHub README",
        "url": "https://github.com/GenerTeam/GENERanno",
        "local_use": "Confirms release lineage and eukaryotic annotation focus; no full training manifest.",
    },
    {
        "name": "GENERanno eukaryote pretraining dataset card",
        "url": "https://huggingface.co/datasets/GenerTeam/pretrain_data_eukaryote",
        "local_use": "RefSeq-derived broad eukaryotic corpus family; public card is not a species exclusion manifest.",
    },
    {
        "name": "M19 local provenance audit",
        "url": "refs/dossiers/m19_generanno_provenance_audit.md",
        "local_use": "Prior project audit already marked Arabidopsis/rice overlap as unknown.",
    },
]

PANEL = [
    {
        "species": "arabidopsis_thaliana",
        "data_path": "data/m1_screen/arabidopsis_thaliana",
        "project_role": "current clean-plant train/val/test screen species",
        "generanno_overlap_status": "overlap_unknown",
        "claim_use": "adaptation/comparability only",
    },
    {
        "species": "oryza_sativa",
        "data_path": "data/m1_screen/oryza_sativa",
        "project_role": "current clean-plant train/val/test screen species",
        "generanno_overlap_status": "overlap_unknown",
        "claim_use": "adaptation/comparability only",
    },
    {
        "species": "arabidopsis_lyrata",
        "data_path": "data/m13_distance_screen/arabidopsis_lyrata",
        "project_role": "close-plant diagnostic candidate",
        "generanno_overlap_status": "overlap_unknown",
        "claim_use": "diagnostic only unless manifest exclusion is obtained",
    },
    {
        "species": "gallus_gallus",
        "data_path": "data/m1_screen/gallus_gallus",
        "project_role": "animal negative-control diagnostic",
        "generanno_overlap_status": "overlap_unknown",
        "claim_use": "diagnostic only",
    },
    {
        "species": "drosophila_melanogaster",
        "data_path": "data/m1_screen/drosophila_melanogaster",
        "project_role": "animal negative-control diagnostic",
        "generanno_overlap_status": "overlap_unknown",
        "claim_use": "diagnostic only",
    },
]


def _exists(path: str) -> bool:
    if path.startswith("refs/") or path.startswith("data/"):
        return (ROOT / path).exists()
    return False


def write_outputs() -> None:
    today = date.today().isoformat()
    status = "BLOCKED_FOR_GENERANNO_CLEAN_HELDOUT_CLAIM"
    recommendation = (
        "Do not spend claim-grade GENERanno GPU on the current Arabidopsis/rice clean-plant panel "
        "until a public or author-provided species/accession manifest excludes the claim species. "
        "Use M18/M19/M20 as adaptation/comparability evidence, or move the clean held-out claim to a "
        "backbone/training protocol with certifiable provenance."
    )
    rows = []
    for item in PANEL:
        rows.append({**item, "path_exists": _exists(item["data_path"])})

    payload = {
        "exp_id": "M20-CLAIM-CLEAN-PANEL-FREEZE",
        "date": today,
        "status": status,
        "decision": "freeze_current_generanno_claim_panel_as_not_claim_clean",
        "recommendation": recommendation,
        "panel": rows,
        "sources": SOURCES,
        "claim_paths": [
            {
                "path": "A",
                "label": "GENERanno adaptation/comparability paper evidence",
                "allowed_claim": "A CDS-specialized pretrained backbone can be adapted with an intron-aware head and FP-aware objective.",
                "forbidden_claim": "No-overlap held-out cross-species SOTA on Arabidopsis/rice.",
            },
            {
                "path": "B",
                "label": "Clean held-out claim",
                "requirement": "Use species/accessions provably absent from the backbone pretraining/fine-tuning data, or use a model/training route whose data provenance we control.",
            },
        ],
    }

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    lines = [
        "# M20-CLAIM-CLEAN-PANEL-FREEZE",
        "",
        f"- Date: {today}",
        f"- Status: `{status}`",
        "- Scope: GENERanno 1.2B CDS-preview / 0.5B base as candidate backbones for clean held-out gene-annotation claims.",
        "",
        "## Verdict",
        "",
        "The current Arabidopsis/rice GENERanno panel is **not claim-clean**. Public documentation still does not expose a complete species/accession manifest for eukaryote pretraining plus CDS-preview annotation tuning, so absence from training cannot be certified.",
        "",
        recommendation,
        "",
        "## Panel Freeze",
        "",
        "| Species | Data path | Path exists | Project role | GENERanno overlap status | Claim use |",
        "|---|---:|---:|---|---|---|",
    ]
    for item in rows:
        lines.append(
            f"| `{item['species']}` | `{item['data_path']}` | {item['path_exists']} | "
            f"{item['project_role']} | `{item['generanno_overlap_status']}` | {item['claim_use']} |"
        )
    lines.extend(
        [
            "",
            "## Claim Consequence",
            "",
            "- `allowed`: present GENERanno results as pretrained-backbone adaptation, same-panel comparability, and mechanism evidence.",
            "- `blocked`: present Arabidopsis/rice M18/M19/M20 as clean no-overlap held-out SOTA evidence.",
            "- `unlock condition`: a complete public/author-provided training and fine-tuning manifest excludes the exact claim species/accessions, or the claim is moved to a backbone/protocol with controlled provenance.",
            "",
            "## Sources Used",
            "",
        ]
    )
    for source in SOURCES:
        lines.append(f"- {source['name']}: {source['url']} — {source['local_use']}")
    lines.append("")
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(lines))
    print(f"wrote {OUT_MD.relative_to(ROOT)}")
    print(f"wrote {OUT_JSON.relative_to(ROOT)}")


if __name__ == "__main__":
    write_outputs()
