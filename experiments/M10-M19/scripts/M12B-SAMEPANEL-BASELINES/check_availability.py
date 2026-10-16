#!/usr/bin/env python3
"""Check local readiness for M12B same-panel baseline runners.

This is intentionally a no-inference preflight: it records which external
baselines are runnable locally before we spend GPU or download large weights.
"""
from __future__ import annotations

import json
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[3]
CONFIG = ROOT / "configs" / "M12B-SAMEPANEL-BASELINES.yaml"


def exists(path: str) -> bool:
    return (ROOT / path).exists()


def nonempty(path: str) -> bool:
    p = ROOT / path
    return p.is_file() and p.stat().st_size > 0


def species_ready(species_paths: dict[str, str]) -> dict[str, dict]:
    out = {}
    for name, rel in species_paths.items():
        out[name] = {
            "path": rel,
            "genome_fasta": nonempty(f"{rel}/genome.fa"),
            "reference_gff3": nonempty(f"{rel}/reference.gff3"),
            "prep_report_present": any((ROOT / rel).glob("*report*.json")),
        }
        out[name]["ready"] = out[name]["genome_fasta"] and out[name]["reference_gff3"]
    return out


def main() -> int:
    cfg = yaml.safe_load(CONFIG.read_text())
    baselines = cfg["baselines"]
    checks = {
        "config": str(CONFIG.relative_to(ROOT)),
        "eval_species": cfg["eval_species"],
        "species": species_ready(cfg["species_paths"]),
        "baselines": {},
    }

    checks["baselines"]["annevo_magnoliopsida"] = {
        "weight": baselines["annevo_magnoliopsida"]["weight"],
        "weight_present": nonempty(baselines["annevo_magnoliopsida"]["weight"]),
        "repo_present": exists(baselines["annevo_magnoliopsida"]["repo"]),
    }
    checks["baselines"]["tiberius_angiosperms"] = {
        "sif": baselines["tiberius_angiosperms"]["sif"],
        "sif_present": nonempty(baselines["tiberius_angiosperms"]["sif"]),
        "model_cfg": baselines["tiberius_angiosperms"]["model_cfg"],
        "model_cfg_present": nonempty(baselines["tiberius_angiosperms"]["model_cfg"]),
        "status_note": "weights_url in model_cfg may download/cache on first run; record artifact before claim",
    }
    checks["baselines"]["helixer_land_plant"] = {
        "sif": baselines["helixer_land_plant"]["sif"],
        "sif_present": nonempty(baselines["helixer_land_plant"]["sif"]),
        "preferred_weight": baselines["helixer_land_plant"]["preferred_weight"],
        "preferred_weight_present": nonempty(baselines["helixer_land_plant"]["preferred_weight"]),
        "download_url": baselines["helixer_land_plant"]["download_url"],
    }
    checks["baselines"]["m9_l12_fixed"] = {
        "dependency": "M12A-FIXEDMODEL-CROSSSPECIES outputs",
        "seed_outputs_present": {
            f"s{seed}": exists(f"outputs/M12A-FIXEDMODEL-CROSSSPECIES-A2R-s{seed}/metrics/metrics.json")
            for seed in (0, 1, 2)
        },
    }

    for name, row in checks["baselines"].items():
        if name == "annevo_magnoliopsida":
            row["ready"] = row["weight_present"] and row["repo_present"]
        elif name == "tiberius_angiosperms":
            row["ready"] = row["sif_present"] and row["model_cfg_present"]
        elif name == "helixer_land_plant":
            row["ready"] = row["sif_present"] and row["preferred_weight_present"]
        elif name == "m9_l12_fixed":
            row["ready"] = all(row["seed_outputs_present"].values())

    checks["ready_baselines"] = [k for k, v in checks["baselines"].items() if v["ready"]]
    checks["blocked_baselines"] = [k for k, v in checks["baselines"].items() if not v["ready"]]
    checks["all_species_ready"] = all(v["ready"] for v in checks["species"].values())
    checks["overall_ready"] = checks["all_species_ready"] and not checks["blocked_baselines"]

    out_path = ROOT / cfg["outputs"]["availability"]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(checks, indent=2, sort_keys=True) + "\n")
    print(json.dumps(checks, indent=2, sort_keys=True))
    return 0 if checks["all_species_ready"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
