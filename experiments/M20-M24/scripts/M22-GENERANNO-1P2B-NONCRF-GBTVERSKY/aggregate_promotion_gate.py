#!/usr/bin/env python3
"""Aggregate per-seed M22 promotion gates after the Slurm array completes."""

import argparse
import json
from pathlib import Path


def load_seed(root: Path, exp_prefix: str, seed: int):
    exp_id = f"{exp_prefix}-s{seed}"
    gate_path = root / "outputs" / exp_id / "metrics" / "m22_promotion_gate.json"
    metrics_path = root / "outputs" / exp_id / "metrics" / "metrics.json"
    status_path = root / "outputs" / exp_id / "STATUS"
    row = {
        "seed": seed,
        "exp_id": exp_id,
        "status_path": str(status_path),
        "gate_path": str(gate_path),
        "metrics_path": str(metrics_path),
        "status": "MISSING",
        "gate": None,
        "gate_error": None,
        "metrics_available": metrics_path.exists(),
    }
    if status_path.exists():
        row["status"] = status_path.read_text().strip()
    if gate_path.exists():
        try:
            with gate_path.open() as handle:
                row["gate"] = json.load(handle)
        except (OSError, json.JSONDecodeError) as exc:
            row["gate_error"] = str(exc)
    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".")
    parser.add_argument("--exp-prefix", default="M22-GENERANNO-1P2B-NONCRF-GBTVERSKY")
    parser.add_argument("--seeds", nargs="+", type=int, default=[0, 1])
    parser.add_argument("--output-json", required=True)
    args = parser.parse_args()

    root = Path(args.root).resolve()
    seeds = [load_seed(root, args.exp_prefix, seed) for seed in args.seeds]
    completed = [row["status"] == "COMPLETED" for row in seeds]
    gates_present = [row["gate"] is not None for row in seeds]
    required_checks = ("hard_fpr_le_0p01", "gbf1_gt_m19_s1", "gene_count_le_1p25")
    promote_flags = [
        bool(
            row["gate"]
            and all(bool(row["gate"].get("checks", {}).get(name)) for name in required_checks)
        )
        for row in seeds
    ]
    hard_fpr_flags = [
        bool(row["gate"] and row["gate"].get("checks", {}).get("hard_fpr_le_0p01"))
        for row in seeds
    ]
    summary = {
        "exp_prefix": args.exp_prefix,
        "seeds": seeds,
        "all_seeds_completed": all(completed),
        "all_seed_gates_present": all(gates_present),
        "any_seed_promote": any(promote_flags),
        "all_seed_fpr_hard_valid": all(hard_fpr_flags) if hard_fpr_flags else False,
        "decision_rule": (
            "continue M22 objective only if any_seed_promote is true; "
            "otherwise abandon/treat as negative without tuning this objective"
        ),
        "continue_route": all(completed) and all(gates_present) and any(promote_flags),
    }
    out_path = Path(args.output_json)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
