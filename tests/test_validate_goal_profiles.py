import json
import subprocess
import sys


def run_validate(goal, metrics, profile, status_file):
    proc = subprocess.run(
        [
            sys.executable,
            "scripts/validate_goal.py",
            "--goal",
            str(goal),
            "--metrics",
            str(metrics),
            "--profile",
            profile,
            "--run-status",
            str(status_file),
        ],
        check=False,
        text=True,
        capture_output=True,
    )
    return proc.returncode, json.loads(proc.stdout)


def test_profile_specific_intergenic_fpr_threshold(tmp_path):
    goal = tmp_path / "goal.json"
    metrics = tmp_path / "metrics.json"
    status = tmp_path / "STATUS"

    goal.write_text(
        json.dumps(
            {
                "primary_metric": "constrained_gene_body_F1",
                "direction": "higher",
                "success_criteria": [
                    {"metric": "constrained_gene_body_F1", "op": ">", "threshold": 0.0}
                ],
                "guardrails": [
                    {
                        "metric": "intergenic_FPR",
                        "op": "<=",
                        "threshold": 0.01,
                        "threshold_by_profile": {
                            "smoke": 0.02,
                            "screen": 0.02,
                            "full": 0.01,
                            "scale": 0.01,
                        },
                    }
                ],
                "screen_anchor": {
                    "metric": "constrained_gene_body_F1",
                    "value": 0.0,
                    "direction": "higher",
                },
                "sota_benchmark": {
                    "metric": "constrained_gene_body_F1",
                    "value": 0.0,
                    "direction": "higher",
                },
                "status": "draft",
            }
        ),
        encoding="utf-8",
    )
    metrics.write_text(
        json.dumps(
            {
                "primary_metric": "constrained_gene_body_F1",
                "constrained_gene_body_F1": 0.5,
                "intergenic_FPR": 0.015,
            }
        ),
        encoding="utf-8",
    )
    status.write_text("COMPLETED\n", encoding="utf-8")

    smoke_code, smoke = run_validate(goal, metrics, "smoke", status)
    full_code, full = run_validate(goal, metrics, "full", status)

    assert smoke_code == 1
    assert smoke["status"] == "progress"
    assert smoke["guardrails_gate"]["ok"] is True
    assert smoke["guardrails_gate"]["rules"][0]["threshold"] == 0.02

    assert full_code == 1
    assert full["status"] == "not_yet"
    assert full["guardrails_gate"]["ok"] is False
    assert full["guardrails_gate"]["rules"][0]["threshold"] == 0.01


def _audit_goal(tmp_path):
    """Goal mirroring ACTIVE_GOAL after M1-AGGREGATION-GATE-AUDIT."""
    goal = tmp_path / "goal.json"
    goal.write_text(
        json.dumps(
            {
                "primary_metric": "constrained_gene_body_F1",
                "direction": "higher",
                "success_criteria": [
                    {"metric": "constrained_gene_body_F1", "op": ">", "threshold": 0.0}
                ],
                "guardrails": [
                    {"metric": "intergenic_FPR", "op": "<=", "threshold": 0.01,
                     "threshold_by_profile": {"smoke": 0.02, "screen": 0.02,
                                              "full": 0.01, "scale": 0.01}},
                    {"metric": "nucleotide_gene_body_F1_drop_vs_anchor", "op": "<=",
                     "threshold": 0.03, "profiles": ["screen", "full", "scale"]},
                    {"metric": "predicted_gene_count_ratio_vs_reference", "op": "<=", "threshold": 1.25},
                ],
                "semantic_success": {
                    "degenerate_bound_check": True,
                    "degenerate_exempt_if_all": [
                        {"metric": "gene_body_F1_unconstrained", "op": ">=", "threshold": 0.05},
                    ],
                },
                "screen_anchor": {"metric": "constrained_gene_body_F1", "value": 0.0,
                                  "direction": "higher"},
                "sota_benchmark": {"metric": "constrained_gene_body_F1", "value": 0.0,
                                   "direction": "higher"},
                "status": "draft",
            }
        ),
        encoding="utf-8",
    )
    return goal


def test_completed_poor_not_failed_run(tmp_path):
    """A finite baseline whose constrained primary is guardrail-hard-zeroed must be
    completed_poor (keep iterating, exit 1), NOT a failed_run (stop+notify, exit 3)."""
    goal = _audit_goal(tmp_path)
    metrics = tmp_path / "metrics.json"
    status = tmp_path / "STATUS"
    metrics.write_text(
        json.dumps(
            {
                "primary_metric": "constrained_gene_body_F1",
                "constrained_gene_body_F1": 0.0,
                "gene_body_F1_unconstrained": 0.7087,
                "intergenic_FPR": 0.0287,
                "predicted_gene_count_ratio_vs_reference": 0.72,
                "semantic_success": True,
            }
        ),
        encoding="utf-8",
    )
    status.write_text("COMPLETED\n", encoding="utf-8")

    code, out = run_validate(goal, metrics, "screen", status)
    assert code == 1
    assert out["status"] == "not_yet"
    assert out["disposition"] == "completed_poor"
    assert out["semantic_ok"] is True


def test_below_floor_signal_is_failed_run(tmp_path):
    """Tightened gate (M1-CONTRACT-REVIEW): a 0.0 primary whose underlying signal is BELOW the
    non-trivial floor (unconstrained 0.04 < 0.05) is degenerate -> failed_run, even though the
    evaluator set semantic_success=true (a constant, no longer accepted as evidence)."""
    goal = _audit_goal(tmp_path)
    metrics = tmp_path / "metrics.json"
    status = tmp_path / "STATUS"
    metrics.write_text(
        json.dumps(
            {
                "primary_metric": "constrained_gene_body_F1",
                "constrained_gene_body_F1": 0.0,
                "gene_body_F1_unconstrained": 0.04,
                "intergenic_FPR": 0.3,
                "predicted_gene_count_ratio_vs_reference": 0.72,
                "semantic_success": True,
            }
        ),
        encoding="utf-8",
    )
    status.write_text("COMPLETED\n", encoding="utf-8")

    code, out = run_validate(goal, metrics, "screen", status)
    assert code == 3
    assert out["status"] == "failed_run"
    assert out["disposition"] is None


def test_gene_explosion_is_guardrail_not_yet_not_failed_run(tmp_path):
    """Gene-count explosion (fragmentation) is a QUALITY problem, NOT an infra/degenerate failure:
    a 0.0 primary with real base signal (F1 0.2 >= floor) + inflated gene-count ratio (12x) is
    completed_poor (exempt from failed_run) but fails the gene-count GUARDRAIL -> not_yet, exit 1.
    It must NOT stop the run as failed_run."""
    goal = _audit_goal(tmp_path)
    metrics = tmp_path / "metrics.json"
    status = tmp_path / "STATUS"
    metrics.write_text(
        json.dumps(
            {
                "primary_metric": "constrained_gene_body_F1",
                "constrained_gene_body_F1": 0.0,
                "gene_body_F1_unconstrained": 0.2,
                "intergenic_FPR": 0.9,
                "predicted_gene_count_ratio_vs_reference": 12.0,
                "semantic_success": True,
            }
        ),
        encoding="utf-8",
    )
    status.write_text("COMPLETED\n", encoding="utf-8")

    code, out = run_validate(goal, metrics, "screen", status)
    assert code == 1
    assert out["status"] == "not_yet"
    assert out["disposition"] == "completed_poor"
    gc = next(r for r in out["guardrails_gate"]["rules"]
              if r["metric"] == "predicted_gene_count_ratio_vs_reference")
    assert gc["ok"] is False


def test_genuine_degenerate_still_failed_run(tmp_path):
    """A 0.0 primary with NO underlying finite signal stays a degenerate failed_run:
    the tripwire must remain intact for truly broken/degenerate output."""
    goal = _audit_goal(tmp_path)
    metrics = tmp_path / "metrics.json"
    status = tmp_path / "STATUS"
    metrics.write_text(
        json.dumps(
            {
                "primary_metric": "constrained_gene_body_F1",
                "constrained_gene_body_F1": 0.0,
                "gene_body_F1_unconstrained": 0.0,
                "intergenic_FPR": 0.5,
            }
        ),
        encoding="utf-8",
    )
    status.write_text("COMPLETED\n", encoding="utf-8")

    code, out = run_validate(goal, metrics, "screen", status)
    assert code == 3
    assert out["status"] == "failed_run"
    assert out["disposition"] is None


def test_profile_scoped_guardrail_skipped_for_smoke(tmp_path):
    """A guardrail scoped to [screen, full, scale] must be SKIPPED (not failed) for smoke."""
    goal = _audit_goal(tmp_path)
    metrics = tmp_path / "metrics.json"
    status = tmp_path / "STATUS"
    metrics.write_text(
        json.dumps(
            {
                "primary_metric": "constrained_gene_body_F1",
                "constrained_gene_body_F1": 0.5,
                "gene_body_F1_unconstrained": 0.5,
                "intergenic_FPR": 0.015,
                "predicted_gene_count_ratio_vs_reference": 1.0,
                # nucleotide_gene_body_F1_drop_vs_anchor intentionally ABSENT
            }
        ),
        encoding="utf-8",
    )
    status.write_text("COMPLETED\n", encoding="utf-8")

    code, out = run_validate(goal, metrics, "smoke", status)
    drop_rule = next(r for r in out["guardrails_gate"]["rules"]
                     if r["metric"] == "nucleotide_gene_body_F1_drop_vs_anchor")
    assert drop_rule.get("skipped") is True
    assert out["guardrails_gate"]["ok"] is True
    assert out["status"] == "progress"
