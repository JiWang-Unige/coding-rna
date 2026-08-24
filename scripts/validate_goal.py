#!/usr/bin/env python3
"""Deterministic goal-validation gate for supervised-autonomy /pursue.

Non-skippable tripwire that fixes the autoloop failure mode. The agent CANNOT
self-judge — this script decides mechanically from metrics + ACTIVE_GOAL.json.

Enforces three disciplines as HARD, scriptable rules (not prose):
  1. Two-tier comparability (--profile): screen is judged ONLY against the
     same-budget `screen_anchor` and can NEVER claim vs published SOTA; only
     full/scale are judged against `sota_benchmark`. Kills the unfair
     "small-sample-vs-large-sample-SOTA" comparison.
  2. Anti-marginal-tuning: when gap_to_target >= tuning_gap_threshold (default
     0.05), `tuning_allowed=false` — /pursue & /pivot MUST pick an architecture
     move, not parameter tuning.
  3. Run/semantic-success: a failed or degenerate run can never pass.

Status (JSON on stdout, exit mirrors severity):
  failed_run (3): run/semantic failure → /pursue STOPS and notifies.
  not_yet (1)   : ran fine, success_criteria not met → continue.
  progress (1)  : progress gate met, claim gate not (or screen profile) → continue.
  success (0)   : criteria + guardrails met AND claim gate ok (full/scale only).
Plus an optional `stale_benchmark` warning when --challenger-sota beats the
current sota_benchmark (a newer SOTA exists → run /revise-goal).

Usage:
  python3 scripts/validate_goal.py --goal ACTIVE_GOAL.json --metrics <m.json> \
      [--profile smoke|screen|full|scale] [--run-status <file>] [--challenger-sota <value>]
"""
import argparse, json, math, sys

OPS = {">=": lambda a, b: a >= b, "<=": lambda a, b: a <= b,
       ">": lambda a, b: a > b, "<": lambda a, b: a < b, "==": lambda a, b: a == b}


def load_json(path):
    try:
        with open(path) as f:
            return json.load(f), None
    except FileNotFoundError:
        return None, f"file not found: {path}"
    except json.JSONDecodeError as e:
        return None, f"unparseable JSON in {path}: {e}"


def check(metrics, rules):
    results, all_ok = [], True
    for r in rules:
        m, op, thr = r["metric"], r.get("op", ">="), r["threshold"]
        val = metrics.get(m)
        if val is None:
            ok, note = False, "metric absent"
        elif not isinstance(val, (int, float)) or not math.isfinite(val):
            ok, note = False, f"non-finite value {val!r}"
        else:
            ok, note = OPS[op](val, thr), f"{val} {op} {thr}"
        results.append({"metric": m, "op": op, "threshold": thr, "value": val, "ok": ok, "note": note})
        all_ok = all_ok and ok
    return all_ok, results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--goal", required=True)
    ap.add_argument("--metrics", required=True)
    ap.add_argument("--profile", default="screen", choices=["smoke", "screen", "full", "scale"])
    ap.add_argument("--run-status", default=None)
    ap.add_argument("--challenger-sota", type=float, default=None,
                    help="a newly-found verified SOTA value; if it beats sota_benchmark -> stale warning")
    ap.add_argument("--prior-screen", type=float, default=None,
                    help="this candidate's own Track-A screen value; full/scale below it -> regression warning (G8)")
    args = ap.parse_args()

    out = {"status": None, "profile": args.profile, "run_ok": False, "semantic_ok": False,
           "comparison_anchor": None, "primary_progress_gate": {}, "guardrails_gate": {},
           "claim_gate": {}, "tuning_allowed": None, "gap_to_target": None,
           "recommended_axis": None, "observed_metrics": {}, "warnings": [], "failures": []}

    goal, gerr = load_json(args.goal)
    if gerr:
        out["status"] = "failed_run"; out["failures"].append(f"goal contract: {gerr}")
        print(json.dumps(out, indent=2)); sys.exit(3)

    # --- Gate 0a: run status ---
    if args.run_status:
        try:
            st = open(args.run_status).read().strip().upper()
        except OSError:
            st = "MISSING"
        if any(k in st for k in ("FAIL", "CANCEL", "TIMEOUT", "OOM", "ERROR", "NODE_FAIL", "MISSING", "STALE", "UNKNOWN")):
            out["status"] = "failed_run"; out["failures"].append(f"run status not OK: {st}")
            print(json.dumps(out, indent=2)); sys.exit(3)
    out["run_ok"] = True

    # --- Gate 0b: semantic success ---
    metrics, merr = load_json(args.metrics)
    if merr:
        out["status"] = "failed_run"; out["failures"].append(f"semantic-success: {merr}")
        print(json.dumps(out, indent=2)); sys.exit(3)
    pm = goal.get("primary_metric")
    pv = metrics.get(pm) if pm else None
    direction = goal.get("direction", "higher")
    if pm and (pv is None or not isinstance(pv, (int, float)) or not math.isfinite(pv)):
        out["status"] = "failed_run"; out["failures"].append(f"semantic-success: primary_metric '{pm}' missing/non-finite ({pv!r})")
        print(json.dumps(out, indent=2)); sys.exit(3)
    # Degenerate/leakage heuristic is only meaningful for higher-is-better, [0,1]-scaled
    # metrics (0.0 = all-wrong, 1.0 = suspiciously perfect / label leakage). For
    # lower-is-better metrics (e.g. loss), 0.0 is a legitimate optimum, so we never
    # auto-fail it here (m1 fix). Range-aware refinement left to project goal.
    if pm and direction == "higher" and pv in (0.0, 1.0):
        out["status"] = "failed_run"; out["failures"].append(f"semantic-success: '{pm}'=={pv} (degenerate/leakage; higher-is-better metric at bound)")
        print(json.dumps(out, indent=2)); sys.exit(3)
    # G9: range-aware "suspiciously good" advisory (does NOT hard-fail — surfaced for
    # pivot/tri-review). Covers leakage that lands just below 1.0 (e.g. AUPRC=0.985) and
    # out-of-range metrics (MCC/Pearson can be negative; the 0/1 degenerate check misses these).
    # Opt-in via goal contract: "sane_upper": <value> and/or "sane_range": [lo, hi].
    if pm and isinstance(pv, (int, float)):
        su = goal.get("sane_upper"); sr = goal.get("sane_range")
        if isinstance(su, (int, float)) and pv > su:
            out["suspicious_high"] = True
            out["warnings"].append(f"suspicious_high: '{pm}'={pv} exceeds sane_upper {su} — possible leakage/eval bug; verify before claim (advisory).")
        if isinstance(sr, (list, tuple)) and len(sr) == 2 and all(isinstance(x, (int, float)) for x in sr):
            if pv < sr[0] or pv > sr[1]:
                out["suspicious_high"] = True
                out["warnings"].append(f"out_of_sane_range: '{pm}'={pv} outside {list(sr)} — likely metric/eval bug; verify (advisory).")

    out["semantic_ok"] = True
    out["observed_metrics"] = {k: v for k, v in metrics.items() if isinstance(v, (int, float))}

    # --- Gate 1: success_criteria + guardrails ---
    crit_ok, crit_res = check(metrics, goal.get("success_criteria", []))
    out["primary_progress_gate"] = {"ok": crit_ok, "rules": crit_res}
    g_ok, g_res = check(metrics, goal.get("guardrails", []))
    out["guardrails_gate"] = {"ok": g_ok, "rules": g_res}

    # --- Two-tier comparison anchor + claim gate (A3/A1 fix) ---
    sota = goal.get("sota_benchmark") or {}
    anchor = goal.get("screen_anchor") or {}
    is_screen = args.profile in ("smoke", "screen")
    # screen/smoke -> ONLY screen_anchor (same-budget, fair). It must NEVER fall back to
    # published SOTA. full/scale -> ONLY sota_benchmark. A missing anchor yields 'none'
    # (claim is withheld, not silently judged against the wrong reference).
    if is_screen:
        cmp = anchor if anchor else None
    else:
        cmp = sota if sota else None
    out["comparison_anchor"] = ("screen_anchor" if (is_screen and anchor) else
                                ("sota_benchmark" if (not is_screen and sota) else "none"))
    if is_screen and not anchor:
        out["warnings"].append("screen/smoke profile but no screen_anchor in goal — comparison_anchor='none'; NOT falling back to published SOTA. Build a screen_anchor first (benchmark-roadmap M1).")
    claim_ok = None
    if cmp:
        m = cmp.get("metric", pm); bench = cmp.get("value"); d = cmp.get("direction", direction)
        val = metrics.get(m)
        if isinstance(val, (int, float)) and isinstance(bench, (int, float)):
            claim_ok = (val > bench) if d == "higher" else (val < bench)
            out["claim_gate"] = {"ok": claim_ok, "metric": m, "observed": val, "anchor_value": bench,
                                 "anchor": out["comparison_anchor"], "note": f"{val} {'>' if d=='higher' else '<'} {bench} (strict)"}

    # --- Anti-marginal-tuning gap ALWAYS references published sota_benchmark when present ---
    # (A3 fix) so a large gap to the REAL SOTA forbids tuning even during screen — exactly
    # the phase where marginal tuning is most tempting. Falls back to the claim anchor only
    # when no sota_benchmark exists at all.
    anti = sota if isinstance(sota.get("value"), (int, float)) else (cmp or {})
    av = anti.get("value"); am = anti.get("metric", pm)
    aval = metrics.get(am)
    if isinstance(av, (int, float)) and isinstance(aval, (int, float)):
        out["gap_to_target"] = round(abs(av - aval), 6)
        out["gap_reference"] = "sota_benchmark" if anti is sota else out["comparison_anchor"]

    # --- Anti-marginal-tuning hard rule (B) ---
    thr = goal.get("tuning_gap_threshold", 0.05)
    if out["gap_to_target"] is not None:
        if out["gap_to_target"] >= thr:
            out["tuning_allowed"] = False
            out["recommended_axis"] = "architecture (gap large — parameter tuning FORBIDDEN; pick head/backbone/objective/decoder/data_view axis)"
        else:
            out["tuning_allowed"] = True
            out["recommended_axis"] = "near target — systematic tuning/scaling now reasonable"

    # --- Staleness (A) ---
    if args.challenger_sota is not None and sota.get("value") is not None:
        sd = sota.get("direction", direction); sv = sota["value"]
        beats = (args.challenger_sota > sv) if sd == "higher" else (args.challenger_sota < sv)
        if beats:
            out["warnings"].append(f"stale_benchmark: a verified SOTA ({args.challenger_sota}) beats sota_benchmark ({sv}) — run /revise-goal before claiming.")

    # --- G8: Track-B scale regression vs the candidate's OWN screen value ---
    # validate normally only compares to fixed anchors, so a scale run that drops BELOW
    # its own Track-A screen (architecture not scaling) is invisible. Surface it (advisory,
    # does not change 4-state): /pivot must then consider backbone/abandon, not tuning.
    if args.prior_screen is not None and not is_screen and isinstance(pv, (int, float)):
        regressed = (pv < args.prior_screen) if direction == "higher" else (pv > args.prior_screen)
        if regressed:
            out["regression"] = True
            out["warnings"].append(
                f"regression: {args.profile} '{pm}'={pv} is worse than this candidate's own screen "
                f"({args.prior_screen}) — architecture may NOT scale. /pivot should weigh "
                f"backbone-change/abandon over tuning, not just continue scaling.")

    # --- Pre-decision guards: a draft/placeholder or criteria-less contract can NEVER pass (m5/m2) ---
    draft = str(goal.get("status", "")).lower() == "draft"
    if draft:
        out["warnings"].append("ACTIVE_GOAL.status=='draft': placeholder contract — success is DISABLED until you fill real thresholds (sota_benchmark/screen_anchor/success_criteria).")
    no_criteria = not goal.get("success_criteria")
    if no_criteria:
        out["warnings"].append("no success_criteria defined — progress gate is vacuous; success withheld until at least one criterion exists.")

    # --- Decide ---
    if not crit_ok:
        out["status"] = "not_yet"; out["failures"] = [r["note"] for r in crit_res if not r["ok"]]
        print(json.dumps(out, indent=2)); sys.exit(1)
    if not g_ok:
        out["status"] = "not_yet"; out["failures"] = ["guardrail: " + r["note"] for r in g_res if not r["ok"]]
        print(json.dumps(out, indent=2)); sys.exit(1)
    # criteria + guardrails met:
    if is_screen:
        # HARD: screen/smoke can never claim vs published SOTA
        out["status"] = "progress"
        out["warnings"].append("screen/smoke profile: cannot claim SOTA; passed screen_anchor progress only. Promote to full/scale to test the published-SOTA claim.")
        print(json.dumps(out, indent=2)); sys.exit(1)
    # full/scale: success requires a FILLED contract AND a usable SOTA anchor that is strictly beaten.
    if draft or no_criteria:
        out["status"] = "progress"
        print(json.dumps(out, indent=2)); sys.exit(1)
    if out["comparison_anchor"] == "none" or claim_ok is None:
        # A1 fix: never treat "no comparison possible" as success. Missing/invalid
        # sota_benchmark means we cannot verify strict SOTA exceedance — withhold success.
        out["status"] = "progress"
        out["warnings"].append("full/scale but no usable sota_benchmark anchor — cannot verify strict SOTA exceedance; success WITHHELD (fill sota_benchmark or run /revise-goal).")
        print(json.dumps(out, indent=2)); sys.exit(1)
    if claim_ok is True:
        out["status"] = "success"
        print(json.dumps(out, indent=2)); sys.exit(0)
    out["status"] = "progress"
    print(json.dumps(out, indent=2)); sys.exit(1)


if __name__ == "__main__":
    main()
