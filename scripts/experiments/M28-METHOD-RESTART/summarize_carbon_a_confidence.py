"""Read-only decomposition of the two fixed Carbon-A primary-exact operating points."""
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
FOLDER = ROOT / "outputs/M28-CARBON-A-TRAIN-EVAL"
FILTERED = "CarbonA_confidence_ge_0.1_or_NA"
UNFILTERED = "CarbonA_unfiltered"


def explain(unfiltered, filtered):
    a, b = unfiltered["exact_chain"], filtered["exact_chain"]
    tu, pu, ref = a["tp"], a["predicted"], a["reference"]
    tf, pf = b["tp"], b["predicted"]
    assert b["reference"] == ref and ref > 0
    assert 0 <= tf <= tu <= min(pu, ref) and tf <= pf <= pu
    removed, lost = pu - pf, tu - tf
    assert 0 <= lost <= removed
    ru, rf = unfiltered["prediction_records"], filtered["prediction_records"]
    du, df = unfiltered["duplicate_chain_records"], filtered["duplicate_chain_records"]
    assert ru - pu == du and rf - pf == df
    assert 0 <= rf <= ru and 0 <= df <= du
    record_removed, duplicate_removed = ru - rf, du - df
    assert record_removed == removed + duplicate_removed
    before, after = 2 * tu / (pu + ref), 2 * tf / (pf + ref)
    assert math.isclose(before, a["f1"], abs_tol=1e-12)
    assert math.isclose(after, b["f1"], abs_tol=1e-12)
    # Ordered algebra, not a feasible intermediate prediction set or a causal effect.
    denominator_term = 200 * tu * (1 / (pf + ref) - 1 / (pu + ref))
    hit_loss_term = -200 * lost / (pf + ref)
    change = 100 * (after - before)
    assert math.isclose(denominator_term + hit_loss_term, change, abs_tol=1e-10)
    return {
        "reference_primary_keys": ref,
        "unfiltered_unique_prediction_chains": pu, "retained_unique_prediction_chains": pf,
        "unfiltered_primary_hits": tu, "retained_primary_hits": tf,
        "prediction_records_removed": record_removed, "unique_prediction_chains_removed": removed,
        "duplicate_chain_records_removed": duplicate_removed, "primary_hits_lost": lost,
        "non_primary_exact_prediction_chains_removed": removed - lost,
        "primary_recall_loss_pp": 100 * lost / ref,
        "exact_chain_F1_change_pp": change,
        "ordered_prediction_denominator_term_pp": denominator_term,
        "ordered_primary_hit_loss_term_pp": hit_loss_term,
    }


def main():
    assert (FOLDER / "STATUS").read_text().strip() == "COMPLETED"
    source = json.loads((FOLDER / "summary.json").read_text())
    assert source["completed"] and source["reference_primary_chains"] == 22688
    assert source["status"] == "separate_CarbonA_original_TRAIN_comparison_not_DEV_or_generalization"
    assert not source["DEV_used"] and not source["sealed_data_used"]
    assert not source["frozen_DEV_family_changed"] and source["optimizer_steps_this_run"] == 0
    assert set(source["methods"]) == {UNFILTERED, FILTERED}
    u, f = source["methods"][UNFILTERED], source["methods"][FILTERED]
    expected = {"NC_003070.9", "NC_003071.7", "NC_003076.8",
                "NC_089043.1", "NC_089044.1", "NC_089046.1"}
    assert set(u["chromosomes"]) == set(f["chromosomes"]) == expected
    assert set(u["per_species"]) == set(f["per_species"]) == {"arabidopsis_thaliana", "oryza_sativa"}
    result = {
        "source_evaluation_job": source["job_id"],
        "source_summary": str(FOLDER / "summary.json"),
        "comparison": "fixed confidence0.1_or_NA minus fixed unfiltered; no threshold selection",
        "chromosomes": {key: explain(u["chromosomes"][key], f["chromosomes"][key]) for key in sorted(expected)},
        "per_species": {key: explain(u["per_species"][key], f["per_species"][key]) for key in sorted(u["per_species"])},
        "pooled": explain(u["pooled"], f["pooled"]),
        "limitations": [
            "Exact-chain F1 uses unique coordinate/strand keys, not GFF records; duplicate-record removal is reported separately.",
            "Non-primary-exact chains are a benchmark residual, not biological false genes; other eligible isoforms can be present.",
            "The ordered F1 terms reconcile endpoints algebraically; they are not independent causal effects or an intermediate model.",
            "Same original TRAIN and raw calls; no new operating point, fitting, independent generalization or decoder-only claim.",
        ],
    }
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
