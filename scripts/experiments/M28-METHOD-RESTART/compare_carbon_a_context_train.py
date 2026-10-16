"""Fixed TRAIN comparison of native Carbon-A window policies, no new inference."""
import json
import math
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "outputs/M28-CARBON-A-CONTEXT24-TRAIN-COMPARISON"
DEFAULT = ROOT / "outputs/M28-CARBON-A-TRAIN-EVAL"
SHORT = ROOT / "outputs/M28-CARBON-A-CONTEXT24-TRAIN-EVAL"
JOINT = ROOT / "outputs/M28-FULL-GRAMMAR-JOINT24-TRAIN-DIAGNOSTIC"
POINTS = ("CarbonA_confidence_ge_0.1_or_NA", "CarbonA_unfiltered")


def completed(folder):
    assert (folder / "STATUS").read_text().strip() == "COMPLETED"
    value = json.loads((folder / "summary.json").read_text())
    assert value["completed"]
    assert str(value["job_id"]) == (folder / "JOBID").read_text().strip()
    return value


def exact(value):
    result = dict(value["exact_chain"])
    expected = 2 * result["tp"] / (result["predicted"] + result["reference"])
    assert math.isclose(result["f1"], expected, rel_tol=0, abs_tol=1e-12)
    return result


def compare(short, default, joint):
    assert set(short["chromosomes"]) == set(default["chromosomes"]) == set(joint["chromosomes"])
    assert set(short["per_species"]) == set(default["per_species"]) == set(joint["per_species"])
    result = {"per_species": {}, "chromosomes": {}}
    for level in ("per_species", "chromosomes"):
        for name in sorted(short[level]):
            a, b = exact(short[level][name]), exact(default[level][name])
            c = exact(joint[level][name])
            assert a["reference"] == b["reference"] == c["reference"]
            result[level][name] = dict(context24576=a, context98304=b, joint24576=c,
                context24576_minus_98304_F1_pp=100 * (a["f1"] - b["f1"]),
                context24576_minus_joint_F1_pp=100 * (a["f1"] - c["f1"]))
    a, b, c = exact(short["pooled"]), exact(default["pooled"]), exact(joint["pooled"])
    assert a["reference"] == b["reference"] == c["reference"] == 22688
    result["pooled"] = dict(context24576=a, context98304=b, joint24576=c,
        context24576_minus_98304_F1_pp=100 * (a["f1"] - b["f1"]),
        context24576_minus_joint_F1_pp=100 * (a["f1"] - c["f1"]))
    for key in ("macro_species_exact_chain_F1", "worst_species_exact_chain_F1"):
        result[key] = dict(context24576=short[key], context98304=default[key], joint24576=joint[key],
            context24576_minus_98304_F1_pp=100 * (short[key] - default[key]),
            context24576_minus_joint_F1_pp=100 * (short[key] - joint[key]))
    return result


def main():
    default, short, joint = completed(DEFAULT), completed(SHORT), completed(JOINT)
    assert str(default["job_id"]) == "13503419" and str(joint["job_id"]) == "13453780"
    assert default.get("context_bp", 98304) == 98304
    assert short["context_bp"] == 24576 and short["context_tokens"] == 4096
    assert short["overlap_bp"] == 6144 and short["window_policy_modified"]
    for value in (default, short):
        assert value["reference_policy"] == joint["reference_policy"]
        assert value["reference_primary_chains"] == 22688
        assert value["optimizer_steps_this_run"] == 0
        assert not any(value[key] for key in ("DEV_used", "sealed_data_used", "reference_used_in_inference", "frozen_DEV_family_changed"))
    assert default["fixed_revision"] == short["fixed_revision"]
    assert default["shared_backbone_source_revision"] == short["shared_backbone_source_revision"]
    report = dict(completed=True, job_id=os.environ["SLURM_JOB_ID"],
        status="original_TRAIN_window_policy_diagnostic_not_causal_attribution_or_generalization",
        reference_policy=joint["reference_policy"], reference_primary_chains=22688,
        sources={"context98304": str(DEFAULT.relative_to(ROOT)), "context24576": str(SHORT.relative_to(ROOT)),
                 "joint24576": str(JOINT.relative_to(ROOT))},
        context_bp=[98304, 24576], absolute_overlap_bp=6144, overlap_fraction=[0.0625, 0.25],
        inference_stride_bp=[92160, 18432], optimizer_steps_this_run=0,
        DEV_used=False, sealed_data_used=False, frozen_DEV_family_changed=False,
        both_fixed_operating_points_retained=True,
        limitations="Window length and relative overlap/averaging geometry change together. Within-Carbon contrast does not isolate the entire Carbon-versus-grammar gap; TRAIN is not independent generalization or a speed benchmark.",
        operating_points={})
    for point in POINTS:
        report["operating_points"][point] = compare(short["methods"][point], default["methods"][point], joint)
    with (OUT / "summary.json").open("x") as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({"completed": True, "reference_primary_chains": 22688}), flush=True)


if __name__ == "__main__":
    main()
