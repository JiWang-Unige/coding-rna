"""Recover a completed artificial native stage after the missing-SciPy ruler failure."""
import importlib.metadata as metadata
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

started = time.monotonic()
ROOT = Path(__file__).resolve().parents[3]
RUNTIME = ROOT / "refs/carbon-a/runtimes/torch271-cu126-tf4562-fa283-r1"
OUT = Path(os.environ["OUT"]).resolve()
source_job = os.environ["SOURCE_NATIVE_JOB"]
assert source_job == "13502103"
assert Path(sys.prefix).resolve() == RUNTIME.resolve()
source_folder = ROOT / "outputs/M28-CARBON-A-PREFLIGHT" / source_job
accounting = subprocess.check_output(
    ["sacct", "-j", source_job, "--format=JobIDRaw,State,ExitCode", "-n", "-P"],
    text=True, env={**os.environ, "TZ": "UTC"},
)
states = {row[0]: row[1:] for line in accounting.splitlines() if (row := line.strip().split("|")) and len(row) == 3}
assert states[source_job] == ["FAILED", "1:0"]
assert states[source_job + ".batch"] == ["FAILED", "1:0"]
assert states[source_job + ".extern"] == ["COMPLETED", "0:0"]
error = (ROOT / "outputs/M28-PILOT-R5/logs" / f"M28CARBONAPRE_{source_job}.err").read_text()
assert "ModuleNotFoundError: No module named 'scipy'" in error
source = json.loads((source_folder / "summary.json").read_text())
assert source["fixed_revision"] == "f6789ee79e6ef397de333577ffec1cb734496620"
assert source["flash_attention_bf16_kernel_passed"]
assert source["effective_attention"] == "flash_attention_2" and source["viterbi_backend"] == "gpu"
assert source["coordinate_fixture"]["coordinate_cds_protein_roundtrip"]
assert source["native_output"]["coordinate_cds_protein_roundtrip"]
assert source["node_local_transformers"]["source_archive_job"] == os.environ["ARCHIVE_JOB"]
assert not source["research_data_read"] and source["training_updates"] == 0
assert "[annotation] complete:" in (source_folder / "native.log").read_text()

def package_versions():
    return {d.metadata["Name"].lower().replace("_", "-"): d.version for d in metadata.distributions()}

before = package_versions()
assert before["torch"] == "2.7.1+cu126" and before["transformers"] == "4.56.2"
print("[qualification] add only scipy1.15.3 to the isolated runtime", flush=True)
with (OUT / "pip-install.log").open("x") as handle:
    subprocess.run([sys.executable, "-m", "pip", "install", "--no-deps", "--only-binary=:all:",
                    "scipy==1.15.3"], check=True, stdout=handle, stderr=subprocess.STDOUT)
subprocess.run([sys.executable, "-m", "pip", "check"], check=True)
after = package_versions()
assert after["scipy"] == "1.15.3"
assert all(after.get(name) == version for name, version in before.items())
assert set(after) - set(before) <= {"scipy"}
scipy_source = Path(metadata.distribution("scipy").locate_file("scipy")).resolve()
assert RUNTIME.resolve() in scipy_source.parents
with (OUT / "pip-freeze.txt").open("x") as handle:
    subprocess.run([sys.executable, "-m", "pip", "freeze"], check=True, stdout=handle)
(OUT / "dependency_repair.json").write_text(json.dumps({
    "scipy_version": after["scipy"], "scipy_source": str(scipy_source),
    "existing_package_versions_preserved": True, "original_generanno_modified": False,
    "source_native_job": source_job, "source_native_job_state": "FAILED",
}, indent=2) + "\n")

# Only the already completed artificial evidence is copied. The failed source is immutable.
for name in ("artificial.fna", "known_coordinate_fixture.gff", "known_coordinate_fixture.fna",
             "known_coordinate_fixture.faa", "native.log"):
    shutil.copy2(source_folder / name, OUT / name)
native_source = Path(source["native_directory"]).resolve()
assert native_source.parent == source_folder.resolve()
assert native_source.name in ("native", "native_WARNING")
native = OUT / native_source.name
shutil.copytree(native_source, native)
summary = dict(source)
summary.update(
    completed=False, qualification_stage="running_common_ruler_after_dependency_repair",
    job_id=os.environ["SLURM_JOB_ID"], source_native_job=source_job,
    source_native_job_state="FAILED", source_native_failure_stage="common_ruler_import_missing_scipy",
    native_GPU_inference_reexecuted=False, native_GPU_evidence_source_job=source_job,
    source_native_elapsed_seconds=source["elapsed_seconds"],
    native_directory=str(native), parquet=str(native / Path(source["parquet"]).name),
    scipy_version=after["scipy"], existing_package_versions_preserved=True,
    original_generanno_modified=False,
)
summary["qualification_node_local_transformers"] = {
    "version": metadata.version("transformers"),
    "source_file": str(Path(importlib.util.find_spec("transformers").origin).resolve()),
    "package_root": os.environ["TASK_TF_STAGE"], "source_archive_job": os.environ["ARCHIVE_JOB"],
}
(OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print("[qualification] execute the actual common ruler on copied artificial fixtures", flush=True)
subprocess.run([sys.executable, str(ROOT / "scripts/experiments/M28-METHOD-RESTART/check_carbon_a_train_integration.py"),
                str(OUT)], check=True)
ruler = json.loads((OUT / "common_ruler_checks.json").read_text())
assert ruler["completed"] and ruler["artificial_only"] and not ruler["real_reference_read"]
summary.update(
    completed=True, qualification_stage="complete_native_evidence_and_common_ruler",
    common_ruler_checks_passed=True, dependency_repair_completed=True,
    elapsed_seconds=time.monotonic() - started,
    elapsed_scope="CPU dependency repair and artificial common-ruler recovery; archive extraction excluded; native GPU inference not repeated",
)
(OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps({"completed": True, "source_native_job": source_job,
                  "source_native_job_state": "FAILED", "native_GPU_inference_reexecuted": False}), flush=True)
