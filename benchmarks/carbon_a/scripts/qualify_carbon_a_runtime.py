"""Recover complete import evidence from 13500426, without importing ML again."""
import importlib.metadata
import json
import os
import sys
from pathlib import Path

root = Path("/home/users/j/jwang/coding-rna")
runtime = (root / "refs/carbon-a/runtimes/torch271-cu126-tf4562-fa283-r1").resolve()
source_folder = root / "outputs/M28-CARBON-A-RUNTIME/13500426"
source = json.loads((source_folder / "summary.json").read_text())
assert Path(sys.prefix).resolve() == runtime
assert Path(source["runtime"]) == runtime
assert source["installation_job_id"] == "13500092"
assert source["isolated_supported_imports_passed"]
assert source["training_updates"] == 0
assert not any(source[key] for key in ("gpu_forward_executed", "raw_annotation_executed", "research_data_read"))
expected = ("torch", "transformers", "flash_attn", "huggingface_hub", "safetensors",
            "numpy", "Bio", "pandas", "pyarrow", "tqdm", "numba", "ninja", "einops")
assert set(source["packages"]) == set(expected)
for name in expected:
    record = source["packages"][name]
    path = Path(record["path"])
    assert path.is_file() and path.resolve().is_relative_to(runtime), (name, path)
    distribution = {"Bio": "biopython", "flash_attn": "flash-attn"}.get(name, name)
    assert importlib.metadata.version(distribution) == record["version"], name
for name, version in {"torch": "2.7.1+cu126", "transformers": "4.56.2", "flash_attn": "2.8.3"}.items():
    assert source["packages"][name]["version"] == version, name

summary = {
    **source, "qualification_job_id": os.environ["SLURM_JOB_ID"],
    "source_import_job_id": "13500426", "source_import_job_state": "TIMEOUT",
    "source_import_job_completed": False, "source_import_evidence": str(source_folder / "summary.json"),
    "source_complete_import_evidence_validated": True,
    "current_package_metadata_matches_import_evidence": True,
    "package_import_executed_in_qualification": False,
    "limitation": "recovered complete CPU import evidence; GPU/native inference still unverified",
}
output = Path(os.environ["OUT"])
(output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary), flush=True)
