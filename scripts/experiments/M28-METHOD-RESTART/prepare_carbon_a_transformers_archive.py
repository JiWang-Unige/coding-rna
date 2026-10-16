"""Archive unchanged Transformers software for node-local import I/O, no ML imports."""
import importlib.metadata
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUNTIME = ROOT / "refs/carbon-a/runtimes/torch271-cu126-tf4562-fa283-r1"
OUT = Path(os.environ["OUT"]).resolve()


def main():
    started = time.monotonic()
    assert Path(sys.prefix).resolve() == RUNTIME.resolve()
    assert importlib.metadata.version("transformers") == "4.56.2"
    site = RUNTIME / "lib/python3.10/site-packages"
    package = site / "transformers"
    metadata = site / "transformers-4.56.2.dist-info"
    assert package.is_dir() and metadata.is_dir()
    archive = OUT / "transformers-4.56.2.tar.gz"
    assert not archive.exists()
    subprocess.run(["tar", "--exclude=__pycache__", "-czf", str(archive),
                    "-C", str(site), "transformers", metadata.name], check=True)
    assert archive.stat().st_size > 0
    summary = dict(completed=True, job_id=os.environ["SLURM_JOB_ID"],
                   package="transformers", version="4.56.2", source_runtime=str(RUNTIME),
                   source_package=str(package), archive=str(archive),
                   archive_bytes=archive.stat().st_size,
                   software_only=True, original_packages_modified=False,
                   package_upgrade=False, ML_import_executed=False,
                   model_weights_read=False, research_data_read=False, training_updates=0,
                   elapsed_seconds=time.monotonic() - started,
                   limitation="Package preparation only; node-local import and native annotation unverified")
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
