"""Read six small software-library samples on gpu035; no ML or research data."""
import json
import os
import socket
import time
from pathlib import Path

root = Path("/home/users/j/jwang/coding-rna")
runtime = root / "refs/carbon-a/runtimes/torch271-cu126-tf4562-fa283-r1"
output = Path(os.environ["OUT"])
chunk = 1024 * 1024
samples = []
for name in ("libtorch_cpu.so", "libtorch_cuda.so"):
    path = runtime / "lib/python3.10/site-packages/torch/lib" / name
    size = path.stat().st_size
    with path.open("rb") as handle:
        for offset in (0, size // 2, size - chunk):
            print(f"[library-io] reading {name} offset={offset}", flush=True)
            started = time.monotonic()
            handle.seek(offset)
            data = handle.read(chunk)
            seconds = time.monotonic() - started
            if len(data) != chunk:
                raise ValueError(f"Short software-library read: {name} offset={offset}")
            row = {"library": name, "file_bytes": size, "offset": offset,
                   "bytes_read": len(data), "seconds": seconds}
            samples.append(row)
            print(json.dumps(row), flush=True)

summary = {
    "completed": True, "job_id": os.environ["SLURM_JOB_ID"], "node": socket.gethostname(),
    "runtime": str(runtime.resolve()), "samples": samples,
    "total_sampled_bytes": sum(row["bytes_read"] for row in samples),
    "software_libraries_only": True, "package_import_executed": False,
    "model_weights_read": False, "research_data_read": False, "training_updates": 0,
    "limitation": "small random-read diagnostic; not full import validation or inference speed",
}
(output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary), flush=True)
