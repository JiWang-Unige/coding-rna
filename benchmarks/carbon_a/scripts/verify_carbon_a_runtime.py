"""Check the isolated Carbon-A runtime; no DNA, model load, or GPU inference."""
import importlib
import importlib.metadata
import json
import os
import platform
import sys
from pathlib import Path

runtime, output = (Path(arg).resolve() for arg in sys.argv[1:])
assert Path(sys.prefix).resolve() == runtime, (sys.prefix, runtime)
print("[runtime] importing torch", flush=True)
import torch
print("[runtime] importing transformers", flush=True)
import transformers
print("[runtime] importing flash_attn", flush=True)
import flash_attn
print("[runtime] core imports complete", flush=True)

assert torch.__version__ == "2.7.1+cu126", torch.__version__
assert transformers.__version__ == "4.56.2", transformers.__version__
assert flash_attn.__version__ == "2.8.3", flash_attn.__version__
packages = {}
for name in ("torch", "transformers", "flash_attn", "huggingface_hub", "safetensors",
             "numpy", "Bio", "pandas", "pyarrow", "tqdm", "numba", "ninja", "einops"):
    print(f"[runtime] checking {name}", flush=True)
    module = importlib.import_module(name)
    path = Path(module.__file__).resolve()
    assert path.is_relative_to(runtime), (name, path)
    distribution = {"Bio": "biopython", "flash_attn": "flash-attn"}.get(name, name)
    packages[name] = {"version": importlib.metadata.version(distribution), "path": str(path)}
summary = {
    "runtime": str(runtime), "python": platform.python_version(),
    "base_prefix": sys.base_prefix, "packages": packages,
    "installation_job_id": os.environ.get("CARBON_A_INSTALL_JOB", os.environ["SLURM_JOB_ID"]),
    "cuda_build": torch.version.cuda, "cxx11_abi": torch._C._GLIBCXX_USE_CXX11_ABI,
    "isolated_supported_imports_passed": True,
    "gpu_forward_executed": False, "raw_annotation_executed": False,
    "research_data_read": False, "training_updates": 0,
}
(output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
