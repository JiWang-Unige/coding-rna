"""Test the native pipeline's four Transformers APIs from an unchanged local copy."""
import importlib.metadata
import inspect
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUNTIME = ROOT / "refs/carbon-a/runtimes/torch271-cu126-tf4562-fa283-r1"
OUT = Path(os.environ["OUT"]).resolve()
STAGE = Path(os.environ["TASK_TF_STAGE"]).resolve()


def main():
    started = time.monotonic()
    assert Path(sys.prefix).resolve() == RUNTIME.resolve()
    assert importlib.metadata.version("transformers") == "4.56.2"
    print("[local Transformers] import Torch", flush=True)
    import torch
    assert torch.__version__ == "2.7.1+cu126"
    assert RUNTIME.resolve() in Path(torch.__file__).resolve().parents
    print("[local Transformers] import four actual native pipeline APIs", flush=True)
    import transformers
    from transformers import AutoModelForTokenClassification, AutoTokenizer, PreTrainedModel, PreTrainedTokenizer
    assert transformers.__version__ == "4.56.2"
    assert STAGE in Path(transformers.__file__).resolve().parents
    classes = (AutoModelForTokenClassification, AutoTokenizer, PreTrainedModel, PreTrainedTokenizer)
    paths = {cls.__name__: str(Path(inspect.getfile(cls)).resolve()) for cls in classes}
    assert all(STAGE in Path(path).parents for path in paths.values())
    assert issubclass(PreTrainedModel, torch.nn.Module)
    summary = dict(completed=True, job_id=os.environ["SLURM_JOB_ID"],
                   source_archive_job=os.environ["ARCHIVE_JOB"],
                   source_runtime=str(RUNTIME), node_local_package_root=str(STAGE),
                   package_version="4.56.2", torch_version=torch.__version__,
                   native_four_API_imports_passed=True, class_module_paths=paths,
                   torch_remains_in_qualified_runtime=True,
                   original_packages_modified=False, package_upgrade=False,
                   model_weights_read=False, gpu_forward_executed=False,
                   research_data_read=False, training_updates=0,
                   elapsed_seconds=time.monotonic() - started,
                   limitation="Actual API import only; full-context GPU native annotation still required")
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
