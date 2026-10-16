"""Artificial-only Carbon-A GPU/raw entry and coordinate preflight."""
import importlib.metadata
import importlib.util
import json
import os
import subprocess
import sys
import time
from pathlib import Path

started = time.monotonic()
ROOT = Path(__file__).resolve().parents[3]
RUNTIME = ROOT / "refs/carbon-a/runtimes/torch271-cu126-tf4562-fa283-r1"
STAGE = Path(os.environ["TASK_TF_STAGE"]).resolve()
probe_job = os.environ["IMPORT_PROBE_JOB"]
probe = json.loads((ROOT / "outputs/M28-CARBON-A-LOCAL-TF-PROBE" / probe_job / "summary.json").read_text())
assert probe["completed"] and probe["native_four_API_imports_passed"]
assert probe["source_archive_job"] == os.environ["ARCHIVE_JOB"]
assert Path(sys.prefix).resolve() == RUNTIME.resolve()
assert importlib.metadata.version("transformers") == "4.56.2"
transformers_source = Path(importlib.util.find_spec("transformers").origin).resolve()
assert STAGE in transformers_source.parents
print("[preflight] import qualified Torch and native dependencies", flush=True)
import numpy as np
import pyarrow.parquet as pq
import torch
from flash_attn import flash_attn_func
from carbon_a_output import read_outputs, rc

ASSETS = ROOT / "refs/carbon-a/f6789ee79e6ef397de333577ffec1cb734496620"
OUT = Path(os.environ["OUT"]).resolve()
PIPELINE = ASSETS / "annotation_pipeline"
sys.path.insert(0, str(PIPELINE / "src"))
import constants
import utils




assert torch.__version__ == "2.7.1+cu126"
assert RUNTIME.resolve() in Path(torch.__file__).resolve().parents
print("[preflight] BF16 FlashAttention GPU kernel", flush=True)
assert torch.cuda.is_available() and torch.cuda.device_count() == 1
device = torch.cuda.get_device_properties(0)
assert device.major >= 8 and torch.cuda.is_bf16_supported()
q = torch.randn(1, 128, 4, 64, device="cuda", dtype=torch.bfloat16)
attention = flash_attn_func(q, q, q)
assert torch.isfinite(attention).all().item()
torch.cuda.synchronize()
del q, attention
torch.cuda.empty_cache()
driver = subprocess.check_output(["nvidia-smi", "--query-gpu=name,driver_version", "--format=csv,noheader"], text=True).strip()

# Eight full-context windows exercise batch8 memory, overlap and the 6-mer tail.
length = 740017
rng = np.random.default_rng(20261011)
bases = list("".join(rng.choice(list("ACGT"), size=length)))
cds = "ATG" + "GCC" * 40 + "TAA"
block = cds[:61] + "GT" + "C" * 96 + "AG" + cds[61:]
bases[1000:1000 + len(block)] = block
bases[6000:6000 + len(block)] = rc(block)
partial_plus, partial_minus = "C" + "GCC" * 4, "CC" + "GCC" * 4
bases[:len(partial_plus)] = partial_plus
bases[-len(partial_minus):] = rc(partial_minus)
sequence = "".join(bases)
assert len(sequence) == length
fasta = OUT / "artificial.fna"
fasta.write_text(">artificial\n" + "\n".join(sequence[i:i + 80] for i in range(0, length, 80)) + "\n")

# These are known writer/assembly fixtures, not model predictions or labels for fitting.
fixture_specs = [
    ("+", [(1000, 1060), (1161, 1225)], utils.START_A, utils.STOP, 0.25),
    ("-", [(length - 1 - 6225, length - 1 - 6165),
            (length - 1 - 6064, length - 1 - 6000)], utils.START_A, utils.STOP, 0.05),
    ("+", [(0, 12)], utils.EXON_PREFIX_TYPE_0, utils.EXON_PREFIX_TYPE_0, None),
    ("-", [(0, 13)], utils.EXON_PREFIX_TYPE_5, utils.EXON_PREFIX_TYPE_0, None),
]
fixture_genes = []
for strand, exons, first, last, confidence in fixture_specs:
    gene = utils.assemble_gene_model("artificial", sequence, exons, strand, length, first, last)
    gene["confidence"] = confidence
    fixture_genes.append(gene)
expected = {
    "gene000001": {"strand": "+", "segments": [(1000, 1061, 0, "+"), (1161, 1226, 2, "+")], "cds": cds},
    "gene000002": {"strand": "-", "segments": [(6165, 6226, 0, "-"), (6000, 6065, 2, "-")], "cds": cds},
    "gene000003": {"strand": "+", "segments": [(0, 13, 1, "+")], "cds": partial_plus},
    "gene000004": {"strand": "-", "segments": [(length - 14, length, 2, "-")], "cds": partial_minus},
}
fixture = OUT / "known_coordinate_fixture"
utils.write_gff3(fixture_genes, str(fixture.with_suffix(".gff")))
utils.write_cds_fasta(fixture_genes, str(fixture.with_suffix(".fna")))
utils.write_protein_fasta(fixture_genes, str(fixture.with_suffix(".faa")), utils.load_codon_table(str(constants.CODON_TABLE_PATH)))
coordinate_check = read_outputs(fixture, sequence, "artificial", expected)
print("[preflight] coordinate fixture passed; begin full native annotation", flush=True)
native_dir = OUT / "native"
command = ["bash", str(PIPELINE / "run_annotation.sh"), "--input", str(fasta),
           "--accession", "artificial", "--model", str(ASSETS), "--output-dir", str(native_dir),
           "--gpu-count", "1", "--batch-size", "8", "--dispatch-chunk-size", "4",
           "--attn-implementation", "flash_attention_2", "--viterbi-backend", "gpu",
           "--min-record-length", "100000", "--min-cds-length", "96",
           "--max-boundary-intron-length", "20000", "--chunk-mbp", "500",
           "--codon-table", str(constants.CODON_TABLE_PATH)]
log = OUT / "native.log"
with log.open("w") as handle:
    subprocess.run(command, cwd=OUT, check=True, stdout=handle, stderr=subprocess.STDOUT)
if not native_dir.is_dir():
    native_dir = Path(str(native_dir) + "_WARNING")
log_text = log.read_text()
assert "effective=flash_attention_2" in log_text
parquets = list(native_dir.glob("*.raw_probs.parquet"))
assert len(parquets) == 1, parquets
table = pq.read_table(parquets[0])
assert table.num_rows == 1
row = table.to_pylist()[0]
assert row[constants.RECORD_ID_COLUMN] == "artificial"
assert row[constants.SEQUENCE_COLUMN] == sequence
for column in (constants.POSITIVE_STRAND_PROBABILITY_COLUMN, constants.NEGATIVE_STRAND_PROBABILITY_COLUMN):
    probabilities = np.asarray(row[column], dtype=np.float32)
    assert len(probabilities) == length
    assert np.isfinite(probabilities).all()
    assert ((probabilities >= 0) & (probabilities <= 1)).all()
native_check = read_outputs(native_dir / "prediction", sequence, "artificial")
summary = {
    "fixed_revision": ASSETS.name, "command": command, "gpu": device.name,
    "driver_inventory": driver, "cuda_build": torch.version.cuda,
    "flash_attention_bf16_kernel_passed": True,
    "effective_attention": "flash_attention_2", "viterbi_backend": "gpu",
    "artificial_record_bp": length, "full_context_bp": 98304,
    "coordinate_fixture": coordinate_check, "native_output": native_check,
    "native_directory": str(native_dir), "parquet": str(parquets[0]),
    "codon_warning": native_dir.name.endswith("_WARNING"),
    "elapsed_seconds": time.monotonic() - started,
    "original_train_executed": False, "accuracy_claim": False,
    "research_data_read": False, "training_updates": 0,
    "runtime": str(RUNTIME), "torch_version": torch.__version__,
    "node_local_transformers": {
        "version": "4.56.2", "source_file": str(transformers_source),
        "package_root": str(STAGE), "source_archive_job": os.environ["ARCHIVE_JOB"],
        "qualified_import_probe_job": probe_job,
        "original_packages_modified": False, "package_upgrade": False,
    },
    "elapsed_scope": "Python imports, artificial fixtures and native annotation; archive extraction excluded",
}
(OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
subprocess.run([sys.executable, str(ROOT / "scripts/experiments/M28-METHOD-RESTART/check_carbon_a_train_integration.py"),
                str(OUT)], check=True)
print(json.dumps(summary, indent=2), flush=True)
