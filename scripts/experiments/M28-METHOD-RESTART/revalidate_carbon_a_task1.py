"""CPU-only recovery of task1's completed native output after literal-codon validator repair."""
import json
import os
from pathlib import Path
import subprocess
import time

started = time.monotonic()
ROOT = Path(__file__).resolve().parents[3]
OUT = Path(os.environ["OUT"])
SOURCE = ROOT / "outputs/M28-CARBON-A-TRAIN-RAW/NC_003071.7"
print("[revalidation] import numeric and FASTA readers; no model imports", flush=True)
import numpy as np
import pyarrow.parquet as pq
from Bio import SeqIO
from carbon_a_output import read_outputs, literal_standard_translation

assert literal_standard_translation("ATGGCNTAA") == "MX*"
assert literal_standard_translation("ATGGCCTAA") == "MA*"
assert literal_standard_translation("ATGNNNA") == "MX"
assert (SOURCE / "STATUS").read_text().strip() == "FAILED"
assert (SOURCE / "JOBID").read_text().strip() == "13503402"
rows = subprocess.check_output(
    ["sacct", "-j", "13503402", "--format=JobIDRaw,State,ExitCode,ElapsedRaw", "-n", "-P"],
    text=True, env={**os.environ, "TZ": "UTC"},
)
states = {line.split("|")[0]: line.split("|")[1:] for line in rows.strip().splitlines()}
assert states["13503402"][:2] == ["FAILED", "1:0"]
assert states["13503402.batch"][:2] == ["FAILED", "1:0"]
assert states["13503402.extern"][:2] == ["COMPLETED", "0:0"]
error = (ROOT / "outputs/M28-PILOT-R5/logs/M28CARBONATRAIN_13503400_1.err").read_text()
assert "AssertionError: gene004167" in error and "carbon_a_output.py" in error
seqid, species, length = "NC_003071.7", "arabidopsis_thaliana", 19698289
fasta = ROOT / "outputs/M28-FULL-GRAMMAR-ORIGINAL-TRAIN-RAW" / seqid / "input.fasta"
records = list(SeqIO.parse(fasta, "fasta"))
assert len(records) == 1 and records[0].id == seqid and len(records[0]) == length
sequence = str(records[0].seq)
command = json.loads((SOURCE / "commands.json").read_text())
flags = dict(zip(command[2::2], command[3::2]))
assert flags["--input"] == str(fasta) and flags["--accession"] == seqid
assert Path(flags["--model"]).name == "f6789ee79e6ef397de333577ffec1cb734496620"
for key, value in {"--gpu-count":"1", "--batch-size":"8", "--dispatch-chunk-size":"4",
                   "--attn-implementation":"flash_attention_2", "--viterbi-backend":"gpu",
                   "--min-record-length":"100000", "--min-cds-length":"96",
                   "--max-boundary-intron-length":"20000", "--chunk-mbp":"500"}.items():
    assert flags[key] == value
assert Path(flags["--codon-table"]).name == "standard.json"
log = (SOURCE / "command.log").read_text()
assert "effective=flash_attention_2" in log and "[annotation] complete:" in log
native = SOURCE / "native"
if not native.is_dir():
    native = SOURCE / "native_WARNING"
parquets = list(native.glob("*.raw_probs.parquet"))
assert len(parquets) == 1
print("[revalidation] check full original DNA, both probability strands and every GFF/CDS/protein", flush=True)
table = pq.read_table(parquets[0])
assert table.num_rows == 1 and table["record_name"][0].as_py() == seqid
assert table["sequence"][0].as_py() == sequence
for column in ("pred_prob_positive_strand_cds", "pred_prob_negative_strand_cds"):
    probabilities = table[column][0].values.to_numpy(zero_copy_only=False)
    assert len(probabilities) == length and np.isfinite(probabilities).all()
    assert ((probabilities >= 0) & (probabilities <= 1)).all()
checked = read_outputs(native / "prediction", sequence, seqid)
summary = dict(
    completed=True, job_id=os.environ["SLURM_JOB_ID"], task=1, species=species, seqid=seqid, length=length,
    fixed_revision="f6789ee79e6ef397de333577ffec1cb734496620",
    shared_backbone_source_revision="b0483c23b6b63787b61a6d3a204a9b517d6ba345",
    runtime_job="13500755", preflight_job="13503202", input_fasta=str(fasta),
    raw_GFF=str((native / "prediction.gff").relative_to(ROOT)), parquet=str(parquets[0].relative_to(ROOT)),
    native_directory=str(native.relative_to(ROOT)), output_check=checked,
    codon_warning=native.name.endswith("_WARNING"), input_mode="fixed_original_TRAIN_shared_raw_FASTA",
    reference_used=False, labels_or_feature_cache_used=False, DEV_used=False, sealed_data_used=False,
    test_setaria_used=False, optimizer_steps=0, native_inference_executed=True,
    effective_attention="flash_attention_2", viterbi_backend="gpu",
    source_inference_job="13503402", source_inference_job_state="FAILED",
    source_failure="independent_validator_used_ambiguity_resolving_translation_for_GCN",
    raw_outputs_modified=False, native_GPU_inference_reexecuted=False,
    literal_standard_codon_table_with_non_ACGT_to_X=True,
    seconds=None, source_allocation_elapsed_seconds=int(states["13503402"][2]),
    source_wrapper_seconds_missing_due_to_failed_validation=True,
    CPU_revalidation_seconds=time.monotonic()-started,
    limitations="original TRAIN diagnostic; accepted native artifact via transparent CPU validation recovery, not independent generalization or fair speed benchmark",
)
(OUT / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
print(json.dumps({"completed":True, "source_inference_job":"13503402",
                  "source_inference_job_state":"FAILED", "output_check":checked,
                  "native_GPU_inference_reexecuted":False}), flush=True)
