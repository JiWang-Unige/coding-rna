"""Fixed Carbon-A pipeline on the six shared original TRAIN FASTAs, no labels."""
import argparse
import importlib.metadata
import importlib.util
import json
import os
import subprocess
import time
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
from Bio import SeqIO
from carbon_a_output import read_outputs
import run_native_train_baseline as N

ROOT = N.ROOT
SCOPE = N.SCOPE
REVISION = "f6789ee79e6ef397de333577ffec1cb734496620"
ASSETS = ROOT / "refs/carbon-a" / REVISION
OUT = ROOT / "outputs/M28-CARBON-A-TRAIN-RAW"
CONTEXT24_OUT = ROOT / "outputs/M28-CARBON-A-CONTEXT24-TRAIN-RAW"


def output_root(context_bp):
    return {98304: OUT, 24576: CONTEXT24_OUT}[context_bp]


def pipeline_commands(fasta, native, seqid, context_bp):
    command = ["bash", str(ASSETS / "annotation_pipeline/run_annotation.sh"),
               "--input", str(fasta), "--accession", seqid, "--model", str(ASSETS),
               "--output-dir", str(native), "--gpu-count", "1", "--batch-size", "8",
               "--dispatch-chunk-size", "4", "--attn-implementation", "flash_attention_2",
               "--viterbi-backend", "gpu", "--min-record-length", "100000", "--min-cds-length", "96",
               "--max-boundary-intron-length", "20000", "--chunk-mbp", "500",
               "--codon-table", str(ASSETS / "annotation_pipeline/codon_tables/standard.json")]
    if context_bp == 98304:
        return [command]
    probability_command = ["python", str(ASSETS / "annotation_pipeline/src/model_inference.py"),
                           "--input", str(fasta), "--output_path", str(native),
                           "--model_name", str(ASSETS), "--gpu_count", "1", "--batch_size", "8",
                           "--dispatch_chunk_size", "4", "--attn_implementation", "flash_attention_2",
                           "--context_length", "4096", "--overlap_length", "1024"]
    return [probability_command, command + ["--skip-inference"]]


def run(task, runtime_job, preflight_job, context_bp=98304):
    out = output_root(context_bp)
    if task not in range(len(SCOPE)):
        raise ValueError("Only the original six TRAIN records are allowed")
    runtime_folder = ROOT / "outputs/M28-CARBON-A-RUNTIME" / runtime_job
    preflight_folder = ROOT / "outputs/M28-CARBON-A-PREFLIGHT" / preflight_job
    for folder in (runtime_folder, preflight_folder):
        if (folder / "STATUS").read_text().strip() != "COMPLETED":
            raise ValueError("Supported runtime and artificial native preflight must complete first")
    runtime = json.loads((runtime_folder / "summary.json").read_text())
    preflight = json.loads((preflight_folder / "summary.json").read_text())
    assert runtime["isolated_supported_imports_passed"]
    assert preflight["fixed_revision"] == REVISION
    assert preflight["effective_attention"] == "flash_attention_2"
    assert preflight["viterbi_backend"] == "gpu" and preflight["flash_attention_bf16_kernel_passed"]
    assert preflight["coordinate_fixture"]["coordinate_cds_protein_roundtrip"]
    assert preflight["native_output"]["coordinate_cds_protein_roundtrip"]
    assert preflight["node_local_transformers"]["source_archive_job"] == os.environ["ARCHIVE_JOB"]
    stage = Path(os.environ["TASK_TF_STAGE"]).resolve()
    transformers_source = Path(importlib.util.find_spec("transformers").origin).resolve()
    assert stage in transformers_source.parents
    assert importlib.metadata.version("transformers") == "4.56.2"
    source = N.reconcile_train_input(task)
    species, seqid, length = SCOPE[task]
    fasta = N.input_fasta(seqid)
    N.check_fasta(fasta, seqid, length)
    records = list(SeqIO.parse(fasta, "fasta"))
    assert len(records) == 1 and records[0].id == seqid
    sequence = str(records[0].seq)
    folder = out / seqid
    folder.mkdir(parents=True, exist_ok=False)
    status = folder / "STATUS"
    status.write_text("RUNNING\n")
    (folder / "JOBID").write_text(os.environ["SLURM_JOB_ID"] + "\n")
    started = time.monotonic()
    try:
        native = folder / "native"
        commands = pipeline_commands(fasta, native, seqid, context_bp)
        logged_commands = commands[0] if len(commands) == 1 else commands
        (folder / "commands.json").write_text(json.dumps(logged_commands, indent=2) + "\n")
        with (folder / "command.log").open("x") as handle:
            for command in commands:
                subprocess.run(command, cwd=folder, check=True, stdout=handle, stderr=subprocess.STDOUT)
        if not native.is_dir():
            native = Path(str(native) + "_WARNING")
        log = (folder / "command.log").read_text()
        assert "effective=flash_attention_2" in log and "[annotation] complete:" in log
        parquets = list(native.glob("*.raw_probs.parquet"))
        assert len(parquets) == 1
        table = pq.read_table(parquets[0])
        assert table.num_rows == 1
        assert table["record_name"][0].as_py() == seqid
        assert table["sequence"][0].as_py() == sequence
        for column in ("pred_prob_positive_strand_cds", "pred_prob_negative_strand_cds"):
            probabilities = table[column][0].values.to_numpy(zero_copy_only=False)
            assert len(probabilities) == length
            assert np.isfinite(probabilities).all()
            assert ((probabilities >= 0) & (probabilities <= 1)).all()
        checked = read_outputs(native / "prediction", sequence, seqid)
        summary = {
            "completed": True, "task": task, "species": species, "seqid": seqid, "length": length,
            "job_id": os.environ["SLURM_JOB_ID"], "fixed_revision": REVISION,
            "context_bp": context_bp, "context_tokens": context_bp // 6,
            "overlap_bp": 6144, "overlap_tokens": 1024,
            "window_policy_modified": context_bp != 98304,
            "shared_backbone_source_revision": "b0483c23b6b63787b61a6d3a204a9b517d6ba345",
            "runtime_job": runtime_job, "preflight_job": preflight_job,
            "node_local_transformers": {
                "version": "4.56.2", "source_archive_job": os.environ["ARCHIVE_JOB"],
                "source_file": str(transformers_source),
            },
            "input_fasta": str(fasta), "original_input_receipt_job_id": source["job_id"],
            "raw_GFF": str((native / "prediction.gff").relative_to(ROOT)),
            "parquet": str(parquets[0].relative_to(ROOT)), "output_check": checked,
            "native_directory": str(native.relative_to(ROOT)),
            "codon_warning": native.name.endswith("_WARNING"),
            "input_mode": N.INPUT_MODE, "reference_used": False, "labels_or_feature_cache_used": False,
            "DEV_used": False, "sealed_data_used": False, "test_setaria_used": False,
            "optimizer_steps": 0, "native_inference_executed": True,
            "effective_attention": "flash_attention_2", "viterbi_backend": "gpu",
            "seconds": time.monotonic() - started,
            "limitations": "original TRAIN diagnostic; not independent generalization or fair speed benchmark",
        }
        with (folder / "summary.json").open("x") as handle:
            json.dump(summary, handle, indent=2, allow_nan=False)
            handle.write("\n")
        status.write_text("COMPLETED\n")
        print(json.dumps(summary), flush=True)
    except BaseException:
        status.write_text("FAILED\n")
        raise


def completed_runs(context_bp=98304):
    runs = {}
    out = output_root(context_bp)
    for task, (species, seqid, length) in enumerate(SCOPE):
        folder = out / seqid
        recovered_original = context_bp == 98304 and seqid == "NC_003071.7"
        if recovered_original:
            folder = ROOT / "outputs/M28-CARBON-A-TRAIN-RAW-VALIDATION" / seqid
        assert (folder / "STATUS").read_text().strip() == "COMPLETED"
        row = json.loads((folder / "summary.json").read_text())
        if recovered_original:
            assert row["source_inference_job"] == "13503402" and row["source_inference_job_state"] == "FAILED"
            assert row["literal_standard_codon_table_with_non_ACGT_to_X"]
            assert not row["native_GPU_inference_reexecuted"] and not row["raw_outputs_modified"]
        if context_bp == 24576:
            assert (row["context_bp"], row["context_tokens"], row["overlap_bp"], row["overlap_tokens"]) == (24576, 4096, 6144, 1024)
            assert row["window_policy_modified"]
        assert row["completed"] and row["native_inference_executed"]
        assert (row["task"], row["species"], row["seqid"], row["length"]) == (task, species, seqid, length)
        assert row["fixed_revision"] == REVISION and row["input_fasta"] == str(N.input_fasta(seqid))
        assert row["input_mode"] == N.INPUT_MODE and row["optimizer_steps"] == 0
        assert not any(row[key] for key in ("reference_used", "labels_or_feature_cache_used", "DEV_used", "sealed_data_used", "test_setaria_used"))
        assert row["output_check"]["coordinate_cds_protein_roundtrip"]
        assert row["effective_attention"] == "flash_attention_2" and row["viterbi_backend"] == "gpu"
        assert (ROOT / row["raw_GFF"]).is_file()
        row["accepted_run_receipt"] = str((folder / "summary.json").relative_to(ROOT))
        runs[seqid] = row
    return runs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", type=int, required=True)
    parser.add_argument("--runtime-job", required=True)
    parser.add_argument("--preflight-job", required=True)
    parser.add_argument("--context-bp", type=int, choices=(98304, 24576), default=98304)
    args = parser.parse_args()
    run(args.task, args.runtime_job, args.preflight_job, args.context_bp)
