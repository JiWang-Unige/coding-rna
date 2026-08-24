#!/usr/bin/env bash
#SBATCH --job-name=M25GENSTRUCTs0
#SBATCH --partition=private-teodoro-gpu
#SBATCH --time=6-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --gres=gpu:nvidia_geforce_rtx_3090:1
#SBATCH --output=outputs/M25-GENERANNO-1P2B-STRUCTURAL-HEADS-s0/logs/%x_%j.out
#SBATCH --error=outputs/M25-GENERANNO-1P2B-STRUCTURAL-HEADS-s0/logs/%x_%j.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

set -euo pipefail

ROOT="/home/users/j/jwang/coding-rna"
EXP_ID="M25-GENERANNO-1P2B-STRUCTURAL-HEADS-s0"
CONFIG="${ROOT}/configs/${EXP_ID}.yaml"
OUT_DIR="${ROOT}/outputs/${EXP_ID}"
STATUS_FILE="${OUT_DIR}/STATUS"
SETARIA_FASTA="${ROOT}/data/m25_setaria/GCF_005286985.2/genome.fa"
RESOLVED_CONFIG="${OUT_DIR}/config_resolved.yaml"
ALLOWLIST="${OUT_DIR}/primary_chromosome_allowlist.json"
CHECKPOINT="${OUT_DIR}/checkpoints/selected_checkpoint.pt"
DECODE_PARAMS="${OUT_DIR}/validation_decode_parameters.json"
FULL_GFF="${OUT_DIR}/predictions/setaria_viridis.full.gff3"
ABLATION_GFF="${OUT_DIR}/predictions/setaria_viridis.ablation.gff3"
MARKER="${OUT_DIR}/SETARIA_EMBARGO_RELEASED.json"
METRICS="${OUT_DIR}/metrics/setaria_viridis.structure_metrics.json"
REFERENCE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/005/286/985/GCF_005286985.2_Setaria_viridis_v4.0/GCF_005286985.2_Setaria_viridis_v4.0_genomic.gff.gz"
REFERENCE_DIR="${OUT_DIR}/evaluation_reference"
REFERENCE_GFF_GZ="${REFERENCE_DIR}/GCF_005286985.2_Setaria_viridis_v4.0_genomic.gff.gz"
REFERENCE_GFF="${REFERENCE_DIR}/GCF_005286985.2_Setaria_viridis_v4.0_genomic.gff3"

cd "${ROOT}"

if [[ ! -f "${CONFIG}" ]]; then
  echo "Missing frozen config: ${CONFIG}" >&2
  exit 2
fi
if [[ ! -f "${SETARIA_FASTA}" ]]; then
  echo "Missing frozen Setaria assembly FASTA: ${SETARIA_FASTA}" >&2
  exit 2
fi
python - "${CONFIG}" <<'PY'
import sys
import yaml

with open(sys.argv[1], encoding="utf-8") as handle:
    setaria = yaml.safe_load(handle)["data"]["setaria"]
if setaria["reference_annotation"] is not None or setaria["reference_annotation_source"] is not None:
    raise SystemExit("Setaria annotation path/source must remain null before training and inference")
PY
if [[ -e "${MARKER}" ]]; then
  echo "Refusing to reuse an existing embargo-release marker: ${MARKER}" >&2
  exit 3
fi
if [[ -d "${OUT_DIR}" ]] && find "${OUT_DIR}" -mindepth 1 -type f ! -path "${OUT_DIR}/logs/*" ! -name JOBID -print -quit | grep -q .; then
  echo "Refusing to overwrite non-log output in ${OUT_DIR}" >&2
  exit 3
fi

for annotation in $(find "${ROOT}/data" -type f \( \
  -iname '*.gff' -o -iname '*.gff3' -o -iname '*.gtf' -o -iname '*.gbff' -o \
  -iname '*.gff.gz' -o -iname '*.gff3.gz' -o -iname '*.gtf.gz' -o -iname '*.gbff.gz' \
\)); do
  annotation_lower="${annotation,,}"
  if [[ "${annotation_lower}" == *setaria* || \
        "${annotation_lower}" == *viridis* || \
        "${annotation_lower}" == *gcf_005286985* || \
        "${annotation_lower}" == *setaria_viridis_v4.0* ]]; then
    echo "Setaria annotation embargo violation before training: ${annotation}" >&2
    exit 4
  fi
done

mkdir -p "${OUT_DIR}/logs" "${OUT_DIR}/metrics"
cp "${CONFIG}" "${RESOLVED_CONFIG}"
echo "RUNNING" > "${STATUS_FILE}"
mark_failed() {
  rc=$?
  if [[ ${rc} -ne 0 ]]; then
    echo "FAILED" > "${STATUS_FILE}"
  fi
  exit ${rc}
}
trap mark_failed EXIT

source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export HF_HOME="${HOME}/.cache/huggingface"
export PYTHONUNBUFFERED=1

{
  echo "exp_id=${EXP_ID}"
  echo "model_name=GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview"
  echo "seed=0"
  echo "fits=1"
  echo "slurm_job_id=${SLURM_JOB_ID:-NA}"
  echo "host=$(hostname)"
  date -Is
  python -V
  python - <<'PY'
import torch, transformers, peft
print("torch", torch.__version__, "cuda", torch.cuda.is_available())
print("transformers", transformers.__version__, "peft", peft.__version__)
PY
} > "${OUT_DIR}/env.txt" 2>&1

python -u -m src.foundation_probe.train_generanno_structural_heads \
  --config "${CONFIG}" \
  --species data/m1_screen/arabidopsis_thaliana data/m1_screen/oryza_sativa \
  --out-dir "${OUT_DIR}" \
  --setaria-fasta "${SETARIA_FASTA}" \
  --model-name GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview \
  --model-task token_classification \
  --window 6144 \
  --sample-fraction 0.12 \
  --epochs 3 \
  --batch-size 1 \
  --predict-batch 1 \
  --head-lr 8e-4 \
  --lora-lr 2e-5 \
  --seed 0 \
  --fp-lambda 2.5 \
  --lora-r 8 \
  --lora-alpha 16 \
  --lora-dropout 0.05 \
  --lora-targets q_proj,k_proj,v_proj,o_proj \
  --attn-implementation sdpa \
  --bf16 \
  --limit-train-windows 1536 \
  --limit-val-windows 768

if [[ "$(<"${STATUS_FILE}")" == "STOP_M25_BRANCH" ]]; then
  echo "Validation selection found no admissible tuple; Setaria inference was not run."
  exit 0
fi

for required in "${RESOLVED_CONFIG}" "${ALLOWLIST}" "${CHECKPOINT}" "${DECODE_PARAMS}" "${FULL_GFF}" "${ABLATION_GFF}"; do
  if [[ ! -s "${required}" ]]; then
    echo "Required immutable artifact is missing or empty: ${required}" >&2
    exit 5
  fi
done

python - "${RESOLVED_CONFIG}" "${ALLOWLIST}" "${CHECKPOINT}" "${DECODE_PARAMS}" \
  "${SETARIA_FASTA}" "${FULL_GFF}" "${ABLATION_GFF}" "${MARKER}" <<'PY'
import hashlib
import json
import sys

labels = ("config_resolved", "primary_chromosome_allowlist", "checkpoint", "validation_decode_selection",
          "genome_fasta", "full_prediction_gff3", "ablation_prediction_gff3")
paths = sys.argv[1:8]
marker = sys.argv[8]

def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

payload = {label: {"path": path, "sha256": sha256(path)}
           for label, path in zip(labels, paths)}
with open(marker, "x", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY

mkdir -p "${REFERENCE_DIR}"
curl -fL "${REFERENCE_URL}" -o "${REFERENCE_GFF_GZ}"
gzip -dc "${REFERENCE_GFF_GZ}" > "${REFERENCE_GFF}"

python -u scripts/eval_m25_structure.py \
  --embargo-marker "${MARKER}" \
  --allowlist "${ALLOWLIST}" \
  --genome-fasta "${SETARIA_FASTA}" \
  --reference-gff3 "${REFERENCE_GFF}" \
  --full-gff3 "${FULL_GFF}" \
  --ablation-gff3 "${ABLATION_GFF}" \
  --output-json "${METRICS}"

python - "${METRICS}" "${STATUS_FILE}" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    decision = json.load(handle)["gate"]["decision"]
if decision not in {"PASSED_DISCOVERY_GATE", "STOP_M25_BRANCH"}:
    raise SystemExit(f"Unexpected M25 evaluator decision: {decision}")
with open(sys.argv[2], "w", encoding="utf-8") as handle:
    handle.write(decision + "\n")
PY

case "$(<"${STATUS_FILE}")" in
  PASSED_DISCOVERY_GATE|STOP_M25_BRANCH) ;;
  *)
    echo "Evaluator did not write a terminal M25 status to ${STATUS_FILE}" >&2
    exit 7
    ;;
esac

echo "DONE ${EXP_ID}: $(<"${STATUS_FILE}")"
trap - EXIT
