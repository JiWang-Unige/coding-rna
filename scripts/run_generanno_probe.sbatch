#!/usr/bin/env bash
#SBATCH --job-name=generanno-native-probe
#SBATCH --partition=private-teodoro-gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=05:00:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/GENERANNO-NATIVE-PROBE/logs/probe_%j.out
#SBATCH --error=outputs/GENERANNO-NATIVE-PROBE/logs/probe_%j.err

set -euo pipefail
cd /srv/beegfs/scratch/shares/ds4dh/common/coding-rna

EXP_ID="GENERANNO-NATIVE-PROBE"
OUT_DIR="outputs/${EXP_ID}"
METRIC_DIR="${OUT_DIR}/metrics"
SP="arabidopsis_thaliana"
SP_PATH="data/m1_screen/arabidopsis_thaliana"
mkdir -p "${OUT_DIR}/logs" "${METRIC_DIR}"

set +u
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
conda activate generanno
set -u
echo "env=$(which python) host=$(hostname)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# ---- native inference on the arabidopsis TEST split (official decode, official postprocess) ----
python scripts/probe_generanno_native.py \
  --species-path "${SP_PATH}" --species-name "${SP}" \
  --out-dir "${OUT_DIR}" --exp-id "${EXP_ID}" \
  --context-length 6000 --overlap-length 1024 --micro-bs 1 --bf16 --postproc on

# ---- evaluate on OUR ruler (--span-mode cds, screen profile) ----
python scripts/eval_gene_body_mask.py \
  --reference-gtf "${OUT_DIR}/eval_subsets/${SP}/reference.gff3" \
  --prediction-gtf "${OUT_DIR}/predictions/${SP}.gff" \
  --genome-fasta "${OUT_DIR}/eval_subsets/${SP}/genome.fa" \
  --output-json "${METRIC_DIR}/${SP}.metrics.json" \
  --experiment-id "${EXP_ID}_${SP}" --profile screen --span-mode cds

cp "${METRIC_DIR}/${SP}.metrics.json" "${METRIC_DIR}/metrics.json"
echo "===== METRICS ====="
python -c "import json; d=json.load(open('${METRIC_DIR}/metrics.json')); print(json.dumps({k:d[k] for k in ['primary_metric','intergenic_specificity','intergenic_FPR','gene_body_F1_unconstrained','gene_body_precision','gene_body_recall','predicted_gene_count','predicted_gene_count_ratio_vs_reference'] if k in d}, indent=2))"
echo "COMPLETED" > "${OUT_DIR}/STATUS"
echo "PROBE_PIPELINE_DONE"
