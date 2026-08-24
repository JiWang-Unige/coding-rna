#!/usr/bin/env bash
#SBATCH --job-name=generanno-probe-smoke
#SBATCH --partition=private-teodoro-gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=00:40:00
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --output=outputs/GENERANNO-NATIVE-PROBE/logs/smoke_%j.out
#SBATCH --error=outputs/GENERANNO-NATIVE-PROBE/logs/smoke_%j.err

set -euo pipefail
cd /srv/beegfs/scratch/shares/ds4dh/common/coding-rna

EXP_ID="GENERANNO-NATIVE-PROBE-SMOKE"
OUT_DIR="outputs/GENERANNO-NATIVE-PROBE/smoke"
METRIC_DIR="${OUT_DIR}/metrics"
SP="arabidopsis_thaliana"
SP_PATH="data/m1_screen/arabidopsis_thaliana"
mkdir -p "outputs/GENERANNO-NATIVE-PROBE/logs" "${METRIC_DIR}"

set +u
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
conda activate generanno
set -u
echo "env=$(which python) host=$(hostname)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

# smoke: ONE small non-test seqid (chloroplast 154kb) -> validate load+inference+decode+eval path
python scripts/probe_generanno_native.py \
  --species-path "${SP_PATH}" --species-name "${SP}" \
  --out-dir "${OUT_DIR}" --exp-id "${EXP_ID}" \
  --context-length 6000 --overlap-length 1024 --micro-bs 1 --bf16 --postproc on \
  --seqids NC_000932.1

python scripts/eval_gene_body_mask.py \
  --reference-gtf "${OUT_DIR}/eval_subsets/${SP}/reference.gff3" \
  --prediction-gtf "${OUT_DIR}/predictions/${SP}.gff" \
  --genome-fasta "${OUT_DIR}/eval_subsets/${SP}/genome.fa" \
  --output-json "${METRIC_DIR}/${SP}.metrics.json" \
  --experiment-id "${EXP_ID}_${SP}" --profile screen --span-mode cds

echo "===== SMOKE METRICS ====="
python -c "import json; d=json.load(open('${METRIC_DIR}/${SP}.metrics.json')); print(json.dumps({k:d[k] for k in ['primary_metric','intergenic_specificity','gene_body_F1_unconstrained','predicted_gene_count'] if k in d}, indent=2))"
echo "SMOKE_OK" > "${OUT_DIR}/STATUS"
echo "SMOKE_PIPELINE_DONE"
