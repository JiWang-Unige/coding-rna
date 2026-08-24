#!/usr/bin/env bash
#SBATCH --job-name=FPEXTM7
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --partition=private-teodoro-gpu
#SBATCH --gres=gpu:1
#SBATCH --output=outputs/fp_segnt_logs/%x_%j.out
#SBATCH --error=outputs/fp_segnt_logs/%x_%j.err
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043

# REANCHOR-HELDOUT-M7 CK2 — SegmentNT feature cache for held-out/UTR-rich species
# (arabidopsis_thaliana full, gallus_gallus SUBSET <=20Mb chromosomes). Idempotent skip-if-exists.
# Generous --time so the rematerialization-heavy forward never gets TIMEOUT-killed mid-savez.
# Usage: ssh baobab 'cd .. && sbatch scripts/run_fp_extract_M7.sbatch'
set -uo pipefail
ROOT="/srv/beegfs/scratch/shares/ds4dh/common/coding-rna"
cd "${ROOT}"
mkdir -p outputs/fp_segnt_logs
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u; conda activate coding-rna; set -u
export HF_HOME="${HOME}/.cache/huggingface"
CACHE="outputs/FP-SEGMENTNT-FEATCACHE/segment_nt_multi_species"
STATUS="${CACHE}/EXTRACT_M7_STATUS"
mkdir -p "${CACHE}"
echo "RUNNING" > "${STATUS}"
mark_failed() { rc=$?; [[ $rc -ne 0 ]] && echo "FAILED" > "${STATUS}"; exit $rc; }
trap mark_failed EXIT

for SP in arabidopsis_thaliana gallus_gallus; do
  if [[ -s "${CACHE}/${SP}.npz" ]]; then
    echo ">>> ${SP} cache exists, skipping"
    continue
  fi
  echo ">>> extracting ${SP} $(date +%T)"
  python -m src.foundation_probe.extract_segmentnt \
    --species "data/m1_screen/${SP}" --out-dir "${CACHE}" \
    --model segment_nt_multi_species --tile-tokens 1000 2>&1 \
    | grep -ivE 'rematerialization|bfc_allocator|Current allocation|^W[0-9]|captured constants'
  echo ">>> ${SP} done $(date +%T)"
done
echo "ALL_EXTRACT_M7_DONE" > "${STATUS}"
echo "ALL_EXTRACT_M7_DONE"
trap - EXIT
