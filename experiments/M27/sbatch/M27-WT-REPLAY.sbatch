#!/usr/bin/env bash
#SBATCH --job-name=M27WT
#SBATCH --no-requeue
#SBATCH --partition=private-teodoro-gpu
#SBATCH --gres=gpu:nvidia_geforce_rtx_3090:1
#SBATCH --exclude=gpu023,gpu024,gpu036,gpu037,gpu038,gpu039,gpu040,gpu041,gpu042,gpu043
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=outputs/M27-WT-REPLAY-R1/logs/%x_%j.out
#SBATCH --error=outputs/M27-WT-REPLAY-R1/logs/%x_%j.err
set -euo pipefail
cd /home/users/j/jwang/coding-rna
M27_ROOT=$(pwd -P)
M27_OUT="$M27_ROOT/outputs/M27-WT-REPLAY-R1"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" > "$M27_OUT/STATUS"; fi' EXIT
printf 'RUNNING\n' > "$M27_OUT/STATUS"
source /opt/ebsofts/Mamba/23.1.0-4/etc/profile.d/conda.sh
set +u
conda activate generanno
set -u
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=8
python -m pytest -q -p no:cacheprovider tests/test_m27_wt_replay.py --junitxml="$M27_OUT/tests.xml"
mkdir "$M27_OUT/ANNEVO" "$M27_OUT/Tiberius"
test ! -e refs/weights/tiberius-m27-nosm-v2/tiberius_nosm_weights_v2
tar -xzf refs/weights/tiberius-m27-nosm-v2/tiberius_nosm_weights_v2.tar.gz -C refs/weights/tiberius-m27-nosm-v2
set +u
conda activate annevo
set -u
python -c 'import torch; print("Torch", torch.__version__, "CUDA", torch.version.cuda, "devices", torch.cuda.device_count()); assert torch.cuda.is_available() and torch.cuda.device_count() == 1'
(
    cd refs/repos/annevo-2026
    /usr/bin/time -v python -u prediction.py \
        -g "$M27_ROOT/data/m27_grch38_chr22/chr22.unmasked.fa" \
        -m saved_model/ANNEVO_Mammalia.pt -l Mammalia -s 100 \
        -p "$M27_OUT/ANNEVO/WT.h5" --batch_size 2 --num_workers 2
) > "$M27_OUT/ANNEVO/prediction.log" 2>&1
for mode in native replay; do
    (
        cd refs/repos/annevo-2026
        /usr/bin/time -v python -u decoding.py \
            -g "$M27_ROOT/data/m27_grch38_chr22/chr22.unmasked.fa" \
            -p "$M27_OUT/ANNEVO/WT.h5" -o "$M27_OUT/ANNEVO/$mode.gff3" \
            -t 8 --min_intron_length 20 --min_prot_length 100 --show_log
    ) > "$M27_OUT/ANNEVO/$mode.log" 2>&1
    if rg -q 'Process failed:' "$M27_OUT/ANNEVO/$mode.log"; then
        printf 'ANNEVO suppressed a worker exception; partial output is not qualified\n' >&2
        exit 1
    fi
done
for mode in record replay; do
    label=native
    if [[ "$mode" == replay ]]; then label=replay; fi
    /usr/bin/time -v apptainer exec --nv --cleanenv \
        --env "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES" \
        --env PYTHONDONTWRITEBYTECODE=1 --env OMP_NUM_THREADS=8 \
        --bind "$M27_ROOT:/work" --pwd /work \
        refs/repos/tiberius-2024/singularity/tiberius_2.0.5.sif \
        python -u scripts/experiments/M27-ALLELIC-INTEGRITY/wt_tiberius.py \
        --cache-mode "$mode" --cache-dir outputs/M27-WT-REPLAY-R1/Tiberius/cache \
        --trace "outputs/M27-WT-REPLAY-R1/Tiberius/${label}_trace.json" \
        --genome data/m27_grch38_chr22/chr22.unmasked.fa \
        --model refs/weights/tiberius-m27-nosm-v2/tiberius_nosm_weights_v2 \
        --no_softmasking --seq_len 400050 --batch_size 2 \
        --out "outputs/M27-WT-REPLAY-R1/Tiberius/$label.gtf" \
        > "$M27_OUT/Tiberius/$label.log" 2>&1
done
python -u scripts/experiments/M27-ALLELIC-INTEGRITY/wt_compare.py --output-dir "$M27_OUT"
du -sb "$M27_OUT" data/m27_grch38_chr22 refs/weights/tiberius-m27-nosm-v2
printf 'COMPLETED\n' > "$M27_OUT/STATUS"
