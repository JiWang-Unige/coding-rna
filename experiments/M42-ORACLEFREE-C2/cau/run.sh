#!/bin/bash
# M42 gate-1 GPU job on the CAU ChaoBaiHe partition (protocol 35 §4).  All parameters in the config file ($1):
#   JOB    tag of this job (times.tsv)            STAGES  space list of train / infer (in order)
#   ARM    U | L | F | U3 (train)                 STOP    optional --stop-at (smoke)       TAG   output subdir (train / smoke)
#   MODEL  infer: model file (default: the trained arm's model_<STOP or 512>.pt; Z = T2 step 2304)   NAME  infer output name    REGIONS  infer regions
# Code, data and the PlantCAD2 weights live under $G/home (HOME is pointed there so the repository's ~/ paths resolve).  Finished stages are skipped.
umask 077; set -uo pipefail
G=$HOME/.cache/fontconfig-x2/g10; export HOME=$G/home; source "$1"
R=$HOME/coding-rna-reorg; S=$R/experiments/M42-ORACLEFREE-C2/scripts; O=$R/outputs/M42-ORACLEFREE-C2; PY=$G/env/bin/python
export M34_BACKEND=pc2 PYTHONDONTWRITEBYTECODE=1 HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 PYTHONUNBUFFERED=1 \
       PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True TRITON_CACHE_DIR=$G/tc/$JOB
T=$G/times.tsv; st() { echo -e "$JOB\t$1\t$(date +%s)\t$(date)\t${SLURM_JOB_ID:-}" >> $T; }
D=$R/data
for f in $D/m42/ziziphus_jujuba/{positives.json,background_V1.json,background_V3.json,split.json,eval_regions.txt} $D/m1_screen/ziziphus_jujuba/{genome.fa,lengths.json} \
         $R/outputs/M34-NTV2-SCALEUP/pc2_T2/model_002304.pt $HOME/plantcad2_work/models/PlantCAD2-Small-l24-d0768/pytorch_model.bin $S/m42_train.py $S/m42_infer.py; do
  [ -s $f ] || { echo "MISSING $f"; st missing_input; exit 4; }
done
nvidia-smi --query-gpu=name,driver_version,memory.total,memory.used --format=csv,noheader > $G/logs/gpu_$JOB.txt
$PY -c "import torch,mamba_ssm,causal_conv1d,transformers,numpy;assert torch.cuda.is_available();print(torch.__version__,mamba_ssm.__version__,causal_conv1d.__version__,transformers.__version__,numpy.__version__,torch.cuda.get_device_name(0))" >> $G/logs/gpu_$JOB.txt || { st env_fail; exit 3; }
cd $S
CPUS=$(taskset -pc $$ | sed 's/.*: //')
for stage in $STAGES; do
  case $stage in
  train)
    OD=$O/${TAG:-train}/${ARM}_ziziphus_jujuba; mkdir -p $OD; END=${STOP:-512}
    if ls $OD/DONE_$END.json >/dev/null 2>&1; then continue; fi
    st train_$ARM start
    $PY m42_train.py run --arm $ARM --out $OD --workers ${WORKERS:-8} --cpus "$CPUS" ${STOP:+--stop-at $STOP} --resume >> $OD/train.log 2>&1
    ls $OD/DONE_$END.json >/dev/null 2>&1 || { st train_$ARM fail; exit 1; }
    st train_$ARM end ;;
  infer)
    M=${MODEL:-$O/${TAG:-train}/${ARM}_ziziphus_jujuba/model_$(printf %06d ${STOP:-512}).pt}; E=$O/${TAG:-train}/eval/$NAME; mkdir -p $E
    if ls $E/DONE_*.json >/dev/null 2>&1; then continue; fi
    st infer_$NAME start
    $PY m42_infer.py --model $M --regions ${REGIONS:-$(cat $D/m42/ziziphus_jujuba/eval_regions.txt)} --out $E --decode-workers ${DW:-4} --cpus "$CPUS" >> $E/infer.log 2>&1
    ls $E/DONE_*.json >/dev/null 2>&1 || { st infer_$NAME fail; exit 1; }
    rm -rf $E/pot
    st infer_$NAME end ;;
  esac
done
st all done; echo CHAIN_DONE
