#!/bin/bash
# Vipsania (M38 env, Streptophyta model j9m0cdmk pinned, pretrained-only mode = author defaults, no finetuning) on the 8 new-species internal-validation
# pseudo-contigs, scored by the unchanged M34 own_species ruler (m34_val.py).  Runs AFTER the three M40 arms (m40_run_all.sh calls it last).  ~3 GPU min (M38: Arabidopsis 119 Mb pretrained 240 s).
# Overlap with Vipsania's own data (protocol 34 section 6): Dioscorea rotundata GCF_009730915.1 is in the Streptophyta TRAINING set (sequence only);
# Solanum lycopersicum GCF_036512215.1 is in its TEST set; the other six are in neither.  Output: ~/coding-rna/outputs/M40-SEMISUP-MLM/vip_newsp/
source $HOME/coding-rna-reorg/experiments/M40-SEMISUP-MLM/scripts/m40_env.sh
PYPC2=$PY
source $HOME/coding-rna-reorg/experiments/M38-VIPSANIA/scripts/vip_env.sh
D=$OUT40/vip_newsp; mkdir -p $D $HOME/coding-rna/logs
LOG=$HOME/coding-rna/logs/m40-vip-newsp
CPUS=${VIP_CPUS:-2,4,6,8,10,12,14,16-24}
ev() { echo "{\"event\":\"$1\",\"time\":\"$(date -Is)\"$2}" >> $D/events.jsonl; }
[ -e $D/DONE ] && exit 0
exec 9>$D/.lock; flock -n 9 || exit 3
(cd $VIPSANIA_CACHE/models/$MODEL_ID && sha256sum -c --quiet $M38/model_$MODEL_ID.sha256) >/dev/null 2>&1 || { ev refused ",\"reason\":\"model checksum\""; exit 11; }
REG=$(cd $M40S && $PYPC2 m40_vip_newsp.py regions 2>/dev/null)
[ -s $D/newsp_val8M.fa ] || (cd $M40S && $PYPC2 m40_vip_newsp.py fasta $D/newsp_val8M.fa 2>>$LOG.err) || { ev fasta_failed; exit 9; }
export TF_FORCE_GPU_ALLOW_GROWTH=true
export LD_LIBRARY_PATH=$($VPY -c "import site, glob, os; print(':'.join(sorted(glob.glob(os.path.join(site.getsitepackages()[0], 'nvidia', '*', 'lib')))))" 2>/dev/null)${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
if [ ! -s $D/pre/vipsania.gff3 ]; then
  mkdir -p $D/pre
  for B in auto 4 2 1; do
    for a in $(seq 1 60); do f=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1); [ "${f:-0}" -ge 11520 ] && break; sleep 120; done
    ARGS=(); [ $B = auto ] || ARGS=(-B $B)
    ev annotate_start ",\"B\":\"$B\""
    timeout 5h nice -n 10 taskset -c $CPUS $VIP annotate $MODEL_ID $D/newsp_val8M.fa -o $D/pre/vipsania.gff3.part "${ARGS[@]}" >> $LOG.out 2>> $LOG.err < /dev/null && \
      mv $D/pre/vipsania.gff3.part $D/pre/vipsania.gff3 && break
    ev annotate_failed ",\"B\":\"$B\""
  done
  [ -s $D/pre/vipsania.gff3 ] || { ev gave_up; exit 9; }
  ev annotate_done ",\"genes\":$(awk -F'\t' '$3=="gene"{n++} END{print n+0}' $D/pre/vipsania.gff3)"
fi
python3.12 -I $M38_SCRIPTS/vip_score.py convert --gff $D/pre/vipsania.gff3 --out $D/run_pre --contigs $REG > $D/convert_pre.json || { ev convert_failed; exit 9; }
cd $M34S && M34_BACKEND=pc2 taskset -c 18-23 $PYPC2 m34_val.py $D/run_pre 0 $REG --tag own_species --arm vipsania_pre --metrics $D/metrics.jsonl >> $LOG.out 2>> $LOG.err || { ev score_failed; exit 9; }
touch $D/DONE; ev done
