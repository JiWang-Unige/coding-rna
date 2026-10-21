#!/bin/bash
# Detached queue (runs on ont, nohup).  GPU order (user decision 2026-10-21): M33 -> M35 (PlantCAD2 + LoRA) -> backbone decision -> M34.  Starts NOTHING until BOTH
#   ~/coding-rna/outputs/M35-PLANTCAD2-LORA/M35_DONE   and   ~/coding-rna/outputs/M34-NTV2-SCALEUP/GO_STAGE1   exist (GO_STAGE1 is created by the coordinator after the user chose
#   the backbone; its first word selects the backend: `ntv2` (default when empty) or `pc2`) and no M32/M33/M35 process is left (m34_gate.py OK in 3 consecutive checks 5 min apart).
# Then: smoke -> arm T2 -> arm T10 (stage one of protocol 19 §3.5), each exactly once; finally M34_STAGE1_DONE.  Other arms are never started here.
# For pc2 the M35 LoRA module (experiments/M35-PLANTCAD2-LORA/scripts/pc2_lora.py) must exist, otherwise the queue stops (event pc2_lora_module_missing).
O=$HOME/coding-rna/outputs/M34-NTV2-SCALEUP
mkdir -p $O
exec 8>$O/.queue.lock
flock -n 8 || { echo "queue already running"; exit 3; }
cd $HOME/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts
q() { echo "{\"event\":\"$1\",\"time\":\"$(date -Is)\"$2}" >> $O/queue.jsonl; }
q queued
PY0=$HOME/envs/generanno/bin/python
ok=0
while [ $ok -lt 3 ]; do
  r=$($PY0 m34_gate.py 2>&1 | tail -1)
  if [ "$r" = "OK" ]; then ok=$((ok+1)); else ok=0; fi
  [ $ok -lt 3 ] && sleep 300
done
B=$(head -c 40 $O/GO_STAGE1 | awk '{print $1}'); B=${B:-ntv2}
case "$B" in ntv2|pc2) ;; *) q bad_go_flag ",\"content\":\"$B\""; exit 1;; esac
export M34_BACKEND=$B
if [ "$B" = pc2 ]; then
  [ -f $HOME/coding-rna-reorg/experiments/M35-PLANTCAD2-LORA/scripts/pc2_lora.py ] || { q pc2_lora_module_missing; exit 1; }
  export M34_LORA=on
fi
source ./m34_env.sh
q gate_open ",\"backend\":\"$B\""
if [ ! -f $OUTROOT/${ARMPFX}smoke/SMOKE_OK ]; then
  bash m34_smoke.sh || { q smoke_failed; exit 1; }
fi
q smoke_ok
for ARM in T2 T10; do
  bash m34_stage.sh $ARM || { q arm_failed ",\"arm\":\"$ARM\""; exit 1; }
  q arm_done ",\"arm\":\"$ARM\""
done
touch $O/M34_STAGE1_DONE
q stage1_done
