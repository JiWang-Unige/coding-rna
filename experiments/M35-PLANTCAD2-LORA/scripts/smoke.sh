#!/bin/bash
# GPU smoke for M35 (run only after the M33 gate): (1) real-model LoRA unit tests, (2) 20 training steps (batch 8, 12 DP workers) with checkpoint, (3) end-to-end inference +
# internal-validation scoring of the smoke checkpoint on a 400-kb held-out region (checks GFF3 write / parse), (4) SMOKE_DONE summary.  Does NOT start the full run.
O=$HOME/coding-rna/outputs/M35-PLANTCAD2-LORA; S=$O/smoke
PY=$HOME/envs/plantcad2/bin/python
E=$HOME/coding-rna-reorg/experiments/M35-PLANTCAD2-LORA
mkdir -p $S $E/results
cd $E/tests
echo "{\"event\":\"gpu_tests_start\",\"time\":\"$(date -Is)\"}" >> $S/smoke.jsonl
HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 PYTHONPATH=$HOME/envs/m32_extra taskset -c 18-23 $PY -m pytest -q -x test_lora_gpu.py --tb=short > $S/gpu_tests.log 2>&1
echo "{\"event\":\"gpu_tests_done\",\"rc\":$?,\"time\":\"$(date -Is)\"}" >> $S/smoke.jsonl
cd $E/scripts
echo "{\"event\":\"train_smoke_start\",\"time\":\"$(date -Is)\"}" >> $S/smoke.jsonl
bash c1_supervise_lora.sh m35-smoke $S 1536 20 --keep-every 20
echo "{\"event\":\"train_smoke_done\",\"time\":\"$(date -Is)\"}" >> $S/smoke.jsonl
if [ -f $S/model_000020.pt ]; then
  bash c1_infer_retry.sh $S/infer $S/model_000020.pt NC_003071.7:5000000-5400000 4
  taskset -c ${M35_SCORE_CPUS:-18-23} $PY c1_val.py $S/infer 20 NC_003071.7:5000000-5400000 --tag smoke_region >> $S/val.log 2>&1
fi
echo "{\"event\":\"smoke_done\",\"time\":\"$(date -Is)\"}" >> $S/smoke.jsonl
date -Is > $S/SMOKE_DONE
