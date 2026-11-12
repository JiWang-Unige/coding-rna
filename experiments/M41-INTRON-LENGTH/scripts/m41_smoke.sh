#!/bin/bash
# M41 GPU smoke (protocol 33 section 5; ~20-30 min GPU, run by the maintenance session when the queue allows): 20 training steps of T2_ILS on the T2
# schedule into a SEPARATE directory, M41 inference + M34 validation scoring on a 2-Mb slice of the Ath validation chromosome, and checks.
# Never touches DEV or set B.  Output: ~/coding-rna/outputs/M41-INTRON-LENGTH/smoke/ ; SMOKE_OK on success.
EXPS=$HOME/coding-rna-reorg/experiments
S41=$EXPS/M41-INTRON-LENGTH/scripts; S34=$EXPS/M34-NTV2-SCALEUP/scripts
O=$HOME/coding-rna/outputs/M41-INTRON-LENGTH/smoke
LOGS=$HOME/coding-rna/logs
mkdir -p $O $LOGS
export M34_BACKEND=pc2 M34_LORA=on HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONUNBUFFERED=1 \
       PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
cd $S34 && source ./m34_env.sh
CPUS=${M41_CPUS:-2,4,6,10,12,14,16-24}
REG=NC_003071.7:0-2000000
taskset -c $CPUS $PY $S41/m41_train.py run --arm T2 --label T2_ILS_smoke --out $O/train --steps 3072 --stop-at 20 --batch 8 --workers 12 --cpus $CPUS \
  > $O/train.out 2> $O/train.err < /dev/null || { echo "train failed"; tail -5 $O/train.err; exit 1; }
taskset -c $CPUS $PY $S41/m41_infer.py --model $O/train/model_000020.pt --regions $REG --out $O/val --decode-workers 2 --cpus $CPUS \
  > $O/infer.out 2> $O/infer.err < /dev/null || { echo "infer failed"; tail -5 $O/infer.err; exit 1; }
taskset -c 18-23 $PY m34_val.py $O/val 20 $REG --tag smoke --arm T2_ILS_smoke --metrics $O/metrics.jsonl > $O/val.out 2>&1 || { echo "val failed"; exit 1; }
$PY - $O <<'PY' || exit 1
import json, sys, statistics as st
O = sys.argv[1]
rows = [json.loads(l) for l in open(f'{O}/train/training.jsonl')]
assert len(rows) == 20 and all(r['infeasible'] == 0 for r in rows), 'training rows'
t2 = [json.loads(l) for l in open('/home/ji/coding-rna/outputs/M34-NTV2-SCALEUP/pc2_T2/training.jsonl')][:20]
inst = sum('m41_decoder_installed' in l for l in open(f'{O}/infer.err'))
val = json.loads(open(f'{O}/metrics.jsonl').read().splitlines()[-1])
out = dict(steps=len(rows), theta_first=rows[0]['theta'], theta_last=rows[-1]['theta'], pi_short=rows[-1]['pi_short'],
           short_mean_bp=rows[-1]['short_mean_bp'], long_mean_bp=rows[-1]['long_mean_bp'],
           median_step_s=st.median(r['step_s'] for r in rows[3:]), median_step_s_T2_first20=st.median(r['step_s'] for r in t2[3:]),
           median_dp_sum_s=st.median(r['dp_sum_s'] for r in rows[3:]), median_dp_sum_s_T2_first20=st.median(r['dp_sum_s'] for r in t2[3:]),
           loss_first=rows[0]['batch_loss'], loss_T2_first=t2[0]['batch_loss'], decoder_installs=inst, smoke_val=val)
assert inst >= 1, 'M41 decoder was not installed in the decode workers'
assert abs(rows[0]['batch_loss'] - t2[0]['batch_loss']) < .5 * max(1, t2[0]['batch_loss']), 'first-step loss far from T2'
json.dump(out, open(f'{O}/smoke_summary.json', 'w'), indent=1)
print(json.dumps(out))
PY
touch $O/SMOKE_OK
