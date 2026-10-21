#!/bin/bash
# GPU smoke (<= 20 steps) -- runs only after the gate (m34_gate.py OK).  Writes ~/coding-rna/outputs/M34-NTV2-SCALEUP/smoke/{summary.json, SMOKE_OK | SMOKE_FAIL}.
#  (1) T10 on a stress schedule (genes with the most exons / longest CDS span of every species + real draws): 10 steps, stop, RESUME, 10 more steps (checkpoint-resume test)
#  (2) T10_k50U: 3 steps (partial-label treatment U on the GPU path);  T10_k50F: 3 steps (treatment F)
#  (3) end-to-end inference + validation scoring of the 20-step model on a 300-kb slice of a new-species validation pseudo-contig
CPUS=${M34_CPUS:-2,4,6,8,10,12,14,16-24}
mkdir -p $O $HOME/coding-rna/logs
cd $HOME/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts
source ./m34_env.sh
O=$OUTROOT/${ARMPFX}smoke
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONUNBUFFERED=1
rm -f $O/SMOKE_OK $O/SMOKE_FAIL
fail() { echo "{\"event\":\"smoke_fail\",\"why\":\"$1\",\"time\":\"$(date -Is)\"}" >> $O/smoke_log.jsonl; touch $O/SMOKE_FAIL; exit 1; }
LOG=$HOME/coding-rna/logs/m34-smoke-$(date +%Y%m%d-%H%M%S)
$PY m34_smoke_schedule.py T10 $O/schedule_smoke_T10.jsonl >> $LOG.out 2>> $LOG.err || fail schedule
T() { # arm outdir stop [schedule]
  A=$1; D=$2; S=$3; SCH=$4; [ -n "$SCH" ] && SCH="--schedule $SCH"
  taskset -c $CPUS $PY m34_train_lora.py run --arm $A --out $D --steps 3072 --stop-at $S --batch 8 --workers 12 --cpus $CPUS --ckpt-every 5 --keep-every 10 --resume ${M34_LORA:+--lora $M34_LORA} $SCH >> $LOG.out 2>> $LOG.err < /dev/null
}
T T10 $O/t10 10 $O/schedule_smoke_T10.jsonl || fail train_t10_a
T T10 $O/t10 20 $O/schedule_smoke_T10.jsonl || fail train_t10_resume
T T10_k50U $O/k50U 3 || fail train_k50U
T T10_k50F $O/k50F 3 || fail train_k50F
[ -f $O/t10/model_000020.pt ] || fail no_model
bash m34_infer_retry.sh $O/infer_slice $O/t10/model_000020.pt NC_037298.2__val8M:0-300000 2 || fail infer
taskset -c 18-23 $PY m34_val.py $O/infer_slice 20 NC_037298.2__val8M:0-300000 --tag smoke --arm T10 --metrics $O/metrics.jsonl >> $LOG.out 2>> $LOG.err || fail val
$PY - <<PYEOF
import json, statistics, glob, os
O = os.path.expanduser('~/coding-rna/outputs/M34-NTV2-SCALEUP/') + os.environ.get('ARMPFX', '') + 'smoke'
def rows(d): return [json.loads(l) for l in open(f'{O}/{d}/training.jsonl')]
t10 = rows('t10')
steps = [r['step'] for r in t10]
st = [r['step_s'] for r in t10[1:]]
done = json.load(open(f'{O}/t10/DONE_20.json'))
inf = sum(r['infeasible'] for r in t10); skip = sum(r['skipped_all_unknown'] for r in t10)
resumed = steps == list(range(1, 21))
out = dict(t10_steps=len(t10), resumed_continuous=resumed, step_s_median=statistics.median(st), step_s_mean=sum(st) / len(st), step_s_max=max(st), step_s_first=t10[0]['step_s'],
           infeasible=inf, skipped_all_unknown=skip, draws=len(t10) * 8, bwd_s_median=statistics.median(r['bwd_s'] for r in t10), dp_sum_s_median=statistics.median(r['dp_sum_s'] for r in t10),
           loss_first=t10[0]['batch_loss'], loss_last=t10[-1]['batch_loss'], peak_gpu_alloc_GB=done['peak_gpu_alloc_GB'], peak_gpu_reserved_GB=done['peak_gpu_reserved_GB'],
           k50U_steps=len(rows('k50U')), k50F_steps=len(rows('k50F')), k50U_loss=[r['batch_loss'] for r in rows('k50U')], k50F_loss=[r['batch_loss'] for r in rows('k50F')],
           k50U_infeasible=sum(r['infeasible'] for r in rows('k50U')), k50F_infeasible=sum(r['infeasible'] for r in rows('k50F')),
           slice_val=json.loads(open(f'{O}/metrics.jsonl').readlines()[-1]))
est = out['step_s_median'] * 3072 / 3600
out['projected_train_hours_3072_steps_at_median'] = est
ok = resumed and inf <= 16 and out['k50U_steps'] == 3 and out['k50F_steps'] == 3 and all(l == l and l < 5 for l in [out['loss_last']]) and out['step_s_median'] < 60
json.dump(out, open(f'{O}/summary.json', 'w'), indent=1)
open(f'{O}/' + ('SMOKE_OK' if ok else 'SMOKE_FAIL'), 'w').write(json.dumps(out))
print(json.dumps(out))
PYEOF
[ -f $O/SMOKE_OK ] || exit 1
echo "{\"event\":\"smoke_ok\",\"time\":\"$(date -Is)\"}" >> $O/smoke_log.jsonl
