#!/bin/bash
# GPU smoke of M39 (copy of m34_smoke.sh; protocol 29 §5.4; <= 20 steps, NT-v2 backend).  Writes ~/coding-rna/outputs/M39-INSECT-C2/smoke/{summary.json, SMOKE_OK | SMOKE_FAIL}.
#  (1) I10_T on a stress schedule (genes with the most exons / longest CDS span of all 10 species + real draws): 10 steps, stop, RESUME, 10 more steps
#  (2) I10_k50U: 3 steps (treatment U on the GPU path);  I10_k50F: 3 steps (treatment F)
#  (3) end-to-end inference + validation scoring of the 20-step model on a 300-kb slice of the first new species' validation crop
# Gate (unchanged from the M34 NT-v2 setting): resume continuous, infeasible <= 16, k50U / k50F ran 3 steps, final loss finite and < 5, median step time < 60 s.
# Does NOT check the GPU queue: the caller (queue of the coordinating session) must start it only when no other coding-rna GPU job runs.
CPUS=${M39_CPUS:-2,4,6,8,10,12,14,16-24}
cd $HOME/coding-rna-reorg/experiments/M39-INSECT-C2/scripts
source ./m39_env.sh
O=$OUTROOT/smoke
mkdir -p $O $HOME/coding-rna/logs
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONUNBUFFERED=1
rm -f $O/SMOKE_OK $O/SMOKE_FAIL
fail() { echo "{\"event\":\"smoke_fail\",\"why\":\"$1\",\"time\":\"$(date -Is)\"}" >> $O/smoke_log.jsonl; touch $O/SMOKE_FAIL; exit 1; }
LOG=$HOME/coding-rna/logs/m39-smoke-$(date +%Y%m%d-%H%M%S)
$PY m39_smoke_schedule.py I10_T $O/schedule_smoke_I10_T.jsonl >> $LOG.out 2>> $LOG.err || fail schedule
SLICE=$($PY -c "import m39_train_common as T; print(T.val_contig(T.new_keys()[0]) + ':0-300000')" 2>/dev/null)
[ -n "$SLICE" ] || fail slice
T() { # arm outdir stop [schedule]
  A=$1; D=$2; S=$3; SCH=$4; [ -n "$SCH" ] && SCH="--schedule $SCH"
  taskset -c $CPUS $PY m39_train_lora.py run --arm $A --out $D --steps 3072 --stop-at $S --batch 8 --workers 12 --cpus $CPUS --ckpt-every 5 --keep-every 10 --resume $SCH >> $LOG.out 2>> $LOG.err < /dev/null
}
T I10_T $O/t 10 $O/schedule_smoke_I10_T.jsonl || fail train_t_a
T I10_T $O/t 20 $O/schedule_smoke_I10_T.jsonl || fail train_t_resume
T I10_k50U $O/k50U 3 || fail train_k50U
T I10_k50F $O/k50F 3 || fail train_k50F
[ -f $O/t/model_000020.pt ] || fail no_model
bash m39_infer_retry.sh $O/infer_slice $O/t/model_000020.pt $SLICE 2 || fail infer
taskset -c 18-23 $PY m39_val.py $O/infer_slice 20 $SLICE --tag smoke --arm I10_T --metrics $O/metrics.jsonl >> $LOG.out 2>> $LOG.err || fail val
$PY - <<PYEOF
import json, statistics, os
O = os.path.expanduser('~/coding-rna/outputs/M39-INSECT-C2/smoke')
def rows(d): return [json.loads(l) for l in open(f'{O}/{d}/training.jsonl')]
t = rows('t')
steps = [r['step'] for r in t]
st = [r['step_s'] for r in t[1:]]
done = json.load(open(f'{O}/t/DONE_20.json'))
inf = sum(r['infeasible'] for r in t); skip = sum(r['skipped_all_unknown'] for r in t)
resumed = steps == list(range(1, 21))
out = dict(t_steps=len(t), resumed_continuous=resumed, step_s_median=statistics.median(st), step_s_mean=sum(st) / len(st), step_s_max=max(st), step_s_first=t[0]['step_s'],
           infeasible=inf, skipped_all_unknown=skip, draws=len(t) * 8, bwd_s_median=statistics.median(r['bwd_s'] for r in t), dp_sum_s_median=statistics.median(r['dp_sum_s'] for r in t),
           loss_first=t[0]['batch_loss'], loss_last=t[-1]['batch_loss'], peak_gpu_alloc_GB=done['peak_gpu_alloc_GB'], peak_gpu_reserved_GB=done['peak_gpu_reserved_GB'],
           k50U_steps=len(rows('k50U')), k50F_steps=len(rows('k50F')), k50U_loss=[r['batch_loss'] for r in rows('k50U')], k50F_loss=[r['batch_loss'] for r in rows('k50F')],
           k50U_infeasible=sum(r['infeasible'] for r in rows('k50U')), k50F_infeasible=sum(r['infeasible'] for r in rows('k50F')),
           slice_val=json.loads(open(f'{O}/metrics.jsonl').readlines()[-1]))
out['projected_train_hours_3072_steps_at_median'] = out['step_s_median'] * 3072 / 3600
out['step_s_limit'] = 60
ok = resumed and inf <= 16 and out['k50U_steps'] == 3 and out['k50F_steps'] == 3 and all(l == l and l < 5 for l in [out['loss_last']]) and out['step_s_median'] < out['step_s_limit']
json.dump(out, open(f'{O}/summary.json', 'w'), indent=1)
open(f'{O}/' + ('SMOKE_OK' if ok else 'SMOKE_FAIL'), 'w').write(json.dumps(out))
print(json.dumps(out))
PYEOF
[ -f $O/SMOKE_OK ] || exit 1
echo "{\"event\":\"smoke_ok\",\"time\":\"$(date -Is)\"}" >> $O/smoke_log.jsonl
