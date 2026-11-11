#!/bin/bash
# M40 GPU smoke (~1 GPU h; run only when the maintainer's queue gives it the GPU).  Writes ~/coding-rna/outputs/M40-SEMISUP-MLM/smoke/{summary.json, SMOKE_OK | SMOKE_FAIL}.
#  (a) lambda = 0 identity: M40_k25U_MLM with the protocol config, 3 steps -> batch_loss must match the M34 pc2_T10_k25U run's first 3 steps (same schedule,
#      seeds, init; lambda = 0 before step 307 means the exact M34 code path) within 2 % (bf16 / non-deterministic GPU kernels)
#  (b) MLM on from step 0 (lambda = 0.1, no ramp) on M40_k25U_MLM with separate CRF / MLM LoRA gradient norms: 6 steps, stop, RESUME, 4 more steps;
#      step-1 MLM accuracy of the untouched pre-trained encoder (LoRA B = 0) must be >= 0.40 (report 16: Small 0.55-0.56 at 15 % all-[MASK]) -> LM-head path is right
#  (c) M40_k0 with MLM on from step 0: 3 steps (new-species draws are MLM-only, legacy draws CRF)
#  (d) inference + validation scoring of the (b) 10-step model on a 300-kb slice (checkpoint is readable by the unchanged M34 inference)
source $HOME/coding-rna-reorg/experiments/M40-SEMISUP-MLM/scripts/m40_env.sh
CPUS=${M34_CPUS:-2,4,6,8,10,12,14,16-24}
O=$OUT40/smoke
mkdir -p $O $HOME/coding-rna/logs
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONUNBUFFERED=1
rm -f $O/SMOKE_OK $O/SMOKE_FAIL
fail() { echo "{\"event\":\"smoke_fail\",\"why\":\"$1\",\"time\":\"$(date -Is)\"}" >> $O/smoke_log.jsonl; touch $O/SMOKE_FAIL; exit 1; }
LOG=$HOME/coding-rna/logs/m40-smoke-$(date +%Y%m%d-%H%M%S)
cd $M40S
[ -f $OUT40/repeats/summary.json ] && [ -f $OUT40/schedules/schedule_M40_k0.jsonl ] || fail cpu_prep_missing
Tr() { # arm outdir stop [extra args]
  A=$1; D=$2; S=$3; shift 3
  taskset -c $CPUS $PY m40_train.py run --arm $A --out $D --steps 3072 --stop-at $S --batch 8 --workers 12 --cpus $CPUS --ckpt-every 2 --keep-every 10 --resume "$@" >> $LOG.out 2>> $LOG.err < /dev/null
}
Tr M40_k25U_MLM $O/eq 3 || fail train_eq
Tr M40_k25U_MLM $O/on 6 --mlm-start-frac 0 --mlm-ramp-frac 0 --diag-grad || fail train_on_a
Tr M40_k25U_MLM $O/on 10 --mlm-start-frac 0 --mlm-ramp-frac 0 --diag-grad || fail train_on_resume
Tr M40_k0 $O/k0 3 --mlm-start-frac 0 --mlm-ramp-frac 0 || fail train_k0
[ -f $O/on/model_000010.pt ] || fail no_model
cd $M34S
bash m34_infer_retry.sh $O/infer_slice $O/on/model_000010.pt NC_037298.2__val8M:0-300000 2 || fail infer
taskset -c 18-23 $PY m34_val.py $O/infer_slice 10 NC_037298.2__val8M:0-300000 --tag smoke --arm M40_k25U_MLM --metrics $O/metrics.jsonl >> $LOG.out 2>> $LOG.err || fail val
$PY - <<PYEOF
import json, os, statistics
O = os.path.expanduser('~/coding-rna/outputs/M40-SEMISUP-MLM/smoke')
K25 = os.path.expanduser('~/coding-rna/outputs/M34-NTV2-SCALEUP/pc2_T10_k25U/training.jsonl')
def rows(d): return [json.loads(l) for l in open(f'{O}/{d}/training.jsonl')]
eq, on, k0 = rows('eq'), rows('on'), rows('k0')
ref = {}
for l in open(K25):
    r = json.loads(l); ref.setdefault(r['step'], r)          # first occurrence of each step (the M34 run was resumed later; steps 1-3 ran once)
eq_rel = [abs(r['batch_loss'] - ref[r['step']]['batch_loss']) / max(1e-9, abs(ref[r['step']]['batch_loss'])) for r in eq]
st = [r['step_s'] for r in on[1:]]
out = dict(eq_steps=len(eq), eq_batch_loss=[r['batch_loss'] for r in eq], ref_batch_loss=[ref[r['step']]['batch_loss'] for r in eq], eq_max_rel_diff=max(eq_rel),
           eq_lambda=[r['mlm_lambda'] for r in eq],
           on_steps=[r['step'] for r in on], on_resumed_continuous=[r['step'] for r in on] == list(range(1, 11)),
           on_mlm_acc_step1=on[0]['mlm_acc'], on_mlm_ce=[r['mlm_ce'] for r in on], on_mlm_positions=[r['mlm_positions'] for r in on],
           on_batch_loss=[r['batch_loss'] for r in on], on_mlm_loss=[r['mlm_loss'] for r in on],
           diag_grad_norm_crf=[r.get('diag_grad_norm_crf') for r in on], diag_grad_norm_mlm=[r.get('diag_grad_norm_mlm') for r in on],
           k0_crf_draws=[r['crf_draws'] for r in k0], k0_mlm_draws=[r['mlm_draws'] for r in k0], k0_mlm_acc=[r['mlm_acc'] for r in k0],
           infeasible=sum(r['infeasible'] for r in eq + on + k0), step_s_median_on=statistics.median(st), step_s_k0=[r['step_s'] for r in k0],
           step_s_eq=[r['step_s'] for r in eq], slice_val=json.loads(open(f'{O}/metrics.jsonl').readlines()[-1]))
out['grad_ratio_mlm_over_crf_median'] = statistics.median(m / c for m, c in zip(out['diag_grad_norm_mlm'], out['diag_grad_norm_crf']) if m and c)
out['projected_train_hours_3072_steps'] = out['step_s_median_on'] * 3072 / 3600
ok = (out['eq_max_rel_diff'] < 0.02 and out['on_resumed_continuous'] and (out['on_mlm_acc_step1'] or 0) >= 0.40
      and all(x is not None and x == x and x < 5 for x in out['on_mlm_ce']) and all(x > 0 for x in out['k0_mlm_draws'])
      and out['infeasible'] <= 16 and out['step_s_median_on'] < 150)
json.dump(out, open(f'{O}/summary.json', 'w'), indent=1)
open(f'{O}/' + ('SMOKE_OK' if ok else 'SMOKE_FAIL'), 'w').write(json.dumps(out))
print(json.dumps(out))
PYEOF
[ -f $O/SMOKE_OK ] || exit 1
echo "{\"event\":\"smoke_ok\",\"time\":\"$(date -Is)\"}" >> $O/smoke_log.jsonl
