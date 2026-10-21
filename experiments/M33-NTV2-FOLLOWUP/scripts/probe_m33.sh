#!/bin/bash
# usage: probe_m33.sh  -> line1: "<n_infer> <n_postq> <gpu_done 0/1> <n_err>" ; then compact status
O=$HOME/coding-rna/outputs/M33-NTV2-FOLLOWUP
ni=$(pgrep -fc '[m]33_gpu_pipeline|[m]33_infer.py'); np=$(pgrep -fc '[m]33_post_queue'); g=0; [ -e $O/gpu_inference.done ] && g=1
ne=$(cat $HOME/coding-rna/logs/m33-infer-*.err 2>/dev/null | grep -c -E 'Traceback|Killed|CUDA out of memory')
echo "$ni $np $g $ne"
tail -n 2 $O/interruptions.jsonl 2>/dev/null | cut -c1-200
for f in $O/postq_*.jsonl; do [ -e $f ] && { echo "$f: done=$(grep -c '"done"' $f) fail=$(grep -c '"fail"' $f) last=$(tail -n1 $f | cut -c1-150)"; }; done
ls $O/*/chains 2>/dev/null | grep -c json
nvidia-smi --query-gpu=memory.used,utilization.gpu --format=csv,noheader
