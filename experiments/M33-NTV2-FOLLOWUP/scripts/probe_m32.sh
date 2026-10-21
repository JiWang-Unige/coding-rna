#!/bin/bash
# prints: "<n_m32_processes> <m32_done 0/1>" then last stage_log / metrics lines of M32.  Path literals are bracketed so this probe never matches itself.
n1=0
for p in /proc/[0-9]*; do
  case "$(readlink $p/cwd 2>/dev/null)" in *M32-C1-PLANTCAD2*) n1=$((n1+1));; esac
done
n2=$(pgrep -fc '[s]tart_after_m31')
O=$HOME/coding-rna/outputs/M32-C1-PLANTCAD[2]
d=0; ls $O/M32_DONE >/dev/null 2>&1 && d=1
echo "$((n1+n2)) $d"
tail -n1 $O/runA/stage_log.jsonl 2>/dev/null | cut -c1-200
tail -n1 $O/runA/metrics.jsonl 2>/dev/null | cut -c1-200
nvidia-smi --query-gpu=memory.used,utilization.gpu --format=csv,noheader
