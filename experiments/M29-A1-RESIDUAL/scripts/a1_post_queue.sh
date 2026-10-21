#!/bin/bash
# rolling posterior jobs: waits for each contig/strand chain file; P parallel processes
RUN=$1; P=${2:-6}; shift 2
cd $HOME/coding-rna-reorg/experiments/M29-A1-RESIDUAL/scripts
one() { RUN=$1; s=$2; st=$3
  while [ ! -f $RUN/chains/${s}_${st}.json ]; do sleep 20; done
  for a in 1 2 3; do
    taskset -c 18-31 $HOME/envs/generanno/bin/python a1_posterior.py $RUN $s $st >> $HOME/coding-rna/logs/a1-post-$(basename $RUN).out 2>> $HOME/coding-rna/logs/a1-post-$(basename $RUN).err < /dev/null && break
    echo "{\"event\":\"posterior_retry\",\"seq\":\"$s\",\"strand\":\"$st\",\"attempt\":$a,\"time\":\"$(date -Is)\"}" >> $RUN/interruptions.jsonl; sleep 20
  done; }
export -f one
for s in "$@"; do for st in + -; do echo "$RUN $s $st"; done; done | xargs -P $P -L 1 bash -c 'one $0 $1 $2'
