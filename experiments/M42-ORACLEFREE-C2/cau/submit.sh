#!/bin/bash
# submit one M42 job: submit.sh <config name> <job name> <time limit> [dependency job id]   (one node, one 4090, explicit time limit)
umask 077; G=$HOME/.cache/fontconfig-x2/g10; mkdir -p $G/logs
DEP=(); [ -n "${4:-}" ] && DEP=(--dependency=afterany:$4)
sbatch --parsable -J "$2" -p ChaoBaiHe --nodes=1 --gres=gpu:1 --cpus-per-task=${NCPU:-10} --mem=60G --time="$3" "${DEP[@]}" -D $G \
  -o $G/logs/$2.%j.out -e $G/logs/$2.%j.out $G/run.sh $G/cfg/$1.env
