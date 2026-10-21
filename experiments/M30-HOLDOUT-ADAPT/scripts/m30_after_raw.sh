#!/bin/bash
# as each contig's Carbon-A raw probabilities appear: build the fixed offset (all contigs) and run the Carbon-A native post-processing (B contigs)
cd $HOME/coding-rna-reorg/experiments/M30-HOLDOUT-ADAPT/scripts
O=$HOME/coding-rna/outputs/M30-HOLDOUT-ADAPT
declare -A ORDER
for SP in dryas_octopetala eschscholzia_californica; do
  for C in $(ls $HOME/coding-rna/data/m1_screen/$SP/inputs | sed 's/.fasta//'); do echo "$SP $C"; done
done > /tmp/m30_contigs.txt
while read SP C; do
  while ! grep -q "done $C" $O/carbonA_raw/log.txt 2>/dev/null; do sleep 30; done
  taskset -c 18-23 $HOME/envs/generanno/bin/python m30_offsets.py $C >> $HOME/coding-rna/logs/m30-offsets.out 2>> $HOME/coding-rna/logs/m30-offsets.err
  case $C in *__val20M) ;; esac
  if python3 -c "import sys,json;s=json.load(open('$HOME/coding-rna/data/m1_screen/$SP/splits.json'));sys.exit(0 if '$C' in s['B'] else 1)"; then
    bash m30_native.sh $SP $C >> $HOME/coding-rna/logs/m30-native.out 2>&1
  fi
  echo "$(date -Is) prepared $C" >> $O/prep_log.txt
done < /tmp/m30_contigs.txt
echo "$(date -Is) ALL_PREPARED" >> $O/prep_log.txt
