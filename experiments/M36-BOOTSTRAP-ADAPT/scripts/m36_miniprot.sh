#!/bin/bash
# miniprot alignment of the 10-species protein queries to the LABEL-SET A contigs of one target species (protocol 22 §2 S2; §7 item 2).
# Evaluation-set B contigs are never concatenated, indexed or aligned.  CPU: cores 24-31, nice 10, one process (8 threads).
# usage: m36_miniprot.sh <dryas_octopetala|eschscholzia_californica> [threads=8] [query.faa] [output tag, e.g. _pilot]
set -euo pipefail
SP=$1; T=${2:-8}; TAG=${4:-}
D=~/coding-rna/data/m36; Q=${3:-$D/proteins/queries_10sp.faa}
MP=~/envs/miniprot/bin/miniprot
M1=~/coding-rna/data/m1_screen/$SP
case $SP in
  dryas_octopetala) A="NC_141366.1 NC_141368.1 NC_141370.1 NC_141372.1 NC_141374.1__val20M";;
  eschscholzia_californica) A="NW_028238198.1 NW_028238420.1 NW_028238472.1__val20M";;
  *) echo "unknown species $SP"; exit 1;;
esac
mkdir -p $D/genomes_A $D/miniprot_raw
G=$D/genomes_A/${SP}_A.fa
if [ ! -s $G ]; then
  for c in $A; do cat $M1/inputs/$c.fasta; done > $G.tmp && mv $G.tmp $G
fi
grep '^>' $G
IDX=$D/genomes_A/${SP}_A.mpi
if [ ! -s $IDX ]; then
  taskset -c 24-31 nice -n 10 $MP -t $T -d $IDX.tmp $G && mv $IDX.tmp $IDX
fi
OUT=$D/miniprot_raw/${SP}${TAG}.gff.gz
LOG=$D/miniprot_raw/${SP}${TAG}.log
date > $LOG
# -I: max intron 3.6*sqrt(refLen); --outc 0.9 / --outs 0.97 / --outn 5: keep a few near-best alignments per protein (the S2 filter is applied afterwards)
taskset -c 24-31 nice -n 10 $MP -t $T -I --gff --outc=0.9 --outs=0.97 --outn=5 -P ${SP:0:3} $IDX $Q 2>>$LOG | gzip -1 > $OUT.tmp && mv $OUT.tmp $OUT
date >> $LOG
echo DONE >> $LOG
