#!/bin/bash
# M37 protocol §3: minimap2 -cx asm10 --cs  A_chr (target)  B_chr (query);  CPU only, cores 25-31, nice.
O=$HOME/coding-rna/outputs/M37-JOINT-DECODE-SMOKE/align; mkdir -p $O
MM=$HOME/envs/m37/bin/minimap2
A=$HOME/coding-rna/data/m37/align/NC_089044.1.fa
B=$HOME/coding-rna/data/m37/oryza_sativa_9311/inputs/CM000135.2.fasta
echo "minimap2 $($MM --version) $(date -Is) : $MM -cx asm10 --cs -t 4 $A $B" > $O/command.txt
nice -n 10 taskset -c 25-31 $MM -cx asm10 --cs -t 4 $A $B > $O/A_vs_B.paf 2> $O/minimap2.log
echo "rc=$? $(date -Is) lines=$(wc -l < $O/A_vs_B.paf)" >> $O/command.txt
