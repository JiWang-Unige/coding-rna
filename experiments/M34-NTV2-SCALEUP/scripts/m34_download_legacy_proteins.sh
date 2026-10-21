#!/bin/bash
# Proteomes (RefSeq protein.faa.gz) of the two legacy training species, for the homology audit only (genomes/GFFs already in data/m1_screen).
D=$HOME/coding-rna/data/m34_train/raw/legacy_proteins; mkdir -p $D; cd $D
for BASE in https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/140/825/GCF_034140825.1_ASM3414082v1; do
  S=$(basename $BASE); F=${S}_protein.faa.gz
  [ -s $F ] || curl -fsSL --retry 5 -o $F $BASE/$F
  curl -fsSL --retry 5 -o $S.md5checksums.txt $BASE/md5checksums.txt
  L=$(md5sum $F | cut -d' ' -f1); R=$(grep " \./$F\$" $S.md5checksums.txt | cut -d' ' -f1)
  echo "$L  $F  ncbi=$R  $([ "$L" = "$R" ] && echo OK || echo MISMATCH)  $(date -Is)" | tee -a local_md5.txt
done
