#!/bin/bash
# Main-panel proteomes (docs/review_20261016/10 §6.1; protein.faa.gz only, ~100 MB) for the training-vs-panel homology audit.  Dryas/Eschscholzia reuse the M30 copies.
D=$HOME/coding-rna/data/m34_train/raw/panel_proteins; mkdir -p $D; cd $D
B=https://ftp.ncbi.nlm.nih.gov/genomes/all
for P in GCF/017/979/445/GCF_017979445.1_ASM1797944v1 GCF/041/222/555/GCF_041222555.1_UWA_Hpros_2 GCF/948/465/745/GCF_948465745.1_dcCheAlbu1.1 \
         GCA/048/537/475/GCA_048537475.1_ASM4853747v1 GCA/057/662/145/GCA_057662145.1_ThtriWBW_v1.0 GCA/027/726/005/GCA_027726005.1_MPIPZ_junEff1m_1 \
         GCA/027/562/975/GCA_027562975.1_MPIPZ_rhyBre1m_1.0 GCA/030/737/835/GCA_030737835.1_ASM3073783v1 GCF/047/967/505/GCF_047967505.1_Dip_SZ_1.0 \
         GCF/054/954/075/GCF_054954075.1_AtrichopodaSantaCruz75HAP2_v2.0; do
  S=$(basename $P); F=${S}_protein.faa.gz; BASE=$B/$P
  [ -s $F ] || nice -n 10 curl -fsSL --retry 5 -o $F $BASE/$F
  curl -fsSL --retry 5 -o $S.md5checksums.txt $BASE/md5checksums.txt
  L=$(md5sum $F | cut -d' ' -f1); R=$(grep " \./$F\$" $S.md5checksums.txt | cut -d' ' -f1)
  echo "$L  $F  ncbi=$R  $([ "$L" = "$R" ] && echo OK || echo MISMATCH)  $(date -Is)" | tee -a local_md5.txt
done
