#!/bin/bash
# End-to-end dry run of `vip_score.py dev` and `setb` with STAND-IN inputs that are already-scored non-Vipsania predictions (never Vipsania):
#   DEV   ft = GENERanno fresh-joint GFFs, pre = Carbon-A (unfiltered) GFFs ;  set B  ft = the T10 set-B GFFs, pre = the T2 set-B GFFs (concatenated per genome)
# Checks that the result JSONs have the expected structure and numbers.  Output: $HOME/coding-rna/outputs/M38-VIPSANIA/_score_dryrun (throw-away).   Run on ont.
O=$HOME/coding-rna/outputs/M38-VIPSANIA/_score_dryrun
M34=$HOME/coding-rna/outputs/M34-NTV2-SCALEUP
SC=$HOME/coding-rna-reorg/experiments/M38-VIPSANIA/scripts
rm -rf $O; mkdir -p $O/in
FJ=$HOME/coding-rna/outputs/M28-FULL-GRAMMAR-JOINT-CONTROLS-DEV-RAW-R1/FRESH_JOINT24576
CA=$HOME/coding-rna/outputs/M29-CARBON-A-GAP/G0
cat $FJ/NC_003074.8/prediction/full_grammar.gff3 > $O/in/ft_ara.gff3
cat $FJ/NC_089041.1/prediction/full_grammar.gff3 > $O/in/ft_rice.gff3
cp $CA/ara_pred/prediction.gff $O/in/pre_ara.gff3
cp $CA/rice_pred/prediction.gff $O/in/pre_rice.gff3
for arm in T10 T2; do
  cat $M34/p2b/$arm/NC_141367.1/prediction/full_grammar.gff3 $M34/p2b/$arm/NC_141369.1/prediction/full_grammar.gff3 $M34/p2b/$arm/NC_141371.1/prediction/full_grammar.gff3 $M34/p2b/$arm/NC_141373.1/prediction/full_grammar.gff3 > $O/in/${arm}_dryas.gff3
  cat $M34/p2b/$arm/NW_028238309.1/prediction/full_grammar.gff3 $M34/p2b/$arm/NW_028238461.1/prediction/full_grammar.gff3 $M34/p2b/$arm/NW_028238483.1/prediction/full_grammar.gff3 > $O/in/${arm}_esch.gff3
done
export CUDA_VISIBLE_DEVICES=""
cd $SC
nice -n 19 $HOME/envs/generanno/bin/python -I vip_score.py dev --out $O/dev --ft arabidopsis_thaliana=$O/in/ft_ara.gff3 oryza_sativa=$O/in/ft_rice.gff3 --pre arabidopsis_thaliana=$O/in/pre_ara.gff3 oryza_sativa=$O/in/pre_rice.gff3 || exit 1
nice -n 19 $HOME/envs/generanno/bin/python -I vip_score.py setb --out $O/setb --ft dryas_octopetala=$O/in/T10_dryas.gff3 eschscholzia_californica=$O/in/T10_esch.gff3 --pre dryas_octopetala=$O/in/T2_dryas.gff3 eschscholzia_californica=$O/in/T2_esch.gff3 || exit 1
nice -n 19 $HOME/envs/generanno/bin/python -I vip_score.py dev --out $O/dev --ft arabidopsis_thaliana=$O/in/ft_ara.gff3 oryza_sativa=$O/in/ft_rice.gff3 2>&1 | tail -2    # must refuse: result exists
