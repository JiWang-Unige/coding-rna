# sourced by the M40 shell drivers: PlantCAD2 backend only (M34 pc2 env / python / output prefix); arm outputs live next to the M34 pc2 arms
export M34_BACKEND=pc2
source $HOME/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts/m34_env.sh
M40S=$HOME/coding-rna-reorg/experiments/M40-SEMISUP-MLM/scripts
M34S=$HOME/coding-rna-reorg/experiments/M34-NTV2-SCALEUP/scripts
OUT40=$HOME/coding-rna/outputs/M40-SEMISUP-MLM
export M40S M34S OUT40
