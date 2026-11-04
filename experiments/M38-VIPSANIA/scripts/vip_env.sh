# sourced by the M38-VIPSANIA drivers: the pinned Vipsania environment (own env, own model cache; nothing shared with other projects' Vipsania installs)
# (M38_ROOT / VIP_BIN may be overridden for the stub test of vip_devb.sh only)
export M38=${M38_ROOT:-$HOME/coding-rna/outputs/M38-VIPSANIA}
export M38_SCRIPTS=$HOME/coding-rna-reorg/experiments/M38-VIPSANIA/scripts
export VENV=$HOME/envs/vipsania_m38                  # python 3.12.12, vipsania 1.1.0, tensorflow 2.19.1 (CUDA 12.5 wheels), keras 3.15.1  (env_versions.json)
export VIP=${VIP_BIN:-$VENV/bin/vipsania}
export VPY=$VENV/bin/python
export VIPSANIA_CACHE=$M38/vipsania_cache            # models live in $VIPSANIA_CACHE/models/<model_id>
export MODEL_ID=j9m0cdmk                             # Streptophyta (versions.json 2026-08-20); always pinned by ID, never by clade name
export PYG=$HOME/envs/generanno/bin/python           # the project ruler's environment (scoring only, CPU)
export TF_CPP_MIN_LOG_LEVEL=${TF_CPP_MIN_LOG_LEVEL:-2}
export WANDB_MODE=offline WANDB_SILENT=true WANDB_DISABLED=true     # vipsania imports wandb; nothing may leave the machine
