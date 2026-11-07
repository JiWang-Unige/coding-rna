#!/bin/bash
# PANEL-P6 on win3090 (RTX 3090, WSL2 Ubuntu 24.04): the three native tools only (Carbon-A and Vipsania stay on ont).  Same assets (SIF files byte-identical to Baobab's,
# sha256 checked in assets.sha256), same commands as p6_native.sh; outputs go back to ont $P6/native/<tool>/<key>/ for scoring.
# Gate additions here: TE's win3090 queue drivers (Vipsania chain) block while alive; 24 threads -> CPU set 2-21.
# Launch (on win3090): cd ~/coding-rna-reorg/experiments/PANEL-P6/scripts && (setsid nohup bash p6_win3090.sh > ~/coding-rna/logs/p6-queue.log 2>&1 < /dev/null &)
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export P6_APPT=$HOME/envs/apptainer_p6/bin/apptainer
export P6_TE_DRIVERS=${P6_TE_DRIVERS:-'scripts/[e]xperiments/'}            # any TE driver or job started from its repo (cwd ab-initio-TE): user rule, no TE task running
export P6_TE_CLAIM=${P6_TE_CLAIM:-$HOME/te/ab-initio-TE/outputs/TE_GPU_CLAIM}
export P6_CPUS=${P6_CPUS:-2-21}
export P6_STEPS=${P6_STEPS:-SMOKE ANNEVO Tiberius Helixer}
export P6_SMOKE_TOOLS=${P6_SMOKE_TOOLS:-ANNEVO Tiberius Helixer}
exec bash $HERE/p6_queue.sh
