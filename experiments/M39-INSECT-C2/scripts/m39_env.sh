# sourced by the M39 shell drivers (copy of m34_env.sh): encoder backend, python of the matching env, encoder experiment dir, arm output prefix, output root.
# Backend switch = env M39_BACKEND (ntv2 [default] | pc2; decision pc2-insect-transfer, 2026-11-11).  Deliberately NOT M34_BACKEND: the ont queue
# m34_after_p6.sh exports M34_BACKEND=pc2 and then calls m39_run_all.sh, which must stay NT-v2.  M34_BACKEND is re-exported here from M39_BACKEND only
# (the python side reads M39_BACKEND, see m39_train_common.BACKEND).  ntv2: protocol 29 arms in <OUTROOT>/<arm>;  pc2: PlantCAD2-Small + M35 LoRA in <OUTROOT>/pc2_<arm>;
# onehot: protocol 32 floor arm L (no language model; six-base one-hot tokens into the M32 head) in <OUTROOT>/oh_<arm>.
export M39_BACKEND=${M39_BACKEND:-ntv2}
export M34_BACKEND=$M39_BACKEND
if [ "$M39_BACKEND" = pc2 ]; then PY=$HOME/envs/plantcad2/bin/python; ENC_EXP=M32-C1-PLANTCAD2; ARMPFX=pc2_
elif [ "$M39_BACKEND" = ntv2 ]; then PY=$HOME/envs/generanno/bin/python; ENC_EXP=M31-C1-NTV2; ARMPFX=
elif [ "$M39_BACKEND" = onehot ]; then PY=$HOME/envs/plantcad2/bin/python; ENC_EXP=M32-C1-PLANTCAD2; ARMPFX=oh_; export M39_LORA=none   # protocol 32 arm L: no encoder
else echo "M39_BACKEND must be ntv2, pc2 or onehot (got $M39_BACKEND)" >&2; exit 2; fi
export PY ENC_EXP ARMPFX
OUTROOT=$HOME/coding-rna/outputs/M39-INSECT-C2
