# sourced by the M34 shell drivers: selects the encoder backend (env M34_BACKEND = ntv2 [default] | pc2), the python of the matching env, and the arm output prefix
export M34_BACKEND=${M34_BACKEND:-ntv2}
if [ "$M34_BACKEND" = pc2 ]; then PY=$HOME/envs/plantcad2/bin/python; ENC_EXP=M32-C1-PLANTCAD2; ARMPFX=pc2_; else PY=$HOME/envs/generanno/bin/python; ENC_EXP=M31-C1-NTV2; ARMPFX=; fi
export PY ENC_EXP ARMPFX
OUTROOT=$HOME/coding-rna/outputs/M34-NTV2-SCALEUP
