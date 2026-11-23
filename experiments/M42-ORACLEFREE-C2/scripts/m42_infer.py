"""M42 inference on the evaluation chromosomes E: `m34_infer.py` (imported, not edited) with the target-species chromosomes registered in the frozen m29_common
reader.  `T.register` (which needs the M34 species table) is replaced by the M42 registration in this process and in the spawned decode workers (module-level patch,
re-executed when the workers re-import this script as `__mp_main__`).  The CRF decoder, stitching and GFF export are the M34 / M31 code.
Usage (identical to m34_infer.py):  M34_BACKEND=pc2 m42_infer.py --model model_N.pt --regions NC_029683.1,... --out DIR [--decode-workers 4] [--cpus ...] [--keep-head]
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'M34-NTV2-SCALEUP' / 'scripts'))
import m34_train_common as T                             # first: loads the backend adapter
if str(HERE) not in sys.path:
    sys.path.append(str(HERE))
import m42_common as Z                                   # noqa: E402


def _register():
    import m29_common as MC
    return Z.register(MC)


T.register = _register
import m34_infer as X                                    # noqa: E402

if __name__ == '__main__':
    X.main()
