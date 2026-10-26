"""M34 P2b inference: `m34_infer.py` (imported, not edited) with the M30 evaluation-set-B contigs (Dryas octopetala / Eschscholzia californica) made readable.
`m34_infer` reads sequences only through `T.register()` (M29 helpers + the M34 validation pseudo-contigs), which does not know the M30 contigs.  M33 solved the same
problem with `M30.register()` called at module level (`m33_infer.py`); here `T.register` is wrapped so that every call (main process and spawned decode workers)
also runs `M30.register()`.  The patch sits at module level, so the spawned workers (which re-import this script as `__mp_main__`) see it as well.

Usage (identical to m34_infer.py):  M34_BACKEND=pc2 m34_p2b_infer.py --model model_N.pt --regions NC_141367.1,... --out DIR [--decode-workers 4] [--cpus ...] [--keep-head]
Probe (CPU only, no model): m34_p2b_infer.py --probe CONTIG[,CONTIG...]   prints length / sequence checksum seen by the main process and by a spawned decode worker.
"""
import hashlib
import json
import multiprocessing as mp
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1]                                    # .../experiments

import m34_train_common as T                             # first: loads the backend adapter (pc2: the M32 `a1_common` shim) before the M29-A1 original can shadow it

_M30 = str(EXP / 'M30-HOLDOUT-ADAPT' / 'scripts')
if _M30 not in sys.path:
    sys.path.insert(0, _M30)
import m30_common as M30                                 # noqa: E402  (imports the already loaded `a1_common`; puts the M29-A1/A2 script dirs at the head of sys.path)

for _p in (str(T.ENC_SCRIPTS), str(HERE)):               # keep the backend's script dirs first, as m34_train_common set them up
    if _p in sys.path:
        sys.path.remove(_p)
    sys.path.insert(0, _p)

_register_m34 = T.register


def _register_all():
    """M34 registration (validation pseudo-contigs; patched reader / reference builder) + M30 registration (SEQ2SPECIES / LENGTH of the B contigs)."""
    MC = _register_m34()
    M30.register()
    return MC


T.register = _register_all                               # m34_infer calls `T.register()` at run time (main, _decode_init, _decode_job) -> sees the wrapper
import m34_infer as X                                    # noqa: E402


def _probe(contig):
    """Runs in a (spawned) process: what that process sees for a contig after the registration used by the decode jobs."""
    MC = T.register()
    s = MC.read_sequence(contig)
    return dict(contig=contig, length=MC.LENGTH[contig], seq_len=len(s), md5=hashlib.md5(s.encode()).hexdigest()[:12], species=MC.SEQ2SPECIES[contig])


def _probe_main(contigs):
    MC = T.register()
    main = [_probe(c) for c in contigs]
    pool = mp.get_context('spawn').Pool(1, initializer=X._decode_init, initargs=([],))
    worker = pool.map(_probe, contigs)
    pool.close(); pool.join()
    print(json.dumps(dict(main=main, worker=worker, same=main == worker, backend=T.BACKEND)))


if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[1] == '--probe':
        _probe_main(sys.argv[2].split(','))
    else:
        X.main()
