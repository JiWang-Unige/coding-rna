"""P5 panel inference: `m34_infer.py` (imported, not edited) with the independent-panel main sequences made readable (m34_panel_common.register), same pattern as
`m34_p2b_infer.py`: `T.register` is wrapped at module level so that every call (main process and spawned decode workers) also registers the panel sequences.

Usage (identical to m34_infer.py):  M34_BACKEND=pc2 m34_panel_infer.py --model model_N.pt --regions SEQ,... --out DIR [--decode-workers 4] [--cpus ...] [--keep-head]
Probe (CPU only, no model): m34_panel_infer.py --probe CONTIG[,CONTIG...]   prints length / sequence checksum seen by the main process and by a spawned decode worker.
"""
import hashlib
import json
import multiprocessing as mp
import sys

import m34_train_common as T                             # first: loads the backend adapter (pc2: the M32 `a1_common` shim)
import m34_panel_common as PC

_register_m34 = T.register


def _register_all():
    return PC.register(_register_m34())


T.register = _register_all                               # m34_infer calls `T.register()` at run time (main, _decode_init, _decode_job) -> sees the wrapper
import m34_infer as X                                    # noqa: E402


def _probe(contig):
    MC = T.register()
    s = MC.read_sequence(contig)
    return dict(contig=contig, length=MC.LENGTH[contig], seq_len=len(s), md5=hashlib.md5(s.encode()).hexdigest()[:12], species=MC.SEQ2SPECIES[contig],
                acgt_fraction=round(sum(s.count(b) for b in 'ACGT') / max(len(s), 1), 4))


def _probe_main(contigs):
    T.register()
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
