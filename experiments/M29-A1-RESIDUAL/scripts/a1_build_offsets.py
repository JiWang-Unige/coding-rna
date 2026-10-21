"""Build compact fixed-offset tracks for TRAIN (6) and DEV (2) contigs -> outputs/M29-A1-RESIDUAL/offset/<seq>_<strand>.npy
TRAIN probabilities: saved Carbon-A raw probs (M28-CARBON-A-TRAIN-RAW); transition tracks: M29 `_cache/tt_*` (same as the parity gate).
DEV probabilities: ont SDPA Carbon-A raw probs (M29 G0); transition tracks computed here with the Carbon-A transition network.
Usage: a1_build_offsets.py seq[,seq...]
"""
import glob, sys, time
import numpy as np
import a1_common as K
sys.path.insert(0, str(K.CA_SRC))
import torch
torch.set_num_threads(4)
import m29_gate_parity as GP                    # noqa: E402  (Carbon-A copy imports; frozen gate construction)
from constants import MIN_PROB, MAX_PROB, GC_AG_INTRON_LOG_PRIOR, TRANSITION_NETWORK_PATH   # noqa: E402
from transition_network import load_transition_network   # noqa: E402
import a1_offset as O                            # noqa: E402

DEV_RAW = {'NC_003074.8': 'G0/ara_raw/NC_003074.8_whole_input_*.raw_probs.parquet',
           'NC_089041.1': 'G0/rice_raw3/NC_089041.1_whole_input_*.raw_probs.parquet'}


def load_dev(seq):
    import pyarrow.parquet as pq
    from constants import RECORD_ID_COLUMN, SEQUENCE_COLUMN, POSITIVE_STRAND_PROBABILITY_COLUMN, NEGATIVE_STRAND_PROBABILITY_COLUMN
    from utils import probability_vector_from_parquet
    p = sorted(glob.glob(str(K.M29 / DEV_RAW[seq])))
    assert len(p) == 1, p
    b = list(pq.ParquetFile(p[0]).iter_batches(batch_size=1, columns=[RECORD_ID_COLUMN, SEQUENCE_COLUMN,
             POSITIVE_STRAND_PROBABILITY_COLUMN, NEGATIVE_STRAND_PROBABILITY_COLUMN]))
    assert len(b) == 1
    b = b[0]
    s = b.column(1)[0].as_py().upper()
    return b.column(0)[0].as_py(), s, probability_vector_from_parquet(b.column(2)[0])[:len(s)], probability_vector_from_parquet(b.column(3)[0])[:len(s)]


def main():
    (K.A1 / 'offset').mkdir(parents=True, exist_ok=True)
    tn = None
    for seq in sys.argv[1].split(','):
        t0 = time.time()
        rid, sequence, pos, neg = load_dev(seq) if seq in DEV_RAW else GP.load_tracks(seq)
        assert rid == seq
        tracks = {'+': (sequence, pos), '-': (O.rc(sequence), np.ascontiguousarray(neg[::-1]))}
        for s in '+-':
            out = O.offset_path(seq, s)
            if out.exists():
                continue
            f = K.CACHE / f'tt_{seq}_{s}.npy'
            if f.exists():
                tt = np.load(f)
            else:
                tn = tn or load_transition_network(TRANSITION_NETWORK_PATH, torch.device('cpu'))
                tt = tn.inference(tracks[s][1]).astype(np.float32)
                np.save(f, tt)
            c = O.compact(tracks[s][0], tracks[s][1], tt, MIN_PROB, MAX_PROB, GC_AG_INTRON_LOG_PRIOR)
            np.save(str(out) + '.tmp.npy', c)
            Path = __import__('pathlib').Path
            Path(str(out) + '.tmp.npy').rename(out)
            print(seq, s, c.shape, f'{time.time()-t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
