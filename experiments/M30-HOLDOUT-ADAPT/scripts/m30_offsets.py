"""Compact fixed Carbon-A offset tracks for M30 contigs (same construction as the parity gate / A1).  Usage: m30_offsets.py contig[,contig...]"""
import sys, time
from pathlib import Path
import numpy as np
import m30_common as M
import a1_common as K
sys.path.insert(0, str(K.CA_SRC))
import torch
torch.set_num_threads(4)
import m29_gate_parity as GP                    # noqa: E402
from constants import MIN_PROB, MAX_PROB, GC_AG_INTRON_LOG_PRIOR, TRANSITION_NETWORK_PATH, RECORD_ID_COLUMN, SEQUENCE_COLUMN, POSITIVE_STRAND_PROBABILITY_COLUMN, NEGATIVE_STRAND_PROBABILITY_COLUMN   # noqa: E402
from transition_network import load_transition_network   # noqa: E402
from utils import probability_vector_from_parquet       # noqa: E402
import a1_offset as O                            # noqa: E402


def load(contig):
    import pyarrow.parquet as pq
    b = list(pq.ParquetFile(M.raw_parquet(contig)).iter_batches(batch_size=1, columns=[RECORD_ID_COLUMN, SEQUENCE_COLUMN, POSITIVE_STRAND_PROBABILITY_COLUMN, NEGATIVE_STRAND_PROBABILITY_COLUMN]))
    assert len(b) == 1
    b = b[0]
    s = b.column(1)[0].as_py().upper()
    return b.column(0)[0].as_py(), s, probability_vector_from_parquet(b.column(2)[0])[:len(s)], probability_vector_from_parquet(b.column(3)[0])[:len(s)]


def main():
    (M.OUT / 'offset').mkdir(parents=True, exist_ok=True)
    (M.OUT / 'tt').mkdir(parents=True, exist_ok=True)
    tn = None
    for contig in sys.argv[1].split(','):
        t0 = time.time()
        rid, seq, pos, neg = load(contig)
        assert rid == contig and len(seq) == M.lengths(M.species_of(contig))[contig]
        tr = {'+': (seq, pos), '-': (O.rc(seq), np.ascontiguousarray(neg[::-1]))}
        for s in '+-':
            out = M.offset_path(contig, s)
            if out.exists():
                continue
            f = M.OUT / 'tt' / f'{contig}_{s}.npy'
            if f.exists():
                tt = np.load(f)
            else:
                tn = tn or load_transition_network(TRANSITION_NETWORK_PATH, torch.device('cpu'))
                tt = tn.inference(tr[s][1]).astype(np.float32)
                np.save(f, tt)
            c = O.compact(tr[s][0], tr[s][1], tt, MIN_PROB, MAX_PROB, GC_AG_INTRON_LOG_PRIOR)
            np.save(str(out) + '.tmp.npy', c)
            Path(str(out) + '.tmp.npy').rename(out)
            print(contig, s, c.shape, f'{time.time()-t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
