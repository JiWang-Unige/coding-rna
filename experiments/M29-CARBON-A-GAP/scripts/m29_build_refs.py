"""Cache the frozen-ruler reference objects of the six TRAIN contigs (run one contig per process)."""
import sys, time
import m29_common as M
t = time.time()
for seq in sys.argv[1:]:
    o = M.build_reference(seq)
    print(seq, len(o['refs']), len(o['isoforms']), round(time.time() - t, 1), flush=True)
