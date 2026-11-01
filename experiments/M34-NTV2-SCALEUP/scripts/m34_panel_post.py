"""P5 panel posterior: M33 `m33_posterior.py` (imported, not edited; gene-set posterior log P, tau = 1 as M33/M36) for a panel main sequence.  m33_posterior reads
the sequence through the M29 helpers registered for the M30 contigs only; the panel sequences are registered on top (m34_panel_common.register) before its main() runs.
Usage (identical to m33_posterior.py): m34_panel_post.py <run_dir> <seq> <strand> <taus> <out_dir>
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'M33-NTV2-FOLLOWUP' / 'scripts'))
import m33_posterior as P                                # noqa: E402  (module level: M30.register(), K.repoint(); P.M29 is the m29_common module)
sys.path.insert(0, str(HERE))
import m34_panel_common as PC                            # noqa: E402

PC.register(P.M29)

if __name__ == '__main__':
    P.main()
