"""M29-A2 shared setup: imports the A1 modules (unchanged) and defines the A2 output root."""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
A1_SCRIPTS = HERE.parents[1] / 'M29-A1-RESIDUAL' / 'scripts'
sys.path.insert(0, str(A1_SCRIPTS))
import a1_common as K            # noqa: E402
A2 = K.OUTS / 'M29-A2-TEMPERATURE'
A1 = K.A1
HELDOUT = ('NC_003071.7', 'NC_089044.1')          # internal validation chromosomes (one per species): never trained on, never DEV
SUBSTITUTE = {'NC_003071.7': 'NC_003070.9', 'NC_089044.1': 'NC_089043.1'}   # stratum remap so the frozen schedule keeps its structure
CPUS = [2, 4, 6, 8, 10, 12, 14] + list(range(16, 25))
