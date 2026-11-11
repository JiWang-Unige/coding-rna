import os
import sys
from pathlib import Path
os.environ.setdefault('M34_BACKEND', 'pc2')
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
sys.path.append(str(Path(__file__).resolve().parents[2] / 'M35-PLANTCAD2-LORA' / 'scripts'))   # pc2_lora for the toy-Caduceus tests (last: M32 adapter wins)
