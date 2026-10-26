import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
sys.path.append(str(Path(__file__).resolve().parents[2] / 'M34-NTV2-SCALEUP' / 'scripts'))      # m34_train_common etc. (appended: a backend shim already on the path stays first)
