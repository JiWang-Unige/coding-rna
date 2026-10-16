import importlib.util
from pathlib import Path
import numpy as np

PATH=Path(__file__).resolve().parents[1]/"scripts/experiments/M27-ALLELIC-INTEGRITY/wt_tiberius.py"
spec=importlib.util.spec_from_file_location("native_tiberius_trace",PATH)
M=importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)


def test_native_trace_returns_original_scores_and_creates_no_cache(tmp_path):
    x=np.zeros((1,3,5),dtype=np.float32)
    y=np.ones((1,3,15),dtype=np.float32)
    calls=[]
    def forward():
        calls.append(1)
        return y
    trace=M.NeuralCache(None,"native")
    assert trace.call(x,{"stage":"repredict","strand":"-"},forward) is y
    trace.finish()
    assert calls==[1] and trace.bytes==0
    assert trace.calls[0]["score_shape"]==[1,3,15]
    assert trace.calls[0]["input_shape"]==[1,3,5]
    assert trace.calls[0]["stage"]=="repredict"
    assert not list(tmp_path.iterdir())
