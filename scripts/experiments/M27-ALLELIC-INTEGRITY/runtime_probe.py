#!/usr/bin/env python3
"""Inspect the existing Tiberius container; no model load or inference."""
import importlib.metadata
import inspect
import json
import sys
from pathlib import Path

result={"python":sys.version,"distributions":{},"model_loaded":False,"inference_run":False}
for name in ("tiberius","tensorflow","bricks2marble","hidten"):
    try:
        result["distributions"][name]=importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        result["distributions"][name]=None
try:
    from tiberius.eval_model_class import PredictionGTF
    result["PredictionGTF_file"]=inspect.getfile(PredictionGTF)
    result["interfaces"]={}
    for name in ("load_model","lstm_prediction","hmm_prediction","predict_vit",
                 "predict_function","repredict_function","adapt_batch_size"):
        fun=getattr(PredictionGTF,name,None)
        result["interfaces"][name]=None if fun is None else str(inspect.signature(fun))
    result["hmm_prediction_source"]=inspect.getsource(PredictionGTF.hmm_prediction)
    result["predict_function_source"]=inspect.getsource(PredictionGTF.predict_function)
    result["repredict_function_source"]=inspect.getsource(PredictionGTF.repredict_function)
    result["import_status"]="success"
except Exception as exc:
    result["import_status"]="failed"
    result["error_type"]=type(exc).__name__
    result["error"]=str(exc)
with Path(sys.argv[1]).open("x") as handle:
    json.dump(result,handle,indent=2)
print(json.dumps({k:v for k,v in result.items() if not k.endswith("_source")},indent=2))
