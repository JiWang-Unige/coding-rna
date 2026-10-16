#!/usr/bin/env python3
"""Trace native Tiberius 2.0.5; optionally record/replay WT neural calls."""
import argparse
import json
import sys
import time
from pathlib import Path
import numpy as np


def dump(path, obj):
    with path.open("x") as handle:
        json.dump(obj,handle,indent=2)


class NeuralCache:
    def __init__(self, root, mode, limit=20*1024**3):
        self.root = Path(root) if root is not None else None
        self.mode, self.limit = mode, limit
        self.index, self.bytes, self.calls = 0, 0, []
        if mode == "native":
            self.expected = None
        elif self.root is None:
            raise ValueError("Record/replay requires --cache-dir")
        elif mode == "record":
            self.root.mkdir(exist_ok=False)
            self.expected = None
        else:
            self.expected = json.loads((self.root/"requests.json").read_text())

    def call(self, x, metadata, forward):
        i=self.index
        self.index += 1
        input_path=self.root/f"{i:04d}.input.npy" if self.root is not None else None
        score_path=self.root/f"{i:04d}.scores.npy" if self.root is not None else None
        meta={**metadata,"input_shape":list(x.shape),"input_dtype":str(x.dtype)}
        started=time.monotonic()
        if self.mode == "native":
            y=forward()
            meta.update(score_shape=list(y.shape),score_dtype=str(y.dtype))
        elif self.mode == "record":
            y=forward()
            # Include ample header space for these ordinary numeric arrays.
            if self.bytes+x.nbytes+y.nbytes+4096 > self.limit:
                raise RuntimeError("Neural cache exceeds configured storage limit")
            with input_path.open("xb") as f:
                np.save(f,x,allow_pickle=False)
            with score_path.open("xb") as f:
                np.save(f,y,allow_pickle=False)
            self.bytes += input_path.stat().st_size+score_path.stat().st_size
            meta.update(score_shape=list(y.shape),score_dtype=str(y.dtype))
        else:
            if i >= len(self.expected):
                raise RuntimeError("Replay requested an uncached neural call")
            expected=self.expected[i]
            if any(meta[k] != expected[k] for k in meta):
                raise RuntimeError(f"Replay request metadata differs at call {i}")
            old=np.load(input_path,mmap_mode="r",allow_pickle=False)
            if not np.array_equal(x,old):
                raise RuntimeError(f"Replay neural input differs at call {i}")
            y=np.load(score_path,mmap_mode="r",allow_pickle=False)
            if list(y.shape) != expected["score_shape"] or str(y.dtype) != expected["score_dtype"]:
                raise RuntimeError("Cached score shape/dtype differs")
            meta.update(score_shape=list(y.shape),score_dtype=str(y.dtype))
        self.calls.append({**meta,"seconds":time.monotonic()-started})
        return y

    def finish(self):
        if self.mode == "record":
            dump(self.root/"requests.json",self.calls)
        elif self.mode == "replay" and self.index != len(self.expected):
            raise RuntimeError("Replay consumed fewer requests than native WT")


def snapshot(annotation):
    result=[]
    for sequence in annotation:
        for tx in sequence:
            blocks=sorted((c.start,c.end) for c in tx.cds)
            order=blocks if tx.strand=="+" else blocks[::-1]
            cumulative, phase = 0, {}
            for s,e in order:
                phase[(s,e)]=(3-cumulative%3)%3
                cumulative += e-s
            result.append([tx.sequence,tx.strand,[[s,e,phase[(s,e)]] for s,e in blocks]])
    return sorted(result)


def main():
    parser=argparse.ArgumentParser(add_help=False)
    parser.add_argument("--cache-mode",choices=["record","replay","native"],required=True)
    parser.add_argument("--cache-dir",type=Path)
    parser.add_argument("--cache-limit-bytes",type=int,default=20*1024**3)
    parser.add_argument("--trace",type=Path,required=True)
    ours,rest=parser.parse_known_args()
    sys.argv=[sys.argv[0]]+rest
    from tiberius.tiberius_args import parseCmd
    from tiberius.main import run_tiberius
    from tiberius.eval_model_class import PredictionGTF
    import bricks2marble as b2m
    import tensorflow as tf
    args=parseCmd()
    if len(tf.config.list_physical_devices("GPU")) != 1:
        raise RuntimeError("Exactly one allocated GPU must be visible")
    model_config=json.loads((Path(args.model)/"model_config.json").read_text())
    if model_config.get("inp_size") != 5 or model_config.get("clamsa") or model_config.get("hmm"):
        raise RuntimeError("Frozen no-softmask neural-only weight configuration differs")
    cache=NeuralCache(ours.cache_dir,ours.cache_mode,limit=ours.cache_limit_bytes)
    traces={"argv":rest,"model_config":model_config,"mode":ours.cache_mode,"requests":[],"filters":[],"neural_forward_calls":0}
    context={}
    native_neural=PredictionGTF.lstm_prediction

    def neural(self,x,clamsa_inp=None,batch_size=None):
        if clamsa_inp is not None or x.shape[-1] != 5:
            raise RuntimeError("Frozen sequence-only/no-softmask input violated")
        if not context or context["next"] >= len(context["directions"]):
            raise RuntimeError("Neural call has no native initial/reprediction context")
        direction=context["directions"][context["next"]]
        context["next"] += 1
        metadata={"stage":context["stage"],"strand":direction,
                  "native_sequences":context["sequences"],
                  "chunk_length":int(x.shape[1]),
                  "batch_size":int(batch_size or self.adapted_batch_size),
                  "parallel_factor":int(self.parallel_factor)}
        def forward():
            traces["neural_forward_calls"] += 1
            return native_neural(self,x,clamsa_inp=clamsa_inp,batch_size=batch_size)
        return cache.call(x,metadata,forward)

    def callback(original,stage):
        def run(self,fasta):
            directions=["+","-"]
            if stage=="repredict":
                evidence=fasta.evidence[:,0]
                directions=(["+"] if np.isin(evidence,[0,2]).any() else []) + (["-"] if np.isin(evidence,[1,2]).any() else [])
            context.clear()
            context.update(stage=stage,directions=directions,next=0,
                           sequences=[[s.name,int(s.start),int(s.end),int(s.N),int(s.T)] for s in fasta])
            value=original(self,fasta)
            if context["next"] != len(directions):
                raise RuntimeError("Native callback neural-call count differs")
            context.clear()
            return value
        return run

    PredictionGTF.lstm_prediction=neural
    PredictionGTF.predict_function=callback(PredictionGTF.predict_function,"initial")
    PredictionGTF.repredict_function=callback(PredictionGTF.repredict_function,"repredict")

    def filter_wrapper(original,name):
        def run(annotation,*args,**kwargs):
            before=snapshot(annotation)
            value=original(annotation,*args,**kwargs)
            traces["filters"].append({"name":name,"before":before,"after":snapshot(annotation)})
            return value
        return run
    for name in ("check_min_coding_length","check_inframe_stop_codons"):
        setattr(b2m.tools,name,filter_wrapper(getattr(b2m.tools,name),name))
    started=time.monotonic()
    try:
        run_tiberius(args)
        cache.finish()
        traces["status"]="COMPLETED"
    except BaseException as exc:
        traces["status"]="FAILED"
        traces["error"]=str(exc)
        raise
    finally:
        traces["seconds"]=time.monotonic()-started
        traces["requests"]=cache.calls
        traces["cache_bytes_written"]=cache.bytes
        try:
            traces["tensorflow_GPU_memory"]=tf.config.experimental.get_memory_info("GPU:0")
        except (ValueError,RuntimeError):
            traces["tensorflow_GPU_memory"]=None
        dump(ours.trace,traces)


if __name__=="__main__":
    main()
