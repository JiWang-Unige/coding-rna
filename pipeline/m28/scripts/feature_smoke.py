#!/usr/bin/env python3
"""Twelve preselected training windows; frozen official GENERanno, no fitting."""
import argparse
import json
import re
import time
from pathlib import Path
import torch
from transformers import AutoModelForTokenClassification,AutoTokenizer
import manifest as M

MODEL="GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview"
REVISION="b0483c23b6b63787b61a6d3a204a9b517d6ba345"
BLOCK=6144

def rows(path):
    with path.open() as f: return [json.loads(line) for line in f if line.strip()]

def select_windows():
    base=M.C.ROOT/"outputs/M28-DATA-MANIFEST-R2"
    paired=M.C.ROOT/"outputs/M28-CORE-LABEL-SMOKE-R2/result/paired_training_draws.jsonl"
    by_id={}
    for species in ("arabidopsis_thaliana","oryza_sativa"):
        by_id.update({r["id"]:r for r in rows(base/species/"windows_train.jsonl")})
    selected=[];seen=set()
    for d in rows(paired):
        key=(d["species"],*d["stratum"])
        if key in seen: continue
        seen.add(key)
        selected.append({"window":by_id[d["window_id"]],"draw":d})
    if len(selected)!=12: raise ValueError("Expected exactly two species x three train chromosomes x two strands")
    return selected

def main(out):
    out.mkdir(parents=True,exist_ok=True)
    selection=select_windows()
    with (out/"selection.json").open("x") as f: json.dump(selection,f,indent=2)
    meta=M.C.metadata(M.C.yaml.safe_load(M.C.CONFIG.read_text()))
    sequences={}
    for species,record in meta.items():
        lengths={q:L for q,L in record["lengths"].items() if record["splits"][q]=="train"}
        sequences[species]=M.read_allowed(record["path"]/"genome.fa",lengths)
    tokenizer=AutoTokenizer.from_pretrained(MODEL,revision=REVISION,trust_remote_code=True,local_files_only=True)
    t0=time.perf_counter()
    full=AutoModelForTokenClassification.from_pretrained(
        MODEL,revision=REVISION,trust_remote_code=True,local_files_only=True,
        attn_implementation="sdpa",torch_dtype=torch.bfloat16)
    if int(getattr(full,"k",getattr(full.config,"k",0)))!=6: raise ValueError("Six-base token contract violated")
    backbone=full.model
    del full
    backbone.eval().requires_grad_(False).to("cuda")
    torch.cuda.synchronize()
    load_seconds=time.perf_counter()-t0
    torch.cuda.reset_peak_memory_stats()
    records=[];total_bytes=0
    for ordinal,item in enumerate(selection):
        w=item["window"]
        raw=sequences[w["species"]][w["seqid"]][w["start"]:w["end"]]
        oriented=raw if w["strand"]=="+" else M.rc(raw)
        sequence=re.sub("[^ACGT]","N",oriented)+"N"*w["pad_right_after_orientation"]
        if len(sequence)!=24576: raise ValueError("Window width changed")
        whole_ids=tokenizer(sequence,add_special_tokens=False,return_tensors="pt",truncation=False)["input_ids"]
        blocks=[tokenizer(sequence[a:a+BLOCK],add_special_tokens=False,return_tensors="pt",truncation=False)["input_ids"]
                for a in range(0,len(sequence),BLOCK)]
        if any(tuple(x.shape)!=(1,1024) for x in blocks) or not torch.equal(torch.cat(blocks,dim=1),whole_ids):
            raise ValueError("Whole-window vs four-block tokenization is not identical")
        features=[];times=[]
        for ids in blocks:
            torch.cuda.synchronize();start=time.perf_counter()
            with torch.inference_mode():
                hidden=backbone(input_ids=ids.to("cuda"),
                                attention_mask=torch.ones_like(ids,device="cuda")).last_hidden_state
            torch.cuda.synchronize()
            times.append(time.perf_counter()-start)
            if hidden.shape[1]!=1024 or not torch.isfinite(hidden).all(): raise ValueError("Invalid hidden features")
            features.append(hidden[0].detach().to("cpu",dtype=torch.bfloat16))
        feature=torch.cat(features,dim=0).contiguous()
        path=out/("window_"+str(ordinal).zfill(2)+".pt")
        if path.exists(): raise FileExistsError(path)
        torch.save({"features":feature,"window":w,"draw":item["draw"],
                    "sequence":oriented,"model":MODEL,"revision":REVISION},path)
        size=path.stat().st_size;total_bytes+=size
        if total_bytes>1024**3: raise RuntimeError("Feature artifact cap exceeded")
        records.append({"ordinal":ordinal,"window_id":w["id"],"species":w["species"],"seqid":w["seqid"],
                        "strand":w["strand"],"shape":list(feature.shape),"dtype":str(feature.dtype),
                        "block_forward_seconds":times,"artifact":path.name,"bytes":size})
        print(json.dumps(records[-1]),flush=True)
    report={"experiment":"M28-GLM-FEATURE-SMOKE-R3","model":MODEL,"revision":REVISION,
            "weights":"original released CDS annotator backbone, no old LoRA","window_bp":24576,"local_block_bp":BLOCK,
            "windows":len(records),"block_forwards":len(records)*4,"hidden_size":int(backbone.config.hidden_size),
            "load_seconds":load_seconds,"forward_seconds":sum(sum(r["block_forward_seconds"]) for r in records),
            "peak_gpu_allocated_bytes":torch.cuda.max_memory_allocated(),"peak_gpu_reserved_bytes":torch.cuda.max_memory_reserved(),
            "gpu":torch.cuda.get_device_name(),"cache_bytes":total_bytes,
            "trainable_backbone_parameters":sum(p.numel() for p in backbone.parameters() if p.requires_grad),
            "backbone_parameter_gradients":sum(p.grad is not None for p in backbone.parameters()),
            "fitting":False,"val_test_setaria_used":False,"records":records}
    with (out/"summary.json").open("x") as f: json.dump(report,f,indent=2);f.write("\n")
    print(json.dumps(report,indent=2),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output-dir",type=Path,required=True)
    main(p.parse_args().output_dir)
