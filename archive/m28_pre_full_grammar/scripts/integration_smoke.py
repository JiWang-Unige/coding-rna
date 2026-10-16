#!/usr/bin/env python3
"""Real frozen features + untrained M28 heads: reference-free generation and gradient wiring."""
import argparse
import json
import time
from collections import Counter
from pathlib import Path
import numpy as np
import torch
import manifest as M
from src.m28.core import M28Core,select_nonoverlapping
from src.m28.candidates import generate,Budget
from src.m28.c0_decode import decode_c0
from src.m28.labels import build_targets,candidate_targets,CDS_COLUMNS
from src.m28.training import joint_losses,one_hot

def rows(path):
    with path.open() as f: return [json.loads(line) for line in f if line.strip()]

def norm_groups(model,loss):
    named=list(model.named_parameters())
    gradients=torch.autograd.grad(loss,[p for _,p in named],retain_graph=True,allow_unused=True) if loss.requires_grad else [None]*len(named)
    sums=Counter()
    for (name,_),g in zip(named,gradients):
        group="link" if name.startswith("chain.link.") else name.split(".")[0]
        if g is not None: sums[group]+=float(g.detach().square().sum())
    return {k:float(sums[k]**.5) for k in ("shared","segmentation","endpoint","link","chain")}

def main(species,feature_dir,out):
    out.mkdir(parents=True,exist_ok=True)
    summary=json.loads((feature_dir/"summary.json").read_text())
    files=[r for r in summary["records"] if r["species"]==species]
    if len(files)!=6: raise ValueError("Expected exactly six preselected windows per species")
    catalog={r["id"]:r for r in rows(M.C.ROOT/"outputs/M28-DATA-MANIFEST-R2"/species/"transcripts.jsonl")}
    torch.manual_seed(0);torch.set_num_threads(2)
    model=M28Core("B1",summary["hidden_size"]).eval()
    records=[];gradient_record=None
    for record in files:
        cached=torch.load(feature_dir/record["artifact"],map_location="cpu",weights_only=True)
        w=cached["window"];sequence=cached["sequence"];n=w["valid_bases"]
        feature=cached["features"].float()[None]
        dna=one_hot(sequence,w["width"])[None]
        start=time.perf_counter()
        with torch.no_grad():
            outputs=model(feature,dna,[n])
            h=outputs["features"][0,:n]
            prob=outputs["segmentation_logits"][0,:n].softmax(-1)
            ep=outputs["endpoint_logits"][0,:n].numpy()
            free=generate(sequence,ep,prob[:,list(CDS_COLUMNS)].sum(-1).numpy(),
                          link_score_fn=lambda pairs:model.chain.link_logits(h,pairs).numpy(),budget=Budget())
            scores=model.chain(h,free["chains"])
        free_seconds=time.perf_counter()-start
        # Persist the free inference output BEFORE any reference-based labels/injection.
        free_file=out/("free_"+str(record["ordinal"]).zfill(2)+".json")
        with free_file.open("x") as f:
            json.dump(dict(free,window_id=w["id"],score=scores["logits"].tolist(),
                           additive_only_score=scores["additive_only_logits"].tolist()),f)
        targets=build_targets(w,catalog,sequence)
        labels=candidate_targets(free["chains"],targets)
        matches=len({tuple(c) for c in free["chains"]}&{tuple(c) for c in targets["positive_chains"]})
        c0_start=time.perf_counter()
        c0=decode_c0(sequence,prob.numpy())
        c0_seconds=time.perf_counter()-c0_start
        item={"window_id":w["id"],"ordinal":record["ordinal"],
              "free_chains":len(free["chains"]),"free_reference_labels":dict(Counter(map(str,labels))),
              "supported_training_positives":len(targets["positive_chains"]),"free_exact_training_positives":matches,
              "untrained_B1_selected":len(select_nonoverlapping(free["chains"],scores["logits"])),
              "untrained_additive_selected":len(select_nonoverlapping(free["chains"],scores["additive_only_logits"])),
              "untrained_C0_chains":len(c0["chains"]),"C0_partial_paths":c0["partial_paths"],
              "free_forward_and_candidate_seconds":free_seconds,"C0_decode_seconds":c0_seconds,
              "counts":free["counts"]}
        if gradient_record is None and targets["positive_chains"] and 0 in labels:
            start=time.perf_counter()
            live=model(feature,dna,[n])
            joint=joint_losses(model,live,free,w,catalog,targets,
                               window_weight=cached["draw"]["importance_weight_uniform_within_stratum"])
            norms={name:norm_groups(model,loss) for name,loss in joint["losses"].items()}
            for loss_name,head in (("local","segmentation"),("endpoint","endpoint"),("link","link"),("chain","chain")):
                if norms[loss_name][head]<=0 or norms[loss_name]["shared"]<=0:
                    raise RuntimeError("Missing actual gradient for "+loss_name)
            cg=torch.autograd.grad(joint["total"],joint["chain_scores"]["logits"],retain_graph=True)[0]
            cl=np.array(joint["chain_labels"])
            if not (cg[cl==0]>0).all() or not (cg[cl==1]<0).all() or not (cg[cl<0]==0).all():
                raise RuntimeError("Candidate/null gradient direction or unknown-mask failure")
            eg=torch.autograd.grad(joint["total"],live["endpoint_logits"],retain_graph=True)[0][0]
            masked=~torch.as_tensor(targets["endpoint_known"])
            if eg[masked].abs().sum()!=0: raise RuntimeError("Masked endpoint contributes gradient")
            gradient_record={"window_id":w["id"],"losses":{k:float(v.detach()) for k,v in joint["losses"].items()},
                             "gradient_norms":norms,"chain_labels":dict(Counter(map(str,joint["chain_labels"]))),
                             "link_labels":dict(Counter(map(str,joint["link_labels"]))),
                             "origins":dict(Counter(joint["origins"])),"frozen_features_require_grad":feature.requires_grad,
                             "masked_endpoint_gradient_abs_sum":float(eg[masked].abs().sum()),
                             "positive_negative_unknown_gradient_check":True,
                             "seconds":time.perf_counter()-start}
            with (out/"training_pool_debug.json").open("x") as f:
                json.dump({"chains":joint["training_chains"],"origins":joint["origins"],
                           "chain_labels":joint["chain_labels"],"links":joint["link_pairs"],"link_labels":joint["link_labels"]},f)
            del joint,live,cg,eg
        records.append(item);print(json.dumps(item),flush=True)
    report={"experiment":"M28-B1-INTEGRATION-R3","species":species,"head_seed":0,"head_fitted":False,
            "backbone_revision":summary["revision"],"head_hidden":128,"budget":free["budget"],
            "loss":"per-window class-balanced endpoint/link/chain BCE + local fine/grouped likelihood; coefficients all1 for smoke",
            "gradient_check":gradient_record,"integration_pass":gradient_record is not None,
            "optimizer_steps":0,"reference_used_for_free_generation":False,"val_test_setaria_used":False,
            "scope":"six isolated TRAIN windows, not chromosome-level assessment; no performance/architecture conclusion",
            "records":records}
    with (out/"summary.json").open("x") as f: json.dump(report,f,indent=2);f.write("\n")
    print(json.dumps(report,indent=2),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--species",choices=("arabidopsis_thaliana","oryza_sativa"),required=True)
    p.add_argument("--feature-dir",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    a=p.parse_args();main(a.species,a.feature_dir,a.output_dir)
