#!/usr/bin/env python3
"""Finalize shared draws and inspect actual labels, not model accuracy."""
import argparse
import json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
import manifest as M
from src.m28.labels import build_targets,candidate_targets
from src.m28.core import M28Core

def rows(path):
    with path.open() as f: return [json.loads(line) for line in f if line.strip()]

def save(path,items):
    with path.open("x") as f:
        for item in items: f.write(json.dumps(item)+"\n")

def main(out):
    out.mkdir(parents=True,exist_ok=True)
    base=M.C.ROOT/"outputs/M28-DATA-MANIFEST-R2"
    meta=M.C.metadata(M.C.yaml.safe_load(M.C.CONFIG.read_text()))
    all_draws,all_probs,all_exposure,summaries=[],[],[],{}
    for species in ("arabidopsis_thaliana","oryza_sativa"):
        directory=base/species
        catalog={r["id"]:r for r in rows(directory/"transcripts.jsonl")}
        windows=rows(directory/"windows_train.jsonl");by_id={r["id"]:r for r in windows}
        draws=rows(directory/"training_draws.jsonl")
        for d in draws: all_draws.append(dict(d,species=species))
        strata=defaultdict(list)
        for w in windows: strata[w["seqid"],w["strand"]].append(w)
        counts=M.allocate({q:L for q,L in meta[species]["lengths"].items() if meta[species]["splits"][q]=="train"})
        for (q,strand),group in strata.items():
            p=M.probabilities(group)
            for w,x in zip(group,p):
                all_probs.append({"window_id":w["id"],"species":species,"stratum":[q,strand],
                                  "probability_per_paired_species_draw":float(counts[q]/1536*x),
                                  "importance_weight":float(1/(len(group)*x))})
        exposure=Counter(i for d in draws for i in by_id[d["window_id"]]["positive_ids"])
        for tx in catalog.values():
            if tx["split"]=="train" and tx["training_primary"]:
                all_exposure.append({"species":species,"transcript":tx["id"],
                                     "eligible_positive":tx["eligible_positive"],
                                     "sampled_positive_exposures":exposure[tx["id"]]})
        # Diagnostic subset is a uniform random subset of draw ordinals, never a prefix cap.
        chosen=np.random.default_rng(23).choice(len(draws),size=min(64,len(draws)),replace=False)
        lengths={q:L for q,L in meta[species]["lengths"].items() if meta[species]["splits"][q]=="train"}
        dna=M.read_allowed(meta[species]["path"]/"genome.fa",lengths)
        fine=Counter();region=Counter();ep_positive=np.zeros(4,dtype=np.int64);ep_known=ep_positive.copy()
        positives=0;background=0;seen_strata=set()
        for i in chosen:
            w=by_id[draws[int(i)]["window_id"]]
            sequence=dna[w["seqid"]][w["start"]:w["end"]]
            if w["strand"]=="-": sequence=M.rc(sequence)
            targets=build_targets(w,catalog,sequence)
            fine.update(dict(zip(*np.unique(targets["seg"],return_counts=True))))
            region.update(dict(zip(*np.unique(targets["region"],return_counts=True))))
            ep_positive+=(targets["endpoint"]*targets["endpoint_known"]).sum(0).astype(np.int64)
            ep_known+=targets["endpoint_known"].sum(0)
            positives+=len(targets["positive_chains"]);background+=targets["reference_empty_window"]
            seen_strata.add((w["seqid"],w["strand"]))
            assert all(x==1 for x in candidate_targets(targets["positive_chains"],targets))
        summaries[species]={"diagnostic_draws":len(chosen),"strata_seen":[list(x) for x in sorted(seen_strata)],
                            "fine_label_counts":{str(k):int(v) for k,v in sorted(fine.items())},
                            "region_label_counts":{str(k):int(v) for k,v in sorted(region.items())},
                            "known_endpoint_positive_counts":ep_positive.tolist(),"known_endpoint_base_counts":ep_known.tolist(),
                            "positive_chain_instances":positives,"callable_reference_background_draws":int(background)}
    rng=np.random.default_rng(0)
    paired=[dict(all_draws[int(i)],shared_ordinal=j) for j,i in enumerate(rng.permutation(len(all_draws)))]
    save(out/"paired_training_draws.jsonl",paired)
    save(out/"window_sampling_probabilities.jsonl",all_probs)
    save(out/"training_gene_exposure.jsonl",all_exposure)
    report={"experiment":"M28-CORE-LABEL-SMOKE-R2","no_backbone_or_training":True,
            "parameter_counts_mock_feature_dim64":{arm:sum(p.numel() for p in M28Core(arm,64).parameters()) for arm in ("C0","B1")},
            "shared_draws":len(paired),"species_draw_counts":dict(Counter(d["species"] for d in paired)),
            "paired_probability_sum":sum(d["probability_per_paired_species_draw"] for d in all_probs),
            "label_subset_seed":23,"shared_order_seed":0,"label_diagnostics":summaries,
            "limitations":["Not a trained model or native decoder test","No actual free-inference candidate coverage yet",
                           "Natural background may be rare, especially in Arabidopsis",
                           "CDS-complete training positives and historical evaluation primary policy differ"]}
    with (out/"summary.json").open("x") as f: json.dump(report,f,indent=2);f.write("\n")
    print(json.dumps(report,indent=2))

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output-dir",type=Path,required=True)
    main(p.parse_args().output_dir)
