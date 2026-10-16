#!/usr/bin/env python3
"""Four preregistered native confirmations; recording is not score intervention."""
import argparse
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path
import paired_run as P
import paired_response as R
from wt_compare import parsed, key
from prepare_reference import LENGTH

ROOT=P.ROOT
OLD=ROOT/"outputs/M27-PAIRED-NATIVE-R1"
SELECTED=(("ANNEVO",1,"ENSG00000198062.16"),("Tiberius",9,"ENSG00000196431.4"))


def counter_records(counts):
    return [{"key":k,"count":v} for k,v in sorted(counts.items())]


def chain_difference(old,new):
    a,b=Counter(map(key,old)),Counter(map(key,new))
    return {"equal":a==b,"old_count":sum(a.values()),"new_count":sum(b.values()),
            "lost":counter_records(a-b),"added":counter_records(b-a)}


def target_signature(row,txs,response):
    span=(row["CDS_chain_0based_halfopen"][0][0],row["CDS_chain_0based_halfopen"][-1][1])
    related=[t for t in txs if t["seqid"]=="chr22" and t["strand"]==row["strand"] and t["CDS"]
             and R.intersects((min(s for s,_,_ in t["CDS"]),max(e for _,e,_ in t["CDS"])),span)]
    # Compare identities/multiplicity of all related chains, not generated IDs.
    fields=("category","bypass_interval","assigned_keys","WT_first_CDS_anchor","WT_last_CDS_anchor",
            "edited_codon_positions_0based","coding_status","reason","foreign_coding_genes")
    data={k:response[k] for k in fields if k in response}
    data["related_chains"]=counter_records(Counter(map(key,related)))
    return json.loads(json.dumps(data))


def compare_target(row,old,new,genome,wt,foreign):
    a=target_signature(row,old,R.classify(row,old,genome,wt,foreign))
    b=target_signature(row,new,R.classify(row,new,genome,wt,foreign))
    return {"equal":a==b,"old":a,"new":b}


def compare_filters(row,old,new,genome,wt,foreign):
    names=["check_min_coding_length","check_inframe_stop_codons"]
    if [f["name"] for f in old["filters"]]!=names or [f["name"] for f in new["filters"]]!=names:
        raise RuntimeError("Frozen native filter sequence differs")
    result=[]
    for a,b in zip(old["filters"],new["filters"]):
        item={"name":a["name"]}
        for point in ("before","after"):
            at,bt=R.snapshot_txs(a[point]),R.snapshot_txs(b[point])
            item[point]={"target":compare_target(row,at,bt,genome,wt,foreign),
                         "whole_chr":chain_difference(at,bt)}
        result.append(item)
    return result


def request_metadata(trace):
    return [{k:v for k,v in r.items() if k!="seconds"} for r in trace["requests"]]


def run_case(method,case,row,out,wt,foreign):
    base=out/method/case["case_id"]
    base.mkdir(parents=True,exist_ok=False)
    receipt={"method":method,"case_id":case["case_id"],"allele":case["allele"],"gene_id":row["gene_id"],
             "status":"RUNNING","commands":[],"slurm_job_id":os.environ.get("SLURM_JOB_ID")}
    started=time.monotonic()
    try:
        fasta=ROOT/case["fasta"]
        with fasta.open() as handle:
            if next(handle).strip()!=">chr22":
                raise RuntimeError("Frozen input header differs")
            genome="".join(line.strip() for line in handle)
        if len(genome)!=LENGTH or genome[case["position_1based"]-1]!=case["alternate_base"]:
            raise RuntimeError("Frozen mutant input length/base differs")
        if method=="ANNEVO":
            P.storage(4*1024**3)
            repo=ROOT/"refs/repos/annevo-2026"
            h5=base/"native_scores.h5"
            output=base/"prediction.gff3"
            P.run_command([sys.executable,"-u","prediction.py","-g",fasta,"-m","saved_model/ANNEVO_Mammalia.pt",
                           "-l","Mammalia","-s","100","-p",h5,"--batch_size","2","--num_workers","2"],
                          base/"prediction.log",receipt,repo)
            P.run_command([sys.executable,"-u","decoding.py","-g",fasta,"-p",h5,"-o",output,
                           "-t","8","--min_intron_length","20","--min_prot_length","100","--show_log"],
                          base/"decode.log",receipt,repo)
            if "Process failed:" in (base/"decode.log").read_text(errors="replace"):
                raise RuntimeError("ANNEVO suppressed worker exception")
            import h5py
            with h5py.File(h5,"r") as handle:
                if set(handle.keys())!={"chr22"} or set(handle["chr22"].keys())!={"predictions_forward","predictions_reverse"}:
                    raise RuntimeError("ANNEVO cache chromosome/dataset scope differs")
                cache={n:{"shape":list(handle["chr22"][n].shape),"dtype":str(handle["chr22"][n].dtype)}
                       for n in handle["chr22"]}
                if any(v!={"shape":[LENGTH,15],"dtype":"float16"} for v in cache.values()):
                    raise RuntimeError("ANNEVO native cache shape/precision differs")
            receipt["retained_cache"]={"bytes":h5.stat().st_size,"datasets":cache}
        else:
            P.storage(1024**3)
            # Reserve 1 GiB for logs/traces/metadata outside array files.
            cache_limit=min(20*1024**3,P.BUDGET-P.storage()-1024**3)
            relative=base.relative_to(ROOT)
            output=base/"prediction.gtf"
            P.run_command(["apptainer","exec","--nv","--cleanenv",
                           "--env","CUDA_VISIBLE_DEVICES="+os.environ["CUDA_VISIBLE_DEVICES"],
                           "--env","PYTHONDONTWRITEBYTECODE=1","--env","OMP_NUM_THREADS=8",
                           "--bind",str(ROOT)+":/work","--pwd","/work",
                           ROOT/"refs/repos/tiberius-2024/singularity/tiberius_2.0.5.sif",
                           "python","-u","scripts/experiments/M27-ALLELIC-INTEGRITY/wt_tiberius.py",
                           "--cache-mode","record","--cache-dir",relative/"cache",
                           "--cache-limit-bytes",str(cache_limit),"--trace",relative/"trace.json",
                           "--genome",case["fasta"],"--model","refs/weights/tiberius-m27-nosm-v2/tiberius_nosm_weights_v2",
                           "--no_softmasking","--seq_len","400050","--batch_size","2","--out",relative/"prediction.gtf"],
                          base/"native.log",receipt)
            trace=json.loads((base/"trace.json").read_text())
            if not (trace["status"]=="COMPLETED" and trace["mode"]=="record"
                    and trace["neural_forward_calls"]==len(trace["requests"])>0
                    and 0<trace["cache_bytes_written"]<=cache_limit):
                raise RuntimeError("Native record trace incomplete or over budget")
            receipt["retained_cache"]={"bytes":trace["cache_bytes_written"],"limit_bytes":cache_limit,
                                       "requests":len(trace["requests"]),"neural_forward_calls":trace["neural_forward_calls"]}
        old_base=OLD/method/case["case_id"]
        old_txs=parsed(old_base/output.name)
        new_txs=parsed(output)
        comparison={"method":method,"case_id":case["case_id"],"gene_id":row["gene_id"],
                    "allele":case["allele"],"input_fasta":case["fasta"],
                    "target":compare_target(row,old_txs,new_txs,genome,wt,foreign),
                    "whole_chr":chain_difference(old_txs,new_txs)}
        wanted="WT_CHAIN_RETAINED" if case["allele"]=="syn" else "BYPASS"
        comparison["target_pass"]=(comparison["target"]["equal"] and comparison["target"]["new"]["category"]==wanted)
        comparison["whole_context_equal"]=comparison["whole_chr"]["equal"]
        if method=="Tiberius":
            old_trace=json.loads((old_base/"trace.json").read_text())
            comparison["filters"]=compare_filters(row,old_trace,trace,genome,wt,foreign)
            comparison["request_metadata_equal"]=request_metadata(old_trace)==request_metadata(trace)
            comparison["target_pass"] &= all(f[p]["target"]["equal"] for f in comparison["filters"] for p in ("before","after"))
            comparison["whole_context_equal"] &= (comparison["request_metadata_equal"] and
                all(f[p]["whole_chr"]["equal"] for f in comparison["filters"] for p in ("before","after")))
        comparison["component_attribution_performed"]=False
        comparison["mutant_cache_replay_verified"]=False
        P.write_json(base/"comparison.json",comparison)
        receipt["M27_bytes_after_case"]=P.storage()
        if not comparison["target_pass"]:
            raise RuntimeError("Target structure/competition/filter identity not reproduced; stop without retry")
        receipt["status"]="COMPLETED"
        print(json.dumps({k:comparison[k] for k in ("method","case_id","target_pass","whole_context_equal")}),flush=True)
    except BaseException as error:
        receipt.update(status="FAILED",error=repr(error))
        raise
    finally:
        receipt["seconds"]=time.monotonic()-started
        P.write_json(base/"execution.json",receipt)


def summarize(out):
    cases=[]
    for method,ordinal,gene in SELECTED:
        for allele in ("syn","PTC"):
            base=out/method/f"{ordinal:02d}_{allele}"
            item={"method":method,"gene_id":gene,"allele":allele,"status":"NOT_RUN"}
            if (base/"execution.json").is_file():
                item["execution"]=json.loads((base/"execution.json").read_text())
                item["status"]=item["execution"]["status"]
            if (base/"comparison.json").is_file():
                item["comparison"]=json.loads((base/"comparison.json").read_text())
            cases.append(item)
    complete=all(c["status"]=="COMPLETED" for c in cases)
    result={"stage":"native_reconfirmation_with_score_recording","planned_cases":4,"cases":cases,
            "complete":complete,"all_targets_reproduced":complete and all(c["comparison"]["target_pass"] for c in cases),
            "all_whole_context_equal":complete and all(c["comparison"]["whole_context_equal"] for c in cases),
            "mutant_cache_replay_verified":False,"component_attribution_performed":False,
            "training_runs":0,"Setaria_accessed":False}
    P.write_json(out/"confirmation_result.json",result)
    return result


def main(out):
    out=out.resolve()
    if out.name!="M27-NATIVE-CONFIRM-R1" or (out/"confirmation_result.json").exists():
        raise RuntimeError("New dedicated confirmation directory required")
    P.STORAGE.append(out)
    manifest=json.loads((OLD/"input_manifest.json").read_text())
    results=json.loads((OLD/"paired_result.json").read_text())
    rows={r["ordinal"]:r for r in manifest["registered_loci"]}
    cases={c["case_id"]:c for c in manifest["cases"]}
    genes=R.reference_coding(ROOT/"data/m27_grch38_chr22/gencode.v49.chr22.gtf")
    try:
        for method,ordinal,gene in SELECTED:
            row=rows[ordinal]
            if row["gene_id"]!=gene or results["methods"][method]["first_preregistered_mechanism_candidate"]!=gene:
                raise RuntimeError("Preregistered confirmation target differs")
            ext="gff3" if method=="ANNEVO" else "gtf"
            wt=set(map(key,parsed(ROOT/f"outputs/M27-WT-REPLAY-R1/{method}/native.{ext}")))
            foreign=R.foreign_regions(row,genes)
            for allele in ("syn","PTC"):
                run_case(method,cases[f"{ordinal:02d}_{allele}"],row,out,wt,foreign)
    finally:
        summarize(out)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output-dir",type=Path,required=True)
    main(parser.parse_args().output_dir)
