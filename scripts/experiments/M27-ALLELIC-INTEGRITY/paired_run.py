#!/usr/bin/env python3
"""Execute the finite M27 native paired panel serially; never retry failed cases."""
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path
import paired_response as R
from wt_compare import parsed, key
from prepare_reference import LENGTH

ROOT=Path(__file__).resolve().parents[3]
BUDGET=50*1024**3
STORAGE=[ROOT/p for p in ("data/m27_grch38_chr22","refs/weights/tiberius-m27-nosm-v2",
    "outputs/M27-ALLELIC-PREP-R1","outputs/M27-ALLELIC-PREP-R2","outputs/M27-WT-REPLAY-R1",
    "outputs/M27-PAIRED-NATIVE-R1")]


def write_json(path,data):
    with path.open("x") as handle:
        json.dump(data,handle,indent=2)


def storage(reserve=0):
    output=subprocess.check_output(["du","-sb",*[str(p) for p in STORAGE]],text=True)
    total=sum(int(line.split()[0]) for line in output.splitlines())
    if total+reserve>BUDGET:
        raise RuntimeError("M27 aggregate 50 GiB storage limit would be exceeded")
    return total


def run_command(command,path,receipt,cwd=ROOT):
    actual=["/usr/bin/time","-v",*map(str,command)]
    started=time.monotonic()
    entry={"argv":actual,"cwd":str(cwd)}
    receipt["commands"].append(entry)
    with path.open("x") as log:
        completed=subprocess.run(actual,cwd=cwd,stdout=log,stderr=subprocess.STDOUT)
    entry.update(returncode=completed.returncode,seconds=time.monotonic()-started,log=str(path.relative_to(ROOT)))
    completed.check_returncode()


def run_case(method,case,row,out,wt_keys,foreign):
    storage(4*1024**3 if method=="ANNEVO" else 1024**3)
    base=out/method/case["case_id"]
    base.mkdir(parents=True)
    receipt={"case_id":case["case_id"],"method":method,"gene_id":row["gene_id"],
             "allele":case["allele"],"status":"RUNNING","commands":[],"slurm_job_id":os.environ.get("SLURM_JOB_ID")}
    started=time.monotonic()
    try:
        fasta=ROOT/case["fasta"]
        with fasta.open() as handle:
            if next(handle).strip()!=">chr22":
                raise RuntimeError("Single-chr22 input header differs")
            genome="".join(line.strip() for line in handle)
        if len(genome)!=LENGTH or genome[case["position_1based"]-1]!=case["alternate_base"]:
            raise RuntimeError("Mutant input length or registered alternate base differs")
        if method=="ANNEVO":
            repo=ROOT/"refs/repos/annevo-2026"
            h5=base/"temporary_native_scores.h5"
            output=base/"prediction.gff3"
            run_command([sys.executable,"-u","prediction.py","-g",fasta,"-m","saved_model/ANNEVO_Mammalia.pt",
                         "-l","Mammalia","-s","100","-p",h5,"--batch_size","2","--num_workers","2"],
                        base/"prediction.log",receipt,repo)
            run_command([sys.executable,"-u","decoding.py","-g",fasta,"-p",h5,"-o",output,
                         "-t","8","--min_intron_length","20","--min_prot_length","100","--show_log"],
                        base/"decode.log",receipt,repo)
            if "Process failed:" in (base/"decode.log").read_text(errors="replace"):
                raise RuntimeError("ANNEVO suppressed worker exception; output is incomplete")
        else:
            output=base/"prediction.gtf"
            relative=base.relative_to(ROOT)
            run_command(["apptainer","exec","--nv","--cleanenv",
                         "--env","CUDA_VISIBLE_DEVICES="+os.environ["CUDA_VISIBLE_DEVICES"],
                         "--env","PYTHONDONTWRITEBYTECODE=1","--env","OMP_NUM_THREADS=8",
                         "--bind",str(ROOT)+":/work","--pwd","/work",
                         ROOT/"refs/repos/tiberius-2024/singularity/tiberius_2.0.5.sif",
                         "python","-u","scripts/experiments/M27-ALLELIC-INTEGRITY/wt_tiberius.py",
                         "--cache-mode","native","--trace",relative/"trace.json","--genome",case["fasta"],
                         "--model","refs/weights/tiberius-m27-nosm-v2/tiberius_nosm_weights_v2",
                         "--no_softmasking","--seq_len","400050","--batch_size","2","--out",relative/"prediction.gtf"],
                        base/"native.log",receipt)
        predictions=parsed(output)
        response=R.classify(row,predictions,genome,wt_keys,foreign)
        response.update(method=method,case_id=case["case_id"],gene_id=row["gene_id"],allele=case["allele"],
                        whole_chromosome_CDS_transcripts=len(predictions))
        if method=="Tiberius":
            trace=json.loads((base/"trace.json").read_text())
            if not (trace["status"]=="COMPLETED" and trace["mode"]=="native" and trace["cache_bytes_written"]==0
                    and trace["neural_forward_calls"]==len(trace["requests"])>0
                    and [f["name"] for f in trace["filters"]]==["check_min_coding_length","check_inframe_stop_codons"]):
                raise RuntimeError("Native-only request/filter trace is incomplete")
            response["neural_requests"]=len(trace["requests"])
            response["request_stage_counts"]=dict(R.Counter(r["stage"] for r in trace["requests"]))
            response["filter_process_labels"]=R.filter_labels(row,trace,genome,wt_keys,foreign)
        receipt["M27_bytes_before_temporary_cleanup"]=storage()
        if method=="ANNEVO":
            receipt["temporary_cache_bytes_deleted"]=h5.stat().st_size
            h5.unlink()
        write_json(base/"response.json",response)
        receipt["status"]="COMPLETED"
        print(json.dumps({"case_id":case["case_id"],"method":method,"status":"COMPLETED","category":response["category"]}),flush=True)
    except BaseException as error:
        receipt["status"]="FAILED"
        receipt["error"]=repr(error)
        raise
    finally:
        receipt["seconds"]=time.monotonic()-started
        write_json(base/"execution.json",receipt)


def main(out):
    manifest=json.loads((out/"input_manifest.json").read_text())
    if manifest["planned_method_cases"]!={"ANNEVO":24,"Tiberius":20}:
        raise RuntimeError("Frozen planned case counts differ")
    cases={c["case_id"]:c for c in manifest["cases"]}
    genes=R.reference_coding(ROOT/"data/m27_grch38_chr22/gencode.v49.chr22.gtf")
    wt={}
    for m,ext in (("ANNEVO","gff3"),("Tiberius","gtf")):
        wt[m]=set(map(key,parsed(ROOT/f"outputs/M27-WT-REPLAY-R1/{m}/native.{ext}")))
    try:
        for row in manifest["registered_loci"]:
            foreign=R.foreign_regions(row,genes)
            for method in ("ANNEVO","Tiberius"):
                if row["method_status"][method]!="paired_native_planned":
                    continue
                for allele in ("syn","PTC"):
                    case=cases[f'{row["ordinal"]:02d}_{allele}']
                    run_case(method,case,row,out,wt[method],foreign)
                    R.summarize(out)
    finally:
        R.summarize(out)
    storage()


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    main(args.output_dir.resolve())
