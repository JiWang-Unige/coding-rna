#!/usr/bin/env python3
"""Compare native/replayed WT chains and reference-selected M27 target identities."""
import argparse
import json
import sys
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
import eval_structure_diagnostic as E

LENGTH=50818468


def key(tx):
    return (tx["seqid"],tx["strand"],tuple(sorted((s,e,int(p)) for s,e,p in tx["CDS"])))


def reference_key(row):
    blocks=[tuple(x) for x in row["CDS_chain_0based_halfopen"]]
    phase,cumulative={},0
    for s,e in (blocks if row["strand"]=="+" else blocks[::-1]):
        phase[(s,e)]=(3-cumulative%3)%3
        cumulative+=e-s
    return ("chr22",row["strand"],tuple((s,e,phase[(s,e)]) for s,e in blocks))


def match_target(row, predictions):
    wanted=reference_key(row)
    ids=[tx["id"] for tx in predictions if key(tx)==wanted]
    return {"exact_prediction_ids":ids,"exact_count":len(ids),
            "unique_WT_exact":len(ids)==1,"ambiguous_duplicate_exact":len(ids)>1}


def parsed(path):
    result=E.parse_annotation(path,{"chr22":LENGTH})
    return [tx for tx in result["transcripts"].values() if tx["CDS"]]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    out=args.output_dir
    targets=json.loads((ROOT/"outputs/M27-ALLELIC-PREP-R1/preparation.json").read_text())["selected"]
    summary={"stage":"WT_native_cache_replay_only","target_count":len(targets),
             "training_runs":0,"mutant_inference_runs":0,"methods":{},"targets":[]}
    failures=[]
    native_txs={}
    for method,extension in (("ANNEVO","gff3"),("Tiberius","gtf")):
        base=out/method
        native=parsed(base/f"native.{extension}")
        replay=parsed(base/f"replay.{extension}")
        a,b=Counter(map(key,native)),Counter(map(key,replay))
        same=a==b
        summary["methods"][method]={"native_transcripts":len(native),"replay_transcripts":len(replay),
            "native_unique_chains_phase":len(a),"replay_unique_chains_phase":len(b),
            "whole_chr_chain_strand_phase_equal":same,
            "native_only_records":sum((a-b).values()),"replay_only_records":sum((b-a).values())}
        if not same:
            failures.append(method+"_whole_chromosome_output_mismatch")
        native_txs[method]=native
    native_trace=json.loads((out/"Tiberius/native_trace.json").read_text())
    replay_trace=json.loads((out/"Tiberius/replay_trace.json").read_text())
    trace_ok=(native_trace["status"]==replay_trace["status"]=="COMPLETED"
              and native_trace["filters"]==replay_trace["filters"]
              and {f["name"] for f in native_trace["filters"]}=={"check_min_coding_length","check_inframe_stop_codons"}
              and native_trace["neural_forward_calls"]==len(native_trace["requests"])>0
              and replay_trace["neural_forward_calls"]==0
              and len(native_trace["requests"])==len(replay_trace["requests"]))
    summary["Tiberius_trace"]={"trace_equal":trace_ok,
        "requests":len(native_trace["requests"]),"stage_counts":dict(Counter(r["stage"] for r in native_trace["requests"])),
        "native_neural_forward_calls":native_trace["neural_forward_calls"],
        "replay_neural_forward_calls":replay_trace["neural_forward_calls"],
        "native_seconds":native_trace["seconds"],"replay_seconds":replay_trace["seconds"],
        "cache_bytes":native_trace["cache_bytes_written"],
        "filters":[{"name":f["name"],"before":len(f["before"]),"after":len(f["after"])} for f in native_trace["filters"]]}
    if not trace_ok:
        failures.append("Tiberius_request_or_filter_replay_mismatch")
    import h5py
    with h5py.File(out/"ANNEVO/WT.h5","r") as h:
        shapes={name:list(h["chr22"][name].shape) for name in h["chr22"]}
        if set(h.keys())!={"chr22"} or shapes!={"predictions_forward":[LENGTH,15],"predictions_reverse":[LENGTH,15]}:
            failures.append("ANNEVO_neural_cache_scope_mismatch")
        summary["ANNEVO_cache"]={"shapes":shapes,
            "dtypes":{name:str(h["chr22"][name].dtype) for name in h["chr22"]},
            "bytes":(out/"ANNEVO/WT.h5").stat().st_size}
    for row in targets:
        summary["targets"].append({"gene_id":row["gene_id"],"gene_name":row["gene_name"],
            "reference_transcript":row["transcript_id"],"strand":row["strand"],
            "position_1based":row["position_1based"],
            **{m:match_target(row,txs) for m,txs in native_txs.items()}})
    for m in native_txs:
        records=summary["targets"]
        n=sum(r[m]["unique_WT_exact"] for r in records)
        summary["methods"][m].update(WT_exact_unique=n,WT_not_exact=sum(r[m]["exact_count"]==0 for r in records),
            WT_ambiguous_duplicate_exact=sum(r[m]["ambiguous_duplicate_exact"] for r in records),
            conditional_analysis_has_WT_targets=n>0)
    summary["common_WT_exact_unique"]=sum(all(r[m]["unique_WT_exact"] for m in native_txs) for r in summary["targets"])
    summary["output_bytes"]=sum(p.stat().st_size for p in out.rglob("*") if p.is_file())
    summary["failure_reasons"]=failures
    summary["replay_qualified"]=not failures
    with (out/"wt_result.json").open("x") as handle:
        json.dump(summary,handle,indent=2)
    print(json.dumps({k:v for k,v in summary.items() if k!="targets"},indent=2))
    if failures:
        raise RuntimeError("WT replay not qualified; no mutant inference allowed")


if __name__=="__main__":
    main()
