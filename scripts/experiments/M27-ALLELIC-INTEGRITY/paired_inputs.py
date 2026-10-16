#!/usr/bin/env python3
"""Construct separately edited whole-chr22 inputs for the fixed M27 WT-conditional panel."""
import argparse
import json
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from prepare_reference import LENGTH, bases, ordered_positions

ROOT=Path(__file__).resolve().parents[3]


def mutant(genome,row,allele):
    position=row["position_1based"]-1
    alt=row["synonymous_alt" if allele=="syn" else "stop_alt"]
    positions=ordered_positions(row["CDS_chain_0based_halfopen"],row["strand"])
    start=3*(row["codon_index_1based"]-1)
    codon_positions=positions[start:start+3]
    if genome[position]!=row["reference_base"] or bases(genome,codon_positions,row["strand"])!=row["WT_codon"]:
        raise ValueError("Registered reference base/codon differs")
    if alt==genome[position] or position not in codon_positions:
        raise ValueError("Registered SNV is not a change in its target codon")
    edited=genome[:position]+alt+genome[position+1:]
    changed=np.flatnonzero(np.frombuffer(genome.encode("ascii"),dtype=np.uint8)!=np.frombuffer(edited.encode("ascii"),dtype=np.uint8))
    if changed.tolist()!=[position]:
        raise ValueError("Input is not exactly the single registered SNV")
    expected=row["synonymous_codon" if allele=="syn" else "stop_codon"]
    if bases(edited,codon_positions,row["strand"])!=expected:
        raise ValueError("Strand-aware mutant codon differs from preregistered pair")
    return edited


def build(out):
    preparation=json.loads((ROOT/"reports/M27-ALLELIC-INTEGRITY/preparation.json").read_text())
    wt=json.loads((ROOT/"reports/M27-ALLELIC-INTEGRITY/wt_result.json").read_text())
    if not wt["replay_qualified"]:
        raise ValueError("WT replay is not qualified")
    wt_rows={r["gene_id"]:r for r in wt["targets"]}
    with (ROOT/"data/m27_grch38_chr22/chr22.unmasked.fa").open() as handle:
        if next(handle).strip()!=">chr22":
            raise ValueError("Expected the fixed single chr22 FASTA")
        genome="".join(line.strip() for line in handle)
    if len(genome)!=LENGTH or set(genome)-set("ACGTN"):
        raise ValueError("Frozen chr22 length or uppercase sequence-only policy differs")
    inputs=out/"inputs"
    inputs.mkdir()
    records,cases=[],[]
    for number,row in enumerate(preparation["selected"],1):
        methods=[m for m in ("ANNEVO","Tiberius") if wt_rows[row["gene_id"]][m]["unique_WT_exact"]]
        records.append({**row,"ordinal":number,
            "method_status":{m:("paired_native_planned" if m in methods else "WT_not_exact_conditional_response_not_measured")
                             for m in ("ANNEVO","Tiberius")}})
        if not methods:
            continue
        for allele in ("syn","PTC"):
            case_id=f"{number:02d}_{allele}"
            sequence=mutant(genome,row,allele)
            path=inputs/(case_id+".fa")
            with path.open("x") as handle:
                handle.write(">chr22\n")
                for start in range(0,len(sequence),60):
                    handle.write(sequence[start:start+60]+"\n")
            cases.append({"case_id":case_id,"ordinal":number,"allele":allele,"gene_id":row["gene_id"],
                          "methods":methods,"fasta":str(path.relative_to(ROOT)),"length":len(sequence),
                          "position_1based":row["position_1based"],"reference_base":row["reference_base"],
                          "alternate_base":row["synonymous_alt" if allele=="syn" else "stop_alt"]})
    manifest={"stage":"paired_native_inputs_no_model_result","registered_loci":records,"cases":cases,
              "planned_method_cases":{m:sum(m in c["methods"] for c in cases) for m in ("ANNEVO","Tiberius")},
              "single_SNV_per_whole_chromosome":True,"source_WT_commit":"223902893a27d73a51fccb592277ee23a5657c8b"}
    if manifest["planned_method_cases"]!={"ANNEVO":24,"Tiberius":20} or len(records)!=21:
        raise ValueError("Frozen WT-conditional panel differs; do not reselect loci")
    with (out/"input_manifest.json").open("x") as handle:
        json.dump(manifest,handle,indent=2)
    print(json.dumps({"registered":len(records),"distinct_mutant_FASTAs":len(cases),
                      "planned_method_cases":manifest["planned_method_cases"]},indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    build(args.output_dir.resolve())
