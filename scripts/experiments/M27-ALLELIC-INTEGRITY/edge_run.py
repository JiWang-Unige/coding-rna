#!/usr/bin/env python3
"""Two fixed whole-chromosome CPU decodes, gated on exact native replay."""
import argparse
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path
import numpy as np
import paired_run as P
import paired_response as R
from wt_compare import parsed, key, reference_key
from confirm_run import chain_difference, compare_target, counter_records
from edge_decode import POS, SPAN, LIMIT, TOLERANCE, release_edge, components
from prepare_reference import LENGTH

CONF = P.ROOT/"outputs/M27-NATIVE-CONFIRM-R1"
OLD = P.ROOT/"outputs/M27-PAIRED-NATIVE-R1"


def read(path):
    return json.loads(path.read_text())


def trace_comparison(out):
    traces = {arm: [read(p) for p in sorted((out/arm).glob("pass_*.json"))]
              for arm in ("baseline", "relax")}
    if read(out/"baseline/regions.json") != read(out/"relax/regions.json"):
        raise RuntimeError("Candidate region selection changed")
    arrays = {arm: [dict(np.load(p)) for p in sorted((out/arm).glob("pass_*.npz"))]
              for arm in traces}
    a, b = arrays["baseline"][-1], arrays["relax"][-1]
    ta, tb = traces["baseline"][-1], traces["relax"][-1]
    if ta["region"] != tb["region"] or any(not np.array_equal(x["emissions"], a["emissions"])
            for arm in arrays for x in arrays[arm]):
        raise RuntimeError("Fixed target region/emission evidence changed")
    passes = {arm: [x["min_intron_length"] for x in traces[arm]] for arm in traces}
    same = passes["baseline"] == passes["relax"]
    start = ta["region"][0]
    outside = (np.arange(len(a["path"]))+start < SPAN[0]) | (np.arange(len(a["path"]))+start >= SPAN[1])
    names = []
    for t, x in ((ta, a), (tb, b)):
        inverse = np.array([n for n, i in sorted(t["states"].items(), key=lambda z: z[1])])
        names.append(inverse[x["path"]])
    result = {"native_pass_sequences": passes, "same_pass_sequence": same,
              "emissions_equal_all_passes": True, "regions_equal": True,
              "outside_target_hidden_state_differences": int(np.sum((names[0] != names[1]) & outside)),
              "target_emission_difference_relax_minus_baseline": tb["target_span"]["emission"]-ta["target_span"]["emission"],
              "full_region_emission_difference": tb["full_region"]["emission"]-ta["full_region"]["emission"],
              "final_pass": {"baseline": ta, "relax": tb},
              "all_passes": traces, "common_relaxed_graph_scores": None,
              "interpretation": "single_transition_pipeline_effect" if same else "controller_mediated_pipeline_effect"}
    if same:
        if ta["states"] != tb["states"] or any(not np.array_equal(a[k], b[k]) for k in ("matrices", "codes", "columns")):
            raise RuntimeError("Common native graph differs despite identical pass sequence")
        m, c = release_edge(a["matrices"], a["codes"], POS-start, ta["states"])
        result["common_relaxed_graph_scores"] = {
            arm: components(x["emissions"], x["columns"], x["path"], m, c)
            for arm, x in (("baseline", a), ("relax", b))}
    return result


def main(out):
    out = out.resolve()
    if out.name != "M27-SINGLE-EDGE-R1" or (out/"edge_result.json").exists() or (out/"baseline").exists():
        raise RuntimeError("Fresh dedicated single-edge output required")
    P.STORAGE.extend([CONF, out])
    result = {"stage": "single_edge_implementation_attribution", "status": "RUNNING",
              "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "commands": [],
              "training_runs": 0, "neural_forward_calls": 0, "allocated_GPUs": 0,
              "Setaria_accessed": False, "selective_risk_validated": False,
              "numerical_tolerance_nats": TOLERANCE, "arms": {}}
    began = time.monotonic()
    try:
        P.storage(LIMIT)
        manifest = read(OLD/"input_manifest.json")
        row = next(r for r in manifest["registered_loci"] if r["ordinal"] == 1)
        case = next(c for c in manifest["cases"] if c["case_id"] == "01_PTC")
        refkey = reference_key(row)
        if (row["gene_id"], row["strand"], refkey[2][0][0], refkey[2][-1][1], case["position_1based"]) != (
                "ENSG00000198062.16", "+", *SPAN, POS+1):
            raise RuntimeError("Frozen locus identity differs")
        positions = R.positions({"strand": "+", "CDS": refkey[2]})
        end = 3*row["codon_index_1based"]
        if positions[end-3:end] != [POS-2, POS-1, POS] or not any(s+30 <= POS-2 and POS+31 < e for s, e, _ in refkey[2]):
            raise RuntimeError("PTC is not registered internal in-frame codon; no alternate edge")
        fasta = P.ROOT/case["fasta"]
        with fasta.open() as handle:
            if next(handle).strip() != ">chr22":
                raise RuntimeError("Frozen input header differs")
            genome = "".join(line.strip() for line in handle)
        if len(genome) != LENGTH or genome[POS-2:POS+1] != "TAA":
            raise RuntimeError("Frozen PTC genome differs")
        wt = set(map(key, parsed(P.ROOT/"outputs/M27-WT-REPLAY-R1/ANNEVO/native.gff3")))
        foreign = R.foreign_regions(row, R.reference_coding(P.ROOT/"data/m27_grch38_chr22/gencode.v49.chr22.gtf"))
        original = parsed(CONF/"ANNEVO/01_PTC/prediction.gff3")
        outputs, responses = {}, {}
        for arm in ("baseline", "relax"):
            base = out/arm
            base.mkdir()
            P.run_command([sys.executable, "-u", Path(__file__).with_name("edge_decode.py"),
                           "--arm", arm, "--trace-dir", base, "-g", fasta,
                           "-p", CONF/"ANNEVO/01_PTC/native_scores.h5", "-o", base/"prediction.gff3",
                           "-t", "8", "--min_intron_length", "20", "--min_prot_length", "100", "--show_log"],
                          base/"decode.log", result, P.ROOT/"refs/repos/annevo-2026")
            if "Process failed:" in (base/"decode.log").read_text(errors="replace"):
                raise RuntimeError("Native worker exception; no retry")
            outputs[arm] = parsed(base/"prediction.gff3")
            responses[arm] = R.classify(row, outputs[arm], genome, wt, foreign)
            result["arms"][arm] = responses[arm]
            P.write_json(base/"response.json", responses[arm])
            if arm == "baseline":
                gate = {"whole_chr": chain_difference(original, outputs[arm]),
                        "target": compare_target(row, original, outputs[arm], genome, wt, foreign)}
                gate["passed"] = gate["whole_chr"]["equal"] and gate["target"]["equal"] and responses[arm]["category"] == "BYPASS"
                result["baseline_replay"] = gate
                P.write_json(out/"baseline_check.json", gate)
                if not gate["passed"]:
                    raise RuntimeError("Same-DNA cache replay not exact; intervention prohibited")
        trace = trace_comparison(out)
        result["trace_comparison"] = trace
        result["whole_chr_difference"] = chain_difference(outputs["baseline"], outputs["relax"])
        background = {}
        for arm in outputs:
            background[arm] = Counter(map(key, outputs[arm]))-Counter(responses[arm]["assigned_keys"])
        result["non_target_chains"] = {"equal": background["baseline"] == background["relax"],
            "lost": counter_records(background["baseline"]-background["relax"]),
            "added": counter_records(background["relax"]-background["baseline"])}
        restored = (responses["relax"]["category"] == "WT_CHAIN_RETAINED"
                    and responses["relax"]["assigned_keys"] == [refkey]
                    and responses["relax"].get("coding_status", {}).get("internal_stop_count") == 1)
        result["same_coordinate_chain_with_PTC_restored"] = restored
        candidate = (restored and trace["final_pass"]["relax"]["edge_used"] and trace["same_pass_sequence"]
                     and trace["outside_target_hidden_state_differences"] == 0
                     and result["non_target_chains"]["equal"]
                     and trace["target_emission_difference_relax_minus_baseline"] > TOLERANCE)
        result["clean_local_conflict_candidate"] = candidate
        result["next_decision"] = ("PROPOSE_REFERENCE_FREE_SIGNAL_AND_BIDIRECTIONAL_VALIDATION_NOT_YET_VALIDATED"
                                   if candidate else "CLOSE_FIXED_PILOT_NO_ADDITIONAL_EDGE_RESCUE")
        result["M27_bytes"] = P.storage()
        import subprocess
        result["stage_bytes"] = int(subprocess.check_output(["du", "-sb", str(out)], text=True).split()[0])
        if result["stage_bytes"] > LIMIT:
            raise RuntimeError("Stage storage limit exceeded")
        result["status"] = "COMPLETED"
    except BaseException as error:
        result.update(status="FAILED", error=repr(error))
        raise
    finally:
        result["seconds"] = time.monotonic()-began
        P.write_json(out/"edge_result.json", result)
        print(json.dumps({k: v for k, v in result.items() if k in (
            "status", "error", "seconds", "same_coordinate_chain_with_PTC_restored", "clean_local_conflict_candidate", "next_decision")}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    main(parser.parse_args().output_dir)
