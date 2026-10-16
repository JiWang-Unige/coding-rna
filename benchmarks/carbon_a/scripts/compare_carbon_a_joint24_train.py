#!/usr/bin/env python3
"""Fixed original-TRAIN reference intersections: Carbon-A versus completed joint24.

Descriptive recall accounting only. Both native confidence points are retained;
reference bins are fixed before execution, without model or threshold selection.
"""
import json
import os
import time
from collections import Counter

import common_ruler as C
import diagnose_native_train_reference as D

ROOT = D.ROOT
OUT = ROOT / "outputs/M28-CARBON-A-JOINT24-PAIRED"
JOINT = ROOT / "outputs/M28-FULL-GRAMMAR-JOINT24-TRAIN-DIAGNOSTIC"
CARBON = ROOT / "outputs/M28-CARBON-A-TRAIN-EVAL"
POINTS = {"CarbonA_confidence_ge_0.1_or_NA": 0.1, "CarbonA_unfiltered": 0.0}
FIELDS = ("reference_chains", "both", "joint_only", "CarbonA_only", "neither",
          "joint_hits", "CarbonA_hits", "CarbonA_minus_joint_primary_hits")


def completed(folder, job):
    value = D.completed(folder)
    if (folder / "JOBID").read_text().strip() != job or value["job_id"] != job:
        raise ValueError("Completed source job identity differs: " + str(folder))
    return value


def pair(rows, method):
    counts = Counter()
    for row in rows:
        a, b = row["joint_hit"], row["carbon_hits"][method]
        counts["both" if a and b else "CarbonA_only" if b else "joint_only" if a else "neither"] += 1
    n = len(rows)
    joint, carbon = counts["both"] + counts["joint_only"], counts["both"] + counts["CarbonA_only"]
    return dict(reference_chains=n, **{key: counts[key] for key in ("both", "joint_only", "CarbonA_only", "neither")},
                joint_hits=joint, CarbonA_hits=carbon, CarbonA_minus_joint_primary_hits=carbon-joint,
                joint_reference_recall=joint/n if n else None, CarbonA_reference_recall=carbon/n if n else None,
                CarbonA_minus_joint_recall_pp=100*(carbon-joint)/n if n else None)


def describe(rows, method):
    value = pair(rows, method)
    groups = {
        "CDS_segments": {"one": [r for r in rows if r["CDS_segments"] == 1],
                         "two_to_five": [r for r in rows if 2 <= r["CDS_segments"] <= 5],
                         "six_or_more": [r for r in rows if r["CDS_segments"] >= 6]},
        "CDS_genomic_span_bp": {"at_most_24576": [r for r in rows if r["span_bp"] <= 24576],
                                "over_24576": [r for r in rows if r["span_bp"] > 24576]},
    }
    value["reference_defined_strata"] = {name: {key: pair(part, method) for key, part in bins.items()}
                                         for name, bins in groups.items()}
    for bins in value["reference_defined_strata"].values():
        for field in FIELDS:
            if sum(part[field] for part in bins.values()) != value[field]:
                raise ValueError("Reference-defined partitions must retain every primary chain")
    return value


def confidence_by_gene(path):
    result = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            row = line.rstrip("\n").split("\t")
            if row[2] == "gene":
                attrs = dict(item.split("=", 1) for item in row[8].split(";"))
                value = attrs["confidence"]
                result[attrs["ID"]] = None if value == "NA" else float(value)
    return result


def checked_predictions(path, species, seqid, length, refs, official, threshold=None):
    txs = C.S.E.primary_transcripts(C.S.E.parse_annotation(path, {seqid: length}), complete_only=False)
    if threshold is not None:
        confidence = confidence_by_gene(path)
        if {t["gene_id"] for t in txs} != set(confidence):
            raise ValueError("Native confidence and parsed gene identities differ")
        txs = [t for t in txs if confidence[t["gene_id"]] is None or confidence[t["gene_id"]] >= threshold]
    keys = {(species,) + C.S.key(t) for t in txs}
    if (D.exact_counts(keys, refs) != {k: official["exact_chain"][k] for k in ("tp", "predicted", "reference")}
            or len(txs) != official["prediction_records"]
            or len(txs)-len(keys) != official["duplicate_chain_records"]):
        raise ValueError("GFF unique-chain replay differs from completed evaluation: " + str(path))
    return keys


def main():
    began = time.perf_counter()
    joint, carbon = completed(JOINT, "13453780"), completed(CARBON, "13503419")
    if (joint["reference_policy"] != carbon["reference_policy"]
            or any(s["DEV_used"] or s["sealed_data_used"] or s["optimizer_steps_this_run"]
                   or s["reference_used_in_inference"] or not s["unsupported_reference_strata_retained"]
                   or not s["training_eligibility_does_not_filter_reference"] for s in (joint, carbon))
            or set(joint["chromosomes"]) != set(D.EXPECTED) or set(carbon["methods"]) != set(POINTS)
            or carbon["fixed_revision"] != "f6789ee79e6ef397de333577ffec1cb734496620"
            or carbon["reference_primary_chains"] != 22688 or carbon["frozen_DEV_family_changed"]):
        raise ValueError("Two completed fixed original-TRAIN evaluations required")
    reference = completed(D.OUT, "13461637")
    if (reference["reference_policy"] != joint["reference_policy"] or reference["reference_primary_chains"] != 22688
            or reference["source_original_evaluation_job"] != "13443329" or reference["reference_filters_added"]
            or not reference["exact_hit_counts_reproduce_completed_evaluations"]
            or reference["DEV_used"] or reference["sealed_data_used"] or reference["optimizer_steps_this_run"]):
        raise ValueError("Complete original primary-reference profile required")
    rows = D.read_rows(D.OUT / "references.jsonl")
    scope = {(r["species"], seqid): D.EXPECTED[seqid] for seqid, r in joint["chromosomes"].items()}
    if (Counter((r["species"], r["seqid"]) for r in rows) != Counter(scope)
            or len({D.reference_key(r) for r in rows}) != 22688
            or any(r["CDS_segments"] != len(r["CDS"]) or not r["CDS"]
                   or r["span_bp"] != r["CDS"][-1][1]-r["CDS"][0][0] for r in rows)):
        raise ValueError("Unique primary identities or saved reference geometry differ")
    for seqid, official_joint in joint["chromosomes"].items():
        species, length = official_joint["species"], official_joint["length"]
        part = [r for r in rows if (r["species"], r["seqid"]) == (species, seqid)]
        refs = {D.reference_key(r) for r in part}
        joint_path = ROOT / official_joint["raw_prediction_path"]
        joint_keys = checked_predictions(joint_path, species, seqid, length, refs, official_joint)
        carbon_keys = {}
        for method, threshold in POINTS.items():
            official = carbon["methods"][method]["chromosomes"][seqid]
            if ((official["species"], official["seqid"], official["length"]) != (species, seqid, length)
                    or official["confidence_threshold"] != threshold or not official["NA_confidence_retained"]
                    or any(official[k] != official_joint[k] for k in ("reference_primary_records",
                           "reference_assessable_isoform_records", "retained_unsupported_reference_strata"))):
                raise ValueError("Paired population, geometry or fixed operating point differs")
            carbon_keys[method] = checked_predictions(ROOT / official["raw_prediction_path"], species, seqid,
                                                     length, refs, official, threshold)
        for row in part:
            key = D.reference_key(row)
            row["joint_hit"] = key in joint_keys
            row["carbon_hits"] = {method: key in keys for method, keys in carbon_keys.items()}
        print(json.dumps(dict(replayed_TRAIN=seqid, primary_references=len(part))), flush=True)
    comparisons = {}
    for method in POINTS:
        views = {"pooled": rows, **{s: [r for r in rows if r["species"] == s] for s in joint["per_species"]},
                 **{q: [r for r in rows if r["seqid"] == q] for q in D.EXPECTED}}
        measured = {scope: describe(part, method) for scope, part in views.items()}
        for scope, value in measured.items():
            a = joint["pooled"] if scope == "pooled" else joint["per_species"][scope] if scope in joint["per_species"] else joint["chromosomes"][scope]
            b = carbon["methods"][method]
            b = b["pooled"] if scope == "pooled" else b["per_species"][scope] if scope in b["per_species"] else b["chromosomes"][scope]
            if (value["joint_hits"] != a["exact_chain"]["tp"] or value["CarbonA_hits"] != b["exact_chain"]["tp"]
                    or value["reference_chains"] != a["exact_chain"]["reference"]
                    or value["reference_chains"] != b["exact_chain"]["reference"]):
                raise ValueError("Paired counts must reconstruct chromosome, species and pooled endpoints")
        comparisons[method] = measured
    report = dict(completed=True, job_id=os.environ["SLURM_JOB_ID"],
                  status="descriptive_original_TRAIN_CarbonA_joint24_reference_pairs_not_DEV_or_generalization",
                  source_evaluation_jobs=dict(joint="13453780", CarbonA="13503419"), source_reference_job="13461637",
                  reference_policy=joint["reference_policy"], reference_primary_chains=22688,
                  operating_points=POINTS, comparisons=comparisons,
                  exact_hits_unique_prediction_counts_and_duplicates_reproduce_completed_evaluations=True,
                  every_reference_partition_reconciled=True, reference_filters_added=False,
                  DNA_read=False, real_reference_GFF_read=False, saved_primary_reference_profile_read=True,
                  prediction_GFF_read=True, checkpoint_read=False, native_forward=False, predictions_modified=False,
                  optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False, frozen_DEV_family_changed=False,
                  scientific_winner_or_threshold_selected=False, seconds=time.perf_counter()-began,
                  limitations=["Reference intersections and bins are recall-only; no subset prediction precision or F1.",
                               "Primary hits are annotated references, not novel genes; auxiliary exposure remains unresolved.",
                               "CDS-span24576 is an encoding-context comparison, not a hard whole-contig grammar limit.",
                               "CDS-segment and span views each partition all references; their contributions must not be added together.",
                               "Shared main weights with differing contexts, heads, learned auxiliary scores and filtering do not isolate decoder causality.",
                               "Original TRAIN only; no independent generalization, significance, ensemble construction or new frozen DEV candidate."])
    with (OUT / "summary.json").open("x") as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps(dict(completed=True, reference_primary_chains=22688)), flush=True)


if __name__ == "__main__":
    main()
