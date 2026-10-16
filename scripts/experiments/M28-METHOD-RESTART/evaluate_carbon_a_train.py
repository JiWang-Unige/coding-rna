"""Separate, fixed Carbon-A TRAIN comparison using the existing coordinate ruler."""
import argparse
import json
import os
import time
import run_carbon_a_train as N
import evaluate_native_train_baselines as B
import evaluate_full_grammar_original_train as E

C = E.C
ROOT = N.ROOT
OUT = ROOT / "outputs/M28-CARBON-A-TRAIN-EVAL"
CONTEXT24_OUT = ROOT / "outputs/M28-CARBON-A-CONTEXT24-TRAIN-EVAL"
OPERATING_POINTS = {"CarbonA_confidence_ge_0.1_or_NA": 0.1, "CarbonA_unfiltered": 0.0}


def gene_confidences(path):
    result = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            row = line.rstrip("\n").split("\t")
            if row[2] != "gene":
                continue
            attrs = dict(item.split("=", 1) for item in row[8].split(";"))
            value = attrs["confidence"]
            result[attrs["ID"]] = None if value == "NA" else float(value)
    return result


def confidence_selection(txs, confidences, threshold):
    return [tx for tx in txs if confidences[tx["gene_id"]] is None
            or confidences[tx["gene_id"]] >= threshold]


def main(context_bp=98304):
    started = time.monotonic()
    out = {98304: OUT, 24576: CONTEXT24_OUT}[context_bp]
    runs = N.completed_runs(context_bp)  # All six raw predictions finish before reference/metrics access.
    assert (E.OUT / "STATUS").read_text().strip() == "COMPLETED"
    control = json.loads((E.OUT / "summary.json").read_text())
    assert control["completed"] and set(control["chromosomes"]) == set(B.EXPECTED)
    assert control["reference_policy"] == "M26_CDS_assessable_longest_eligible_CDS_tie_transcript_ID"
    assert (B.OUT / "STATUS").read_text().strip() == "COMPLETED"
    baselines = json.loads((B.OUT / "summary.json").read_text())
    assert baselines["completed"] and baselines["reference_primary_chains"] == 22688
    assert baselines["input_mode"] == N.N.INPUT_MODE
    rows = {name: [] for name in OPERATING_POINTS}
    for species, seqid, length in N.SCOPE:
        sequences = C.M.read_allowed(N.N.input_fasta(seqid), {seqid: length})
        data = ROOT / "data/m1_screen" / species
        coding = C.S.E.parse_annotation(data / "reference.gff3", {seqid: length}, protein_coding_only=True)
        annotation = C.S.E.parse_annotation(data / "reference.gff3", {seqid: length})
        with (data / "reference.gff3").open() as handle:
            _, _, raw_cds, _ = C.P.raw_features(handle, seqid)
        refs, isoforms = C.select_reference(coding, raw_cds, sequences)
        assert len(refs) == B.EXPECTED[seqid] and len({C.S.key(tx) for tx in refs}) == B.EXPECTED[seqid]
        assert all(C.R7.phase_chain_consistent(tx) for tx in refs)
        unsupported = C.Counter()
        overlaps = C.M.C.overlapping_ids(refs)
        for tx in refs:
            reasons, _ = C.M.chain_checks(tx, sequences[seqid])
            if tx["id"] in overlaps:
                reasons.append("same_strand_primary_overlap")
            if len(tx["CDS"]) > 128:
                reasons.append("over_128_exons")
            if tx["CDS"][-1][1] - tx["CDS"][0][0] > C.M.WIDTH:
                reasons.append("longer_than_native_window")
            unsupported.update(set(reasons))
        native = runs[seqid]
        gff = ROOT / native["raw_GFF"]
        predicted = C.S.E.parse_annotation(gff, {seqid: length})
        confidences = gene_confidences(gff)
        assert set(confidences) == set(predicted["genes"])
        txs = C.S.E.primary_transcripts(predicted, complete_only=False)
        assert len(txs) == native["output_check"]["genes"]
        for name, threshold in OPERATING_POINTS.items():
            selected = confidence_selection(txs, confidences, threshold)
            row = C.score(selected, refs, isoforms, length, annotation)
            row["prediction_nonexact_primary"] = row.pop("reference_nonexact_primary")
            row.update(species=species, seqid=seqid, length=length, reference_primary_records=len(refs),
                       reference_assessable_isoform_records=len(isoforms), retained_unsupported_reference_strata=dict(unsupported),
                       raw_prediction_path=native["raw_GFF"], native_run_receipt=native["accepted_run_receipt"],
                       confidence_threshold=threshold, NA_confidence_retained=True,
                       unfiltered_prediction_records=len(txs), filtered_out_records=len(txs) - len(selected),
                       phase_chain_consistent_prediction_records=sum(C.R7.phase_chain_consistent(tx) for tx in selected),
                       native_runtime_seconds=native["seconds"], codon_table_warning=native["codon_warning"],
                       source_raw_wrapper_seconds_missing=native.get("source_wrapper_seconds_missing_due_to_failed_validation", False))
            assert B.reference_metadata(row) == B.reference_metadata(control["chromosomes"][seqid])
            rows[name].append(row)
        print(json.dumps({"scored_TRAIN": seqid, "reference_primary_chains": len(refs)}), flush=True)
    report = dict(completed=True, status="separate_CarbonA_original_TRAIN_comparison_not_DEV_or_generalization",
                  job_id=os.environ["SLURM_JOB_ID"], fixed_revision=N.REVISION,
                  context_bp=context_bp, context_tokens=context_bp // 6, overlap_bp=6144,
                  window_policy_modified=context_bp != 98304,
                  shared_backbone_source_revision="b0483c23b6b63787b61a6d3a204a9b517d6ba345",
                  reference_policy=control["reference_policy"], reference_primary_chains=22688,
                  reference_used_in_inference=False, unsupported_reference_strata_retained=True,
                  training_eligibility_does_not_filter_reference=True, DEV_used=False, sealed_data_used=False,
                  optimizer_steps_this_run=0, input_mode=N.N.INPUT_MODE, frozen_DEV_family_changed=False,
                  control_receipt=str((E.OUT / "summary.json").relative_to(ROOT)),
                  native_baseline_receipt=str((B.OUT / "summary.json").relative_to(ROOT)), methods={})
    comparators = {"M28_final4608": control, **baselines["methods"]}
    for name, method_rows in rows.items():
        species_names = sorted({row["species"] for row in method_rows})
        per_species = {species: E.aggregate([row for row in method_rows if row["species"] == species]) for species in species_names}
        pooled = E.aggregate(method_rows)
        f1s = [per_species[species]["exact_chain"]["f1"] for species in species_names]
        macro, worst = sum(f1s) / len(f1s), min(f1s)
        differences = {}
        for comparator, result in comparators.items():
            differences[comparator] = dict(
                per_species={species: 100 * (per_species[species]["exact_chain"]["f1"] - result["per_species"][species]["exact_chain"]["f1"]) for species in species_names},
                pooled=100 * (pooled["exact_chain"]["f1"] - result["pooled"]["exact_chain"]["f1"]),
                macro_species=100 * (macro - result["macro_species_exact_chain_F1"]),
                worst_species=100 * (worst - result["worst_species_exact_chain_F1"]))
        report["methods"][name] = dict(chromosomes={row["seqid"]: row for row in method_rows}, per_species=per_species,
            pooled=pooled, macro_species_exact_chain_F1=macro, worst_species_exact_chain_F1=worst,
            CarbonA_minus_comparator_exact_chain_F1_pp=differences)
    report["evaluation_seconds"] = time.monotonic() - started
    with (out / "summary.json").open("x") as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({"completed": True, "reference_primary_chains": 22688}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--context-bp", type=int, choices=(98304, 24576), default=98304)
    main(parser.parse_args().context_bp)
