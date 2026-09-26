#!/usr/bin/env python3
"""One frozen R2 panel, epoch 1, 158 original windows; phase gate only.

No training, candidate regeneration, threshold search, or Setaria access.
Original-path reconciliation for the entire panel precedes B evaluation.
"""
import json
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import yaml

import phase_panel_preflight as P
import redecode_error_decomposition as D

ROOT = P.ROOT
PANEL = ROOT / "outputs/M25R-E1-PHASE-GATE-PANEL-R2"
OUT = ROOT / "outputs/M25R-E1-PHASE-GATE-BYPASS-R2"
CHECKPOINT = ROOT / "outputs/M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0/checkpoints/epoch_1.pt"
CONFIG = ROOT / "configs/M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0.yaml"
NAMES = D.m25.BOUNDARY_NAMES


def read_json(path):
    with path.open() as handle:
        return json.load(handle)


def orf_checks(sequence, cds):
    coding = "".join(sequence[a:b] for a, b in cds)
    return {"minimum_length": len(coding) >= 6, "frame_length": len(coding) % 3 == 0,
            "start_ATG": coding.startswith("ATG"), "terminal_stop": coding[-3:] in D.m25.STOP_CODONS,
            "no_internal_stop": not any(coding[i:i+3] in D.m25.STOP_CODONS
                                        for i in range(3, len(coding)-3, 3))}


def evaluate_fixed(trace, sequence, boundary_scores, phase_classes, thresholds, bypass=False, reverse=False):
    """Use frozen runs; original motif picker and unchanged downstream grammar.

    boundary_scores: position->4-vector of original float16-logit sigmoid.
    phase_classes: position->argmax of same float16 phase logits.
    Dict indexing intentionally rejects any missing, unbudgeted score position.
    """
    sequence = sequence.upper()
    runs = trace["runs"]
    anchors = {"start": [runs[0][0]], "stop": [runs[-1][1]-3],
               "donor": [left[1] for left, right in zip(runs, runs[1:])],
               "acceptor": [right[0]-2 for left, right in zip(runs, runs[1:])]}
    candidates, picks = {}, {}
    for channel, name in enumerate(NAMES):
        candidates[name] = [D.motif_positions(sequence, name, a) for a in anchors[name]]
        probabilities = {p: values[channel] for p, values in boundary_scores.items()}
        picks[name] = [D.m25._select_motif(sequence, name, a, probabilities, D.RADIUS_BP, False, reverse)
                       for a in anchors[name]]
        if any(pick is None for pick in picks[name]):
            raise AssertionError("frozen candidate no longer has a motif")
    if candidates != trace["candidate_positions"]:
        raise AssertionError("candidate enumeration differs from frozen trace")
    chosen = {name: [int(p) for p, _score in values] for name, values in picks.items()}
    if chosen != trace["chosen_positions"]:
        raise AssertionError("motif choice differs from frozen trace")
    if any(d+2 > a for d, a in zip(chosen["donor"], chosen["acceptor"])):
        raise AssertionError("frozen chosen splice order invalid")
    starts = chosen["start"] + [a+2 for a in chosen["acceptor"]]
    ends = chosen["donor"] + [chosen["stop"][0]+3]
    cds = list(zip(starts, ends))
    if any(a >= b for a, b in cds):
        raise AssertionError("frozen CDS interval nonpositive")
    scores = {name: min((float(s) for _p, s in values), default=1.0) for name, values in picks.items()}
    result = {"chosen_positions": chosen, "phase_checks": [], "boundary_scores": scores,
              "terminal": None, "model": None, "ORF_checks": None, "threshold_checks": None}
    phases, coding_length = [], 0
    for index, (a, b) in enumerate(cds):
        expected = 0 if index == 0 else (3 - coding_length % 3) % 3
        predicted = int(phase_classes[a])
        result["phase_checks"].append({"position": a, "expected": expected+1, "predicted": predicted})
        phases.append(expected)
        if not bypass and predicted != expected+1:
            result["terminal"] = "phase_check"
            return result
        coding_length += b-a
    result["ORF_checks"] = orf_checks(sequence, cds)
    # Diagnostic booleans may expose both gates; terminal ordering stays R4.
    result["threshold_checks"] = {name: scores[name] >= float(thresholds[name]) for name in NAMES}
    if not all(result["ORF_checks"].values()):
        result["terminal"] = "complete_ORF_internal_stop_check"
        return result
    if not all(result["threshold_checks"].values()):
        result["terminal"] = "boundary_threshold_filter"
        return result
    result["terminal"] = "emitted"
    result["model"] = {"cds": cds, "phase": phases, "start_codon": (starts[0], starts[0]+3),
                       "stop_codon": (ends[-1]-3, ends[-1]), "boundary_scores": scores,
                       "region_span": (runs[0][0], runs[-1][1])}
    return result


def assert_original(trace, result):
    expected = trace["failure_stage"] or "emitted"
    if result["terminal"] != expected or result["phase_checks"] != trace["phase_checks"]:
        raise AssertionError("original terminal or saved phase prefix mismatch")
    if trace.get("emitted_model"):
        if D.model_key(result["model"]) != D.model_key(trace["emitted_model"]):
            raise AssertionError("emitted control differs, including saved boundary scores")
        if tuple(result["model"]["region_span"]) != tuple(trace["emitted_model"]["region_span"]):
            raise AssertionError("emitted control region span mismatch")


def validate_manifest(manifest, coverage):
    panel = manifest["panel"]
    if len(panel) != 128 or len({P.identity(r) for r in panel}) != 128:
        raise AssertionError("frozen 128-lineage manifest not intact")
    if dict(Counter(r["panel_group"] for r in panel)) != P.GROUPS:
        raise AssertionError("panel group counts changed")
    if coverage["windows"] != json.loads(json.dumps(P.coverage(panel))):
        raise AssertionError("coverage no longer matches frozen score positions")
    if len(coverage["windows"]) != 158:
        raise AssertionError("approved 158-window scope changed")
    if any((r["species"], r["seqid"]) not in P.SCOPE for r in panel):
        raise AssertionError("panel outside A/rice development scope")
    original_tuple = read_json(P.SOURCE / "epoch_1/diagnostic.json")["frozen_tuple"]
    if manifest["frozen_tuple"] != original_tuple or original_tuple["epoch"] != 1 or original_tuple["enumeration_order"] != 601:
        raise AssertionError("frozen epoch/tuple changed")
    source_traces = {P.identity(t): t for t in P.read_jsonl(P.SOURCE / "epoch_1/candidate_lineages.jsonl")}
    for row in panel:
        if source_traces[P.identity(row)] != row["trace"]:
            raise AssertionError("panel source trace changed")
    return panel, original_tuple


def collect_scores(panel, windows, sequences, tokenizer, forward, device, k):
    import torch
    needed = defaultdict(set)
    for row in panel:
        key = row["species"], row["seqid"], row["strand"]
        needed[key].update(row["required_score_positions"])
    boundary, phase = defaultdict(dict), defaultdict(dict)
    (OUT / "scores").mkdir()
    start_time = time.monotonic()
    with torch.no_grad(), (OUT / "window_runs.jsonl").open("x") as log:
        for index, window in enumerate(windows, 1):
            key = window["species"], window["seqid"], window["strand"]
            sequence = sequences[key]
            start = window["start"]
            end = min(start + P.WINDOW, len(sequence))
            if not 0 <= start < end:
                raise AssertionError("window out of original sequence")
            padded = sequence[start:end] + "A" * (P.WINDOW - (end-start))
            ids = D._tokenize_window(tokenizer, D._clean(padded), P.WINDOW, k).unsqueeze(0).to(device)
            attention = torch.ones_like(ids)
            nucleotide = torch.from_numpy(D.m25._one_hot(padded)).unsqueeze(0).to(device)
            with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
                outputs = forward(ids, attention, nucleotide)
            arrays = [v[0, :end-start].float().cpu().numpy().astype(np.float16) for v in outputs]
            if any(not np.isfinite(a).all() for a in arrays):
                raise AssertionError("nonfinite model scores")
            score_name = f"{index:03d}_{key[1]}_{'plus' if key[2]=='+' else 'minus'}_{start}.npz"
            np.savez(OUT / "scores" / score_name, region=arrays[0], boundary=arrays[1], phase=arrays[2])
            probabilities = D.m25._sigmoid(arrays[1])
            classes = arrays[2].argmax(axis=-1)
            for position in needed[key]:
                if start <= position < end:
                    boundary[key][position] = probabilities[position-start]
                    phase[key][position] = int(classes[position-start])
            record = {**window, "score_file": score_name, "real_end": end,
                      "elapsed_since_first_window": time.monotonic()-start_time}
            log.write(json.dumps(record) + "\n")
            log.flush()
            print(f"window={index}/{len(windows)} elapsed={record['elapsed_since_first_window']:.1f}s", flush=True)
    for key, positions in needed.items():
        if positions != set(boundary[key]) or positions != set(phase[key]):
            raise AssertionError("missing necessary model scores")
    return boundary, phase


def summarize(rows):
    strata = defaultdict(Counter)
    for row in rows:
        labels = [row["panel_group"], row["panel_group"]+"|"+row["species"],
                  row["panel_group"]+"|"+row["strand"],
                  row["panel_group"]+"|"+("single_CDS" if row["CDS_count"]==1 else "multi_CDS")]
        if row["panel_group"] == "exact_phase_failure":
            labels.append(row["panel_group"]+"|"+("conflict" if row["phase_failure_token_label_conflict"] else "no_conflict"))
        for label in labels:
            strata[label]["n"] += 1
            strata[label][row["B"]["terminal"]] += 1
    return dict(strata)


def main():
    import torch
    import transformers
    import peft
    started = time.monotonic()
    if (OUT / "resolved_run.json").exists() or (OUT / "scores").exists():
        raise FileExistsError("isolated run already has artifacts; no overwrite/retry")
    manifest = read_json(PANEL / "frozen_panel.json")
    coverage = read_json(PANEL / "window_coverage.json")
    panel, thresholds = validate_manifest(manifest, coverage)
    config = yaml.safe_load(CONFIG.read_text())
    if config["seed"] != 0 or config["model"]["window_bp"] != 6144:
        raise AssertionError("frozen config changed")
    if (torch.__version__, transformers.__version__, peft.__version__) != ("2.5.1", "4.49.0", "0.19.1"):
        raise AssertionError("original R4 inference environment changed")
    sequences, forward_sequences = {}, {}
    for species, seqid in sorted(P.SCOPE):
        fasta = ROOT / "data/m1_screen" / species / "genome.fa"
        sequence = D.screen_data.read_fasta(str(fasta))[seqid].upper()
        forward_sequences[seqid] = sequence
        sequences[(species, seqid, "+")] = sequence
        sequences[(species, seqid, "-")] = D.m25.reverse_complement(sequence)
    OUT.mkdir(exist_ok=True)
    D.save_json_atomic(OUT / "resolved_run.json", {"checkpoint": str(CHECKPOINT), "config": str(CONFIG),
        "frozen_panel": str(PANEL / "frozen_panel.json"), "thresholds": thresholds, "windows": 158,
        "training": False, "setaria_files_read": False, "candidate_generation": False,
        "gate_change": "phase_rejection_only", "input_versions": {"torch": torch.__version__,
        "transformers": transformers.__version__, "peft": peft.__version__}})
    tokenizer, backbone, heads, forward, device, k = D.load_inference_model(config)
    D.load_checkpoint(CHECKPOINT, backbone, heads)
    for module in (backbone, heads):
        module.eval()
        for parameter in module.parameters():
            parameter.requires_grad_(False)
    print("frozen_model_loaded", flush=True)
    boundary, phase = collect_scores(panel, coverage["windows"], sequences, tokenizer, forward, device, k)
    originals = {}
    with (OUT / "A_replay.jsonl").open("x") as log:
        for row in panel:
            key = row["species"], row["seqid"], row["strand"]
            result = evaluate_fixed(row["trace"], sequences[key], boundary[key], phase[key], thresholds,
                                    reverse=key[2]=="-")
            assert_original(row["trace"], result)
            originals[P.identity(row)] = result
            log.write(json.dumps({"identity": P.identity(row), "A": result})+"\n")
            log.flush()
    D.save_json_atomic(OUT / "replay_check.json", {"all_128_reconciled": True, "emitted_controls": 32,
                                                   "B_evaluated_before_A_complete": False})
    print("A_REPLAY_128_PASS", flush=True)
    results = []
    predictions = {"A": defaultdict(list), "B": defaultdict(list)}
    with (OUT / "paired_results.jsonl").open("x") as log:
        for row in panel:
            key = row["species"], row["seqid"], row["strand"]
            b = evaluate_fixed(row["trace"], sequences[key], boundary[key], phase[key], thresholds,
                               bypass=True, reverse=key[2]=="-")
            a = originals[P.identity(row)]
            if row["panel_group"] == "emitted_control" and D.model_key(a["model"]) != D.model_key(b["model"]):
                raise AssertionError("phase bypass changed an emitted control")
            record = {k: v for k, v in row.items() if k not in {"trace", "required_score_positions", "phase_evidence"}}
            record.update(A=a, B=b)
            results.append(record)
            log.write(json.dumps(record)+"\n")
            for arm, outcome in (("A", a), ("B", b)):
                if outcome["model"] is not None:
                    predictions[arm][key[1]].append(D.m25._map_model_from_orientation(outcome["model"], len(sequences[key]), key[2]))
    audits = {}
    for arm in ("A", "B"):
        path = OUT / f"{arm}_panel_predictions.gff3"
        D.m25._write_gff3(str(path), predictions[arm], "M25R_phase_panel")
        audits[arm] = D.audit_gff3(path, forward_sequences, {k: len(v) for k,v in forward_sequences.items()},
                                 sum(len(v) for v in predictions[arm].values()))
        if audits[arm]["invalid_transcripts"] or audits[arm]["audit_coverage"] != 1.0:
            raise AssertionError("panel GFF3 contradicts passed frozen structural checks")
        D.save_json_atomic(OUT / f"{arm}_structural_audit.json", audits[arm])
    counts = summarize(results)
    D.save_json_atomic(OUT / "summary.json", {"status": "COMPLETED_FROZEN_PANEL_DIAGNOSTIC", "counts": counts,
        "elapsed_seconds": time.monotonic()-started, "windows_forwarded": 158, "A_replay_pass": 128,
        "controls_unchanged": 32, "setaria_files_read": False, "weights_updated": False,
        "threshold_search": False, "full_pipeline_replay": False, "population_estimate": False,
        "exploratory_flags_not_gates": {"exact_at_least_32": counts["exact_phase_failure"].get("emitted",0)>=32,
                "nonexact_at_most_8": counts["nonexact_phase_failure"].get("emitted",0)<=8},
        "gpu_max_allocated_bytes": torch.cuda.max_memory_allocated()})
    (OUT / "STATUS").write_text("COMPLETED_FROZEN_PANEL_DIAGNOSTIC\n")
    print(json.dumps(counts, indent=2), flush=True)


if __name__ == "__main__":
    main()
