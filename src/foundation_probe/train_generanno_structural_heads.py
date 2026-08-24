"""M25 GENERanno structural heads: one frozen raw-FASTA-to-coding-GFF3 experiment.

This is deliberately an experiment path, not a reusable annotation framework.  It keeps the
M19 GENERanno/LoRA and FP-aware region contracts, adds nucleotide identity plus boundary/phase
heads, calibrates one deterministic decoder on pooled plant validation data, and only then runs
the unchanged full/ablation decoders on Setaria FASTA.
"""

import argparse
import hashlib
import itertools
import json
import os
import re
import shutil
import sys
import time
from collections import defaultdict

import numpy as np


ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCRIPTS = os.path.join(ROOT, "scripts")
for _path in (ROOT, SCRIPTS):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from eval_structure_diagnostic import parse_annotation, primary_transcripts  # noqa: E402
from src.foundation_probe.train_generanno_lora_3class import (  # noqa: E402
    MODEL,
    _clean,
    _tokenize_window,
)
from src.foundation_probe.train_probe_head import _ConvLSTMHead, fp_penalty  # noqa: E402
from src.screen_anchor import data as D  # noqa: E402


I, C, G = 0, 1, 2
BOUNDARY_NAMES = ("start", "stop", "donor", "acceptor")
STOP_CODONS = {"TAA", "TAG", "TGA"}
MOTIFS = {"start": {"ATG"}, "stop": STOP_CODONS, "donor": {"GT"}, "acceptor": {"AG"}}
TRANSITIONS = {
    (I, C): ("start",),
    (G, C): ("start", "acceptor"),
    (C, G): ("stop", "donor"),
    (C, I): ("stop",),
}
_NON_ID = re.compile(r"[^A-Za-z0-9_.-]")


def reverse_complement(sequence):
    """Return the DNA reverse complement, preserving N and mapping other IUPAC symbols."""
    table = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn")
    return sequence.translate(table)[::-1]


def map_rc_interval(start, end, sequence_length):
    """Map one 0-based half-open interval between a sequence and its reverse complement."""
    return sequence_length - end, sequence_length - start


def select_primary_chromosomes(sequences, seqids):
    missing = [seqid for seqid in seqids if seqid not in sequences]
    if missing:
        raise ValueError(f"primary chromosome seqids absent from FASTA: {missing}")
    return {seqid: sequences[seqid] for seqid in seqids}


def region_state_path(region_logits, genic_threshold=0.5):
    """Decode M19-order region logits into I/C/G states with the frozen M25 rule."""
    logits = np.asarray(region_logits, dtype=np.float32)
    shifted = logits - logits.max(axis=-1, keepdims=True)
    prob = np.exp(shifted)
    prob /= prob.sum(axis=-1, keepdims=True)
    state = np.zeros(logits.shape[0], dtype=np.int8)
    genic = 1.0 - prob[:, I]
    on = genic >= genic_threshold
    state[on] = np.where(prob[on, C] >= prob[on, G], C, G)
    return state


def transition_candidates(states):
    """Return allowed coding-boundary events at every adjacent region-state transition.

    Each row is ``(boundary_coordinate, left_state, right_state, allowed_events)``.  The
    coordinate is the 0-based boundary immediately before the right state.
    """
    states = np.asarray(states)
    return [
        (position, int(states[position - 1]), int(states[position]), TRANSITIONS[pair])
        for position in range(1, len(states))
        if (pair := (int(states[position - 1]), int(states[position]))) in TRANSITIONS
    ]


def _oriented_transcript(transcript, sequence_length, strand):
    tx = {
        "id": transcript.get("id", "tx"),
        "strand": "+",
        "partial": bool(transcript.get("partial", False)),
        "CDS": [],
        "exon": [],
    }
    if strand == "+":
        tx["CDS"] = [(int(a), int(b), str(p)) for a, b, p in transcript.get("CDS", [])]
        tx["exon"] = [(int(a), int(b), str(p)) for a, b, p in transcript.get("exon", [])]
    else:
        tx["CDS"] = [(*map_rc_interval(int(a), int(b), sequence_length), str(p))
                     for a, b, p in transcript.get("CDS", [])]
        tx["exon"] = [(*map_rc_interval(int(a), int(b), sequence_length), str(p))
                      for a, b, p in transcript.get("exon", [])]
    tx["CDS"].sort()
    tx["exon"].sort()
    return tx


def build_strand_targets(sequence_length, transcripts, strand):
    """Build orientation-specific region/boundary/phase targets.

    ``transcripts`` are M24-parser transcript dictionaries.  Only the requested genomic strand
    contributes.  Minus-strand intervals are returned in reverse-complement coordinates.  Region
    supervision includes partial/overlapping primary models; their boundary and phase spans are
    masked.  Boundary channels label motif starts: ATG, stop codon, GT donor and AG acceptor.
    """
    selected = [
        _oriented_transcript(tx, sequence_length, strand)
        for tx in transcripts
        if tx.get("strand") == strand and tx.get("CDS")
    ]
    spans = [(tx["CDS"][0][0], tx["CDS"][-1][1]) for tx in selected]
    ambiguous = set()
    for left in range(len(selected)):
        for right in range(left + 1, len(selected)):
            if max(spans[left][0], spans[right][0]) < min(spans[left][1], spans[right][1]):
                ambiguous.update((left, right))

    region = np.zeros(sequence_length, dtype=np.int8)
    boundary = np.zeros((sequence_length, 4), dtype=np.uint8)
    phase = np.zeros(sequence_length, dtype=np.int8)
    structural_mask = np.ones(sequence_length, dtype=bool)
    for tx in selected:
        exons = tx["exon"] or tx["CDS"]
        span_start = min(a for a, _b, _p in exons)
        span_end = max(b for _a, b, _p in exons)
        region[span_start:span_end] = G
    for tx in selected:
        for start, end, _gff_phase in tx["CDS"]:
            region[start:end] = C
    for index, tx in enumerate(selected):
        cds = tx["CDS"]
        exons = tx["exon"] or cds
        span_start = min(a for a, _b, _p in exons)
        span_end = max(b for _a, b, _p in exons)
        if tx["partial"] or index in ambiguous:
            structural_mask[span_start:span_end] = False
            continue
        boundary[cds[0][0], 0] = 1
        boundary[cds[-1][1] - 3, 1] = 1
        coding_offset = 0
        for cds_index, (start, end, gff_phase) in enumerate(cds):
            segment_phase = int(gff_phase) if gff_phase in ("0", "1", "2") else (3 - coding_offset % 3) % 3
            positions = np.arange(end - start, dtype=np.int64)
            phase[start:end] = 1 + ((positions + segment_phase) % 3).astype(np.int8)
            coding_offset += end - start
            if cds_index + 1 < len(cds):
                next_start = cds[cds_index + 1][0]
                boundary[end, 2] = 1
                boundary[next_start - 2, 3] = 1
    return {
        "region": region,
        "boundary": boundary,
        "phase": phase,
        "structural_mask": structural_mask,
    }


def _sigmoid(values):
    values = np.asarray(values, dtype=np.float32)
    return 1.0 / (1.0 + np.exp(-np.clip(values, -30.0, 30.0)))


def _select_motif(sequence, event, anchor, probabilities, radius, ablation, reverse_mapped):
    motifs = MOTIFS[event]
    motif_length = len(next(iter(motifs)))
    positions = [
        position
        for position in range(max(0, anchor - radius), min(len(sequence) - motif_length, anchor + radius) + 1)
        if sequence[position:position + motif_length] in motifs
    ]
    if not positions:
        return None
    coordinate_key = (lambda value: -value) if reverse_mapped else (lambda value: value)
    if ablation:
        position = min(positions, key=lambda value: (abs(value - anchor), coordinate_key(value)))
    else:
        position = min(
            positions,
            key=lambda value: (-float(probabilities[value]), abs(value - anchor), coordinate_key(value)),
        )
    return position, float(probabilities[position])


def _coding_blocks(states):
    blocks = []
    position = 0
    while position < len(states):
        while position < len(states) and states[position] == I:
            position += 1
        block_start = position
        while position < len(states) and states[position] != I:
            position += 1
        if block_start == position:
            continue
        cds_runs = []
        cursor = block_start
        while cursor < position:
            while cursor < position and states[cursor] != C:
                cursor += 1
            start = cursor
            while cursor < position and states[cursor] == C:
                cursor += 1
            if start < cursor:
                cds_runs.append((start, cursor))
        if cds_runs:
            blocks.append(cds_runs)
    return blocks


def decode_orientation(
    sequence,
    region_logits,
    boundary_logits,
    phase_logits,
    thresholds,
    radius=6,
    ablation=False,
    reverse_mapped=False,
):
    """Decode complete canonical coding models in one plus-oriented sequence.

    The same function handles a genomic forward sequence and its reverse complement.  Returned
    coordinates are 0-based half-open in that orientation.  The ablation uses identical region
    states and motif grammar but ignores boundary/phase scores and their thresholds.
    """
    sequence = sequence.upper()
    states = region_state_path(region_logits, float(thresholds["region"]))
    boundary_probability = _sigmoid(boundary_logits)
    phase_class = np.asarray(phase_logits).argmax(axis=-1)
    models = []
    for runs in _coding_blocks(states):
        first_start, final_end = runs[0][0], runs[-1][1]
        if first_start == 0 or final_end == len(states):
            continue
        if "start" not in TRANSITIONS.get((int(states[first_start - 1]), int(states[first_start])), ()):
            continue
        if "stop" not in TRANSITIONS.get((int(states[final_end - 1]), int(states[final_end])), ()):
            continue
        chosen = {}
        start_pick = _select_motif(
            sequence, "start", first_start, boundary_probability[:, 0], radius, ablation, reverse_mapped
        )
        stop_pick = _select_motif(
            sequence, "stop", final_end - 3, boundary_probability[:, 1], radius, ablation, reverse_mapped
        )
        if start_pick is None or stop_pick is None:
            continue
        chosen["start"] = [start_pick]
        chosen["stop"] = [stop_pick]
        chosen["donor"], chosen["acceptor"] = [], []
        failed = False
        for left, right in zip(runs, runs[1:]):
            donor = _select_motif(
                sequence, "donor", left[1], boundary_probability[:, 2], radius, ablation, reverse_mapped
            )
            acceptor = _select_motif(
                sequence, "acceptor", right[0] - 2, boundary_probability[:, 3], radius, ablation, reverse_mapped
            )
            if donor is None or acceptor is None or donor[0] + 2 > acceptor[0]:
                failed = True
                break
            chosen["donor"].append(donor)
            chosen["acceptor"].append(acceptor)
        if failed:
            continue
        if not ablation and any(
            score < float(thresholds[name])
            for name, picks in chosen.items()
            for _position, score in picks
        ):
            continue

        cds = []
        for index in range(len(runs)):
            start = start_pick[0] if index == 0 else chosen["acceptor"][index - 1][0] + 2
            end = stop_pick[0] + 3 if index == len(runs) - 1 else chosen["donor"][index][0]
            if end <= start:
                failed = True
                break
            cds.append((start, end))
        if failed:
            continue
        phases = []
        coding_length = 0
        for index, (start, end) in enumerate(cds):
            expected = 0 if index == 0 else (3 - coding_length % 3) % 3
            phases.append(expected)
            if not ablation and int(phase_class[start]) != expected + 1:
                failed = True
                break
            coding_length += end - start
        if failed:
            continue
        coding_sequence = "".join(sequence[start:end] for start, end in cds)
        if (
            len(coding_sequence) < 6
            or len(coding_sequence) % 3
            or not coding_sequence.startswith("ATG")
            or coding_sequence[-3:] not in STOP_CODONS
            or any(coding_sequence[offset:offset + 3] in STOP_CODONS
                   for offset in range(3, len(coding_sequence) - 3, 3))
        ):
            continue
        boundary_scores = {
            name: min((score for _position, score in picks), default=1.0)
            for name, picks in chosen.items()
        }
        models.append({
            "cds": cds,
            "phase": phases,
            "start_codon": (start_pick[0], start_pick[0] + 3),
            "stop_codon": (stop_pick[0], stop_pick[0] + 3),
            "boundary_scores": boundary_scores,
            "region_span": (runs[0][0], runs[-1][1]),
        })
    return models


def _one_hot(sequence):
    encoded = np.frombuffer(sequence.upper().encode("ascii", "replace"), dtype=np.uint8)
    lut = np.full(256, 4, dtype=np.int64)
    for base, index in (("A", 0), ("C", 1), ("G", 2), ("T", 3)):
        lut[ord(base)] = index
    indices = lut[encoded]
    result = np.zeros((len(sequence), 5), dtype=np.float32)
    result[np.arange(len(sequence)), indices] = 1.0
    return result


def _primary_with_partial(annotation):
    selected = []
    for tx in primary_transcripts(annotation, complete_only=False):
        copied = dict(tx)
        copied["partial"] = bool(tx.get("partial") or annotation["genes"][tx["gene_id"]].get("partial"))
        selected.append(copied)
    return selected


class OrientationWindowDataset:
    def __init__(self, species, split, tokenizer, window, sample_fraction, seed, limit):
        self.tokenizer = tokenizer
        self.window = window
        descriptors = []
        for species_name, record in species.items():
            for seqid, split_name in record["splits"].items():
                if split_name != split:
                    continue
                for start in range(0, len(record["seqs"][seqid]) - window + 1, window):
                    descriptors.extend(((species_name, seqid, start, "+"), (species_name, seqid, start, "-")))
        rng = np.random.default_rng(seed)
        if sample_fraction < 1.0 and descriptors:
            count = max(1, int(round(len(descriptors) * sample_fraction)))
            indices = sorted(rng.choice(len(descriptors), size=count, replace=False))
            descriptors = [descriptors[index] for index in indices]
        if limit is not None and limit > 0:
            descriptors = descriptors[:limit]

        by_seqid = defaultdict(list)
        for descriptor in descriptors:
            by_seqid[descriptor[:2]].append(descriptor)
        self.examples = []
        for (species_name, seqid), rows in by_seqid.items():
            record = species[species_name]
            sequence = record["seqs"][seqid]
            transcripts = record["transcripts_by_seqid"][seqid]
            for strand in ("+", "-"):
                strand_rows = [row for row in rows if row[3] == strand]
                if not strand_rows:
                    continue
                targets = build_strand_targets(len(sequence), transcripts, strand)
                oriented_sequence = sequence if strand == "+" else reverse_complement(sequence)
                for _species_name, _seqid, genomic_start, _strand in strand_rows:
                    orient_start = genomic_start if strand == "+" else len(sequence) - genomic_start - window
                    orient_end = orient_start + window
                    self.examples.append((
                        oriented_sequence[orient_start:orient_end],
                        targets["region"][orient_start:orient_end].copy(),
                        targets["boundary"][orient_start:orient_end].copy(),
                        targets["phase"][orient_start:orient_end].copy(),
                        targets["structural_mask"][orient_start:orient_end].copy(),
                    ))

    def __len__(self):
        return len(self.examples)

    def __getitem__(self, index):
        import torch

        sequence, region, boundary, phase, structural_mask = self.examples[index]
        ids = _tokenize_window(self.tokenizer, _clean(sequence), self.window, 6)
        return (
            ids,
            torch.ones_like(ids),
            torch.from_numpy(_one_hot(sequence)),
            torch.from_numpy(region.astype(np.int64)),
            torch.from_numpy(boundary.astype(np.float32)),
            torch.from_numpy(phase.astype(np.int64)),
            torch.from_numpy(structural_mask),
        )


def _collate(batch):
    import torch

    return tuple(torch.stack([row[column] for row in batch]) for column in range(len(batch[0])))


def _f1(predicted, reference):
    matched = len(predicted & reference)
    precision = matched / len(predicted) if predicted else 0.0
    recall = matched / len(reference) if reference else 0.0
    return 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0


def _merge(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def _overlap_length(left, right):
    total = 0
    i = j = 0
    while i < len(left) and j < len(right):
        total += max(0, min(left[i][1], right[j][1]) - max(left[i][0], right[j][0]))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return total


def _validation_metrics(predictions, references, lengths):
    pred_intervals, ref_intervals, pred_chains, ref_chains = set(), set(), set(), set()
    pred_spans, ref_full = defaultdict(list), defaultdict(list)
    for key, models in predictions.items():
        species, seqid = key
        for model in models:
            strand = model["strand"]
            chain = tuple(model["cds"])
            pred_chains.add((species, seqid, strand, chain))
            pred_spans[(species, seqid)].append((chain[0][0], chain[-1][1]))
            pred_intervals.update((species, seqid, strand, start, end) for start, end in chain)
    for key, transcripts in references.items():
        species, seqid = key
        for tx in transcripts:
            chain = tuple((start, end) for start, end, _phase in tx["CDS"])
            ref_chains.add((species, seqid, tx["strand"], chain))
            ref_intervals.update((species, seqid, tx["strand"], start, end) for start, end in chain)
            parts = tx["exon"] or tx["CDS"]
            ref_full[(species, seqid)].append((parts[0][0], parts[-1][1]))
    false_positive_bases = 0
    intergenic_bases = 0
    for key, length in lengths.items():
        pred = _merge(pred_spans[key])
        ref = _merge(ref_full[key])
        predicted_bases = sum(end - start for start, end in pred)
        false_positive_bases += predicted_bases - _overlap_length(pred, ref)
        intergenic_bases += length - sum(end - start for start, end in ref)
    return {
        "exact_CDS_chain_F1": _f1(pred_chains, ref_chains),
        "exact_CDS_interval_F1": _f1(pred_intervals, ref_intervals),
        "intergenic_FPR": false_positive_bases / intergenic_bases if intergenic_bases else 0.0,
        "gene_count_ratio": len(predictions_flat(predictions)) / len(ref_chains) if ref_chains else 0.0,
        "structurally_valid_complete_fraction": 1.0 if predictions_flat(predictions) else 0.0,
    }


def predictions_flat(predictions):
    return [model for models in predictions.values() for model in models]


def _map_model_from_orientation(model, sequence_length, strand):
    mapped = dict(model)
    if strand == "+":
        mapped["cds"] = list(model["cds"])
        mapped["start_codon"] = model["start_codon"]
        mapped["stop_codon"] = model["stop_codon"]
    else:
        paired = [(*map_rc_interval(start, end, sequence_length), phase)
                  for (start, end), phase in zip(model["cds"], model["phase"])]
        paired.sort()
        mapped["cds"] = [(start, end) for start, end, _phase in paired]
        mapped["phase"] = [phase for _start, _end, phase in paired]
        mapped["start_codon"] = map_rc_interval(*model["start_codon"], sequence_length)
        mapped["stop_codon"] = map_rc_interval(*model["stop_codon"], sequence_length)
    mapped["strand"] = strand
    return mapped


def _filter_models(candidate_models, boundary_thresholds):
    return [
        model for model in candidate_models
        if all(model["boundary_scores"][name] >= boundary_thresholds[name] for name in BOUNDARY_NAMES)
    ]


def _write_gff3(path, models_by_seqid, source):
    with open(path, "w") as handle:
        handle.write("##gff-version 3\n")
        serial = 0
        for seqid in sorted(models_by_seqid):
            for model in sorted(models_by_seqid[seqid], key=lambda row: (row["cds"][0][0], row["strand"], row["cds"])):
                serial += 1
                gene_id = f"M25_gene_{serial:07d}"
                transcript_id = f"{gene_id}.t1"
                gene_start, gene_end = model["cds"][0][0], model["cds"][-1][1]
                common = (seqid, source)
                handle.write(f"{common[0]}\t{common[1]}\tgene\t{gene_start + 1}\t{gene_end}\t.\t{model['strand']}\t.\tID={gene_id}\n")
                handle.write(f"{common[0]}\t{common[1]}\tmRNA\t{gene_start + 1}\t{gene_end}\t.\t{model['strand']}\t.\tID={transcript_id};Parent={gene_id}\n")
                for (start, end), phase in zip(model["cds"], model["phase"]):
                    handle.write(f"{common[0]}\t{common[1]}\tCDS\t{start + 1}\t{end}\t.\t{model['strand']}\t{phase}\tParent={transcript_id}\n")
                for feature in ("start_codon", "stop_codon"):
                    start, end = model[feature]
                    handle.write(f"{common[0]}\t{common[1]}\t{feature}\t{start + 1}\t{end}\t.\t{model['strand']}\t0\tParent={transcript_id}\n")


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _save_json(path, payload):
    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--species", nargs=2, required=True)
    parser.add_argument("--setaria-fasta", required=True)
    parser.add_argument("--exp-id", default="M25-GENERANNO-1P2B-STRUCTURAL-HEADS-s0")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--model-name", default=MODEL)
    parser.add_argument("--model-task", default="token_classification", choices=("token_classification", "masked_lm"))
    parser.add_argument("--window", type=int, default=6144)
    parser.add_argument("--sample-fraction", type=float, default=0.12)
    parser.add_argument("--epochs", type=int, default=3)
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--predict-batch", type=int, default=1)
    parser.add_argument("--head-lr", type=float, default=8e-4)
    parser.add_argument("--lora-lr", type=float, default=2e-5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--fp-lambda", type=float, default=2.5)
    parser.add_argument("--lora-r", type=int, default=8)
    parser.add_argument("--lora-alpha", type=int, default=16)
    parser.add_argument("--lora-dropout", type=float, default=0.05)
    parser.add_argument("--lora-targets", default="q_proj,k_proj,v_proj,o_proj")
    parser.add_argument("--attn-implementation", default="sdpa", choices=("eager", "sdpa", "flash_attention_2"))
    parser.add_argument("--bf16", action="store_true")
    parser.add_argument("--limit-train-windows", type=int, default=1536)
    parser.add_argument("--limit-val-windows", type=int, default=768)
    parser.add_argument("--development-smoke", action="store_true")
    args = parser.parse_args()

    if args.epochs != 3:
        raise ValueError("M25 requires exactly three epochs")
    if args.window != 6144:
        raise ValueError("M25 freezes the GENERanno window at 6144 bp")

    import torch
    import torch.nn.functional as F
    import yaml
    from peft import LoraConfig, get_peft_model
    from torch.utils.data import DataLoader
    from transformers import AutoModelForMaskedLM, AutoModelForTokenClassification, AutoTokenizer

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    dtype = torch.bfloat16 if args.bf16 else torch.float32
    use_amp = device == "cuda" and args.bf16
    os.makedirs(args.out_dir, exist_ok=True)
    checkpoint_dir = os.path.join(args.out_dir, "checkpoints")
    os.makedirs(checkpoint_dir, exist_ok=True)

    with open(args.config) as handle:
        frozen_config = yaml.safe_load(handle)
    primary_seqids = frozen_config["data"]["primary_chromosome_seqids"]
    _save_json(os.path.join(args.out_dir, "primary_chromosome_allowlist.json"), primary_seqids)

    config = {
        "experiment_id": args.exp_id,
        "model": args.model_name,
        "species": args.species,
        "setaria_fasta": args.setaria_fasta,
        "seed": args.seed,
        "window": args.window,
        "epochs": 3,
        "early_stopping": False,
        "train_window_cap": args.limit_train_windows,
        "validation_window_cap": args.limit_val_windows,
        "sample_fraction": args.sample_fraction,
        "head_lr": args.head_lr,
        "lora_lr": args.lora_lr,
        "lora": {
            "r": args.lora_r,
            "alpha": args.lora_alpha,
            "dropout": args.lora_dropout,
            "targets": args.lora_targets.split(","),
        },
        "batch_size": args.batch_size,
        "bf16": args.bf16,
        "loss": "M19_region_FP + 2.0*boundary_focal_gamma2 + 0.5*phase_CE",
        "fp_lambda": args.fp_lambda,
        "decode": {
            "region_thresholds": [0.4, 0.5, 0.6],
            "boundary_thresholds": [0.1, 0.2, 0.3, 0.4, 0.5],
            "snap_radius": 6,
            "selection": "chain_F1; interval_F1; lower_FPR; enumeration_order",
            "constraints": {"intergenic_FPR_max": 0.02, "gene_count_ratio": [0.8, 1.2], "valid_fraction_min": 0.99},
        },
    }
    tokenizer = AutoTokenizer.from_pretrained(args.model_name, trust_remote_code=True)
    loader = AutoModelForMaskedLM if args.model_task == "masked_lm" else AutoModelForTokenClassification
    full_model = loader.from_pretrained(
        args.model_name,
        trust_remote_code=True,
        attn_implementation=args.attn_implementation,
    )
    k = int(getattr(full_model, "k", getattr(full_model.config, "k", 6)))
    if k != 6:
        raise ValueError(f"M25 requires 6-mer upsampling, model reports k={k}")
    backbone = full_model.model
    del full_model
    backbone.gradient_checkpointing_enable()
    backbone.to(device=device, dtype=dtype)
    for parameter in backbone.parameters():
        parameter.requires_grad_(False)
    backbone = get_peft_model(backbone, LoraConfig(
        r=args.lora_r,
        lora_alpha=args.lora_alpha,
        target_modules=args.lora_targets.split(","),
        lora_dropout=args.lora_dropout,
        bias="none",
    ))
    if hasattr(backbone, "enable_input_require_grads"):
        backbone.enable_input_require_grads()

    class StructuralHeads(torch.nn.Module):
        def __init__(self, hidden_size):
            super().__init__()
            channels = hidden_size + 5
            self.region = _ConvLSTMHead(channels, n_classes=3)
            self.boundary = torch.nn.Conv1d(channels, 4, kernel_size=9, padding=4)
            self.phase = torch.nn.Conv1d(channels, 4, kernel_size=1)

        def forward(self, per_base, nucleotide):
            features = torch.cat((per_base, nucleotide.to(per_base.dtype)), dim=-1).transpose(1, 2)
            return self.region(features), self.boundary(features).transpose(1, 2), self.phase(features).transpose(1, 2)

    heads = StructuralHeads(int(backbone.config.hidden_size)).to(device)

    species = {}
    for species_path in args.species:
        name = os.path.basename(species_path.rstrip("/"))
        all_sequences = D.read_fasta(os.path.join(species_path, "genome.fa"))
        all_splits = D.assign_splits(list(all_sequences))
        sequences = select_primary_chromosomes(all_sequences, primary_seqids[name])
        splits = {seqid: all_splits[seqid] for seqid in sequences}
        lengths = {seqid: len(sequence) for seqid, sequence in sequences.items()}
        annotation = parse_annotation(os.path.join(species_path, "reference.gff3"), lengths, protein_coding_only=True)
        transcripts = _primary_with_partial(annotation)
        by_seqid = defaultdict(list)
        for transcript in transcripts:
            by_seqid[transcript["seqid"]].append(transcript)
        species[name] = {
            "seqs": sequences,
            "lengths": lengths,
            "annotation": annotation,
            "transcripts": transcripts,
            "transcripts_by_seqid": by_seqid,
            "splits": splits,
        }

    train_dataset = OrientationWindowDataset(
        species, "train", tokenizer, args.window, args.sample_fraction, args.seed, args.limit_train_windows
    )
    validation_dataset = OrientationWindowDataset(
        species, "val", tokenizer, args.window, 1.0, args.seed, args.limit_val_windows
    )
    train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True, collate_fn=_collate)
    validation_loader = DataLoader(validation_dataset, batch_size=args.batch_size, shuffle=False, collate_fn=_collate)

    region_counts = np.zeros(3, dtype=np.int64)
    boundary_positive = np.zeros(4, dtype=np.int64)
    boundary_total = np.zeros(4, dtype=np.int64)
    phase_counts = np.zeros(3, dtype=np.int64)
    for _sequence, region, boundary, phase, structural_mask in train_dataset.examples:
        region_counts += np.bincount(region, minlength=3)
        boundary_positive += boundary[structural_mask].sum(axis=0, dtype=np.int64)
        boundary_total += int(structural_mask.sum())
        valid_phase = phase[(phase > 0) & structural_mask] - 1
        phase_counts += np.bincount(valid_phase, minlength=3)
    if not region_counts.all() or not boundary_positive.all() or not phase_counts.all():
        raise ValueError(
            f"training sample lacks supervised classes: region={region_counts.tolist()} "
            f"boundary={boundary_positive.tolist()} phase={phase_counts.tolist()}"
        )
    region_raw = 1.0 / np.sqrt(region_counts)
    region_weight = torch.tensor(region_raw / region_raw.mean(), dtype=torch.float32, device=device)
    boundary_positive_weight = torch.tensor(
        (boundary_total - boundary_positive) / boundary_total, dtype=torch.float32, device=device
    )
    boundary_negative_weight = torch.tensor(
        boundary_positive / boundary_total, dtype=torch.float32, device=device
    )
    phase_raw = 1.0 / np.sqrt(phase_counts)
    phase_weight = torch.tensor(phase_raw / phase_raw.mean(), dtype=torch.float32, device=device)
    phase_ce_weight = torch.cat((torch.zeros(1, device=device), phase_weight))

    optimizer = torch.optim.Adam([
        {"params": heads.parameters(), "lr": args.head_lr},
        {"params": [parameter for parameter in backbone.parameters() if parameter.requires_grad], "lr": args.lora_lr},
    ])

    def forward(ids, attention, nucleotide):
        hidden = backbone(input_ids=ids, attention_mask=attention).last_hidden_state
        per_base = hidden.repeat_interleave(k, dim=1)
        if per_base.shape[1] != args.window:
            raise ValueError(f"hidden/base alignment failed: {per_base.shape[1]} != {args.window}")
        return heads(per_base, nucleotide)

    def loss_value(region_logits, boundary_logits, phase_logits, region, boundary, phase, structural_mask):
        region_loss = F.cross_entropy(region_logits.reshape(-1, 3), region.reshape(-1), weight=region_weight)
        region_loss = region_loss + args.fp_lambda * fp_penalty(region_logits, region)
        probability = torch.sigmoid(boundary_logits)
        positive = -boundary_positive_weight * boundary * (1.0 - probability).pow(2) * F.logsigmoid(boundary_logits)
        negative = -boundary_negative_weight * (1.0 - boundary) * probability.pow(2) * F.logsigmoid(-boundary_logits)
        mask = structural_mask.unsqueeze(-1).expand_as(boundary)
        boundary_loss = (positive + negative)[mask].mean()
        phase_mask = structural_mask & (phase > 0)
        phase_loss = (
            F.cross_entropy(phase_logits[phase_mask], phase[phase_mask], weight=phase_ce_weight)
            if phase_mask.any()
            else phase_logits.sum() * 0.0
        )
        return region_loss + 2.0 * boundary_loss + 0.5 * phase_loss

    if args.development_smoke:
        smoke_index = next(
            index for index, example in enumerate(train_dataset.examples)
            if example[2].sum() and (example[3] > 0).any()
        )
        batch = _collate([train_dataset[smoke_index]])
        ids, attention, nucleotide, region, boundary, phase, structural_mask = [value.to(device) for value in batch]
        optimizer.zero_grad(set_to_none=True)
        with torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=use_amp):
            outputs = forward(ids, attention, nucleotide)
            loss = loss_value(*outputs, region, boundary, phase, structural_mask)
        if not torch.isfinite(loss):
            raise ValueError(f"development smoke produced non-finite loss: {loss.item()}")
        loss.backward()
        optimizer.step()
        checkpoint_path = os.path.join(checkpoint_dir, "development_smoke.pt")
        torch.save({"heads": heads.state_dict()}, checkpoint_path)
        _save_json(os.path.join(args.out_dir, "development_smoke.json"), {
            "loss": float(loss.item()),
            "train_window": args.window,
            "setaria_inference_run": False,
        })
        print("M25_DEVELOPMENT_SMOKE_DONE", flush=True)
        return

    validation_references = {}
    validation_lengths = {}
    for species_name, record in species.items():
        complete = primary_transcripts(record["annotation"])
        by_seqid = defaultdict(list)
        for transcript in complete:
            if record["splits"][transcript["seqid"]] == "val":
                by_seqid[transcript["seqid"]].append(transcript)
        for seqid, split in record["splits"].items():
            if split == "val":
                validation_references[(species_name, seqid)] = by_seqid[seqid]
                validation_lengths[(species_name, seqid)] = record["lengths"][seqid]

    def predict_sequence(sequence):
        score = [np.zeros((len(sequence), width), dtype=np.float16) for width in (3, 4, 4)]
        heads.eval()
        backbone.eval()
        with torch.no_grad():
            for start in range(0, len(sequence), args.window):
                real_end = min(start + args.window, len(sequence))
                padded = sequence[start:real_end] + "A" * (args.window - (real_end - start))
                ids = _tokenize_window(tokenizer, _clean(padded), args.window, k).unsqueeze(0).to(device)
                attention = torch.ones_like(ids)
                nucleotide = torch.from_numpy(_one_hot(padded)).unsqueeze(0).to(device)
                with torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=use_amp):
                    outputs = forward(ids, attention, nucleotide)
                width = real_end - start
                for target, output in zip(score, outputs):
                    target[start:real_end] = output[0, :width].float().cpu().numpy().astype(np.float16)
        return score

    def validation_scores():
        scores = {}
        for species_name, record in species.items():
            for seqid, split in record["splits"].items():
                if split != "val":
                    continue
                sequence = record["seqs"][seqid]
                scores[(species_name, seqid)] = {}
                for strand, oriented in (("+", sequence), ("-", reverse_complement(sequence))):
                    scores[(species_name, seqid)][strand] = predict_sequence(oriented)
        return scores

    def validation_candidates(scores, region_threshold):
        candidates = {}
        for key, orientation_scores in scores.items():
            species_name, seqid = key
            sequence = species[species_name]["seqs"][seqid]
            models = []
            for strand, oriented in (("+", sequence), ("-", reverse_complement(sequence))):
                region_score, boundary_score, phase_score = orientation_scores[strand]
                decoded = decode_orientation(
                    oriented, region_score, boundary_score, phase_score,
                    {"region": region_threshold, **{name: 0.0 for name in BOUNDARY_NAMES}},
                    reverse_mapped=strand == "-",
                )
                models.extend(_map_model_from_orientation(model, len(sequence), strand) for model in decoded)
            candidates[key] = models
        return candidates

    best = None
    epoch_rows = []
    enumeration_order = 0
    for epoch in range(1, 4):
        heads.train()
        backbone.train()
        total_loss = 0.0
        started = time.time()
        for batch in train_loader:
            ids, attention, nucleotide, region, boundary, phase, structural_mask = [value.to(device) for value in batch]
            optimizer.zero_grad(set_to_none=True)
            with torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=use_amp):
                outputs = forward(ids, attention, nucleotide)
                loss = loss_value(*outputs, region, boundary, phase, structural_mask)
            loss.backward()
            optimizer.step()
            total_loss += float(loss.item())
        checkpoint_path = os.path.join(checkpoint_dir, f"epoch_{epoch}.pt")
        torch.save({
            "epoch": epoch,
            "heads": heads.state_dict(),
            "lora": {key: value.detach().cpu() for key, value in backbone.state_dict().items() if "lora_" in key},
        }, checkpoint_path)
        epoch_row = {
            "epoch": epoch,
            "train_loss": total_loss / max(len(train_loader), 1),
            "seconds": time.time() - started,
            "checkpoint": checkpoint_path,
        }
        epoch_rows.append(epoch_row)

        pooled_validation_scores = validation_scores()
        for region_threshold in config["decode"]["region_thresholds"]:
            candidates = validation_candidates(pooled_validation_scores, region_threshold)
            for boundary_values in itertools.product(config["decode"]["boundary_thresholds"], repeat=4):
                enumeration_order += 1
                boundary_thresholds = dict(zip(BOUNDARY_NAMES, boundary_values))
                predictions = {
                    key: _filter_models(models, boundary_thresholds) for key, models in candidates.items()
                }
                metrics = _validation_metrics(predictions, validation_references, validation_lengths)
                valid = (
                    metrics["intergenic_FPR"] <= 0.020
                    and 0.80 <= metrics["gene_count_ratio"] <= 1.20
                    and metrics["structurally_valid_complete_fraction"] >= 0.99
                )
                if not valid:
                    continue
                row = {
                    "epoch": epoch,
                    "region": region_threshold,
                    **boundary_thresholds,
                    "metrics": metrics,
                    "enumeration_order": enumeration_order,
                    "checkpoint": checkpoint_path,
                }
                rank = (
                    metrics["exact_CDS_chain_F1"],
                    metrics["exact_CDS_interval_F1"],
                    -metrics["intergenic_FPR"],
                    -enumeration_order,
                )
                if best is None or rank > best[0]:
                    best = (rank, row)

    summary = {
        "experiment_id": args.exp_id,
        "train_windows": len(train_dataset),
        "validation_windows_for_loss": len(validation_dataset),
        "region_counts": region_counts.tolist(),
        "boundary_positive": dict(zip(BOUNDARY_NAMES, boundary_positive.tolist())),
        "phase_counts": phase_counts.tolist(),
        "epochs": epoch_rows,
    }
    if best is None:
        summary["selection"] = None
        _save_json(os.path.join(args.out_dir, "train_summary.json"), summary)
        with open(os.path.join(args.out_dir, "STATUS"), "w") as handle:
            handle.write("STOP_M25_BRANCH\n")
        print("STOP_M25_BRANCH", flush=True)
        return

    selection = best[1]
    summary["selection"] = selection
    _save_json(os.path.join(args.out_dir, "validation_decode_parameters.json"), selection)
    checkpoint = torch.load(selection["checkpoint"], map_location="cpu")
    selected_checkpoint = os.path.join(checkpoint_dir, "selected_checkpoint.pt")
    shutil.copyfile(selection["checkpoint"], selected_checkpoint)
    heads.load_state_dict(checkpoint["heads"])
    current = backbone.state_dict()
    current.update(checkpoint["lora"])
    backbone.load_state_dict(current)

    setaria_sequences = select_primary_chromosomes(
        D.read_fasta(args.setaria_fasta), primary_seqids["setaria_viridis"]
    )
    full_models = defaultdict(list)
    ablation_models = defaultdict(list)
    raw_dir = os.path.join(args.out_dir, "raw_scores")
    os.makedirs(raw_dir, exist_ok=True)
    raw_manifest = {}
    thresholds = {name: selection[name] for name in ("region",) + BOUNDARY_NAMES}
    for seqid, sequence in setaria_sequences.items():
        payload = {}
        for strand, oriented in (("+", sequence), ("-", reverse_complement(sequence))):
            region_score, boundary_score, phase_score = predict_sequence(oriented)
            prefix = "forward" if strand == "+" else "reverse_complement"
            payload[f"{prefix}_region"] = region_score
            payload[f"{prefix}_boundary"] = boundary_score
            payload[f"{prefix}_phase"] = phase_score
            full = decode_orientation(
                oriented, region_score, boundary_score, phase_score, thresholds,
                reverse_mapped=strand == "-",
            )
            ablation = decode_orientation(
                oriented, region_score, boundary_score, phase_score, thresholds,
                ablation=True, reverse_mapped=strand == "-",
            )
            full_models[seqid].extend(_map_model_from_orientation(model, len(sequence), strand) for model in full)
            ablation_models[seqid].extend(
                _map_model_from_orientation(model, len(sequence), strand) for model in ablation
            )
        raw_path = os.path.join(raw_dir, f"{_NON_ID.sub('_', seqid)}.npz")
        np.savez_compressed(raw_path, **payload)
        raw_manifest[seqid] = {"path": raw_path, "sha256": _sha256(raw_path)}

    prediction_dir = os.path.join(args.out_dir, "predictions")
    os.makedirs(prediction_dir, exist_ok=True)
    full_path = os.path.join(prediction_dir, "setaria_viridis.full.gff3")
    ablation_path = os.path.join(prediction_dir, "setaria_viridis.ablation.gff3")
    _write_gff3(full_path, full_models, "M25_GENERanno_structural_heads")
    _write_gff3(ablation_path, ablation_models, "M25_GENERanno_region_motif_ablation")
    hashes = {
        "setaria_fasta": _sha256(args.setaria_fasta),
        "selected_checkpoint": _sha256(selected_checkpoint),
        "full_prediction": _sha256(full_path),
        "ablation_prediction": _sha256(ablation_path),
        "raw_scores": raw_manifest,
    }
    summary.update({
        "selection": selection,
        "selected_checkpoint_sha256": hashes["selected_checkpoint"],
        "setaria_full_models": sum(map(len, full_models.values())),
        "setaria_ablation_models": sum(map(len, ablation_models.values())),
        "artifact_hashes": hashes,
    })
    _save_json(os.path.join(args.out_dir, "artifact_hashes.json"), hashes)
    _save_json(os.path.join(args.out_dir, "train_summary.json"), summary)
    resolved_config = frozen_config
    resolved_config["selected_decode"] = selection
    resolved_config["artifact_hashes"] = hashes
    with open(os.path.join(args.out_dir, "config_resolved.yaml"), "w") as handle:
        yaml.safe_dump(resolved_config, handle, sort_keys=False)
    with open(os.path.join(args.out_dir, "STATUS"), "w") as handle:
        handle.write("COMPLETED\n")
    print("M25_TRAIN_AND_SETARIA_INFERENCE_DONE", flush=True)


if __name__ == "__main__":
    main()
