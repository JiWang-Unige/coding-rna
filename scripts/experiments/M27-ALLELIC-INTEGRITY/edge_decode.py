#!/usr/bin/env python3
"""One registered local transition intervention; native ANNEVO source is untouched."""
import argparse
import json
import sys
import subprocess
from pathlib import Path
import numpy as np
import paired_run as P

POS = 15695751
SPAN = (15690077, 15719777)
LIMIT = 2 * 1024**3
TOLERANCE = 1e-3
CLASSES = ("INTERGENIC", "CODING_EXON_0", "CODING_EXON_2", "CODING_EXON_1",
           "INTRON_0", "INTRON_2", "INTRON_1", "DSS_0", "DSS_2", "DSS_1",
           "ASS_0", "ASS_2", "ASS_1", "START", "END")
CONTEXT = {}


def target_index(region):
    a, b, seqid, strand = region[:4]
    if a is None or seqid != "chr22" or strand != 1 or not a <= POS < b:
        return None
    if not a <= SPAN[0] < SPAN[1] <= b:
        raise RuntimeError("Target region does not include both frozen anchors")
    return POS - a


def release_edge(matrices, codes, index, states):
    source, dest = states["CDS1_TA"], states["CDS2"]
    ordinary = matrices[0, states["CDS1"], dest]
    if codes[index] != 0 or not np.isneginf(matrices[0, source, dest]) or not np.isfinite(ordinary):
        raise RuntimeError("Registered A-specific forbidden edge mapping differs")
    special = matrices[0].copy()
    special[source, dest] = ordinary
    changed = np.argwhere(special != matrices[0])
    if changed.tolist() != [[source, dest]]:
        raise RuntimeError("Intervention changed more than one matrix entry")
    new_codes = codes.copy()
    new_codes[index] = len(matrices)
    return np.concatenate((matrices, special[None]), axis=0), new_codes


def components(emissions, columns, path, matrices, codes, start=1, end=None):
    # Native DP starts with dp[0,intergenic]=0; no emission is scored at t=0.
    sites = np.arange(max(1, start), len(path) if end is None else end)
    e = float(emissions[sites, columns[path[sites]]].sum(dtype=np.float64))
    terms = matrices[codes[sites], path[sites-1], path[sites]]
    bad = sites[~np.isfinite(terms)].tolist()
    transition = None if bad else float(terms.sum(dtype=np.float64))
    return {"positions": len(sites), "emission": e, "transition": transition,
            "total": None if bad else e+transition, "impossible_transition_indices": bad}


def budget(reserve=0):
    used = int(subprocess.check_output(["du", "-sb", str(TRACE.parent)], text=True).split()[0])
    if used+reserve > LIMIT:
        raise RuntimeError("Single-edge stage exceeds 2 GiB limit")
    P.storage(reserve)
    return used


def observed_batch(regions, *args):
    result = []
    for region in regions:
        CONTEXT.clear()
        CONTEXT.update(region=region[:4], index=target_index(region), passes=0)
        result.extend(ORIGINAL_BATCH([region], *args))
    return result


def observed_viterbi(predictions, sequence, states, count, columns, min_intron_length, **kwargs):
    if CONTEXT.get("index") is not None:
        i = CONTEXT["index"]
        if sequence[i-2:i+1] != "TAA":
            raise RuntimeError("Registered codon is not TAA in actual decoding sequence")
        CONTEXT.update(states=states, columns=columns, minimum=min_intron_length)
        CONTEXT["passes"] += 1
    return ORIGINAL_VITERBI(predictions, sequence, states, count, columns, min_intron_length, **kwargs)


def observed_core(emissions, matrices, codes, original, precision):
    if CONTEXT.get("index") is None:
        return original(emissions, matrices, codes)
    states, groups, i = CONTEXT["states"], CONTEXT["columns"], CONTEXT["index"]
    released, new_codes = release_edge(matrices, codes, i, states)
    active_m, active_c = (released, new_codes) if ARM == "relax" else (matrices, codes)
    path = original(emissions, active_m, active_c)
    state_columns = np.full(emissions.shape[1], -1, dtype=np.int32)
    reduced = np.empty((len(path), 15), dtype=emissions.dtype)
    for c, name in enumerate(CLASSES):
        ids = groups[name]
        reduced[:, c] = emissions[:, ids[0]]
        for state in ids:
            if not np.array_equal(emissions[:, state], reduced[:, c]):
                raise RuntimeError("Native emission class mapping differs")
            state_columns[state] = c
    if np.any(state_columns < 0):
        raise RuntimeError("Unmapped native states")
    arrays = dict(emissions=reduced, columns=state_columns, path=path,
                  matrices=matrices, codes=codes)
    budget(sum(x.nbytes for x in arrays.values())+8192)
    prefix = TRACE/f"pass_{CONTEXT['passes']:02d}"
    with prefix.with_suffix(".npz").open("xb") as handle:
        np.savez(handle, **arrays)
    a = CONTEXT["region"][0]
    item = {"region": CONTEXT["region"], "pass": CONTEXT["passes"],
            "min_intron_length": CONTEXT["minimum"], "core_precision": precision,
            "states": states, "edited_position_0based": POS, "region_index": i,
            "edge": ["CDS1_TA", "CDS2"], "ordinary_score": float(released[-1, states["CDS1_TA"], states["CDS2"]]),
            "edge_used": bool(path[i-1] == states["CDS1_TA"] and path[i] == states["CDS2"]),
            "full_region": components(reduced, state_columns, path, active_m, active_c),
            "target_span": components(reduced, state_columns, path, active_m, active_c, SPAN[0]-a, SPAN[1]-a),
            "path_under_native_constraints": components(reduced, state_columns, path, matrices, codes)}
    P.write_json(prefix.with_suffix(".json"), item)
    return path


def observed32(e, m, c):
    return observed_core(e, m, c, ORIGINAL32, "float32")


def observed64(e, m, c):
    return observed_core(e, m, c, ORIGINAL64, "float64")


def observed_write(regions, *args, **kwargs):
    if sum(target_index(r) is not None for r in regions) != 1:
        raise RuntimeError("Expected exactly one registered forward target region")
    P.write_json(TRACE/"regions.json", [r[:4] for r in regions])
    return ORIGINAL_WRITE(regions, *args, **kwargs)


def main():
    global TRACE, ARM, ORIGINAL_BATCH, ORIGINAL_VITERBI, ORIGINAL32, ORIGINAL64, ORIGINAL_WRITE
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--arm", choices=("baseline", "relax"), required=True)
    parser.add_argument("--trace-dir", type=Path, required=True)
    args, remaining = parser.parse_known_args()
    TRACE, ARM = args.trace_dir.resolve(), args.arm
    if not TRACE.is_dir() or list(TRACE.glob("pass_*")):
        raise RuntimeError("Fresh existing trace directory required")
    P.STORAGE.extend([P.ROOT/"outputs/M27-NATIVE-CONFIRM-R1", TRACE.parent])
    budget()
    sys.path.insert(0, str(P.ROOT/"refs/repos/annevo-2026"))
    from src import HMM as H
    from src import gene_decoding as G
    import decoding
    ORIGINAL_BATCH, ORIGINAL_VITERBI, ORIGINAL_WRITE = G.process_gene_segment_batch, G.viterbi_decoding, G.decode_and_write
    ORIGINAL32, ORIGINAL64 = H._viterbi_core_numba_float32, H._viterbi_core_numba_float64
    G.process_gene_segment_batch, G.viterbi_decoding, G.decode_and_write = observed_batch, observed_viterbi, observed_write
    H._viterbi_core_numba_float32, H._viterbi_core_numba_float64 = observed32, observed64
    sys.argv = [sys.argv[0], *remaining]
    decoding.main()
    if not list(TRACE.glob("pass_*.json")):
        raise RuntimeError("No target path observed")
    budget()


if __name__ == "__main__":
    main()
