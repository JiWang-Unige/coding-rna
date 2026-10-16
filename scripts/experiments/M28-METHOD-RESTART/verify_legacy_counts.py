#!/usr/bin/env python3
"""Verify the sampler reconstruction against archived M25R label aggregates."""
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
import census
from src.foundation_probe import train_generanno_structural_heads as M

def main():
    config = census.yaml.safe_load(census.CONFIG.read_text())
    records = {}
    for name, item in census.metadata(config).items():
        annotation = M.parse_annotation(item["path"] / "reference.gff3", item["lengths"], protein_coding_only=True)
        by_q = defaultdict(list)
        for tx in M._primary_with_partial(annotation):
            by_q[tx["seqid"]].append(tx)
        records[name] = {"seqs": {q: "N" * n for q,n in item["lengths"].items()},
                         "splits": item["splits"], "transcripts_by_seqid": by_q}
    # Actual legacy constructor; DNA letters do not enter any counted label field.
    dataset = M.OrientationWindowDataset(records,"train",None,6144,.12,0,1536)
    region = np.zeros(3,dtype=np.int64)
    boundary = np.zeros(4,dtype=np.int64)
    phase = np.zeros(3,dtype=np.int64)
    for _,r,b,p,mask in dataset.examples:
        region += np.bincount(r,minlength=3)
        boundary += b[mask].sum(axis=0,dtype=np.int64)
        phase += np.bincount(p[(p>0)&mask]-1,minlength=3)
    reconstructed = {"region_counts":region.tolist(),
                     "boundary_positive":dict(zip(M.BOUNDARY_NAMES,map(int,boundary))),
                     "phase_counts":phase.tolist(),"train_windows":len(dataset.examples)}
    historical = json.loads((census.ROOT / "reports/M25R-GENERANNO-1P2B-STRUCTURAL-HEADS-s0/train_summary.json").read_text())
    comparison = {key: reconstructed[key] == historical[key] for key in reconstructed}
    payload = {"status":"matched" if all(comparison.values()) else "mismatch",
               "method":"unmodified current legacy dataset constructor; length-matched dummy DNA; label construction unchanged",
               "not_verified":"No original per-window identity trace exists here. Equal aggregates are not an identity proof.",
               "setaria_access":False,"neural_forward_or_training":False,
               "reconstructed":reconstructed,"archived":{key:historical[key] for key in reconstructed},
               "comparison":comparison}
    out = census.ROOT / "outputs/M28-METHOD-RESTART-R1/legacy_counts"
    out.mkdir(parents=True,exist_ok=True)
    with (out / "summary.json").open("x") as f:
        json.dump(payload,f,indent=2)
        f.write("\n")
    print(json.dumps(payload,indent=2))
    if not all(comparison.values()):
        raise RuntimeError("Reconstructed labels disagree with historical aggregates; do not attribute historical failure.")
if __name__ == "__main__":
    main()
