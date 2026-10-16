"""Actual common-parser, exact-key and confidence/NA checks on artificial GFF only."""
import json
import sys
from pathlib import Path
import evaluate_carbon_a_train as E

folder = Path(sys.argv[1])
length = json.loads((folder / "summary.json").read_text())["artificial_record_bp"]
gff = folder / "known_coordinate_fixture.gff"
annotation = E.C.S.E.parse_annotation(gff, {"artificial": length})
txs = E.C.S.E.primary_transcripts(annotation, complete_only=False)
assert len(txs) == 4
expected = {
    "gene000001": ("+", [(1000, 1061, "0"), (1161, 1226, "2")]),
    "gene000002": ("-", [(6000, 6065, "2"), (6165, 6226, "0")]),
    "gene000003": ("+", [(0, 13, "1")]),
    "gene000004": ("-", [(length - 14, length, "2")]),
}
for tx in txs:
    assert (tx["strand"], tx["CDS"]) == expected[tx["gene_id"]]
confidence = E.gene_confidences(gff)
selected = E.confidence_selection(txs, confidence, 0.1)
assert {tx["gene_id"] for tx in selected} == {"gene000001", "gene000003", "gene000004"}
unfiltered = E.confidence_selection(txs, confidence, 0.0)
assert len(unfiltered) == 4
scores = {}
for name, predictions, count in (("confidence0.1_or_NA", selected, 3), ("unfiltered", unfiltered, 4)):
    row = E.C.score(predictions, txs, txs, length, annotation)
    assert (row["exact_chain"]["tp"], row["exact_chain"]["predicted"], row["exact_chain"]["reference"]) == (count, count, 4)
    assert row["duplicate_chain_records"] == 0
    scores[name] = row["exact_chain"]
summary = {"completed": True, "artificial_only": True, "common_parser_records": 4,
           "confidence_below_0.1_removed": True, "NA_edge_records_retained": 2,
           "strand_sensitive_coordinate_key_passed": True, "scores": scores,
           "real_reference_read": False, "DEV_used": False, "sealed_data_used": False}
(folder / "common_ruler_checks.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary), flush=True)
