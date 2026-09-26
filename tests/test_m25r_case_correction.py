import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "case_correction", ROOT / "scripts/experiments/M25R-DEV-REDECODE-ERROR-DECOMPOSITION/posthoc_r4_case_correction.py")
CORR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CORR)
D = CORR.diag
TOY_SPEC = importlib.util.spec_from_file_location("original_cases", ROOT / "tests/test_m25r_redecode_error_decomposition.py")
TOY = importlib.util.module_from_spec(TOY_SPEC)
TOY_SPEC.loader.exec_module(TOY)


def test_diagnostic_entry_normalizes_without_changing_reader(monkeypatch):
    raw = {"chr": "aTgtaa"}
    monkeypatch.setattr(D.screen_data, "read_fasta", lambda _path: dict(raw))
    monkeypatch.setattr(D.screen_data, "assign_splits", lambda _ids: {"chr": "val"})
    monkeypatch.setattr(D, "parse_annotation", lambda *_a, **_kw: {})
    monkeypatch.setattr(D.m25, "_primary_with_partial", lambda _a: [])
    config = {"data": {"development_species": ["arabidopsis_thaliana"],
                       "primary_chromosome_seqids": {"arabidopsis_thaliana": ["chr"]}}}
    result = D.load_species(Path("/unused"), config)
    assert result["arabidopsis_thaliana"]["seqs"]["chr"] == "ATGTAA"
    assert D.screen_data.read_fasta("unused")["chr"] == "aTgtaa"


@pytest.mark.parametrize("reverse", [False, True])
def test_case_validity_decoder_and_reconstruction(tmp_path, reverse):
    sequence, region, boundary, phase = TOY.coding_example()
    thresholds = dict.fromkeys(("region", "start", "stop", "donor", "acceptor"), 0.5)
    states, traces, prefilter, emitted = D.trace_decode_orientation(sequence, region, boundary, phase, thresholds)
    lower = D.trace_decode_orientation(sequence.lower(), region, boundary, phase, thresholds)
    np.testing.assert_array_equal(states, lower[0])
    assert (traces, prefilter, emitted) == lower[1:]
    np.testing.assert_array_equal(states, CORR.reconstruct_states(len(sequence), traces))
    # trace_decode_orientation itself checks exact equality with the production decoder.
    model = D.m25._map_model_from_orientation(emitted[0], len(sequence), "-" if reverse else "+")
    genome = D.m25.reverse_complement(sequence) if reverse else sequence
    gff = tmp_path / "prediction.gff3"
    D.m25._write_gff3(str(gff), {"chr": [model]}, "test")
    raw_audit = D.audit_gff3(gff, {"chr": genome.lower()}, {"chr": len(genome)}, 1)
    corrected = D.audit_gff3(gff, {"chr": genome.upper()}, {"chr": len(genome)}, 1)
    assert raw_audit["invalid_transcripts"] == 1
    assert corrected["valid_transcripts"] == 1
    ref = {"events": D.truth_events([(5, 12), (20, 28)])}
    assert CORR.indexed_reachability(ref, CORR.reachability_index(states, sequence)) == (True, True)
    assert CORR.indexed_reachability(ref, CORR.reachability_index(states, sequence.lower())) == (True, False)
    assert D.canonical_truth(sequence, [(5, 12), (20, 28)])
    assert not D.canonical_truth(sequence.lower(), [(5, 12), (20, 28)])


def test_index_matches_original_reachability_for_supported_events():
    sequence, region, boundary, phase = TOY.coding_example()
    states, _traces, _pre, _emit = D.trace_decode_orientation(
        sequence, region, boundary, phase, dict.fromkeys(("region", "start", "stop", "donor", "acceptor"), 0.5))
    for seq in (sequence, sequence.lower()):
        index = CORR.reachability_index(states, seq)
        for donor_positions in ([12], [11], [12, 13], [0], []):
            ref = {"events": {"start": [5], "stop": [25], "donor": donor_positions, "acceptor": [18]}}
            assert CORR.indexed_reachability(ref, index) == D.transition_upper_bounds(states, seq, ref)


def test_reconstruction_rejects_overlapping_lineages():
    with pytest.raises(AssertionError, match="overlapping"):
        CORR.reconstruct_states(30, [{"block_span": [2, 12], "runs": [[3, 8]]},
                                     {"block_span": [10, 20], "runs": []}])

