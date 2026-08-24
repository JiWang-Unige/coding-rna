import importlib.util
import sys
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts/eval_structure_diagnostic.py"
sys.path.insert(0, str(SCRIPT.parent))
SPEC = importlib.util.spec_from_file_location("eval_structure_diagnostic", SCRIPT)
M24 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M24)


def write_fixture(path, rows):
    path.write_text("##gff-version 3\n" + "\n".join(rows) + "\n")
    return path


def test_primary_transcript_and_exact_candidate_semantics(tmp_path):
    reference = write_fixture(tmp_path / "reference.gff3", [
        "chr1\tRefSeq\tgene\t1\t100\t.\t+\t.\tID=g1;gene_biotype=protein_coding",
        "chr1\tRefSeq\tmRNA\t1\t100\t.\t+\t.\tID=t2;Parent=g1",
        "chr1\tRefSeq\tCDS\t1\t10\t.\t+\t0\tParent=t2",
        "chr1\tRefSeq\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1",
        "chr1\tRefSeq\tCDS\t1\t20\t.\t+\t0\tParent=t1",
        "chr1\tRefSeq\tCDS\t41\t60\t.\t+\t1\tParent=t1",
    ])
    candidate = write_fixture(tmp_path / "candidate.gff", [
        'chr1\tpred\tCDS\t1\t20\t.\t+\t0\tgene_id "p1"; transcript_id "p1";',
        'chr1\tpred\tCDS\t41\t60\t.\t+\t0\tgene_id "p1"; transcript_id "p1";',
    ])
    ref = M24.parse_annotation(reference, {"chr1": 100}, protein_coding_only=True)
    pred = M24.parse_annotation(candidate, {"chr1": 100})
    assert [tx["id"] for tx in M24.primary_transcripts(ref)] == ["t1"]
    result = M24.exact_structure_metrics(ref, pred, "candidate")
    assert result["exact_CDS_run"]["f1"] == 1.0
    assert result["pseudo_CDS_chain"]["f1"] == 1.0
    assert result["strand_aware_CDS"] == M24.NA
    assert result["phase_accuracy"] == M24.NA


def test_real_baseline_syntaxes_group_features(tmp_path):
    fixtures = {
        "tiberius.gtf": [
            'chr1\tTiberius\tgene\t1\t90\t.\t-\t.\tgene_id "g";',
            'chr1\tTiberius\ttranscript\t1\t90\t.\t-\t.\tgene_id "g"; transcript_id "t";',
            'chr1\tTiberius\tCDS\t10\t30\t.\t-\t2\tgene_id "g"; transcript_id "t";',
        ],
        "annevo.gff": [
            "chr1\tANNEVO\tgene\t1\t90\t.\t+\t.\tID=g",
            "chr1\tANNEVO\tmRNA\t1\t90\t.\t+\t.\tID=t;Parent=g",
            "chr1\tANNEVO\texon\t10\t30\t.\t+\t.\tParent=t;transcript_id=display_t",
            "chr1\tANNEVO\tCDS\t10\t30\t.\t+\t0\tParent=t;transcript_id=display_t",
        ],
        "helixer.gff3": [
            "chr1\tHelixer\tgene\t1\t90\t.\t+\t.\tID=g",
            "chr1\tHelixer\tmRNA\t1\t90\t.\t+\t.\tID=t;Parent=g",
            "chr1\tHelixer\tfive_prime_UTR\t1\t9\t.\t+\t.\tParent=t",
            "chr1\tHelixer\texon\t1\t30\t.\t+\t.\tParent=t",
            "chr1\tHelixer\tCDS\t10\t30\t.\t+\t0\tParent=t",
        ],
    }
    for name, rows in fixtures.items():
        parsed = M24.parse_annotation(write_fixture(tmp_path / name, rows), {"chr1": 100})
        selected = M24.primary_transcripts(parsed)
        assert len(selected) == 1
        assert selected[0]["CDS"] == [(9, 30, "2" if name == "tiberius.gtf" else "0")]

    tiberius = M24.parse_annotation(tmp_path / "tiberius.gtf", {"chr1": 100})
    baseline = M24.exact_structure_metrics(tiberius, tiberius, "baseline_coding")
    assert baseline["strand_accuracy"] == 1.0
    assert baseline["phase_accuracy"] == 1.0
    assert baseline["exact_intron"]["f1"] == 0.0


def test_out_of_bounds_is_not_silently_skipped(tmp_path):
    path = write_fixture(tmp_path / "bad.gff", [
        'chr1\tpred\tCDS\t1\t101\t.\t+\t0\tgene_id "g"; transcript_id "t";'
    ])
    try:
        M24.parse_annotation(path, {"chr1": 100})
    except ValueError as error:
        assert "out-of-bounds" in str(error)
    else:
        raise AssertionError("out-of-bounds input did not fail")
