import json
import subprocess
import sys


def write(path, text):
    path.write_text(text, encoding="utf-8")


def test_gene_count_ratio_uses_reference_gene_count_not_transcripts(tmp_path):
    fasta = tmp_path / "genome.fa"
    ref = tmp_path / "ref.gtf"
    pred = tmp_path / "pred.gtf"
    metrics = tmp_path / "metrics.json"

    write(fasta, ">chr1\n" + "A" * 1000 + "\n")
    write(
        ref,
        "\n".join(
            [
                'chr1\tref\tCDS\t101\t150\t.\t+\t0\tgene_id "g1"; transcript_id "t1";',
                'chr1\tref\tintron\t151\t199\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
                'chr1\tref\tCDS\t200\t250\t.\t+\t0\tgene_id "g1"; transcript_id "t1";',
                'chr1\tref\tCDS\t101\t150\t.\t+\t0\tgene_id "g1"; transcript_id "t2";',
                'chr1\tref\tintron\t151\t199\t.\t+\t.\tgene_id "g1"; transcript_id "t2";',
                'chr1\tref\tCDS\t200\t250\t.\t+\t0\tgene_id "g1"; transcript_id "t2";',
                "",
            ]
        ),
    )
    write(
        pred,
        "\n".join(
            [
                'chr1\tpred\tgene\t101\t250\t.\t+\t.\tgene_id "p1";',
                'chr1\tpred\ttranscript\t101\t250\t.\t+\t.\tgene_id "p1"; transcript_id "p1.t1";',
                'chr1\tpred\tCDS\t101\t150\t.\t+\t0\tgene_id "p1"; transcript_id "p1.t1";',
                'chr1\tpred\tintron\t151\t199\t.\t+\t.\tgene_id "p1"; transcript_id "p1.t1";',
                'chr1\tpred\tCDS\t200\t250\t.\t+\t0\tgene_id "p1"; transcript_id "p1.t1";',
                "",
            ]
        ),
    )

    subprocess.run(
        [
            sys.executable,
            "scripts/eval_gene_body_mask.py",
            "--reference-gtf",
            str(ref),
            "--prediction-gtf",
            str(pred),
            "--genome-fasta",
            str(fasta),
            "--output-json",
            str(metrics),
            "--experiment-id",
            "TEST",
            "--profile",
            "smoke",
        ],
        check=True,
    )
    data = json.loads(metrics.read_text(encoding="utf-8"))

    assert data["reference_gene_count"] == 1
    assert data["reference_transcript_count"] == 2
    assert data["predicted_gene_count"] == 1
    assert data["predicted_gene_count_ratio_vs_reference"] == 1.0
    assert data["predicted_transcript_count_ratio_vs_reference"] == 0.5

