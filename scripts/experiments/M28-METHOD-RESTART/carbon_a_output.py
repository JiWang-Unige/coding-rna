"""Independent GFF/CDS/protein roundtrip used by Carbon-A artificial and TRAIN runs."""
from collections import defaultdict
from Bio import SeqIO
from Bio.Data import CodonTable

_STANDARD = CodonTable.unambiguous_dna_by_id[1]
_AMINO_ACIDS = {**_STANDARD.forward_table, **{codon: "*" for codon in _STANDARD.stop_codons}}


def literal_standard_translation(sequence):
    """Match the fixed author's table: every non-ACGT codon is X, even GCN."""
    return "".join(_AMINO_ACIDS.get(sequence[i:i + 3], "X")
                   for i in range(0, len(sequence) // 3 * 3, 3))


def rc(sequence):
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def read_outputs(prefix, sequence, record_id, expected=None):
    """Independently reconstruct written CDS from 1-based closed GFF intervals."""
    genes, transcripts, segments = {}, {}, defaultdict(list)
    for line in prefix.with_suffix(".gff").read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        assert len(fields) == 9, line
        record, _, kind, start, end, _, strand, phase, attributes = fields
        assert record == record_id and strand in ("+", "-")
        start, end = int(start) - 1, int(end)
        assert 0 <= start < end <= len(sequence), line
        attrs = dict(item.split("=", 1) for item in attributes.split(";"))
        if kind == "gene":
            assert attrs["ID"] not in genes
            confidence = attrs["confidence"]
            assert confidence == "NA" or 0 <= float(confidence) <= 1
            genes[attrs["ID"]] = (strand, start, end, attrs)
        elif kind == "mRNA":
            transcripts[attrs["ID"]] = attrs["Parent"]
        elif kind == "CDS":
            assert int(phase) in (0, 1, 2)
            segments[attrs["Parent"]].append((start, end, int(phase), strand))
        else:
            raise AssertionError(kind)
    cds = {record.id: str(record.seq) for record in SeqIO.parse(prefix.with_suffix(".fna"), "fasta")}
    proteins = {record.id: str(record.seq) for record in SeqIO.parse(prefix.with_suffix(".faa"), "fasta")}
    assert set(cds) == set(genes), (set(cds), set(genes))
    assert set(transcripts.values()) == set(genes)
    assert set(segments) == set(transcripts)
    for transcript_id, gene_id in transcripts.items():
        strand, gene_start, gene_end, _ = genes[gene_id]
        ordered = sorted(segments[transcript_id], reverse=strand == "-")
        assert all(segment[3] == strand for segment in ordered)
        assert min(segment[0] for segment in ordered) == gene_start
        assert max(segment[1] for segment in ordered) == gene_end
        reconstructed = "".join(sequence[start:end] if strand == "+" else rc(sequence[start:end])
                                for start, end, _, _ in ordered)
        initial_phase = ordered[0][2]
        consumed = 0
        for start, end, phase, _ in ordered:
            assert phase == (initial_phase - consumed) % 3
            consumed += end - start
        assert reconstructed == cds[gene_id], gene_id
        coding = reconstructed[initial_phase:]
        protein = literal_standard_translation(coding)
        if protein.endswith("*"):
            protein = protein[:-1]
        assert "*" not in protein, gene_id
        assert proteins.get(gene_id, "") == protein, gene_id
        if expected is not None:
            wanted = expected[gene_id]
            assert strand == wanted["strand"]
            assert ordered == wanted["segments"], (ordered, wanted)
            assert reconstructed == wanted["cds"]
    if expected is not None:
        assert set(genes) == set(expected)
    return {"genes": len(genes), "proteins": len(proteins), "coordinate_cds_protein_roundtrip": True}
