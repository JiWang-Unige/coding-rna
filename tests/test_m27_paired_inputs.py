import importlib.util
from pathlib import Path
import pytest

PATH=Path(__file__).resolve().parents[1]/"scripts/experiments/M27-ALLELIC-INTEGRITY/paired_inputs.py"
spec=importlib.util.spec_from_file_location("paired_inputs",PATH)
M=importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)


@pytest.mark.parametrize("strand,genome,position,ref,syn,stop",[
    ("+","AAATACGGG",6,"C","T","A"),
    ("-","AAAGTAGGG",4,"G","A","T"),
])
def test_single_SNV_preserves_whole_context_and_strand_codon(strand,genome,position,ref,syn,stop):
    row={"position_1based":position,"reference_base":ref,"synonymous_alt":syn,"stop_alt":stop,
         "CDS_chain_0based_halfopen":[[0,9]],"strand":strand,"codon_index_1based":2,
         "WT_codon":"TAC","synonymous_codon":"TAT","stop_codon":"TAA"}
    for allele,alt in (("syn",syn),("PTC",stop)):
        sequence=M.mutant(genome,row,allele)
        assert len(sequence)==len(genome)
        assert [i for i,(a,b) in enumerate(zip(sequence,genome)) if a!=b]==[position-1]
        assert sequence[position-1]==alt


def test_wrong_registered_reference_stops_before_input_creation():
    row={"position_1based":6,"reference_base":"T","synonymous_alt":"T","stop_alt":"A",
         "CDS_chain_0based_halfopen":[[0,9]],"strand":"+","codon_index_1based":2,
         "WT_codon":"TAC","synonymous_codon":"TAT","stop_codon":"TAA"}
    with pytest.raises(ValueError,match="reference base"):
        M.mutant("AAATACGGG",row,"PTC")
