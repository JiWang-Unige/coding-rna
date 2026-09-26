import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('m26', ROOT / 'scripts/experiments/M26-SAME-SCOPE-MECHANISM/same_scope.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def tx(name, strand, cds):
    return {'id': name, 'seqid': 'chr', 'strand': strand, 'CDS': [(a,b,'0') for a,b in cds], 'exon': []}


def test_strand_and_background_denominators():
    scope = ('sp', 'chr')
    refs = {scope: [tx('r', '+', [(10,20)])]}
    preds = {scope: [tx('p', '-', [(10,20)]), tx('extra', '+', [(40,50)])]}
    all_ann = {scope: {'genes': {'a': {'start': 10, 'end': 20}, 'b': {'start': 40, 'end': 50}}}}
    result = m.compare(preds, refs, {scope: 100}, all_ann)
    assert result['exact_chain']['tp'] == 0
    assert result['CDS_base']['tp'] == 10
    assert result['frozen_metrics']['intergenic_FPR'] == 10/90
    assert result['all_annotated_gene_background_FPR'] == 0
    assert 'structurally_valid_complete_fraction' not in result['frozen_metrics']


def test_nonexact_is_not_automatically_biological_false_positive():
    primary = tx('r', '+', [(10,20), (30,40)])
    iso = tx('iso', '+', [(11,20), (30,42)])
    terminal = tx('terminal', '+', [(12,20), (30,43)])
    antisense = tx('anti', '-', [(10,20)])
    separate = tx('new', '+', [(70,80)])
    rows = m.error_rows([primary, iso, terminal, antisense, separate], [primary], [primary, iso])
    assert [r['category'] for r in rows] == ['exact_primary', 'exact_other_complete_isoform',
        'same_intron_chain_terminal_disagreement', 'opposite_strand_reference_overlap_only',
        'no_complete_primary_CDS_span_overlap']


def test_zero_based_half_open_and_negative_strand_sequence():
    seqs = {'chr': 'TTATTTCAT'}
    result = m.sequence_convention([tx('minus', '-', [(0,9)])], seqs)
    assert result['terminal_stop_in_CDS'] == 1
    assert result['ATG_start'] == 1
