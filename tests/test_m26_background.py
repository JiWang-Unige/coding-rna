import importlib.util
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts/experiments/M26-SAME-SCOPE-MECHANISM'))
import background_attribution as B


def test_region_partition_preserves_partial_and_noncoding_distinction():
    genes = {'eligible': {'start': 10, 'end': 30, 'partial': False},
             'partial': {'start': 25, 'end': 50, 'partial': True},
             'other_coding': {'start': 60, 'end': 70, 'partial': False},
             'nc': {'start': 65, 'end': 80, 'partial': False}}
    primary = [{'gene_id': 'eligible', 'CDS': [(15,20,'0')], 'exon': [(12,23,'.')]}]
    labels = B.regions(100, {'genes': genes}, {'genes': {k:v for k,v in genes.items() if k!='nc'}}, primary)
    assert np.bincount(labels, minlength=6).tolist() == [11, 9, 20, 10, 10, 40]
    stats = B.burden([{'CDS': [(20,45,'0')]}, {'CDS': [(40,65,'0')]}], labels)
    assert sum(stats['predicted_union_span_bp_by_region'].values()) == 45
    assert stats['predicted_union_span_bp_by_region']['partial_flagged_coding_gene'] == 20
    assert stats['frozen_background_FPR'] == 42/89
    assert stats['all_coding_gene_background_FPR'] == 10/50
