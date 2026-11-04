"""Unit tests of the stdlib part of vip_score.py (converter + independent chain derivation).  Run: python3 -m unittest experiments/M38-VIPSANIA/tests/test_vip_convert.py"""
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import vip_score as V          # noqa: E402

# a Vipsania-style GFF3 (bricks2marble writer): gene, mRNA(Parent=gene), CDS+exon children, frames; "+" gene with 2 exons, "-" gene with 2 exons (listed in + order),
# a short single-exon gene (60 bp, dropped by the CDS>=96 rule), and a second contig
GFF = """##gff-version 3
ctg1\tVipsania\tgene\t101\t400\t.\t+\t.\tID=vipsania_g1;Name=vipsania_g1
ctg1\tVipsania\tmRNA\t101\t400\t.\t+\t.\tID=vipsania_g1.t1;Name=vipsania_g1.t1;Parent=vipsania_g1
ctg1\tVipsania\tCDS\t101\t200\t.\t+\t0\tID=vipsania_g1.t1.CDS.1;Parent=vipsania_g1.t1
ctg1\tVipsania\texon\t101\t200\t.\t+\t0\tID=vipsania_g1.t1.exon.1;Parent=vipsania_g1.t1
ctg1\tVipsania\tCDS\t301\t400\t.\t+\t2\tID=vipsania_g1.t1.CDS.2;Parent=vipsania_g1.t1
ctg1\tVipsania\texon\t301\t400\t.\t+\t2\tID=vipsania_g1.t1.exon.2;Parent=vipsania_g1.t1
ctg1\tVipsania\tgene\t1001\t1500\t.\t-\t.\tID=vipsania_g2;Name=vipsania_g2
ctg1\tVipsania\tmRNA\t1001\t1500\t.\t-\t.\tID=vipsania_g2.t1;Name=vipsania_g2.t1;Parent=vipsania_g2
ctg1\tVipsania\tCDS\t1001\t1100\t.\t-\t1\tID=vipsania_g2.t1.CDS.1;Parent=vipsania_g2.t1
ctg1\tVipsania\tCDS\t1401\t1500\t.\t-\t0\tID=vipsania_g2.t1.CDS.2;Parent=vipsania_g2.t1
ctg1\tVipsania\tgene\t2001\t2060\t.\t+\t.\tID=vipsania_g3;Name=vipsania_g3
ctg1\tVipsania\tmRNA\t2001\t2060\t.\t+\t.\tID=vipsania_g3.t1;Parent=vipsania_g3
ctg1\tVipsania\tCDS\t2001\t2060\t.\t+\t0\tID=vipsania_g3.t1.CDS.1;Parent=vipsania_g3.t1
ctg2\tVipsania\tgene\t11\t310\t.\t+\t.\tID=vipsania_g4;Name=vipsania_g4
ctg2\tVipsania\tmRNA\t11\t310\t.\t+\t.\tID=vipsania_g4.t1;Parent=vipsania_g4
ctg2\tVipsania\tCDS\t11\t310\t.\t+\t0\tID=vipsania_g4.t1.CDS.1;Parent=vipsania_g4.t1
"""


class ConvertTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.gff = self.dir / 'in.gff3'
        self.gff.write_text(GFF)

    def tearDown(self):
        self.tmp.cleanup()

    def test_chains_and_cds96_rule(self):
        d = V.read_gff(self.gff, keep={'ctg1'})
        self.assertEqual(sorted(d), ['ctg1'])
        all_ = V.chains_of(d['ctg1'], 'ctg1', cds96=False)
        self.assertEqual(len(all_), 3)
        kept = V.chains_of(d['ctg1'], 'ctg1', cds96=True)
        self.assertEqual(sorted(kept), [('ctg1', '+', ((100, 200), (300, 400))), ('ctg1', '-', ((1000, 1100), (1400, 1500)))])   # 0-based half-open, g3 (60 bp) dropped

    def test_written_layout_roundtrip(self):
        rep = V.convert(self.gff, self.dir / 'run', ['ctg1', 'ctg2', 'ctg_empty'])
        self.assertEqual(rep['ctg1']['genes'], 3)
        self.assertEqual(rep['ctg2']['chains'], 1)
        self.assertEqual(rep['ctg_empty']['genes'], 0)
        p = self.dir / 'run' / 'ctg1' / 'prediction' / 'full_grammar.gff3'
        d2 = V.read_gff(p)
        d1 = V.read_gff(self.gff, keep={'ctg1'})
        self.assertEqual(sorted(V.chains_of(d2['ctg1'], 'ctg1', False)), sorted(V.chains_of(d1['ctg1'], 'ctg1', False)))
        txt = p.read_text().splitlines()
        self.assertEqual(txt[0], '##gff-version 3')
        self.assertTrue(any('\tgene\t101\t400\t' in x and 'gene_biotype=protein_coding' in x for x in txt))
        self.assertTrue(all(len(x.split('\t')) == 9 for x in txt if not x.startswith('#')))
        cds = [x.split('\t') for x in txt if '\tCDS\t' in x]
        self.assertEqual([c[7] for c in cds], ['0', '2', '1', '0', '0'])                 # phases are carried over

    def test_gene_filter(self):
        V.convert(self.gff, self.dir / 'run', ['ctg1'], gene_filter=lambda g: g != 'vipsania_g1')
        d = V.read_gff(self.dir / 'run' / 'ctg1' / 'prediction' / 'full_grammar.gff3')
        self.assertEqual(sorted(d['ctg1']['genes']), ['vipsania_g2', 'vipsania_g3'])

    def test_slice_shift_and_rename(self):
        sl = dict(name='ctg1', contig='chrA', start1=1001)
        d = V.read_gff(self.gff, slice_info=sl)
        self.assertEqual(sorted(d), ['chrA'])
        self.assertEqual(d['chrA']['tx']['vipsania_g1.t1']['cds'][0][:2], (1101, 1200))

    def test_children_before_parents_and_gene_parent_cds(self):
        g = self.dir / 'odd.gff3'
        g.write_text("ctg1\tx\tCDS\t11\t220\t.\t+\t0\tID=c1;Parent=t1\n"
                     "ctg1\tx\tmRNA\t11\t220\t.\t+\t.\tID=t1;Parent=g1\n"
                     "ctg1\tx\tgene\t11\t220\t.\t+\t.\tID=g1\n"
                     "ctg1\tx\tCDS\t301\t520\t.\t-\t0\tID=c2;Parent=g2\n")                   # CDS parented on a gene id, no mRNA row
        d = V.read_gff(g)
        self.assertEqual(len(V.chains_of(d['ctg1'], 'ctg1')), 2)

    def test_bad_column_count(self):
        g = self.dir / 'bad.gff3'
        g.write_text('ctg1\tx\tgene\t1\t10\t.\t+\n')
        with self.assertRaises(ValueError):
            V.read_gff(g)


if __name__ == '__main__':
    unittest.main()
