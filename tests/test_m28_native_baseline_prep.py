"""Synthetic scope/case and official model selection checks; no real data/download."""
import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('native_baseline_prep',
    ROOT/'scripts/experiments/M28-METHOD-RESTART/prepare_native_baselines.py')
P = importlib.util.module_from_spec(spec); spec.loader.exec_module(P)


class NativeBaselinePrepTests(unittest.TestCase):
    def test_one_target_preserves_sequence_and_uses_m28_case(self):
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp); source = parent/'in.fa'; out = parent/'out.fa'
            source.write_text('>excluded\nTTTTTT\n>allowed details\naCgN\nrysw\n>other\nGGGG\n')
            old = P.OUT; P.OUT = parent
            try:
                r = P.extract_fasta(source,'allowed',8,out)
            finally:
                P.OUT = old
            self.assertEqual(out.read_text(),'>allowed\nACGNRYSW\n')
            self.assertEqual(r['length'],8)
            self.assertEqual(r['lowercase_source_bases'],6)
            self.assertNotIn('excluded',out.read_text())
            self.assertNotIn('other',out.read_text())

    def test_missing_allowed_target_creates_no_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp); source = parent/'in.fa'; out = parent/'out.fa'
            source.write_text('>excluded\nATGC\n')
            with self.assertRaises(ValueError):
                P.extract_fasta(source,'allowed',4,out)
            self.assertFalse(out.exists())

    def test_length_mismatch_creates_no_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            parent = Path(tmp); source = parent/'in.fa'; out = parent/'out.fa'
            source.write_text('>allowed\nATGC\n')
            with self.assertRaises(ValueError):
                P.extract_fasta(source,'allowed',5,out)
            self.assertFalse(out.exists())

    def test_official_priority_not_filename_or_other_lineage(self):
        text = ('# comment\nlineage,priority,model_file_name\n'
                'land_plant,1100,old.h5\nfungi,1,fungi.h5\n'
                'land_plant,80,current.h5\nland_plant,90,alternative.h5\n')
        self.assertEqual(P.best_plant_model(text)['model_file_name'],'current.h5')


if __name__ == '__main__':
    unittest.main()
