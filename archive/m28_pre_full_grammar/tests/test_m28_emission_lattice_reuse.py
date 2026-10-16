"""Actual saved synthetic topology/answer reconciliation after A42 reuse edits."""
import json
from pathlib import Path
import resource
import time
import unittest

from profile_emission_lattice_synthetic import known_fixture, random_fixture, verify_selected
from src.m28.emission_lattice import decode_emission_lattice


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'outputs/M28-A42-EMISSION-LATTICE-LONG-CHECKS'


class EmissionLatticeReuseTests(unittest.TestCase):
    def test_saved_full_window_chains_gains_and_every_traversal_count_unchanged(self):
        for case in (0, 1):
            folder = BASE/('case_'+str(case))
            expected_job = '13421438_'+str(case)
            if ((folder/'STATUS').read_text().strip() != 'COMPLETED' or
                    (folder/'JOBID').read_text().strip() != expected_job):
                raise ValueError('Actual completed original synthetic array task required')
            saved = json.loads((folder/'summary.json').read_text())
            self.assertTrue(saved['completed'])
            self.assertTrue(saved['synthetic_only'])
            self.assertEqual(saved['job_id'], expected_job)
            self.assertEqual(saved['width'], 24576)
            if case == 0:
                sequence, logits, _, _ = known_fixture()
            else:
                sequence, logits = random_fixture()
            self.assertEqual([r['exon_cap'] for r in saved['runs']], [128, 127] if case == 0 else [128])
            for row in saved['runs']:
                began = time.perf_counter()
                result = decode_emission_lattice(sequence, logits, row['exon_cap'])
                seconds = time.perf_counter() - began
                verify_selected(sequence, logits, result, row['exon_cap'])
                self.assertEqual(result['chains'], tuple(tuple(map(tuple, c)) for c in row['selected_chains']))
                self.assertEqual(result['chain_gains'], tuple(row['chain_gains']))
                self.assertEqual(result['gain'], row['gain'])
                self.assertEqual(result['counts'], row['counts'])
                print(json.dumps(dict(source_job=expected_job,case=case,exon_cap=row['exon_cap'],
                    original_seconds=row['seconds'],reuse_seconds=seconds,
                    peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                    exact_chains_gains_and_all_counts=True),allow_nan=False),flush=True)


if __name__ == '__main__':
    unittest.main()
