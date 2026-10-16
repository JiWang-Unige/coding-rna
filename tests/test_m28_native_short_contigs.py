"""Two artificial-contig native interface failures that change readiness."""
import sys
from pathlib import Path
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import preflight_native_short_contigs as A


class NativeShortContigTests(unittest.TestCase):
    def test_only_fasta_is_replaced_and_all_native_operation_points_remain(self):
        run = Path('/artificial-run')
        self.assertEqual(list(A.LENGTHS.values()), [332, 7536])
        for method in A.N.METHODS:
            original, probe = A.N.commands(method, 'artificial_short_preflight', run, '3', cpus=8)
            actual, actual_probe = A.commands(method, run, '3')
            option = {'ANNEVO': '-g', 'Tiberius': '--genome', 'Helixer': '--fasta-path'}[method]
            expected = original[:]
            expected[expected.index(option)+1] = '/artificial-run/artificial.fa' if method == 'ANNEVO' else '/m28-run/artificial.fa'
            self.assertEqual(actual, expected); self.assertEqual(actual_probe, probe)
            self.assertEqual(sum(x != y for x, y in zip(actual, original)), 1)

    def test_empty_output_is_valid_only_after_completion_and_actual_multicontig_CDS_is_checked(self):
        with tempfile.TemporaryDirectory() as root:
            path, log = Path(root)/'prediction.gff3', Path(root)/'command.log'
            log.write_text(A.N.FINISHED['Helixer']); path.write_text('##gff-version 3\n')
            empty = A.check_output(path, log, 'Helixer', A.LENGTHS)
            self.assertTrue(empty['empty_prediction_is_not_filtered_or_rescued'])
            self.assertFalse(empty['nonempty_CDS_structure_observed'])
            path.write_text('artificial_332\tnative\tCDS\t1\t9\t.\t+\t0\tParent=t1\n'
                            'artificial_7536\tnative\tCDS\t6141\t7535\t.\t-\t2\tParent=t2\n')
            self.assertTrue(A.check_output(path, log, 'Helixer', A.LENGTHS)['nonempty_CDS_structure_observed'])
            for row in ('artificial_332\tx\tCDS\t1\t333\t.\t+\t0\tParent=t1\n',
                        'other_contig\tx\tCDS\t1\t9\t.\t+\t0\tParent=t1\n',
                        'artificial_332\tx\tCDS\t1\t9\t.\t+\t.\tParent=t1\n'):
                path.write_text(row)
                with self.assertRaises(ValueError):
                    A.check_output(path, log, 'Helixer', A.LENGTHS)
            log.write_text('postprocessor failed but main returned zero')
            with self.assertRaises(RuntimeError):
                A.check_output(path, log, 'Helixer', A.LENGTHS)


if __name__ == '__main__':
    unittest.main()
