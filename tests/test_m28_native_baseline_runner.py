import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('native_runner',ROOT/'scripts/experiments/M28-METHOD-RESTART/run_native_baseline.py')
N=importlib.util.module_from_spec(spec); spec.loader.exec_module(N)

class NativeRunnerTests(unittest.TestCase):
    def test_exact_six_task_scope_and_invalid_slots(self):
        self.assertEqual(len({(N.task_scope(i)[0],N.task_scope(i)[1][0]) for i in range(6)}),6)
        for i in (-1,6):
            with self.assertRaises(ValueError): N.task_scope(i)

    def test_native_commands_and_evidence_boundary(self):
        for method in N.METHODS:
            cmd,_=N.commands(method,'arabidopsis_thaliana',Path('/run'), 'GPU-assigned')
            text=' '.join(cmd)
            self.assertNotIn('reference',text); self.assertNotIn('clamsa',text)
            self.assertNotIn('overlap_pred',text); self.assertNotIn('no-overlap',text)
            if method=='ANNEVO':
                self.assertIn('annotation.py',text)
                self.assertEqual(cmd[cmd.index('--batch_size')+1],'32')
            elif method=='Tiberius':
                self.assertIn('/m28-runtime/site/model_weights:ro',text)
                self.assertIn('PYTHONPATH=/m28-runtime/site',cmd)
                self.assertIn('tiberius.main',cmd)
                self.assertNotIn('/opt/Tiberius/tiberius.py',cmd)
                self.assertEqual(cmd[cmd.index('--batch_size')+1],'4')
                self.assertEqual(cmd[cmd.index('--seq_len')+1],'200070')
                self.assertIn('/m28-prep/assets/angiosperms.yaml',cmd)
            else:
                self.assertIn('--deterministic',cmd); self.assertIn('--no-multiprocess',cmd)
                self.assertEqual(cmd[cmd.index('--batch-size')+1],'16')
                self.assertIn('land_plant_v0.3_a_0080.h5',text)

    def test_container_preserves_only_assigned_gpu(self):
        for method in ('Tiberius','Helixer'):
            cmd,probe=N.commands(method,'oryza_sativa',Path('/run'),'GPU-assigned')
            for c in (cmd,probe):
                self.assertIn('--cleanenv',c)
                self.assertIn('CUDA_VISIBLE_DEVICES=GPU-assigned',c)
                self.assertIn('/m28-run', ' '.join(c))
            if method=='Tiberius':
                self.assertIn('PYTHONPATH=/m28-runtime/site',probe)
                self.assertIn("m.version('tiberius')=='2.0.8'",probe[-1])

    def fixture(self,tmp,gff,log):
        path=Path(tmp)/'prediction.gff3'; path.write_text(gff)
        logfile=Path(tmp)/'command.log'; logfile.write_text(log)
        return path,logfile

    def test_valid_gff_and_genuine_zero_prediction(self):
        with tempfile.TemporaryDirectory() as tmp:
            path,log=self.fixture(tmp,'chr\tt\tCDS\t1\t8\t.\t-\t2\tParent=tx\n',
                                  N.FINISHED['Helixer']+' x\n')
            self.assertEqual(N.check_output(path,log,'Helixer','chr',8)['feature_counts'],{'CDS':1})
            path.write_text('# no predictions\n')
            self.assertTrue(N.check_output(path,log,'Helixer','chr',8)['empty_prediction_is_not_filtered_or_rescued'])

    def test_wrong_coordinate_scope_and_phase_rejected(self):
        for row in ('other\tt\tCDS\t1\t8\t.\t+\t0\tParent=tx\n',
                    'chr\tt\tCDS\t0\t8\t.\t+\t0\tParent=tx\n',
                    'chr\tt\tCDS\t1\t9\t.\t+\t0\tParent=tx\n',
                    'chr\tt\tCDS\t1\t8\t.\t+\t.\tParent=tx\n',
                    'chr\tt\tCDS\t1\t8\t.\t+\t0\tID=orphan\n'):
            with tempfile.TemporaryDirectory() as tmp:
                path,log=self.fixture(tmp,row,N.FINISHED['ANNEVO']+'1 seconds\n')
                with self.assertRaises(ValueError): N.check_output(path,log,'ANNEVO','chr',8)

    def test_helixer_post_failure_not_hidden_by_output_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path,log=self.fixture(tmp,'# existing partial file\n',
                                  'An error occurred during post processing. Exiting.\n')
            with self.assertRaises(RuntimeError): N.check_output(path,log,'Helixer','chr',8)

if __name__=='__main__': unittest.main()
