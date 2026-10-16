from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import evaluate_native_baselines as E

def gff(seqid,starts):
    lines=['##gff-version 3']
    for i,start in enumerate(starts):
        attrs='ID=g'+str(i)+';gene_biotype=protein_coding'
        lines.extend([f'{seqid}\tt\tgene\t{start}\t{start+8}\t.\t+\t.\t{attrs}',
                      f'{seqid}\tt\tmRNA\t{start}\t{start+8}\t.\t+\t.\tID=t{i};Parent=g{i}',
                      f'{seqid}\tt\tCDS\t{start}\t{start+8}\t.\t+\t0\tParent=t{i}'])
    return '\n'.join(lines)+'\n'

class NativeEvaluationTests(unittest.TestCase):
    def fixture(self,root):
        scope=(('arabidopsis_thaliana','chrA',100),('oryza_sativa','chrR',100))
        for species,seqid,length in scope:
            data=root/'data/m1_screen'/species; data.mkdir(parents=True)
            starts=[1,31] if species=='arabidopsis_thaliana' else [1]
            (data/'reference.gff3').write_text(gff(seqid,starts))
            inputs=root/'outputs/M28-NATIVE-BASELINE-PREP-R2/inputs'; inputs.mkdir(parents=True,exist_ok=True)
            sequence='ATGAAATAA'+'A'*21+'ATGAAATAA'+'A'*61
            (inputs/(species+'.fa')).write_text('>'+seqid+'\n'+sequence+'\n')
            for method in E.METHODS:
                folder=root/'outputs/M28-NATIVE-BASELINES'/method/species; folder.mkdir(parents=True)
                predstarts=[] if method=='Helixer' else [11] if method=='ANNEVO' and species=='oryza_sativa' else starts
                (folder/'prediction.gff3').write_text(gff(seqid,predstarts))
                (folder/'STATUS').write_text('COMPLETED\n')
                row={'completed':True,'method':method,'species':species,'seqid':seqid,'length':length,
                     'input_mode':'fixed_DEV_raw_FASTA','optimizer_steps':0,'reference_used':False,
                     'labels_or_feature_cache_used':False,'test_setaria_used':False,'raw_GFF':'prediction.gff3',
                     'output_check':{'feature_counts':{'CDS':len(predstarts)}},'seconds':1,
                     'baseline_asset_receipt':{'fixture_only':True}}
                (folder/'summary.json').write_text(json.dumps(row))
        out=root/'outputs/M28-NATIVE-BASELINE-EVAL'; out.mkdir()
        return scope,out

    def patched_main(self,root,scope,out,expected):
        with patch.object(E,'ROOT',root),patch.object(E,'PRED',root/'outputs/M28-NATIVE-BASELINES'),\
             patch.object(E,'ASSET_PREP',root/'outputs/M28-NATIVE-BASELINE-PREP-R2'),\
             patch.object(E,'OUT',out),patch.object(E,'SCOPE',scope),patch.object(E,'EXPECTED',expected),\
             redirect_stdout(io.StringIO()):
            E.main()

    def test_real_common_parser_coordinate_ruler_micro_macro_and_background(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); scope,out=self.fixture(root)
            self.patched_main(root,scope,out,{'arabidopsis_thaliana':2,'oryza_sativa':1})
            result=json.loads((out/'summary.json').read_text())
            self.assertEqual(result['reference_primary_chains'],3)
            self.assertEqual(set(result['pooled']),{'ANNEVO','Tiberius','Helixer'})
            a=result['pooled']['ANNEVO']
            self.assertEqual((a['exact_chain']['tp'],a['exact_chain']['predicted'],a['exact_chain']['reference']),(2,3,3))
            self.assertAlmostEqual(a['exact_chain']['f1'],2/3)
            self.assertEqual(a['macro_exact_chain_F1'],0.5)
            self.assertEqual(a['background']['predicted_span_bp'],9)
            self.assertEqual(result['pooled']['Tiberius']['exact_chain']['f1'],1)
            self.assertEqual(result['pooled']['Helixer']['exact_chain']['predicted'],0)
            self.assertEqual(result['pooled']['Helixer']['exact_chain']['f1'],0)
            self.assertFalse(result['test_setaria_used'])

    def test_incomplete_native_run_prevents_any_reference_open(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); scope,out=self.fixture(root)
            (root/'outputs/M28-NATIVE-BASELINES/Helixer/oryza_sativa/STATUS').write_text('RUNNING\n')
            with patch.object(E.C.M,'read_allowed') as reference_read:
                with self.assertRaisesRegex(ValueError,'Incomplete native'):
                    self.patched_main(root,scope,out,{'arabidopsis_thaliana':2,'oryza_sativa':1})
                reference_read.assert_not_called()
            self.assertFalse((out/'summary.json').exists())

    def test_wrong_reference_denominator_not_silently_dropped(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); scope,out=self.fixture(root)
            with self.assertRaisesRegex(ValueError,'reference denominator'):
                self.patched_main(root,scope,out,{'arabidopsis_thaliana':3,'oryza_sativa':1})
            self.assertFalse((out/'summary.json').exists())

if __name__=='__main__': unittest.main()
