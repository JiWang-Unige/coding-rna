from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import run_native_train_baseline as N
import evaluate_native_train_baselines as E


def gff(seqid, starts):
    rows = ['##gff-version 3']
    for i,start in enumerate(starts):
        rows += [f'{seqid}\tt\tgene\t{start}\t{start+8}\t.\t+\t.\tID=g{i};gene_biotype=protein_coding',
                 f'{seqid}\tt\tmRNA\t{start}\t{start+8}\t.\t+\t.\tID=t{i};Parent=g{i}',
                 f'{seqid}\tt\tCDS\t{start}\t{start+8}\t.\t+\t0\tParent=t{i}']
    return '\n'.join(rows)+'\n'


class NativeTrainTests(unittest.TestCase):
    def test_eighteen_unique_original_TRAIN_pairs(self):
        pairs = [N.task_scope(i) for i in range(18)]
        self.assertEqual(len({(m,s[1]) for m,s in pairs}),18)
        self.assertEqual(N.SCOPE,E.E.I.SCOPE)
        for i in (-1,18):
            with self.assertRaises(ValueError): N.task_scope(i)

    def test_shared_FASTA_changes_only_input_not_native_operating_point(self):
        fasta = Path('/train/chr/input.fasta')
        for method in N.METHODS:
            old,old_probe = N.N.commands(method,'arabidopsis_thaliana',Path('/run'),'GPU-assigned')
            cmd,probe = N.N.commands(method,'arabidopsis_thaliana',Path('/run'),'GPU-assigned',fasta=fasta)
            input_arg = {'ANNEVO':'-g','Tiberius':'--genome','Helixer':'--fasta-path'}[method]
            self.assertEqual(cmd[cmd.index(input_arg)+1],str(fasta) if method=='ANNEVO' else '/m28-input.fa')
            cmd[cmd.index(input_arg)+1] = old[old.index(input_arg)+1]
            if method != 'ANNEVO':
                for command in (cmd,probe):
                    index = command.index(str(fasta)+':/m28-input.fa:ro')
                    self.assertEqual(command[index-1],'--bind')
                    del command[index-1:index+1]
            self.assertEqual(cmd,old)
            self.assertEqual(probe,old_probe)

    def test_shared_FASTA_scope_length_and_encoding(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fasta = root/'input.fasta'
            with patch.object(N,'ROOT',root):
                fasta.write_text('>chr\nACGTNN\n')
                self.assertEqual(N.check_fasta(fasta,'chr',6)['N_bases'],2)
                for text,length in [('>other\nACGTNN\n',6),('>chr\nACGTNN\n',7),
                                    ('>chr\nACGTRN\n',6),('>chr\nacgtNN\n',6),
                                    ('>chr\nACGTNN\n>other\nA\n',7)]:
                    fasta.write_text(text)
                    with self.assertRaises(ValueError): N.check_fasta(fasta,'chr',length)

    def fixture(self,root):
        scope = (('arabidopsis_thaliana','chrA',100),('arabidopsis_thaliana','chrA2',100),
                 ('oryza_sativa','chrR',100))
        pred = root/'pred'; original = root/'original'; control = root/'control'; out = root/'evaluation'
        control.mkdir(); out.mkdir()
        assets = {'baselines':{m:{'fixture_only':True} for m in N.METHODS}}
        control_rows = {}
        for species in ('arabidopsis_thaliana','oryza_sativa'):
            data = root/'data/m1_screen'/species; data.mkdir(parents=True)
            (data/'reference.gff3').write_text(''.join(gff(seqid,[1,31] if seqid=='chrA' else [1])
                for s,seqid,_ in scope if s==species))
        for index,(species,seqid,length) in enumerate(scope):
            folder = original/seqid; folder.mkdir(parents=True)
            sequence = 'ATGAAATAA'+'A'*21+'ATGAAATAA'+'A'*61
            (folder/'input.fasta').write_text('>'+seqid+'\n'+sequence+'\n')
            annotation = E.C.S.E.parse_annotation(root/'data/m1_screen'/species/'reference.gff3',{seqid:length})
            refs = E.C.S.E.primary_transcripts(annotation)
            row = E.C.score(refs,refs,refs,length,annotation)
            row.update(species=species,seqid=seqid,length=length,reference_primary_records=len(refs),
                       reference_assessable_isoform_records=len(refs),retained_unsupported_reference_strata={})
            control_rows[seqid] = row
            for method_index,method in enumerate(N.METHODS):
                run = pred/method/seqid; run.mkdir(parents=True)
                starts = [1,31] if seqid=='chrA' else [1]
                if method=='Helixer' or (method=='ANNEVO' and seqid=='chrA2'): starts=[]
                elif method=='ANNEVO' and seqid=='chrR': starts=[11]
                (run/'prediction.gff3').write_text(gff(seqid,starts))
                (run/'STATUS').write_text('COMPLETED\n')
                summary = dict(completed=True,task=method_index*len(scope)+index,method=method,species=species,
                    seqid=seqid,length=length,input_mode=N.INPUT_MODE,input_fasta=str(folder/'input.fasta'),
                    optimizer_steps=0,native_inference_executed=True,reference_used=False,
                    labels_or_feature_cache_used=False,DEV_used=False,sealed_data_used=False,test_setaria_used=False,
                    raw_GFF='prediction.gff3',output_check={'feature_counts':{'CDS':len(starts)}},
                    seconds=1,baseline_asset_receipt=assets['baselines'][method])
                (run/'summary.json').write_text(json.dumps(summary))
        per_species = {s:E.E.aggregate([r for r in control_rows.values() if r['species']==s])
                       for s in ('arabidopsis_thaliana','oryza_sativa')}
        receipt = dict(completed=True,status='original_TRAIN_end_to_end_diagnostic_not_DEV_or_generalization',
            reference_policy='M26_CDS_assessable_longest_eligible_CDS_tie_transcript_ID',chromosomes=control_rows,
            per_species=per_species,pooled=E.E.aggregate(list(control_rows.values())),
            macro_species_exact_chain_F1=1,worst_species_exact_chain_F1=1)
        (control/'STATUS').write_text('COMPLETED\n'); (control/'summary.json').write_text(json.dumps(receipt))
        return scope,pred,original,control,out,assets

    def execute(self,root,fixture,expected=None):
        scope,pred,original,control,out,assets = fixture
        with patch.object(N,'SCOPE',scope),patch.object(N,'OUT',pred),patch.object(N,'ORIGINAL_TRAIN',original),\
             patch.object(N,'runtime_receipt',return_value=assets),patch.object(E,'ROOT',root),\
             patch.object(E,'CONTROL',control),patch.object(E,'OUT',out),\
             patch.object(E,'EXPECTED',expected or {'chrA':2,'chrA2':1,'chrR':1}),\
             patch.dict('os.environ',{'SLURM_JOB_ID':'synthetic-check'}),redirect_stdout(io.StringIO()):
            E.main()

    def test_real_parser_six_scope_logic_species_macro_and_zero_predictions(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fixture = self.fixture(root)
            self.execute(root,fixture)
            result = json.loads((fixture[4]/'summary.json').read_text())
            a = result['methods']['ANNEVO']
            self.assertEqual(result['reference_primary_chains'],4)
            self.assertAlmostEqual(a['pooled']['exact_chain']['f1'],4/7)
            self.assertAlmostEqual(a['macro_species_exact_chain_F1'],0.4)
            self.assertEqual(a['worst_species_exact_chain_F1'],0)
            self.assertEqual(a['pooled']['background']['predicted_span_bp'],9)
            self.assertAlmostEqual(a['native_minus_original4608_F1_pp']['pooled'],100*(4/7-1))
            self.assertEqual(result['methods']['Tiberius']['pooled']['exact_chain']['f1'],1)
            self.assertEqual(result['methods']['Helixer']['pooled']['exact_chain']['predicted'],0)

    def test_incomplete_native_blocks_reference_access(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fixture = self.fixture(root)
            (fixture[1]/'Helixer/chrR/STATUS').write_text('RUNNING\n')
            with patch.object(E.C.M,'read_allowed') as read_reference:
                with self.assertRaisesRegex(ValueError,'eighteen'):
                    self.execute(root,fixture)
                read_reference.assert_not_called()

    def test_wrong_reference_denominator_is_not_filtered(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fixture = self.fixture(root)
            with self.assertRaisesRegex(ValueError,'reference denominator'):
                self.execute(root,fixture,{'chrA':3,'chrA2':1,'chrR':1})
            self.assertFalse((fixture[4]/'summary.json').exists())


if __name__ == '__main__': unittest.main()
