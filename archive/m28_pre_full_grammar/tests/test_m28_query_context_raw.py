import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import torch
from test_m28_query_decode import fixture

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import infer_query_context as I
import evaluate_query_context as E
from src.m28.export import write_gff3


class SyntheticGenerator:
    def __init__(self,model):
        self.query = model; self.calls = 0

    def __call__(self,features,dna,valid):
        self.calls += 1
        if tuple(dna.shape) != (1,I.M.WIDTH,5) or len(valid) != 1:
            raise AssertionError('Original raw-window DNA contract changed')
        return dict(features=features)


def fitted_receipt():
    return dict(arm=I.ARM,configuration=I.R.CONFIG,feature_revision=I.F.REVISION,
        completed=True,optimizer_steps=4608,trainable_parameters=649457,
        final_checkpoint=str((I.R.FIT/'step_004608.pt').relative_to(I.ROOT)),
        main_evaluation='last_step_004608_only',epochs=[dict(steps=1536)]*3,
        DEV_used_for_selection=False,test_setaria_used=False)


def decoded_receipt(fitted):
    scopes = []
    for ordinal,(species,lengths) in enumerate(I.DEV_LENGTHS.items()):
        for seqid,length in lengths.items():
            for strand in ('+','-'):
                scopes.append(dict(species=species,seqid=seqid,strand=strand,ordinal=ordinal,
                    windows=[[a,min(a+I.M.WIDTH,length)] for a in I.M.starts(length)],
                    file='contig_'+str(ordinal)+'_'+strand+'.jsonl',
                    reference_used=False,labels_or_feature_cache_used=False))
    return dict(inference_complete=True,arm=I.ARM,configuration=I.R.CONFIG,
        optimizer_steps=4608,checkpoint=fitted['final_checkpoint'],encoder_model=I.F.MODEL,
        encoder_revision=I.F.REVISION,feature_spec=I.C.FULL_SPEC,input_mode='fixed_DEV_raw_FASTA',
        reference_used=False,labels_or_feature_cache_used=False,optimizer_steps_this_run=0,
        test_setaria_used=False,DEV_used_for_selection=False,windows=8690,scopes=scopes,
        lengths={q:n for qs in I.DEV_LENGTHS.values() for q,n in qs.items()})


class ContextRawTests(unittest.TestCase):
    def test_encoder_bridge_uses_train_bf16_materialization_once_before_fp32_head(self):
        frozen = torch.tensor([[.125,1.5],[-2.,3.]],dtype=torch.bfloat16)
        with patch.object(I.C,'encode_full_window',return_value=frozen) as native:
            features = I.encode_window('backbone','tokenizer','ATG','cpu')
        native.assert_called_once_with('backbone','tokenizer','ATG')
        self.assertEqual(features.dtype,torch.float32)
        self.assertTrue(torch.equal(features,frozen.float()[None]))
        self.assertFalse(features.requires_grad)

    def test_raw_orientation_normalization_single_window_forwards_and_gff_phase(self):
        oriented = 'ATGCGTTAGTTTAA'; model,features = fixture(oriented,[(0,4),(9,14)])
        generator = SyntheticGenerator(model)
        raw = I.M.rc(oriented).lower()+'r'; seen = []
        def encode(backbone,tokenizer,sequence,device):
            seen.append(sequence)
            # The ambiguous flank is outside this scripted complete CDS chain.
            values = torch.cat([features,torch.zeros(1,4)]) if sequence.startswith(oriented) else torch.zeros(15,4)
            if sequence.startswith('N'+oriented):
                values = torch.cat([torch.zeros(1,4),features])
            return values[None]
        with tempfile.TemporaryDirectory() as directory,patch.object(I,'encode_window',side_effect=encode):
            out = Path(directory)
            txs,scopes = I.annotate_contig(generator,None,None,'synthetic','chr',raw,0,out,'cpu')
            self.assertEqual(generator.calls,2)
            self.assertEqual(seen,[raw[:-1].upper()+'N','N'+oriented])
            self.assertEqual([s['windows'] for s in scopes],[[(0,15)],[(0,15)]])
            minus = [t for t in txs if t['strand'] == '-']
            self.assertEqual(len(minus),1)
            self.assertEqual(minus[0]['CDS'],[(0,5,'2'),(10,14,'0')])
            gff = out/'synthetic.gff3'; write_gff3(gff,txs,{'chr':15},'M28_A10')
            parsed = E.E.S.E.primary_transcripts(E.E.S.E.parse_annotation(gff,{'chr':15}),complete_only=False)
            E.A3.validate_gff_roundtrip(parsed,txs)
            for scope in scopes:
                rows = list(I.I.saved_records(out/scope['file']))
                self.assertEqual(len(rows),1)
                self.assertFalse(rows[0]['reference_used'])
                self.assertFalse(scope['labels_or_feature_cache_used'])

    def test_zero_gain_completed_chains_are_recorded_but_emit_no_genes(self):
        sequence = 'ATGAAATAA'; model,features = fixture(sequence,[(0,9)])
        with torch.no_grad():
            model.exists.bias.zero_()
        with tempfile.TemporaryDirectory() as directory,patch.object(I,'encode_window',return_value=features[None]):
            out = Path(directory)
            txs,scopes = I.annotate_contig(SyntheticGenerator(model),None,None,'synthetic','chr',sequence,0,out,'cpu')
            self.assertEqual(txs,[])
            row = next(I.I.saved_records(out/scopes[0]['file']))
            self.assertTrue(row['chains'])
            self.assertTrue(all(value == 0 for value in row['score']))
            self.assertEqual(scopes[0]['assembly_counts']['selected_unique'],0)

    def test_fit_rejects_early_checkpoint_and_old_context_despite_same_shape_revision(self):
        good = fitted_receipt(); I.validate_fit(good)
        for key,value in (('optimizer_steps',1536),('main_evaluation','best_DEV'),
                          ('configuration',{**I.R.CONFIG,'encoder_context':I.C.C.FEATURE_SPEC})):
            with self.subTest(key=key),self.assertRaises(ValueError):
                I.validate_fit({**good,key:value})

    def test_complete_receipt_rejects_old_features_partial_grid_and_explicit_input(self):
        fitted = fitted_receipt(); good = decoded_receipt(fitted); E.validate_complete(fitted,good)
        for key,value in (('feature_spec',I.C.C.FEATURE_SPEC),('scopes',good['scopes'][:-1]),
                          ('input_mode','explicit_raw_FASTA'),('labels_or_feature_cache_used',True)):
            with self.subTest(key=key),self.assertRaises(ValueError):
                E.validate_complete(fitted,{**good,key:value})
        truncated = copy.deepcopy(good); truncated['scopes'][0]['windows'].pop()
        with self.assertRaises(ValueError):
            E.validate_complete(fitted,truncated)

    def test_partial_actual_records_fail_before_reference_context_is_opened(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); fit = root/'fit'; infer = root/'infer'
            fit.mkdir(); infer.mkdir()
            with patch.object(I,'ROOT',root),patch.object(I.R,'FIT',fit),patch.object(I,'OUT',infer):
                fitted = fitted_receipt(); decoded = decoded_receipt(fitted)
                for folder,receipt in ((fit,fitted),(infer,decoded)):
                    (folder/'STATUS').write_text('COMPLETED\n')
                    (folder/'summary.json').write_text(json.dumps(receipt))
                (infer/decoded['scopes'][0]['file']).write_text('')
                with patch.object(E.A3,'reference_context') as reference:
                    with self.assertRaisesRegex(ValueError,'Partial inference scope'):
                        E.main()
                    reference.assert_not_called()


if __name__ == '__main__':
    unittest.main()
