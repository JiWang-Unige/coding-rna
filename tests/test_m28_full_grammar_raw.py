"""New raw caller wiring with artificial emissions, not a trained GLM run."""
import copy
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import torch
import annotate_full_grammar_coarse_channel as R
from test_m28_emission_lattice import prefer_chain
from test_m28_emission_lattice_owner import sequence_with_genes


class Head:
    def __init__(self,matrices): self.matrices=matrices; self.calls=0
    def __call__(self,features,dna,lengths):
        assert int(features)==self.calls
        value=self.matrices[self.calls]; self.calls+=1
        assert lengths==[len(value)]
        return {'segmentation_logits':value[None]}


class FullGrammarRawTests(unittest.TestCase):
    def invoke(self,sequence,logits,target):
        windows=[(a,min(a+R.M.WIDTH,len(sequence))) for a in R.M.starts(len(sequence))]
        matrices=[]; strings=[]
        for strand in ('+','-'):
            for a,b in windows:
                strings.append(sequence[a:b] if strand=='+' else R.M.rc(sequence[a:b]))
                if strand!=target: value=torch.zeros(b-a,15)
                elif strand=='+': value=logits[a:b].clone()
                else: value=logits[len(sequence)-b:len(sequence)-a].clone()
                matrices.append(value)
        head=Head(matrices); calls=[]
        def encode(backbone,tokenizer,dna,device):
            self.assertEqual(dna,strings[len(calls)]); self.assertEqual(device,'cpu')
            value=torch.tensor(len(calls)); calls.append(dna); return value
        with patch.object(R,'encode_window',encode):
            txs,scopes=R.annotate_contig(head,None,None,'toy',sequence,0,'cpu')
        self.assertEqual(head.calls,len(windows)*2); self.assertEqual(len(calls),head.calls)
        other=next(r for r in scopes if r['strand']!=target)
        self.assertEqual((other['selected_complete_CDS_chains'],other['gain']),(0,0.))
        self.assertTrue(all(r['joint_whole_oriented_contig_decode'] and not r['reference_used'] for r in scopes))
        return txs,scopes

    def exported_CDS(self,txs,length):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'toy.gff3'; R.write_gff3(path,txs,{'toy':length},'M28_FULL')
            lines=path.read_text().splitlines()
        self.assertEqual(lines[0],'##gff-version 3')
        self.assertEqual(lines[1],'##sequence-region toy 1 '+str(length))
        self.assertEqual(sum('\tgene\t' in x for x in lines),len(txs))
        self.assertEqual(sum('\tmRNA\t' in x for x in lines),len(txs))
        return [x.split('\t') for x in lines if '\tCDS\t' in x]

    def test_long_gene_crosses_native_window_and_exports_both_strands(self):
        chain=((3000,33000),); oriented=sequence_with_genes(36864,[chain])
        logits=prefer_chain(oriented,chain).float()
        for strand in ('+','-'):
            sequence=oriented if strand=='+' else R.M.rc(oriented)
            txs,scopes=self.invoke(sequence,logits,strand)
            self.assertEqual(len(txs),1); tx=txs[0]
            expected=[(3000,33000,'0')] if strand=='+' else [(3864,33864,'0')]
            self.assertEqual(tx['CDS'],expected); self.assertEqual(tx['strand'],strand)
            self.assertEqual(tx['score'],300000.)
            row=self.exported_CDS(txs,len(sequence))[0]
            self.assertEqual((int(row[3]),int(row[4]),row[6],row[7]),(expected[0][0]+1,expected[0][1],strand,'0'))
            self.assertEqual(next(s for s in scopes if s['strand']==strand)['selected_complete_CDS_chains'],1)

    def test_actual_codon_carry_minus_GFF_phase_and_intronic_N(self):
        oriented='ATGCGTNAGTTTAA'; chain=((0,4),(9,14))
        txs,_=self.invoke(R.M.rc(oriented),prefer_chain(oriented,chain).float(),'-')
        self.assertEqual(len(txs),1)
        self.assertEqual(txs[0]['CDS'],[(0,5,'2'),(10,14,'0')])
        rows=self.exported_CDS(txs,14)
        self.assertEqual([(int(r[3]),int(r[4]),r[6],r[7]) for r in rows],[(1,5,'-','2'),(11,14,'-','0')])

    def test_shifted_final_grid_owner_rows_without_role_column_permutation(self):
        length=40000; sequence='N'*length
        windows=[(a,min(a+R.M.WIDTH,length)) for a in R.M.starts(length)]
        self.assertEqual(windows,[(0,24576),(12288,36864),(15424,40000)])
        for strand in ('+','-'):
            matrices=[torch.arange((b-a)*15).reshape(b-a,15).float()+i*1000000 for i,(a,b) in enumerate(windows)]
            head=Head(matrices); encoded=[]
            def encode(backbone,tokenizer,dna,device):
                value=torch.tensor(len(encoded)); encoded.append(dna); return value
            with patch.object(R,'encode_window',encode):
                joined=R.contig_emissions(head,None,None,sequence,windows,strand,'cpu')
            expected=torch.empty(length,15)
            for q in range(length):
                p=q if strand=='+' else length-1-q
                owner=min((i for i,(a,b) in enumerate(windows) if a<=p<b),
                    key=lambda i:(abs(2*p+1-sum(windows[i])),windows[i][0]))
                a,b=windows[owner]; local=p-a if strand=='+' else b-1-p
                expected[q]=matrices[owner][local]
            self.assertTrue(torch.equal(joined['local_potentials'],expected))
            self.assertEqual(head.calls,3)

    def test_checkpoint_identity_rejects_old_loss_intermediate_step_and_cap(self):
        state=dict(configuration=dict(R.IDENTITY),step=4608,epoch=3,protocol=R.PROTOCOL,encoder_revision=R.F.REVISION)
        self.assertEqual(R.validate_checkpoint(state),R.IDENTITY)
        for key,value in [('step',1536),('epoch',2),('protocol','old_protocol')]:
            bad=copy.deepcopy(state); bad[key]=value
            with self.assertRaises(ValueError): R.validate_checkpoint(bad)
        for key,value in [('loss','local_CE'),('max_exons',64),('epsilon',0.2),('encoder_revision','other_revision')]:
            bad=copy.deepcopy(state); bad['configuration'][key]=value
            with self.assertRaises(ValueError): R.validate_checkpoint(bad)


if __name__=='__main__':
    torch.set_num_threads(2); unittest.main()
