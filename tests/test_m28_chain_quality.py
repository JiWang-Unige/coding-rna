import sys
import unittest
from pathlib import Path
import torch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.m28.query_beam_pool import completed_union
from src.m28.training import inject_training
from src.m28.core import chain_loss


class ChainQualityTests(unittest.TestCase):
    def test_union_preserves_greedy_only_completed_and_traceable_duplicates(self):
        g = [{'terminated':True,'chain':[(10,16)]},
             {'terminated':True,'chain':[(0,6)]},
             {'terminated':True,'chain':[(0,6)]},
             {'terminated':False,'chain':[(30,36)]}]
        b = [{'completed':[{'chain':[(0,6)]},{'chain':[(20,26)]}]},
             {'completed':[{'chain':[(20,26)]}]}]
        chains,sources = completed_union(g,b)
        self.assertEqual(chains,[[(0,6)],[(10,16)],[(20,26)]])
        self.assertEqual(sources,[['beam4','greedy'],['greedy'],['beam4']])
        self.assertEqual(completed_union(list(reversed(g)),list(reversed(b))),(chains,sources))

    def test_training_injection_never_relabels_free_origin(self):
        free,sources = completed_union([{'terminated':True,'chain':[(0,6)]}],[])
        pool,origins = inject_training(free,[[(0,6)],[(10,16)]])
        self.assertEqual(pool,[[(0,6)],[(10,16)]])
        self.assertEqual(origins,['free','reference_injected'])
        self.assertEqual(free,[[(0,6)]])
        self.assertEqual(sources,[['greedy']])

    def test_known_mean_excludes_protected_gradient_and_all_unknown_is_zero(self):
        logits = torch.tensor([0.,1.,.5],requires_grad=True)
        null = torch.tensor(.2,requires_grad=True)
        scores = {'logits':logits,'null':null}
        loss = chain_loss(scores,[1,0,-1])
        expected = torch.nn.functional.binary_cross_entropy_with_logits(logits[:2],torch.tensor([1.,0.]))
        self.assertTrue(torch.equal(loss,expected))
        loss.backward()
        self.assertAlmostEqual(float(logits.grad[0]),-.25)
        self.assertAlmostEqual(float(logits.grad[1]),float(torch.sigmoid(logits[1]))/2,places=7)
        self.assertEqual(float(logits.grad[2]),0.)
        logits.grad = None
        zero = chain_loss(scores,[-1,-1,-1])
        zero.backward()
        self.assertEqual(float(zero),0.)
        self.assertTrue(torch.equal(logits.grad,torch.zeros(3)))
        self.assertEqual(float(null.grad),0.)


if __name__=='__main__':
    unittest.main()
