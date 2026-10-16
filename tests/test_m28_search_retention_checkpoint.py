"""A30 continuation preparation: synthetic CPU math, no actual fit or source writes."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/experiments/M28-METHOD-RESTART'))
import inspect_search_retention_checkpoint as C


def toy_update(model, optimizer, ordinal):
    optimizer.zero_grad(set_to_none=True)
    inputs = torch.tensor([[1., float(ordinal+1)]])
    (model(inputs)-float(ordinal)/3).square().mean().backward()
    torch.nn.utils.clip_grad_norm_(model.parameters(), 1.)
    optimizer.step()


class CheckpointPreparationTests(unittest.TestCase):
    def test_AdamW_resume_matches_uninterrupted_next_update_not_fresh_optimizer(self):
        torch.manual_seed(0); original = torch.nn.Linear(2, 1)
        optimizer = torch.optim.AdamW(original.parameters(), lr=3e-4, weight_decay=.01, betas=(.9, .999), eps=1e-8)
        for i in range(2):
            toy_update(original, optimizer, i)
        weights, moments = copy.deepcopy(original.state_dict()), copy.deepcopy(optimizer.state_dict())
        resumed, reset = torch.nn.Linear(2, 1), torch.nn.Linear(2, 1)
        resumed.load_state_dict(weights); reset.load_state_dict(weights)
        loaded = C.load_optimizer(resumed, moments, 2)
        fresh = torch.optim.AdamW(reset.parameters(), lr=3e-4, weight_decay=.01, betas=(.9, .999), eps=1e-8)
        for model, opt in ((original, optimizer), (resumed, loaded), (reset, fresh)):
            toy_update(model, opt, 2)
        self.assertTrue(all(torch.equal(v, resumed.state_dict()[k]) for k, v in original.state_dict().items()))
        self.assertTrue(any(not torch.equal(v, reset.state_dict()[k]) for k, v in original.state_dict().items()))

    def test_invalid_optimizer_moments_or_hyperparameters_do_not_resume(self):
        torch.manual_seed(0); model = torch.nn.Linear(2, 1)
        optimizer = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=.01, betas=(.9, .999), eps=1e-8)
        toy_update(model, optimizer, 0); good = optimizer.state_dict()
        for kind in ('lr', 'nan', 'shape', 'step', 'empty'):
            bad = copy.deepcopy(good)
            if kind == 'lr':
                bad['param_groups'][0]['lr'] = .001
            elif kind == 'empty':
                bad['state'] = {}
            else:
                row = next(iter(bad['state'].values()))
                if kind == 'nan': row['exp_avg'].flatten()[0] = float('nan')
                if kind == 'shape': row['exp_avg'] = torch.zeros(8)
                if kind == 'step': row['step'] = torch.tensor(3.)
            with self.subTest(kind=kind), self.assertRaises(ValueError):
                C.load_optimizer(model, bad, 1)

    def test_checkpoint_prefix_keeps_original_order_weight_and_epoch_boundary(self):
        species = list(C.R.A.ALLOW)
        windows = {s: dict(id=s, species=s, split='train', seqid=sorted(C.R.A.ALLOW[s])[0]) for s in species}
        draws = [dict(window_id=species[i % 2], importance_weight_uniform_within_stratum=1.+i/1536)
                 for i in range(1536)]
        rows = []
        for i in range(1537):
            epoch, ordinal = divmod(i, 1536); draw = draws[ordinal]; s = draw['window_id']
            rows.append(dict(step=i+1, epoch=epoch+1, draw_ordinal=ordinal, window_id=s, species=s,
                importance_weight=draw['importance_weight_uniform_within_stratum'], reference_injected_candidates=0,
                losses=dict(local=1., existence_bce=1., boundary_ce=1., termination_bce=1.),
                weighted_original_total=4.*draw['importance_weight_uniform_within_stratum'],
                weighted_total=4.*draw['importance_weight_uniform_within_stratum'],
                search_retention_margin=0., gradient_norm_before_clip=1.))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'training.jsonl'
            def save(values):
                path.write_text(''.join(json.dumps(row)+'\n' for row in values)+'{discarded partial tail')
            save(rows); self.assertEqual(list(C.validated_prefix(path, 1537, draws, windows)), rows)
            self.assertEqual((rows[-1]['epoch'], rows[-1]['draw_ordinal']), (2, 0))
            for key, value in (('draw_ordinal', 1), ('importance_weight', 7.), ('epoch', 1),
                               ('window_id', species[1]), ('reference_injected_candidates', 1),
                               ('weighted_total', float('nan'))):
                bad = copy.deepcopy(rows); bad[-1][key] = value; save(bad)
                with self.subTest(key=key), self.assertRaises(ValueError):
                    list(C.validated_prefix(path, 1537, draws, windows))

    def test_latest_checkpoint_requires_actual_timeout_before_selection(self):
        for state in ('RUNNING', 'PENDING', 'COMPLETED', 'FAILED'):
            with patch.object(C, 'source_receipt', return_value=dict(job=C.SOURCE_JOB, state=state, exit_code='0:0')):
                with self.subTest(state=state), self.assertRaisesRegex(ValueError, 'TIMEOUT'):
                    C.latest_timeout_checkpoint()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'step_000256.pt').touch(); (root/'step_000512.pt').touch()
            receipt = C.root_receipt(C.SOURCE_JOB+'|TIMEOUT|0:15|\n'+C.SOURCE_JOB+'.batch|CANCELLED|0:15|\n')
            with patch.object(C, 'SOURCE', root), patch.object(C, 'source_receipt', return_value=receipt):
                self.assertEqual(C.latest_timeout_checkpoint(), root/'step_000512.pt')
            with self.assertRaises(ValueError):
                C.root_receipt(C.SOURCE_JOB+'.batch|TIMEOUT|0:15|\n')


if __name__ == '__main__':
    unittest.main()
