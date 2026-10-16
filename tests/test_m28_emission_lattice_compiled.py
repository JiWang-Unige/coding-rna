"""CPU-JIT scalar recurrence vs original answers/counters; artificial inputs only."""
import json
import time
import unittest

import numba
import numpy as np
import torch

import test_m28_emission_lattice_batched as B
from src.m28.emission_lattice import _Margins, decode_emission_lattice
from src.m28.emission_lattice_compiled import _solve, decode_compiled_emission_lattice, emission_arrays
from test_m28_emission_lattice import oracle_gain, prefer_chain


class CompiledEmissionLatticeTests(B.BatchedEmissionLatticeTests):
    @classmethod
    def setUpClass(cls):
        began = time.perf_counter()
        sequence = 'ATGAAATAA'
        result = decode_compiled_emission_lattice(sequence, torch.zeros(9, 15), 1)
        if result != decode_emission_lattice(sequence, torch.zeros(9, 15), 1):
            raise AssertionError('Cold compiled zero-null result/counters differ')
        print(json.dumps(dict(stage='artificial_nopython_kernel_ready', numba_version=numba.__version__,
            cold_call_including_compile_seconds=time.perf_counter()-began,
            nopython_signatures=len(_solve.nopython_signatures), fastmath=_solve.targetoptions['fastmath'],
            research_inputs_opened=False, optimizer_steps=0)), flush=True)

    def compare(self, sequence, potentials, cap, oracle=False):
        actual = decode_compiled_emission_lattice(sequence, potentials, cap)
        old = decode_emission_lattice(sequence, potentials, cap)
        self.assertEqual(actual, old)  # Every chain/order/gain and all eight counters.
        self.assertTrue(_solve.nopython_signatures); self.assertIs(_solve.targetoptions['fastmath'], False)
        if oracle:
            expected, legal = oracle_gain(sequence, potentials, cap)
            self.assertEqual(actual['gain'], expected)
            self.assertTrue(all(chain in legal for chain in actual['chains']))
        return actual

    def test_existing_torch_margin_arrays_match_all_FP64_values(self):
        generator = torch.Generator().manual_seed(83)
        for potentials in (torch.randn(128, 15, generator=generator),
                torch.randn(128, 15, generator=generator)*torch.logspace(-30, 30, 128)[:, None]):
            old = _Margins(potentials)
            for actual, expected in zip(emission_arrays(potentials), (old.gain, old.prefix, old.coding)):
                self.assertTrue(np.array_equal(actual, np.asarray(expected)))

    def test_None_long_phase_keeps_actual_short_carry_and_cap(self):
        sequence = 'ATGT'+'GTTAG'+'A'+'GTTAG'+'CAAATAA'
        chain = ((0, 4), (9, 10), (15, 22)); potentials = prefer_chain(sequence, chain).float()
        self.assertEqual(self.compare(sequence, potentials, 3, oracle=True)['chains'], (chain,))
        self.assertNotIn(chain, self.compare(sequence, potentials, 2, oracle=True)['chains'])

    def test_dense_small_random_FP32_integer_and_noninteger_identity(self):
        generator = torch.Generator().manual_seed(71)
        for _ in range(12):
            sequence = 'ATG' + ''.join('ACGT'[i] for i in torch.randint(0, 4, (768,), generator=generator).tolist()) + 'TAA'
            for cap in (2, 4, 8):
                for potentials in (torch.randint(-3, 4, (len(sequence), 15), generator=generator).float(),
                                   torch.randn(len(sequence), 15, generator=generator)):
                    self.compare(sequence, potentials, cap)


if __name__ == '__main__':
    torch.set_num_threads(2); unittest.main()
