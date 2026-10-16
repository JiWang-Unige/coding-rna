"""One fixed CPU synthetic teacher-cost measurement; not a training benchmark."""
import argparse
import json
from pathlib import Path
from time import perf_counter
import torch
from src.m28.query_pointer import GeneQueryPointer
from src.m28.query_training import teacher_chain_components


# TRAIN metadata informs size only. No biological sequence or reference is read.
LENGTH, QUERIES, HIDDEN = 24576, 11, 128
PARITY_ATOL = 1e-5


@torch.no_grad()
def measure():
    torch.set_num_threads(1)
    torch.manual_seed(0)
    model = GeneQueryPointer(num_queries=QUERIES, hidden=HIDDEN)
    window = model.for_window(torch.randn(LENGTH, HIDDEN))
    chain = tuple((j*384, j*384+96) for j in range(60)) + ((23040, LENGTH),)

    # Fixed scalar-then-batch order, one measurement each; initialization shared.
    start = perf_counter()
    scalar = []
    for i in range(QUERIES):
        labels = [-1]*QUERIES
        labels[i] = 1
        scalar.append(window.teacher_components(labels, {i: chain}))
    scalar_seconds = perf_counter()-start
    start = perf_counter()
    batched = teacher_chain_components(window, chain)
    batched_seconds = perf_counter()-start

    errors = {}
    for name in ('boundary_ce', 'termination_bce'):
        expected = torch.stack([part[name] for part in scalar])
        if not torch.isfinite(expected).all() or not torch.isfinite(batched[name]).all():
            raise ValueError('Nonfinite synthetic teacher component')
        errors[name] = float((batched[name]-expected).abs().max())
        if errors[name] > PARITY_ATOL:
            raise ValueError(f'{name} parity exceeded fixed tolerance: {errors[name]}')
    return {'scope': 'single_fixed_order_cpu_synthetic_teacher_cost_only',
            'seed': 0, 'cpu_threads': 1, 'length': LENGTH, 'queries': QUERIES,
            'hidden': HIDDEN, 'exons': len(chain), 'parity_atol': PARITY_ATOL,
            'component_max_abs_errors': errors,
            'scalar_seconds': scalar_seconds, 'batched_seconds': batched_seconds,
            'scalar_over_batched': scalar_seconds/batched_seconds,
            'torch_version': torch.__version__,
            'excluded_from_timing': ['window_initialization', 'backbone', 'assignment',
                                     'gradient_readout', 'backward', 'optimizer',
                                     'genome_assembly', 'io'],
            'limitations': ['random_untrained_tensors_and_authored_chain',
                            'single_measurement_not_repeated_or_order_balanced',
                            'preceding_tests_warm_dispatch',
                            'scalar_includes_original_existence_bce_and_label_checks',
                            'not_real_data_or_gpu_training_throughput',
                            'queries_11_is_fixture_not_frozen_model_capacity']}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    result = measure()
    with (args.output_dir / 'timing.json').open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
