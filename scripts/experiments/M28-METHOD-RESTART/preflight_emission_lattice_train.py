#!/usr/bin/env python3
"""Eight metadata-first original TRAIN windows; frozen A3, no labels or fit."""
import argparse
import json
import resource
import time

import torch

import readout_structured_emission_train as B
import run_local_path_caller as R
from src.m28.core import select_nonoverlapping
from src.m28.emission_lattice import decode_emission_lattice
from src.m28.local_path_score import local_path_scores


ROOT = R.ROOT
PROTOCOL = 'reports/M28-METHOD-RESTART/a42_emission_lattice_train_protocol.md'
POTENTIALS = ROOT/'outputs/M28-A42-EMISSION-LATTICE-TRAIN-POTENTIAL-PREFLIGHT'
DECODE = ROOT/'outputs/M28-A42-EMISSION-LATTICE-TRAIN-DECODE-PREFLIGHT'
CONFIG = dict(source_checkpoint=R.C.GENERATOR, frozen_source_parameters=649457,
    source_optimizer_steps=4608, optimizer_steps_this_run=0, new_parameters=0,
    population='original_1491_unique_TRAIN_windows_first_occurrence_order',
    cases='sorted_species_plus_minus_first_two_windows_each', cases_count=8,
    model_precision='FP32', score_accumulation='FP64', max_exons=128,
    threshold=0., null=0., all_role_weights=1.,
    search='all_motif_legal_oriented_window_emission_lattice_no_query_or_beam')
STOP = {'TAA', 'TAG', 'TGA'}


def cases(free):
    selected = []
    for species in sorted(R.C.R.ALLOW):
        for strand in ('+', '-'):
            ordinals = [i for i, row in enumerate(free)
                        if (row['window']['species'], row['window']['strand']) == (species, strand)][:2]
            if len(ordinals) != 2:
                raise ValueError('Two original ordered windows per species/strand required')
            selected.extend(ordinals)
    if len(selected) != 8 or len(set(selected)) != 8:
        raise ValueError('Exactly eight fixed distinct metadata cases required')
    return selected


def check_chain(sequence, chain):
    if (not chain or len(chain) > 128 or chain[0][1] - chain[0][0] < 3 or
            chain[-1][1] - chain[-1][0] < 3 or
            any(not 0 <= a < b <= len(sequence) for a, b in chain)):
        raise ValueError('Complete original supported valid-length CDS required')
    for (_, b), (a, _) in zip(chain, chain[1:]):
        if a - b < 5 or sequence[b:b + 2] not in ('GT', 'GC') or sequence[a - 2:a] != 'AG':
            raise ValueError('Original GT/GC--AG intron support differs')
    coding = ''.join(sequence[a:b] for a, b in chain)
    codons = [coding[i:i + 3] for i in range(0, len(coding), 3)]
    if (any(c not in 'ACGT' for c in coding) or len(coding) % 3 or coding[:3] != 'ATG' or
            codons[-1] not in STOP or any(c in STOP for c in codons[:-1])):
        raise ValueError('Actual concatenated CDS ORF differs from fixed grammar')


@torch.inference_mode()
def potentials():
    B.F.A.require_jobs(['13421592', '13421848'])
    draws, free, prior, index, _ = B.population()
    ordinals = cases(free)
    model = R.C.load_generator().eval().requires_grad_(False)
    if sum(p.numel() for p in model.parameters()) != CONFIG['frozen_source_parameters']:
        raise ValueError('Original frozen A3 parameter support differs')
    before = {k: v.detach().clone() for k, v in model.state_dict().items()}
    began = time.perf_counter(); receipts = []; torch.cuda.reset_peak_memory_stats()
    for case, ordinal in enumerate(ordinals):
        pool, old = free[ordinal], prior[ordinal]['variants']['A31_A3']
        w = pool['window']; length = w['valid_bases']; started = time.perf_counter()
        cached = torch.load(ROOT/index[w['id']]['artifact'], map_location='cpu', weights_only=True)
        features, dna = R.Q.window_inputs(cached, w, 'cuda')
        hidden, _ = model.shared(features, dna, [length])
        local = model.segmentation(hidden[0])[:length]
        if local.dtype != torch.float32 or tuple(local.shape) != (length, 15) or not bool(torch.isfinite(local).all()):
            raise ValueError('Finite original FP32 valid local logits required')
        scores = local_path_scores(local, pool['chains']).cpu()
        expected = torch.tensor(old['scores'], dtype=torch.float64)
        torch.testing.assert_close(scores, expected, rtol=1e-12, atol=1e-8)
        if select_nonoverlapping(pool['chains'], scores.tolist()) != old['selected']:
            raise ValueError('Original stronger baseline selection differs')
        for chain in pool['chains']:
            check_chain(cached['sequence'].upper(), chain)
        path = POTENTIALS/('case_' + str(case) + '.pt')
        with path.open('xb') as handle:
            torch.save(dict(configuration=CONFIG, protocol=PROTOCOL, case=case, original_ordinal=ordinal,
                window=w, sequence=cached['sequence'], local_logits=local.cpu(),
                original_chains=pool['chains'], original_scores=old['scores'],
                original_selected=old['selected'], revision=R.C.I.F.REVISION,
                reference_used_in_scoring=False), handle)
        torch.cuda.synchronize()
        row = dict(case=case, original_ordinal=ordinal, window=w,
            artifact=str(path.relative_to(ROOT)), original_candidates=len(pool['chains']),
            original_A31_score_max_error=float((scores - expected).abs().max()) if len(scores) else 0.,
            original_selected_indices_exact=True, all_original_pool_chains_grammar_supported=True,
            seconds=time.perf_counter() - started)
        receipts.append(row); print(json.dumps(row, allow_nan=False), flush=True)
        del cached, features, dna, hidden, local, scores, expected
    if (any(not torch.equal(v, before[k]) for k, v in model.state_dict().items()) or
            any(p.requires_grad or p.grad is not None for p in model.parameters())):
        raise ValueError('Frozen source changed or received gradients')
    R.C.R.write_summary(POTENTIALS, dict(completed=True, configuration=CONFIG, protocol=PROTOCOL,
        source_preparation_jobs=['13421592', '13421848'], windows=1491, original_draws=len(draws),
        original_natural_candidates=sum(len(x['chains']) for x in free), records=receipts,
        frozen_source_weights_unchanged=True, all_source_parameters_frozen_and_gradients_absent=True,
        all_potential_files_closed=True, feature_revision=R.C.I.F.REVISION,
        optimizer_steps_this_run=0, checkpoints_written=False, query_decoding=False, generator_updated=False,
        new_candidate_search=False, new_backbone_forward=False, reference_labels_opened=False,
        TRAIN_catalog_opened=False, new_R4_labelling=False, DEV_used=False, test_setaria_used=False,
        A33_outputs_opened=False, seconds=time.perf_counter() - began,
        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),
        limitation='Eight metadata-first frozen TRAIN emissions; no A42 decode or performance result.'))


@torch.inference_mode()
def decode(case, potential_job):
    if not 0 <= case < 8:
        raise ValueError('Only fixed cases 0 through 7 are permitted')
    B.F.A.require_jobs(['13421592', '13421848', potential_job])
    report = R.completed(POTENTIALS, potential_job)
    if (report['configuration'] != CONFIG or report['protocol'] != PROTOCOL or
            len(report['records']) != 8 or not report['all_potential_files_closed'] or
            not report['frozen_source_weights_unchanged'] or
            not report['all_source_parameters_frozen_and_gradients_absent']):
        raise ValueError('Actual complete fixed eight-potential source required')
    receipt = report['records'][case]
    value = torch.load(ROOT/receipt['artifact'], map_location='cpu', weights_only=True)
    if (value['configuration'] != CONFIG or value['protocol'] != PROTOCOL or value['case'] != case or
            value['window'] != receipt['window'] or value['original_ordinal'] != receipt['original_ordinal'] or
            value['revision'] != report['feature_revision'] or value['reference_used_in_scoring']):
        raise ValueError('Actual fixed corresponding TRAIN potential case required')
    w = value['window']; sequence = value['sequence'].upper(); local = value['local_logits']
    if (w['split'] != 'train' or w['seqid'] not in R.C.R.ALLOW[w['species']] or
            len(sequence) != w['valid_bases'] or local.dtype != torch.float32 or
            tuple(local.shape) != (len(sequence), 15)):
        raise ValueError('Original TRAIN split and FP32 local geometry required')
    out = DECODE/('case_' + str(case))
    with (out/'started.json').open('x') as handle:
        json.dump(dict(case=case, window=w, potential_job=potential_job,
                       protocol=PROTOCOL, configuration=CONFIG), handle, allow_nan=False)
    baseline = local_path_scores(local, value['original_chains'])
    torch.testing.assert_close(baseline, torch.tensor(value['original_scores'], dtype=torch.float64),
                               rtol=1e-12, atol=1e-8)
    if select_nonoverlapping(value['original_chains'], baseline.tolist()) != value['original_selected']:
        raise ValueError('Actual CPU original baseline selection differs')
    for chain in value['original_chains']:
        check_chain(sequence, chain)
    began = time.perf_counter()
    natural = decode_emission_lattice(sequence, local, CONFIG['max_exons'])
    seconds = time.perf_counter() - began
    previous_end = 0
    for chain in natural['chains']:
        check_chain(sequence, chain)
        if chain[0][0] < previous_end:
            raise ValueError('Selected complete-gene spans overlap')
        previous_end = chain[-1][1]
    direct = local_path_scores(local, natural['chains'])
    torch.testing.assert_close(direct, torch.tensor(natural['chain_gains'], dtype=torch.float64),
                               rtol=1e-12, atol=1e-8)
    torch.testing.assert_close(direct.sum(), torch.tensor(natural['gain'], dtype=torch.float64),
                               rtol=1e-12, atol=1e-8)
    old_gain = float(baseline[value['original_selected']].sum())
    tolerance = 1e-8 + 1e-12 * abs(old_gain)
    if (natural['reference_used'] or any(s <= 0 for s in natural['chain_gains']) or
            natural['gain'] + tolerance < old_gain):
        raise ValueError('Reference-free full grammar must attain original supported pool objective')
    motifs = dict(start=sum(sequence[i:i + 3] == 'ATG' for i in range(len(sequence) - 2)),
        donor=sum(sequence[i:i + 2] in ('GT', 'GC') for i in range(1, len(sequence) - 1)),
        acceptor=sum(sequence[i:i + 2] == 'AG' for i in range(len(sequence) - 2)),
        stop=sum(sequence[i:i + 3] in STOP for i in range(len(sequence) - 2)))
    R.C.R.write_summary(out, dict(completed=True, configuration=CONFIG, protocol=PROTOCOL,
        case=case, original_ordinal=value['original_ordinal'], window=w, potential_job=potential_job,
        feature_revision=value['revision'], motifs=motifs, natural=natural,
        original_baseline_gain=old_gain, original_baseline_selected_indices_exact=True,
        original_A31_score_max_error=float((baseline - torch.tensor(value['original_scores'],
            dtype=torch.float64)).abs().max()) if len(baseline) else 0.,
        selected_ORF_motif_cap_span_and_direct_role_gains_verified=True,
        all_original_pool_chains_grammar_supported=True, full_grammar_objective_not_lower_than_baseline=True,
        decode_seconds=seconds, process_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        optimizer_steps_this_run=0, checkpoints_written=False, model_or_backbone_forward=False,
        new_candidate_search=True, reference_labels_opened=False, TRAIN_catalog_opened=False,
        new_R4_labelling=False, DEV_used=False, test_setaria_used=False, A33_outputs_opened=False,
        limitation='One metadata-first TRAIN engineering case; no labels, performance, raw assembly or generalization.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('potentials', 'decode'), required=True)
    parser.add_argument('--case', type=int); parser.add_argument('--potential-job')
    args = parser.parse_args(); torch.set_num_threads(2)
    if args.mode == 'potentials':
        potentials()
    elif args.case is not None and args.potential_job:
        decode(args.case, args.potential_job)
    else:
        parser.error('Fixed case and actual potential job required for decode')
