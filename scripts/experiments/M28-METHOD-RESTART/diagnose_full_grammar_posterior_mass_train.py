#!/usr/bin/env python3
"""Read-only numeric localization of the observed whole-TRAIN mass failure."""
import json
import math
import os
import time

import numpy as np
import torch

import readout_full_grammar_gene_posterior_train as P
import run_full_grammar_e_m_energy_pool_train as Q
from src.m28.full_grammar_gene_posterior import _prefix_queries, _suffix_queries
from src.m28.full_grammar_boundary_compiled import ALPHABET, _state_table

OUT = P.ROOT/'outputs/M28-FULL-GRAMMAR-POSTERIOR-MASS-DIAGNOSTIC'
PROTOCOL = P.PROTOCOL


def main():
    began = time.perf_counter(); name, task, strand = 'E_RANK4', 0, '+'
    module, fit_job, reports = P.source(name)
    species, seqid, length = P.E.I.SCOPE[task]
    folder = Q.emissions(name, task, reports[seqid])
    P.write(OUT/'started.json', dict(job_id=os.environ['SLURM_JOB_ID'],protocol=PROTOCOL,
        model=name,task=task,strand=strand,seqid=seqid,length=length,source_emissions_job='13548914',
        source_failed_posterior_job='13548915',source_fit_job=fit_job,reference_used=False))
    sequence = P.C.M.read_allowed(module.OUT/seqid/'input.fasta',{seqid:length})[seqid]
    txs = P.predictions(module,seqid,length,reports[seqid])
    chains = [tuple((a,b) for a,b,_ in t['CDS']) for t in txs if t['strand']==strand]
    array = np.load(folder/'plus.npy',mmap_mode='c',allow_pickle=False)
    if array.shape != (length,15) or array.dtype != np.float32:
        raise ValueError('Exact full-chromosome original FP32 roles required')
    values = torch.from_numpy(array)
    if not torch.isfinite(values).all(): raise ValueError('Source roles must be finite')
    # Same unmodified arithmetic and positions as the failed production call.
    gain = (values.double()-values[:,:1].double()).numpy()
    codes = np.array([ALPHABET.find(base) if base in ALPHABET else 4 for base in sequence],dtype=np.int8)
    starts = np.array(sorted({0,*(c[0][0] for c in chains)}),dtype=np.int64)
    ends = np.array(sorted({length,*(c[-1][1] for c in chains)}),dtype=np.int64)
    states,destinations,roles,background = _state_table(128)
    forward_started = time.perf_counter()
    log_z,_ = _prefix_queries(codes,gain,starts,destinations,roles,background)
    forward_seconds = time.perf_counter()-forward_started
    print(json.dumps(dict(stage='original_forward_finished',log_partition=float(log_z))),flush=True)
    reverse_started = time.perf_counter()
    reverse_z,_ = _suffix_queries(codes,gain,ends,destinations,roles,background)
    reverse_seconds = time.perf_counter()-reverse_started
    finite = math.isfinite(log_z) and math.isfinite(reverse_z)
    P.write(OUT/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],protocol=PROTOCOL,
        model=name,task=task,strand=strand,species=species,seqid=seqid,length=length,
        source_emissions_job='13548914',source_failed_posterior_job='13548915',source_fit_job=fit_job,
        log_partition=float(log_z),reverse_log_partition=float(reverse_z),both_partitions_finite=finite,
        signed_difference=float(log_z-reverse_z),absolute_difference=float(abs(log_z-reverse_z)),
        original_rel_tol=1e-12,original_abs_tol=1e-8,
        original_allowed_difference=max(1e-8,1e-12*max(abs(float(log_z)),abs(float(reverse_z)))),
        original_forward_backward_check_passed=math.isclose(log_z,reverse_z,rel_tol=1e-12,abs_tol=1e-8),
        empty_set_energy_zero_retained=log_z>=-1e-8,state_count=len(states),
        original_algorithm_unchanged=True,original_tolerances_unchanged=True,
        source_DNA_equality_receipt_verified=True,whole_original_TRAIN_contig=True,
        reference_used=False,native_forward=False,checkpoint_read=False,predictions_modified=False,
        weights_modified=False,optimizer_steps_this_run=0,DEV_used=False,sealed_data_used=False,
        forward_seconds=forward_seconds,reverse_seconds=reverse_seconds,seconds=time.perf_counter()-began,
        interpretation='numeric_failure_localization_not_successful_posterior_or_accuracy_result'))


if __name__=='__main__':
    torch.set_num_threads(2)
    main()
