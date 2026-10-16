#!/usr/bin/env python3
"""One whole original E/plus native qualification, not a production replacement."""
import math
import os
import time
from unittest.mock import patch

import numpy as np
import torch

import readout_full_grammar_gene_posterior_train as P
import run_full_grammar_e_m_energy_pool_train as Q
import src.m28.full_grammar_gene_posterior as G
from src.m28.full_grammar_scaled_mass import _scaled_prefix_queries, _scaled_suffix_queries

OUT = P.ROOT/'outputs/M28-FULL-GRAMMAR-SCALED-MASS-NATIVE'
CHECKS = P.ROOT/'outputs/M28-FULL-GRAMMAR-SCALED-MASS-CHECKS/13552225'


def main():
    began = time.perf_counter(); qualification = P.read(CHECKS/'summary.json')
    if ((CHECKS/'STATUS').read_text().strip() != 'COMPLETED' or not qualification['completed']
            or qualification['job_id'] != '13552225' or qualification['protocol'] != P.PROTOCOL
            or qualification['tests_run'] != 8 or qualification['failures'] or qualification['errors']
            or not qualification['original_probability_and_mass_tolerances_unchanged']
            or qualification['production_algorithm_modified']):
        raise ValueError('Completed eight-test candidate qualification required')
    old = P.read(P.ROOT/'outputs/M28-FULL-GRAMMAR-POSTERIOR-MASS-DIAGNOSTIC/summary.json')
    if (not old['completed'] or old['job_id'] != '13551206'
            or old['source_emissions_job'] != '13548914' or not old['original_algorithm_unchanged']
            or not old['original_tolerances_unchanged'] or not old['both_partitions_finite']
            or old['original_forward_backward_check_passed']):
        raise ValueError('Completed actual unchanged-algorithm failure localization required')
    module,fit_job,reports = P.source('E_RANK4')
    species,seqid,length = P.E.I.SCOPE[0]
    folder = Q.emissions('E_RANK4',0,reports[seqid])
    P.write(OUT/'started.json',dict(job_id=os.environ['SLURM_JOB_ID'],protocol=P.PROTOCOL,
        checks_job='13552225',source_mass_diagnostic_job='13551206',source_emissions_job='13548914',
        source_fit_job=fit_job,model='E_RANK4',task=0,strand='+',seqid=seqid,length=length,
        production_algorithm_modified=False,reference_used=False,DEV_used=False,sealed_data_used=False))
    sequence = P.C.M.read_allowed(module.OUT/seqid/'input.fasta',{seqid:length})[seqid]
    txs = sorted((t for t in P.predictions(module,seqid,length,reports[seqid]) if t['strand']=='+'),
        key=lambda t:t['CDS'][0][0])
    chains = [tuple((a,b) for a,b,_ in t['CDS']) for t in txs]
    array = np.load(folder/'plus.npy',mmap_mode='c',allow_pickle=False)
    if array.shape != (length,15) or array.dtype != np.float32:
        raise ValueError('Exact original full-chromosome FP32 roles required')
    values = torch.from_numpy(array)
    # A qualification process substitutes only the two tested numerical kernels.
    # The source file and all production consumers stay unchanged.
    with patch.multiple(G,_prefix_queries=_scaled_prefix_queries,_suffix_queries=_scaled_suffix_queries):
        result = G.complete_gene_set_posteriors(sequence,values,chains,128)
    for tx,gain in zip(txs,result['chain_gains']):
        if not math.isclose(gain,float(tx['score']),rel_tol=6e-8,abs_tol=1e-5):
            raise ValueError('Original whole-chromosome MAP GFF gain must remain unchanged')
    z,rz = result['log_partition'],result['reverse_log_partition']
    P.write(OUT/'summary.json',dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],protocol=P.PROTOCOL,
        checks_job='13552225',source_mass_diagnostic_job='13551206',source_emissions_job='13548914',
        source_fit_job=fit_job,model='E_RANK4',task=0,strand='+',species=species,seqid=seqid,length=length,
        candidate=qualification['candidate'],log_partition=z,reverse_log_partition=rz,
        absolute_difference=abs(z-rz),original_allowed_difference=max(1e-8,1e-12*max(abs(z),abs(rz))),
        original_forward_backward_check_passed=math.isclose(z,rz,rel_tol=1e-12,abs_tol=1e-8),
        unscaled_log_partition=old['log_partition'],unscaled_reverse_log_partition=old['reverse_log_partition'],
        unscaled_absolute_difference=old['absolute_difference'],
        all_queried_genes_literal_DNA_legal=True,original_MAP_GFF_gains_preserved=True,
        denominator_all_legal_gene_sets_not_selected_pool=True,
        original_probability_and_mass_tolerances_unchanged=True,production_algorithm_modified=False,
        whole_original_TRAIN_contig_one_plus_orientation=True,other_eleven_orientations_not_qualified=True,
        genes=[dict(id=t['id'],gain=g,log_probability=p) for t,g,p in
            zip(txs,result['chain_gains'],result['gene_log_probabilities'])],
        queried_set_log_probability=result['queried_set_log_probability'],
        empty_set_log_probability=result['empty_set_log_probability'],state_count=result['state_count'],
        source_DNA_equality_receipt_verified=True,reference_used=False,checkpoint_read=False,native_forward=False,
        predictions_modified=False,weights_modified=False,optimizer_steps_this_run=0,DEV_used=False,
        sealed_data_used=False,seconds=time.perf_counter()-began,
        interpretation='one_original_whole_TRAIN_plus_numerical_qualification_not_accuracy_or_full_twelve_readouts'))


if __name__=='__main__':
    torch.set_num_threads(2)
    main()
