#!/usr/bin/env python3
"""Prepared fixed DEV comparison, not permission to access or execute DEV."""
import json
import os
import re
import subprocess
import time
import common_ruler as C
import infer_full_grammar_fixed_dev as I
import run_native_baseline as N
from evaluate_full_grammar_original_train import aggregate

OUT = I.ROOT/'outputs/M28-FULL-GRAMMAR-FIXED-DEV-COMPARISON'
EXPECTED = {'arabidopsis_thaliana':5437,'oryza_sativa':2291}
METHODS = ('M28_final4608',*N.METHODS)


def check_candidate(outer,raw,task,species,seqid,length,*,checkpoint=I.CHECKPOINT,
                    identity=I.R.IDENTITY,protocol=I.R.PROTOCOL):
    if (not outer['inference_complete'] or
            (outer['task'],outer['species'],outer['seqid'],outer['length']) != (task,species,seqid,length) or
            outer['checkpoint'] != str(checkpoint) or outer['raw_scope'] != I.RAW_SCOPE or
            not outer['DEV_used'] or any(outer[k] for k in
                ('sealed_data_used','reference_used','labels_or_feature_cache_used')) or
            outer['optimizer_steps_this_run'] != 0):
        raise ValueError('Fixed final-checkpoint DEV run identity differs')
    if (not raw['inference_complete'] or raw['checkpoint'] != str(checkpoint) or
            raw['protocol'] != protocol or raw['lengths'] != {seqid:length} or
            any(raw['configuration'].get(k) != v for k,v in identity.items()) or
            any(raw[k] for k in ('reference_used','labels_or_feature_cache_used','observation_model_used_at_inference')) or
            raw['optimizer_steps_this_run'] != 0 or raw['backend'] != 'exact_exon_budget_skyline' or
            len(raw['scopes']) != 2 or {s['strand'] for s in raw['scopes']} != {'+','-'}):
        raise ValueError('Raw caller identity, scope or label-free boundary differs')
    expected = [[a,min(a+C.M.WIDTH,length)] for a in C.M.starts(length)]
    for scope in raw['scopes']:
        if (scope['seqid'] != seqid or scope['windows'] != expected or scope['reference_used'] or
                not scope['joint_whole_oriented_contig_decode'] or scope['context_policy'] != I.R.CONTEXT or
                scope['pruning'] != 'exact_same_stage_actual_carry_lower_exon_budget_dominance_max_only' or
                scope['score_accumulation'] != 'FP64_per_base_forward_order' or
                scope['tie_policy'] != 'background_continuation_then_state_table_arrival_lower_count_dominates_equal_score'):
            raise ValueError('Fixed whole-contig coverage, scores or tie policy differs')


def check_native_receipt(row,method,species,seqid,length):
    if (not row['completed'] or (row['method'],row['species'],row['seqid'],row['length']) != (method,species,seqid,length) or
            row['input_mode'] != 'fixed_DEV_raw_FASTA' or row['raw_GFF'] != 'prediction.gff3' or
            row['case_policy'] != 'uppercase_same_as_M28_raw_caller' or not row['native_inference_executed'] or
            row['optimizer_steps'] != 0 or any(row[k] for k in ('reference_used','labels_or_feature_cache_used','test_setaria_used')) or
            row['GPU_runtime']['visible_gpus'] != 1):
        raise ValueError('Native raw-FASTA run identity differs')
    asset = row['baseline_asset_receipt']
    if method == 'ANNEVO':
        valid = (asset['source_commit'] == '42c920f70f52f96c7c06c665a639ec47d9304fff' and
                 asset['version'] == '2.3.3' and asset['lineage'] == 'Magnoliopsida' and
                 asset['weight'] == 'assets/ANNEVO/saved_model/ANNEVO_Magnoliopsida.pt')
    elif method == 'Tiberius':
        valid = (asset['native_metadata']['source_commit'] == '871d4979568ae09e56de7eca5d7489046ee833ea' and
                 asset['native_metadata']['version'] == '2.0.8' and
                 asset['base_container_actual_metadata']['version'] == '2.0.7' and
                 asset['configuration_commit'] == '871d4979568ae09e56de7eca5d7489046ee833ea' and
                 asset['entry_point'] == 'python3 -m tiberius.main' and
                 asset['weights_directory'] == str((N.ASSET_PREP/'assets/Tiberius-model-weights/angiosperms').relative_to(I.ROOT)))
    elif method == 'Helixer':
        valid = (asset['native_metadata']['version'] == '0.3.7' and
                 asset['native_metadata']['source_commit'] == 'd17bb496b8842542590644b2437798711c0125d5' and
                 asset['model_source_commit'] == 'd17bb496b8842542590644b2437798711c0125d5' and
                 asset['selected_model']['model_file_name'] == 'land_plant_v0.3_a_0080.h5' and
                 int(asset['selected_model']['priority']) == 80)
    else:
        raise ValueError('Unknown fixed comparator')
    if not valid:
        raise ValueError('Frozen native version, weight or lineage differs: '+method)


def check_terminal_rows(rows,job_ids):
    actual = {r[0]:(r[1],r[2]) for r in rows}
    for job in job_ids:
        for key in (job,job+'.batch',job+'.extern'):
            if actual.get(key) != ('COMPLETED','0:0'):
                raise ValueError('Missing or unsuccessful Slurm terminal: '+key)


def structure_metrics(txs,refs):
    predicted,reference = ({C.S.key(t) for t in rows} for rows in (txs,refs))
    result = {}
    for label,predicate in (('single_CDS_segment',lambda key: len(key[2]) == 1),
                            ('multiple_CDS_segments',lambda key: len(key[2]) > 1),
                            ('span_longer_than_native_window',lambda key: key[2][-1][1]-key[2][0][0] > C.M.WIDTH)):
        p = {key for key in predicted if predicate(key)}
        r = {key for key in reference if predicate(key)}
        result[label] = C.S.prf(len(p & r),len(p),len(r))
    return result


def reconcile_inputs(infer_root=I.OUT,checkpoint=I.CHECKPOINT,identity=I.R.IDENTITY,protocol=I.R.PROTOCOL,
                     *,candidate_validator=None):
    if tuple(N.SCOPE) != I.SCOPE:
        raise ValueError('Prepared DEV and fixed native scopes differ')
    candidates,baselines,jobs = {},{},set()
    for task,(species,seqid,length) in enumerate(I.SCOPE):
        folder = infer_root/seqid
        if (folder/'STATUS').read_text().strip() != 'COMPLETED':
            raise ValueError('Complete both candidate predictions before reference access')
        outer = json.loads((folder/'summary.json').read_text())
        raw = json.loads((folder/'prediction/summary.json').read_text())
        validator = check_candidate if candidate_validator is None else candidate_validator
        validator(outer,raw,task,species,seqid,length,checkpoint=checkpoint,identity=identity,protocol=protocol)
        if not (folder/'prediction/full_grammar.gff3').stat().st_size:
            raise ValueError('Candidate GFF missing')
        candidates[species] = (outer,raw)
        jobs.add(outer['job_id'])
        for method in N.METHODS:
            folder = N.OUT/method/species
            if (folder/'STATUS').read_text().strip() != 'COMPLETED':
                raise ValueError('Native prediction incomplete')
            row = json.loads((folder/'summary.json').read_text())
            check_native_receipt(row,method,species,seqid,length)
            commands = json.loads((folder/'commands.json').read_text())
            expected,_ = N.commands(method,species,folder,commands['CUDA_VISIBLE_DEVICES'],8)
            # ANNEVO's saved interpreter is its original annevo environment,
            # not this CPU evaluator's generanno interpreter; compare its arguments.
            if commands['native'][1 if method == 'ANNEVO' else 0:] != expected[1 if method == 'ANNEVO' else 0:]:
                raise ValueError('Frozen native operation point differs')
            if N.check_output(folder/'prediction.gff3',folder/'command.log',method,seqid,length) != row['output_check']:
                raise ValueError('Saved native completion/GFF differs from receipt')
            baselines[(method,species)] = row
            jobs.add((folder/'JOBID').read_text().strip())
    if any(not job.isdigit() for job in jobs):
        raise ValueError('Numeric actual Slurm job IDs required')
    result = subprocess.run(['sacct','-j',','.join(sorted(jobs)),'-n','-P','-o','JobIDRaw,State,ExitCode'],
                            check=True,text=True,stdout=subprocess.PIPE)
    check_terminal_rows([line.split('|') for line in result.stdout.splitlines() if line.strip()],jobs)
    sequences = {}
    for species,seqid,length in I.SCOPE:
        candidate = C.M.read_allowed(infer_root/seqid/'input.fasta',{seqid:length})
        native = C.M.read_allowed(N.ASSET_PREP/'inputs'/(species+'.fa'),{seqid:length})
        original = C.M.read_allowed(I.ROOT/'data/m1_screen'/species/'genome.fa',{seqid:length})
        if candidate != native or native != original:
            raise ValueError('Actual original/candidate/native FASTA records differ; do not score')
        sequences[species] = original
    return candidates,baselines,sequences


def main(infer_root=I.OUT,checkpoint=I.CHECKPOINT,identity=I.R.IDENTITY,protocol=I.R.PROTOCOL,
         method_label=METHODS[0],out=OUT,*,candidate_validator=None):
    began = time.perf_counter()
    if candidate_validator is None:
        candidates,baselines,sequences = reconcile_inputs(infer_root,checkpoint,identity,protocol)
    else:
        candidates,baselines,sequences = reconcile_inputs(infer_root,checkpoint,identity,protocol,
            candidate_validator=candidate_validator)
    methods = (method_label,*N.METHODS)
    report = dict(completed=True,job_id=os.environ['SLURM_JOB_ID'],
        status='fixed_two_DEV_development_comparison_not_independent_generalization',
        checkpoint=str(checkpoint),reference_policy='M26_CDS_assessable_longest_eligible_CDS_tie_transcript_ID',
        DEV_used=True,sealed_data_used=False,optimizer_steps_this_run=0,
        baseline_predictions_reused_without_new_native_runs=True,
        unsupported_reference_strata_retained=True,input_records_directly_equal=True,
        species={},pooled={},differences_percentage_points={},
        limitations='historically_reused_DEV_and_method_dependent_or_unknown_pretrained_exposure_no_clean_zero_shot_or_publication_PASS')
    pooled = {method:[] for method in methods}
    for species,seqid,length in I.SCOPE:
        sequence = sequences[species]
        data = I.ROOT/'data/m1_screen'/species
        coding = C.S.E.parse_annotation(data/'reference.gff3',{seqid:length},protein_coding_only=True)
        annotation = C.S.E.parse_annotation(data/'reference.gff3',{seqid:length})
        with (data/'reference.gff3').open() as handle:
            _,_,raw_cds,_ = C.P.raw_features(handle,seqid)
        refs,isoforms = C.select_reference(coding,raw_cds,sequence)
        if len(refs) != EXPECTED[species] or len({C.S.key(t) for t in refs}) != EXPECTED[species]:
            raise ValueError('Frozen primary denominator differs; do not filter references')
        if any(not C.R7.phase_chain_consistent(t) for t in refs):
            raise ValueError('Reference phase differs; do not filter references')
        overlaps = C.M.C.overlapping_ids(refs); strata = C.Counter()
        for tx in refs:
            reasons,_ = C.M.chain_checks(tx,sequence[seqid])
            if tx['id'] in overlaps: reasons.append('same_strand_primary_overlap')
            if len(tx['CDS']) > 128: reasons.append('over_128_CDS_segments')
            if tx['CDS'][-1][1]-tx['CDS'][0][0] > C.M.WIDTH: reasons.append('span_longer_than_native_window_descriptive_not_grammar_exclusion')
            strata.update(set(reasons))
        item = dict(seqid=seqid,length=length,reference_primary_chains=len(refs),retained_reference_strata=dict(strata),methods={})
        for method in methods:
            path = (infer_root/seqid/'prediction/full_grammar.gff3' if method == methods[0]
                    else N.OUT/method/species/'prediction.gff3')
            txs = C.S.E.primary_transcripts(C.S.E.parse_annotation(path,{seqid:length}),complete_only=False)
            if method == methods[0]:
                if len(txs) != candidates[species][1]['primary_CDS_chain_predictions']:
                    raise ValueError('Candidate GFF roundtrip count differs')
                normalized = re.sub('[^ACGT]','N',sequence[seqid].upper())
                for tx in txs:
                    reasons,_ = C.M.chain_checks(tx,normalized)
                    if reasons or len(tx['CDS']) > 128:
                        raise ValueError('Candidate actual DNA/phase/cap differs')
                for strand in ('+','-'):
                    ordered = sorted((t['CDS'][0][0],t['CDS'][-1][1]) for t in txs if t['strand'] == strand)
                    if any(b > c for (_,b),(c,_) in zip(ordered,ordered[1:])):
                        raise ValueError('Candidate selected same-strand genes overlap')
            elif baselines[(method,species)]['output_check']['feature_counts'].get('CDS',0) and not txs:
                raise ValueError('Native CDS missing from common parser')
            measured = C.score(txs,refs,isoforms,length,annotation)
            measured['prediction_nonexact_primary'] = measured.pop('reference_nonexact_primary')
            measured['phase_chain_consistent_prediction_records'] = sum(C.R7.phase_chain_consistent(t) for t in txs)
            measured['descriptive_structure_metrics_not_reference_filters'] = structure_metrics(txs,refs)
            measured['raw_prediction_path'] = str(path.relative_to(I.ROOT))
            measured['runtime_receipt'] = (candidates[species][0] if method == methods[0] else baselines[(method,species)])
            item['methods'][method] = measured; pooled[method].append(measured)
        report['species'][species] = item
        report['differences_percentage_points'][species] = {
            method:100*(item['methods'][methods[0]]['exact_chain']['f1']-item['methods'][method]['exact_chain']['f1'])
            for method in N.METHODS}
    for method,rows in pooled.items():
        result = aggregate(rows)
        f1s = [r['exact_chain']['f1'] for r in rows]
        result.update(macro_species_exact_chain_F1=sum(f1s)/len(f1s),worst_species_exact_chain_F1=min(f1s))
        report['pooled'][method] = result
    report['evaluation_seconds'] = time.perf_counter()-began
    with (out/'summary.json').open('x') as handle:
        json.dump(report,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps(dict(completed=True,methods=methods,reference_primary_chains=7728)),flush=True)


if __name__ == '__main__':
    main()
