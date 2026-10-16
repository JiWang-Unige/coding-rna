#!/usr/bin/env python3
"""Fixed original TRAIN primaries: completion witness audit, not a new caller."""
import json
import os
import time

import torch

import run_query as Q
from src.m28.candidates import Budget, ORFIndex, STOP
from src.m28.labels import build_targets, oriented_parts
from src.m28.full_grammar_observed_partition import allowed_roles_from_local_targets
from src.m28.phase_evidence_chain import candidate_roles

ROOT = Q.ROOT
OUT = ROOT / 'outputs/M28-FULL-GRAMMAR-R4-COMPLETION-WITNESS'
SOURCE = ROOT / 'outputs/M28-FULL-GRAMMAR-R4-SUPPORT-PROBE'
SOURCE_JOB = '13427583'
PROTOCOL = 'reports/M28-METHOD-RESTART/full_grammar_partition_preparation.md'
CONFIG = dict(scope='same_original_eight_fixed_TRAIN_witnesses_all_original_primaries',
    context='original_TRAIN_cached_DNA_only_no_unindexed_FASTA_scan',
    bounds='union_of_core_and_original_primary_CDS_spans_diagnostic_only_not_inference',
    missing_context='record_unverified_never_fabricate_DNA_or_delete_observations',
    targets='unchanged_original_R4_fine_coarse_unknown', max_exons=128,
    min_intron_bp=5, optimizer_steps=0, partition_or_NLL_computed=False)


def load_cache(w, index, opened):
    artifact = index[w['id']]['artifact']
    cached = torch.load(ROOT / artifact, map_location='cpu', weights_only=True)
    Q.validate_cached(cached, w)
    opened.add(artifact)
    return cached


def missing_intervals(bases):
    spans, start = [], None
    for i, base in enumerate(bases + ['end']):
        if base is None and start is None:
            start = i
        elif base is not None and start is not None:
            spans.append([start, i]); start = None
    return spans


def cached_context(w, core, a, b, windows, index, opened):
    """Join only original TRAIN cache intervals; compare every shared base."""
    bases = [None] * (b - a)
    fragments = []

    def insert(other, sequence):
        genomic = sequence if other['strand'] == '+' else Q.M.rc(sequence)
        left, right = max(a, other['start']), min(b, other['end'])
        shared = 0
        for x in range(left, right):
            base = genomic[x - other['start']]
            if bases[x - a] is not None:
                if bases[x - a] != base:
                    raise ValueError('Original TRAIN cache DNA disagrees at a shared genomic base')
                shared += 1
            bases[x - a] = base
        fragments.append(dict(window_id=other['id'], artifact=index[other['id']]['artifact'],
            genomic_interval_read=[left, right], strand=other['strand'], shared_bases_compared=shared))

    insert(w, core)
    neighbors = sorted((other for other in windows.values()
        if other['id'] != w['id'] and other['species'] == w['species']
        and other['seqid'] == w['seqid'] and other['start'] < b and other['end'] > a),
        key=lambda other: (other['start'], other['id']))
    for other in neighbors:
        left, right = max(a, other['start']), min(b, other['end'])
        if not any(base is None for base in bases[left-a:right-a]):
            continue
        cached = load_cache(other, index, opened)
        insert(other, cached['sequence'])
        del cached
    missing = missing_intervals(bases)
    if missing:
        return None, fragments, [[a+x, a+y] for x, y in missing]
    genomic = ''.join(bases)
    sequence = genomic if w['strand'] == '+' else Q.M.rc(genomic)
    offset = w['start'] - a if w['strand'] == '+' else b - w['end']
    if sequence[offset:offset + len(core)] != core:
        raise ValueError('Completed cache context does not preserve the original oriented core DNA')
    return sequence, fragments, []


def geometry_failures(parts, length):
    failures = []
    if not 1 <= len(parts) <= 128:
        failures.append('exon_cap')
    if any(not 0 <= x < y <= length for x, y in parts):
        failures.append('exon_bounds')
    if parts[0][1] - parts[0][0] < 3 or parts[-1][1] - parts[-1][0] < 3:
        failures.append('split_terminal_codon')
    if any(z-y < Budget().min_intron_bp for (_, y), (z, _) in zip(parts, parts[1:])):
        failures.append('intron_below_5bp')
    return failures


def DNA_checks(sequence, parts):
    """Independent concatenated CDS verdict versus actual-carry ORFIndex."""
    coding = ''.join(sequence[x:y] for x, y in parts)
    failures = []
    if len(coding) % 3:
        failures.append('length_not_mod3')
    if coding[:3] != 'ATG':
        failures.append('non_ATG_start')
    if coding[-3:] not in STOP:
        failures.append('no_terminal_stop')
    if any(coding[i:i+3] in STOP for i in range(0, max(0, len(coding)-3), 3)):
        failures.append('internal_stop')
    if set(coding) - set('ACGT'):
        failures.append('ambiguous_CDS')
    orf, carry = ORFIndex(sequence), ''
    actual = coding[:3] == 'ATG'
    for ordinal, (x, y) in enumerate(parts):
        carry = orf.extend(x, y, carry, ordinal == len(parts) - 1)
        if carry is None:
            actual = False; break
    if actual != (not failures):
        raise ValueError('Actual-carry and direct concatenated-CDS grammar verdicts disagree')
    for (_, y), (z, _) in zip(parts, parts[1:]):
        if sequence[y:y+2] not in ('GT', 'GC') or sequence[z-2:z] != 'AG':
            failures.append('nonstandard_splice')
    return failures, actual


def main():
    if ((SOURCE / 'STATUS').read_text().strip() != 'COMPLETED'
            or (SOURCE / 'JOBID').read_text().strip() != SOURCE_JOB):
        raise ValueError('Matching closed original R4 observation support diagnostic required')
    source = json.loads((SOURCE / 'summary.json').read_text())
    if (not source['completed'] or source['job_id'] != SOURCE_JOB
            or source['actual_witnesses'] != 8 or source['optimizer_steps_this_run'] != 0
            or not source['all_original_witnesses_checked']
            or source['original_R4_policy_changed'] or source['DEV_used']):
        raise ValueError('Original unchanged eight-witness observation source required')
    draws, index, windows, catalogs = Q.load_train_inputs()
    if len(draws) != 1536 or len(windows) != 1491:
        raise ValueError('Original TRAIN population identity changed')
    witnesses = source['records']; ids = [r['window_id'] for r in witnesses]
    if len(set(ids)) != 8 or [r['ordinal'] for r in witnesses] != list(range(8)):
        raise ValueError('Original eight distinct ordered witness identities required')
    opened, records = set(), []
    began = time.perf_counter()
    with (OUT / 'started.json').open('x') as handle:
        json.dump(dict(job_id=os.environ['SLURM_JOB_ID'], source_job=SOURCE_JOB,
            configuration=CONFIG, witness_ids=ids,
            TRAIN_allowlists={s: sorted(v) for s, v in Q.ALLOW.items()}), handle, indent=2)
    with torch.inference_mode(), (OUT / 'windows.jsonl').open('x') as handle:
        for ordinal, original in enumerate(witnesses):
            w = windows[original['window_id']]; n = w['valid_bases']
            cached = load_cache(w, index, opened)
            core = cached['sequence']
            targets = build_targets(w, catalogs[w['species']], core)
            allowed = allowed_roles_from_local_targets(targets, n)
            fine = int((torch.as_tensor(targets['seg'][:n]) >= 0).sum())
            free = int(allowed.all(dim=1).sum()); coarse = n - fine - free
            if (fine, coarse, free) != (original['fine_positions'], original['coarse_positions'],
                                       original['fully_unknown_positions']):
                raise ValueError('Original fine/coarse/unknown role counts changed')
            primary = [catalogs[w['species']][tid] for tid in w['overlapping_transcript_ids']
                       if catalogs[w['species']][tid]['training_primary']]
            if len(primary) != original['original_overlapping_primaries']:
                raise ValueError('Original full primary set identity changed')
            a = min([w['start']] + [tx['CDS'][0][0] for tx in primary])
            b = max([w['end']] + [tx['CDS'][-1][1] for tx in primary])
            extended = dict(start=a, end=b, strand=w['strand'])
            offset = w['start'] - a if w['strand'] == '+' else b - w['end']
            context, fragments, missing = cached_context(w, core, a, b, windows, index, opened)
            chains = [(tx, oriented_parts(tx, extended)) for tx in primary]
            ordered = sorted(chains, key=lambda pair: (pair[1][0][0], pair[1][-1][1], pair[0]['id']))
            overlaps = [dict(left=tx['id'], right=other['id'],
                oriented_overlap=[max(parts[0][0], other_parts[0][0]), min(parts[-1][1], other_parts[-1][1])])
                for j, (tx, parts) in enumerate(ordered) for other, other_parts in ordered[j+1:]
                if other_parts[0][0] < parts[-1][1]]
            roles = torch.zeros(b-a, dtype=torch.long)
            chain_records, renderable = [], True
            for tx, parts in chains:
                geom = geometry_failures(parts, b-a)
                original_parts = oriented_parts(tx, w)
                if context is not None:
                    DNA_source, DNA, DNA_parts = 'completed_original_TRAIN_cache_context', context, parts
                elif original_parts[0][0] >= 0 and original_parts[-1][1] <= n:
                    DNA_source, DNA, DNA_parts = 'original_core_cache', core, original_parts
                else:
                    DNA_source, DNA, DNA_parts = 'external_DNA_unavailable', None, None
                failures, actual = DNA_checks(DNA, DNA_parts) if DNA is not None else (None, None)
                can_render = (not any(reason in geom for reason in ('exon_bounds', 'split_terminal_codon'))
                    and all(z >= y for (_, y), (z, _) in zip(parts, parts[1:])))
                renderable &= can_render
                mismatches = None
                if can_render:
                    current = candidate_roles(parts, b-a)
                    lo, hi = max(0, parts[0][0]-offset), min(n, parts[-1][1]-offset)
                    positions = torch.arange(lo, hi)
                    mismatches = int((~allowed[positions, current[positions+offset]]).sum())
                    if not overlaps:
                        roles[parts[0][0]:parts[-1][1]] = current[parts[0][0]:parts[-1][1]]
                chain_records.append(dict(transcript_id=tx['id'], gene_id=tx['gene_id'],
                    oriented_span=[parts[0][0], parts[-1][1]], exons=len(parts),
                    historical_CDS_complete_supported=tx['CDS_complete_supported'],
                    historical_sequence_structure_failures=tx['sequence_structure_failures'],
                    geometry_failures=geom, DNA_source=DNA_source, direct_DNA_failures=failures,
                    actual_carry_coding_verdict=actual,
                    actual_full_grammar_legal=(not geom and not failures) if DNA is not None else None,
                    individual_core_role_mismatches=mismatches))
            mismatch_positions = None
            if not overlaps and renderable:
                mismatch_positions = torch.where(~allowed[torch.arange(n), roles[offset:offset+n]])[0].tolist()
            if overlaps:
                status = 'original_primary_set_span_overlap'
            elif any(r['geometry_failures'] or r['actual_full_grammar_legal'] is False for r in chain_records):
                status = 'original_primary_chain_grammar_rejection'
            elif mismatch_positions is None or mismatch_positions:
                status = 'original_core_observation_mismatch'
            elif missing:
                status = 'external_DNA_unavailable_not_a_verified_completion'
            else:
                status = 'DNA_verified_completion_witness'
            row = dict(ordinal=ordinal, window_id=w['id'], species=w['species'], category=original['category'],
                original_closed_window_observation_feasible=original['observation_feasible'],
                fine_positions=fine, coarse_positions=coarse, fully_unknown_positions=free,
                original_primary_transcript_ids=[tx['id'] for tx in primary],
                diagnostic_genomic_context=[a,b], core_oriented_offset=offset,
                cached_fragments=fragments, missing_genomic_DNA_intervals=missing,
                primary_span_overlaps=overlaps, primary_chains=chain_records,
                set_core_role_mismatches=None if mismatch_positions is None else len(mismatch_positions),
                mismatch_core_position_sample=None if mismatch_positions is None else mismatch_positions[:12],
                completion_witness_status=status, DNA_verified_completion_witness=status == 'DNA_verified_completion_witness',
                historical_metadata_geometry_and_core_compatible=(not overlaps and mismatch_positions == []
                    and all(r['historical_CDS_complete_supported'] and not r['geometry_failures'] for r in chain_records)))
            records.append(row); handle.write(json.dumps(row, allow_nan=False)+'\n'); handle.flush()
            print(json.dumps(dict(ordinal=ordinal, window_id=w['id'], status=status,
                primary_span_overlaps=len(overlaps), missing_DNA_intervals=len(missing),
                set_core_role_mismatches=row['set_core_role_mismatches'])), flush=True)
            del cached, targets, allowed
    if len(records) != 8 or [r['window_id'] for r in records] != ids:
        raise ValueError('Original eight witness audit incomplete or reordered')
    summary = dict(completed=True, job_id=os.environ['SLURM_JOB_ID'], source_job=SOURCE_JOB,
        protocol=PROTOCOL, configuration=CONFIG, original_draws=1536, original_train_unique_windows=1491,
        actual_witnesses=8, records=records, seconds_after_original_TRAIN_metadata_load=time.perf_counter()-began,
        cached_artifacts_opened=sorted(opened), DNA_verified_completion_witnesses=sum(r['DNA_verified_completion_witness'] for r in records),
        optimizer_steps_this_run=0, model_constructed=False, optimizer_constructed=False,
        model_or_backbone_forward=False, research_weights_read=False, checkpoints_written=False,
        original_R4_policy_changed=False, observation_rows_deleted=False, primary_transcripts_dropped=False,
        external_DNA_fabricated=False, unindexed_FASTA_scanned=False, raw_FASTA_opened=False,
        partition_or_NLL_computed=False, training_or_inference_caller_changed=False,
        DEV_used=False, independent_species_used=False, A33_outputs_opened=False,
        interpretation='fixed_primary_completion_audit_not_general_support_verdict_padding_rule_fit_or_performance')
    with (OUT / 'summary.json').open('x') as handle:
        json.dump(summary, handle, indent=2, allow_nan=False)


if __name__ == '__main__':
    torch.set_num_threads(2)
    main()
