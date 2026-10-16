#!/usr/bin/env python3
"""Terminal A10 GFF/member/numerical receipt, not new inference or selection."""
import csv
import json
import math
from collections import Counter
import diagnose_query_context as D
import verify_completed_a6_a7_a11 as V

ROOT = D.ROOT
OUT = ROOT/'outputs/M28-A10-CONTEXT-RESULT-CHECK'


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    if OUT.exists():
        raise FileExistsError('Keep prior receipt intact')
    evaluated = D.E.read_report(D.E.OUT/'summary.json')
    diagnostic = D.E.read_report(D.OUT/'summary.json')
    paired = table(D.OUT/'references.tsv')
    native = table(D.BASELINE/'references.tsv')
    assert len(paired) == len(native) == 7728
    assert len({V.ref_key(r) for r in paired}) == 7728
    assert {V.ref_key(r) for r in paired} == {V.ref_key(r) for r in native}
    original = {V.ref_key(r): r for r in native}
    records = []
    for raw in paired:
        prior = original[V.ref_key(raw)]
        for field in prior:
            assert raw[field] == prior[field]
        row = D.parse_baseline(raw)
        for stage in D.STAGES:
            value = raw['A10_'+stage]
            assert value in ('True', 'False')
            row['A10_'+stage] = value == 'True'
        row['A10_state'] = raw['A10_state']
        records.append(row)
    assert D.summarize(records) == diagnostic['pooled']
    lengths = {seqid: length for _, seqid, length, _ in V.F.A.E.S.SCOPE}
    parsed = V.F.A.E.S.E.primary_transcripts(V.F.A.E.S.E.parse_annotation(
        D.E.I.OUT/'A10_CONTEXT.gff3', lengths), complete_only=False)
    geometry = {V.F.A.E.S.key(t) for t in parsed}
    assert len(geometry) == len(parsed)
    assert all(V.F.A.E.R7.phase_chain_consistent(t) for t in parsed)
    metrics, turnover = {}, {}
    current = evaluated['methods']['A10_CONTEXT']
    for species, count in D.N.COUNTS.items():
        part = [r for r in records if r['species'] == species]
        assert len(part) == count
        assert D.summarize(part) == diagnostic['per_species'][species]
        states = Counter()
        for r in part:
            actual = (r['seqid'], r['strand'], V.chain(r['CDS'])) in geometry
            assert r['A10_selected_unique'] == actual
            old = r['A3_selected_unique']
            states['both' if old and actual else 'A10_only' if actual else 'A3_only' if old else 'neither'] += 1
        value = current['per_species'][species]['exact_chain']
        tp = states['both'] + states['A10_only']
        predicted = sum(t['seqid'] == part[0]['seqid'] for t in parsed)
        assert (value['tp'], value['predicted'], value['reference']) == (tp, predicted, count)
        for field, fraction in [('precision', tp/predicted), ('recall', tp/count),
                                ('f1', 2*tp/(predicted+count))]:
            assert math.isclose(value[field], fraction, rel_tol=0, abs_tol=1e-14)
        metrics[species] = value
        turnover[species] = {k: states[k] for k in ('both', 'A10_only', 'A3_only', 'neither')}
        assert turnover[species] == diagnostic['per_species'][species]['stages']['selected_unique']
    pooled = current['pooled']['exact_chain']
    assert (pooled['tp'], pooled['predicted'], pooled['reference']) == (
        sum(v['tp'] for v in metrics.values()), len(parsed), 7728)
    for field, fraction in [('precision', pooled['tp']/len(parsed)), ('recall', pooled['tp']/7728),
                            ('f1', 2*pooled['tp']/(len(parsed)+7728))]:
        assert math.isclose(pooled[field], fraction, rel_tol=0, abs_tol=1e-14)
    fixed_a = V.read('M28-A-EVALUATION')['methods']['A']['per_species']
    control = json.loads((ROOT/'reports/M28-METHOD-RESTART/r5_result.json').read_text())['methods']['C0']['per_species']
    gate = D.E.A3.continuation_gate(current['per_species'], fixed_a, control)
    assert gate == evaluated['predeclared_gate'] == diagnostic['predeclared_gate_unchanged']
    receipt = dict(completed=True, reference_primary_chains=7728,
        accepted_reference_geometry_unsupported_and_A3_native_fields_unchanged=True,
        GFF_predictions=len(parsed), GFF_duplicates=0, GFF_all_phase_consistent=True,
        actual_GFF_membership_matches_all_7728_rows=True, diagnostic_pooled_and_species_recomputed=True,
        per_species_metrics=metrics, pooled_metrics=pooled, selected_turnover=turnover,
        original_frozen_gate=gate, optimizer_steps_this_run=0, new_model_inference=False,
        test_setaria_used=False)
    OUT.mkdir()
    with (OUT/'summary.json').open('x') as handle:
        json.dump(receipt, handle, indent=2, allow_nan=False)
    print(json.dumps(receipt, allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
