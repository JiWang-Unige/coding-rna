#!/usr/bin/env python3
"""Closed A37 TRAIN transitions from saved scores only, never model rescue."""
from collections import Counter, defaultdict
from itertools import zip_longest
import infer_context_chain as I
from diagnose_structured_emission_attrition import transitions

D, R = I.D, I.D.R
OUT = I.ROOT/'outputs/M28-A37-CONTEXT-CHAIN-ATTRITION'


def main():
    D.F.A.require_jobs([D.FIT_JOB,I.READOUT_JOB,I.CHECK_JOB])
    fitted = R.completed(D.F.OUT,D.FIT_JOB); D.validate_fit(fitted)
    report = R.completed(D.OUT,I.READOUT_JOB)
    checked = R.completed(D.CHECK,I.CHECK_JOB)
    baseline = R.completed(R.READOUT,13398521)
    if (report['configuration'] != D.CONFIG or checked['configuration'] != D.CONFIG or
            report['protocol'] != D.PROTOCOL or checked['protocol'] != D.PROTOCOL or
            checked['source_jobs'] != dict(fit=D.FIT_JOB,readout=I.READOUT_JOB) or
            any(checked[k] is not True for k in I.CPU_PROOFS) or
            checked['selected_counts_per_species'] != report['selected_counts_per_species'] or
            checked['eligible_reference_pairs'] != report['eligible_reference_pairs'] or
            report['TRAIN_preparation_supported'] is not False or checked['TRAIN_preparation_supported'] is not False or
            D.gate(report['selected_counts_per_species']) != report['predeclared_TRAIN_gate'] or
            checked['predeclared_TRAIN_gate'] != report['predeclared_TRAIN_gate'] or
            all(report['predeclared_TRAIN_gate'].values())):
        raise ValueError('Actual closed numerically verified A37 FAIL required, not native DEV readiness')
    saved = R.C.R.rows(D.OUT/'natural_scores.jsonl')
    prior = R.C.R.rows(R.READOUT/'natural_scores.jsonl')
    labelled = R.C.R.rows(D.OUT/'windows.jsonl')
    groups = {s:defaultdict(Counter) for s in report['selected_counts_per_species']}
    counts = {s:{a:Counter() for a in ('A31_A3','A37')} for s in groups}
    windows,candidates = 0,0
    for new,old,lab in zip_longest(saved,prior,labelled):
        if new is None or old is None or lab is None:
            raise ValueError('Complete paired original windows required')
        source,labels = old['variants']['A31_A3'],lab['labels']
        if (new['window'] != old['window'] or new['window'] != lab['window'] or
                new['chains'] != source['chains'] or new['reference_used_in_scoring'] or
                new['window']['split'] != 'train'):
            raise ValueError('Original reference-free same-natural-pool TRAIN identity differs')
        species = new['window']['species']
        for arm,value in (('A31_A3',source),('A37',new)):
            counts[species][arm].update(R.selected_counts(value,labels))
        for key,value in transitions(new['chains'],labels,source,new).items():
            groups[species][key].update(value)
        windows += 1; candidates += len(labels)
    actual = {s:{a:dict(c) for a,c in values.items()} for s,values in counts.items()}
    if (windows != 1491 or candidates != 16321 or actual != report['selected_counts_per_species'] or
            any(actual[s]['A31_A3'] != baseline['selected_counts_per_species'][s]['A31_A3'] for s in actual) or
            sum(v['eligible'] for v in report['eligible_reference_pairs'].values()) != 2916 or
            sum(v['five_or_more'] for v in report['eligible_reference_pairs'].values()) != 1086):
        raise ValueError('Full transition census must reproduce closed counts and original denominators')
    pooled = defaultdict(Counter)
    for species in groups.values():
        for key,value in species.items():
            pooled[key].update(value)
    R.C.R.write_summary(OUT,dict(completed=True,experiment='M28-A37-CONTEXT-CHAIN-ATTRITION',
        source_jobs=dict(fit=D.FIT_JOB,readout=I.READOUT_JOB,CPU=I.CHECK_JOB,baseline_readout='13398521'),
        scientific_status='A37_TRAIN_FAIL_closed_no_DEV_or_rescue',
        grain='original_unique_TRAIN_window_natural_chain_pair_not_unique_gene',
        windows=windows,natural_candidates=candidates,
        per_species={s:{k:dict(c) for k,c in v.items()} for s,v in groups.items()},
        pooled={k:dict(c) for k,c in pooled.items()},selected_counts_per_species=actual,
        eligible_reference_pairs=report['eligible_reference_pairs'],
        predeclared_TRAIN_gate=report['predeclared_TRAIN_gate'],
        transition_logic='unchanged_diagnose_structured_emission_attrition.transitions',
        model_or_backbone_forward=False,model_or_feature_tensor_opened=False,optimizer_steps_this_run=0,
        new_candidate_search=False,new_R4_labelling=False,TRAIN_catalog_or_reference_opened=False,
        DEV_used=False,test_setaria_used=False,A33_outputs_opened=False,
        limitation='Descriptive paired transitions of already verified fixed scores and selections; '
            'not new training, causal component ablation, threshold rescue, gene F1 or generalization.'))


if __name__ == '__main__':
    main()
