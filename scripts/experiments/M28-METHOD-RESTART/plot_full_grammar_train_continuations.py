#!/usr/bin/env python3
"""Static paper figures from completed original-TRAIN summaries only.

Chart contract: compare complete-CDS-chain F1 across all completed branches
and fixed native tools (12 methods x three scopes), and distinguish gains
from losses against STRUCTURED_CONTINUE final9216 (seven x three scopes).
Use prespecified method order, identical axes within each figure, count-based
pooled F1, explicit reference denominators, direct labels and no uncertainty.
PNG/SVG scientific article surface: blue model root, neutral native tools;
gain blue/fill versus loss orange/hatch plus signed axis, not color alone.
Keep negative results. TRAIN fitting and paired recall accounting do not
establish independent generalization, causal mechanism or model selection.
No provider branding, hosting, new inference, fitting or DEV data access.
"""
import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path

from plot_initial_train_benchmark import CONTIGS, REFERENCES, SPECIES, recompute
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[3]
SOURCES = {
    'family': ('M28-FULL-GRAMMAR-CONTINUATION-TRAIN-FAMILY', '13493626'),
    'four_cell': ('M28-FULL-GRAMMAR-CONTEXT-PRIOR-TRAIN-COMPARISON', '13460976'),
    'recovery': ('M28-FULL-GRAMMAR-TRAIN-RECOVERY', '13461885'),
    'core': ('M28-FULL-GRAMMAR-COMPLETE-CORE-RECOVERY', '13463353'),
    'native': ('M28-NATIVE-TRAIN-BASELINE-EVAL', '13461307'),
}
METHODS = (
    'M28_final4608', 'M28_STRUCTURED_CONTINUE_final9216',
    'M28_LOCAL_CHANNEL_final9216', 'M28_HALF_STRUCTURED_HALF_LOCAL_final9216',
    'M28_NATIVE_CDS_PRIOR_final9216', 'M28_FRESH_DRAWS_final9216',
    'M28_JOINT24576_final9216', 'M28_JOINT24576_NATIVE_CDS_PRIOR_final9216',
    'M28_COMPLETE_CORE_DRAWS_final9216', 'ANNEVO', 'Tiberius', 'Helixer',
)
LABELS = ('C0 (4,608)', 'Structured (9,216)', 'Local (9,216)', 'Mixed (9,216)',
          'Native CDS prior', 'Fresh uniform', 'Joint 24,576 bp', 'Joint + prior',
          'Complete-core draws', 'ANNEVO', 'Tiberius', 'Helixer')
CONTROL, COMBINED, CORE = METHODS[1], METHODS[7], METHODS[8]
SCOPES = (*SPECIES, 'pooled')
TITLES = ('A  Arabidopsis thaliana', 'B  Oryza sativa', 'C  Pooled six contigs')
BLUE, ORANGE = '#2F6B9A', '#B45F28'


def scope_cell(model, scope):
    return model['pooled'] if scope == 'pooled' else model['per_species'][scope]


def load_data():
    source = {}
    for name, (root, job) in SOURCES.items():
        directory = ROOT/'outputs'/root
        report = json.loads((directory/'summary.json').read_text())
        if ((directory/'STATUS').read_text().strip() != 'COMPLETED'
                or (directory/'JOBID').read_text().strip() != job
                or report['completed'] is not True or report['job_id'] != job
                or report['optimizer_steps_this_run'] != 0
                or report['DEV_used'] or report['sealed_data_used']):
            raise ValueError('Expected completed, zero-update TRAIN source: '+name)
        source[name] = report
    family, four, recovery, core, native = (source[k] for k in SOURCES)
    if (tuple(family['methods']) != METHODS[:7]
            or set(native['methods']) != set(METHODS[9:])
            or recovery['reference_policy'] != native['reference_policy']
            or core['reference_policy'] != native['reference_policy']
            or not native['unsupported_reference_strata_retained']
            or not native['training_eligibility_does_not_filter_reference']):
        raise ValueError('Fixed methods and retained reference policy required')
    models = {**family['measurements'], COMBINED: four['measurements']['F11'],
              CORE: core['measurements'][CORE], **native['methods']}
    jobs = {**family['source_evaluation_jobs'],
            COMBINED: four['source_evaluation_jobs']['F11'],
            CORE: core['source_evaluation_jobs'][CORE],
            **{method: native['job_id'] for method in METHODS[9:]}}
    rows = []
    anchor = models[METHODS[0]]
    for method, label in zip(METHODS, LABELS):
        model = models[method]
        if (set(model['chromosomes']) != {c for cs in CONTIGS.values() for c in cs}
                or set(model['per_species']) != set(SPECIES)):
            raise ValueError('Original two species and six TRAIN contigs required')
        if method not in METHODS[9:]:
            steps = 4608 if method == METHODS[0] else 9216
            if (model['cumulative_optimizer_steps'] != steps
                    or not model['checkpoint'].endswith(f'/step_{steps:06d}.pt')):
                raise ValueError('Use fixed final weights, not selected checkpoints')
        for seqid, cell in model['chromosomes'].items():
            base = anchor['chromosomes'][seqid]
            for field in ('seqid', 'species', 'length', 'reference_primary_records',
                          'reference_assessable_isoform_records', 'retained_unsupported_reference_strata'):
                if cell[field] != base[field]:
                    raise ValueError('Reference metadata differs: '+seqid+' '+field)
            recompute(cell['exact_chain'])
        for scope in SCOPES:
            cell, base = scope_cell(model, scope), scope_cell(anchor, scope)
            counts = recompute(cell['exact_chain'])
            seqids = tuple(c for cs in CONTIGS.values() for c in cs) if scope == 'pooled' else CONTIGS[scope]
            if (counts['reference'] != REFERENCES[scope] or any(
                    counts[k] != sum(model['chromosomes'][c]['exact_chain'][k] for c in seqids)
                    for k in ('tp', 'predicted', 'reference'))
                    or cell['CDS_base_unstranded_union']['reference'] != base['CDS_base_unstranded_union']['reference']
                    or cell['background']['bp'] != base['background']['bp']):
                raise ValueError('Scope count aggregation or denominator differs')
            for counterpart in (recovery['measurements'].get(method), core['measurements'].get(method)):
                if counterpart and scope_cell(counterpart, scope)['exact_chain'] != cell['exact_chain']:
                    raise ValueError('Primary counts differ across completed sources')
            rows.append(dict(method=method, label=label, scope=scope, **counts,
                f1_percent=100*counts['f1'], source_evaluation_job=jobs[method],
                checkpoint=model.get('checkpoint'), cumulative_optimizer_steps=model.get('cumulative_optimizer_steps'),
                CDS_base_F1=cell['CDS_base_unstranded_union']['f1'],
                wholly_background_chains=cell['background']['wholly_background_unique_chains'],
                macro_species_exact_chain_F1=model['macro_species_exact_chain_F1'],
                worst_species_exact_chain_F1=model['worst_species_exact_chain_F1']))
        species_f1 = [scope_cell(model, s)['exact_chain']['f1'] for s in SPECIES]
        if (not math.isclose(model['macro_species_exact_chain_F1'], sum(species_f1)/2, abs_tol=1e-12)
                or not math.isclose(model['worst_species_exact_chain_F1'], min(species_f1), abs_tol=1e-12)):
            raise ValueError('Macro/worst species scores differ')
    pairs = []
    core_pair = core['comparisons'][CONTROL]
    for method in METHODS[2:9]:
        method_pairs = {}
        for scope in SCOPES:
            pair = (core_pair['pooled']['all_references'] if scope == 'pooled'
                    else core_pair['per_species'][scope]['all_references']) if method == CORE else recovery['matched_structured9216_reference_pairs'][scope][method]
            n, retained, gain, loss, missed = (pair[k] for k in
                ('reference_chains', 'retained', 'gained', 'lost', 'still_missed'))
            control_tp = scope_cell(models[CONTROL], scope)['exact_chain']['tp']
            candidate_tp = scope_cell(models[method], scope)['exact_chain']['tp']
            if (min(n, retained, gain, loss, missed) < 0 or n != REFERENCES[scope]
                    or retained+gain+loss+missed != n
                    or retained+loss != pair['control_hits'] or pair['control_hits'] != control_tp
                    or retained+gain != pair['candidate_hits'] or pair['candidate_hits'] != candidate_tp
                    or gain-loss != pair['net_reference_hits']
                    or not math.isclose(pair['reference_recall'], candidate_tp/n, abs_tol=1e-12)):
                raise ValueError('Matched reference-hit partition differs: '+method+' '+scope)
            method_pairs[scope] = pair
            pairs.append(dict(method=method, label=LABELS[METHODS.index(method)],
                scope=scope, control=CONTROL, source_job=core['job_id'] if method == CORE else recovery['job_id'], **pair))
        if any(method_pairs['pooled'][k] != sum(method_pairs[s][k] for s in SPECIES)
               for k in ('reference_chains', 'control_hits', 'candidate_hits', 'retained',
                         'gained', 'lost', 'still_missed', 'net_reference_hits')):
            raise ValueError('Pooled pair counts do not aggregate by species')
    return source, rows, pairs


def style_axis(ax, scope, title, n, labels):
    ax.set_ylim(n-.4, -.6)
    ax.set_yticks(range(n), labels)
    ax.set_title(title, loc='left', fontsize=12, fontweight='bold', pad=29)
    ax.text(0, 1.025, f'Primary references: {REFERENCES[scope]:,}',
            transform=ax.transAxes, fontsize=10, color='#5B626B')
    ax.grid(axis='x', color='#E4E7EB', linewidth=.65, zorder=0)
    ax.tick_params(axis='both', length=0, colors='#424850', pad=7)
    for name in ('top', 'right', 'left'):
        ax.spines[name].set_visible(False)
    ax.spines['bottom'].set_color('#7E8790')
    ax.spines['bottom'].set_linewidth(.65)


def save_figure(fig, out, stem):
    fig.savefig(out/(stem+'.png'), dpi=220, facecolor='white')
    fig.savefig(out/(stem+'.svg'), facecolor='white')
    plt.close(fig)


def main(out):
    source, rows, pairs = load_data()
    if any(out.glob('*.png')) or any(out.glob('*.svg')) or (out/'provenance.json').exists():
        raise FileExistsError('Preserve existing figure outputs')
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'svg.fonttype': 'none'})
    fig, axes = plt.subplots(1, 3, figsize=(15.8, 9.0), sharex=True, sharey=True)
    fig.subplots_adjust(left=.17, right=.985, top=.80, bottom=.22, wspace=.15)
    colors = [BLUE]*9+['#8A949F', '#A7AFB8', '#C4CAD0']
    for ax, scope, title in zip(axes, SCOPES, TITLES):
        for y, (method, color) in enumerate(zip(METHODS, colors)):
            row = next(r for r in rows if r['method'] == method and r['scope'] == scope)
            ax.barh(y, row['f1_percent'], height=.58, color=color, edgecolor='#424850',
                    linewidth=.5, hatch='///' if y == 0 else None, zorder=3)
            ax.text(row['f1_percent']+1.3, y, f"{row['f1_percent']:.1f}", va='center', fontsize=10)
        style_axis(ax, scope, title, len(METHODS), LABELS)
        ax.axhline(8.5, color='#9AA1A8', linewidth=.7, linestyle='--')
        ax.set_xlim(0, 100)
        ax.set_xticks((0, 25, 50, 75, 100))
        ax.set_xlabel('Complete-CDS-chain F1 (%)', labelpad=9)
    fig.suptitle('Complete-CDS-chain accuracy on the original TRAIN panel', x=.17,
                 y=.963, ha='left', fontsize=18, fontweight='bold', color='#242A31')
    fig.text(.17, .907, 'All completed branches and fixed native tools; same six raw FASTA contigs and retained references', fontsize=11)
    fig.text(.17, .144, 'All model continuations: final 9,216 cumulative updates. C0: final 4,608. Order is not a performance ranking.', fontsize=10)
    fig.text(.17, .095, 'Pooled F1 uses pooled counts. Native training budgets and upstream annotation exposure are not matched.', fontsize=10, color='#5B626B')
    fig.text(.17, .047, 'TRAIN diagnosis, not independent generalization, causal attribution or DEV-based selection. No uncertainty interval.', fontsize=10, color='#5B626B')
    save_figure(fig, out, 'train_continuation_complete_chain')

    fig, axes = plt.subplots(1, 3, figsize=(15.8, 7.8), sharex=True, sharey=True)
    fig.subplots_adjust(left=.17, right=.985, top=.74, bottom=.25, wspace=.15)
    for ax, scope, title in zip(axes, SCOPES, TITLES):
        for y, method in enumerate(METHODS[2:9]):
            pair = next(p for p in pairs if p['method'] == method and p['scope'] == scope)
            ax.barh(y, pair['gained'], height=.55, color=BLUE, edgecolor=BLUE, zorder=3)
            ax.barh(y, -pair['lost'], height=.55, color='white', edgecolor=ORANGE, hatch='///', linewidth=.9, zorder=3)
            inside = pair['gained'] >= 2000
            ax.text(pair['gained']-110 if inside else pair['gained']+85, y,
                    f"{pair['gained']:,}", va='center', ha='right' if inside else 'left',
                    color='white' if inside else BLUE, fontsize=9)
            ax.text(-pair['lost']-85, y, f"{pair['lost']:,}", va='center', ha='right', color=ORANGE, fontsize=9)
            ax.text(4850, y, f"{pair['net_reference_hits']:+,}", va='center', ha='right', fontsize=9)
        style_axis(ax, scope, title, 7, LABELS[2:9])
        ax.axvline(0, color='#424850', linewidth=.9, zorder=2)
        ax.set_xlim(-5200, 5000)
        ax.set_xticks((-4000, -2000, 0, 2000, 4000))
        ax.set_xticklabels(('-4,000', '-2,000', '0', '2,000', '4,000'))
        ax.set_xlabel('Reference-hit change (count)', labelpad=9)
        ax.text(4850, -.8, 'Net', ha='right', fontsize=9, color='#424850')
    fig.suptitle('Reference hits gained and lost against structured continuation',
                 x=.17, y=.963, ha='left', fontsize=18, fontweight='bold', color='#242A31')
    fig.text(.17, .900, 'Same final 9,216-update operating points; fixed primary-reference population; no subset precision or F1', fontsize=11)
    fig.legend(handles=[Patch(facecolor=BLUE, edgecolor=BLUE, label='Gained reference hits (right)'),
                        Patch(facecolor='white', edgecolor=ORANGE, hatch='///', label='Lost reference hits (left)')],
               loc='upper left', bbox_to_anchor=(.165, .87), ncol=2, frameon=False, fontsize=10)
    fig.text(.17, .157, 'Control: STRUCTURED_CONTINUE final9216. Net = gained minus lost; all negative branches retained.', fontsize=10)
    fig.text(.17, .101, 'Counts are descriptive recall accounting on annotated TRAIN chains, not newly discovered genes or independent tests.', fontsize=10, color='#5B626B')
    fig.text(.17, .049, 'Single fits; upstream exposure and effective optimization are not controlled. No causal or significance claim.', fontsize=10, color='#5B626B')
    save_figure(fig, out, 'train_continuation_reference_gain_loss')

    provenance = dict(completed=True, experiment='M28-TRAIN-CONTINUATION-FIGURES',
        job_id=os.environ['SLURM_JOB_ID'], generated_UTC=datetime.now(timezone.utc).isoformat(),
        sources={k: dict(path=str(ROOT/'outputs'/root/'summary.json'), job_id=job) for k, (root, job) in SOURCES.items()},
        reference_policy=source['native']['reference_policy'], plotted_metric_rows=rows, plotted_reference_pairs=pairs,
        method_order=list(METHODS), scope_order=list(SCOPES),
        counts_reference_metadata_and_source_overlap_reconciled=True,
        reference_ID_or_geometry_equality_retested=False,
        no_DNA_reference_GFF_predictions_checkpoints_or_features_opened=True,
        model_forward=False, optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
        scientific_winner_or_DEV_candidate_selected=False,
        chart_contract=dict(surface='static scientific article PNG and SVG',
            accuracy_unit='percent, 100 x 2TP / (predicted + reference)',
            pair_unit='annotated primary reference-chain count',
            palette=dict(models=BLUE, gained=BLUE, lost=ORANGE, native='neutral gray'),
            noncolor_identification='method labels, signed pair axis, open hatched losses',
            axis_ranges=dict(accuracy=[0, 100], gain_loss=[-5200, 5000]),
            uncertainty='none; single fits, no independent replicates'),
        limitation='Original TRAIN fitting and paired recall accounting, not independent generalization or causal mechanism.')
    with (out/'provenance.json').open('x') as handle:
        json.dump(provenance, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(completed=True, metric_rows=len(rows), paired_rows=len(pairs),
                         files=sorted(p.name for p in out.iterdir())), allow_nan=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    main(parser.parse_args().output_dir)
