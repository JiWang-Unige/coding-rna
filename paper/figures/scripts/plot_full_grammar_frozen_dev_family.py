#!/usr/bin/env python3
"""Publication figures from the completed fixed DEV family and extension.

Chart contract: compare all ten fixed methods' complete-CDS-chain F1 and
show precision/recall separately so an F1 advantage is not a recall claim.
Scientific manuscript PNG/SVG, not an inline widget or hosted report.
Thirty method/scope rows: ten methods x AT/Rice/pooled two DEV contigs.
F1 uses horizontal bars; precision/recall uses paired horizontal bars.
Keep the original eight-method order, visibly append both predeclared
extensions, retain negative results, use identical 0-100% axes per figure.
Blue root plus neutral native references; extension open/hatch, and recall
open/hatch distinguish meaning without color alone. No provider branding
on this third-party scientific manuscript. Final QA is the exported PNGs.
Takeaway: joint encoding exceeds Helixer's F1 in both DEV species but not
ANNEVO/Tiberius; its pooled recall is below Helixer. Historical DEV reuse,
single fits and upstream exposure preclude an independent-generalization
or causal claim. No model selection, fitting or new inference is performed.
"""
import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path

from plot_initial_train_benchmark import recompute
from plot_full_grammar_train_continuations import BLUE, save_figure
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[3]
SOURCES = {
    'original_eight': ('M28-FULL-GRAMMAR-CONTINUATION-DEV-FAMILY', '13493572'),
    'separate_ten': ('M28-FULL-GRAMMAR-SAMPLING-CONTEXT-DEV-FAMILY', '13493579'),
}
METHODS = (
    'M28_final4608', 'M28_STRUCTURED_CONTINUE_final9216',
    'M28_LOCAL_CHANNEL_final9216', 'M28_HALF_STRUCTURED_HALF_LOCAL_final9216',
    'M28_NATIVE_CDS_PRIOR_final9216', 'ANNEVO', 'Tiberius', 'Helixer',
    'M28_FRESH_DRAWS_final9216', 'M28_JOINT24576_final9216',
)
LABELS = ('C0 (4,608)', 'Structured (9,216)', 'Local (9,216)', 'Mixed (9,216)',
          'Native CDS prior', 'ANNEVO', 'Tiberius', 'Helixer',
          'Fresh uniform (extension)', 'Joint 24,576 (extension)')
SPECIES = ('arabidopsis_thaliana', 'oryza_sativa')
SCOPES = (*SPECIES, 'pooled')
REFERENCES = {SPECIES[0]: 5437, SPECIES[1]: 2291, 'pooled': 7728}
TITLES = ('A  Arabidopsis thaliana', 'B  Oryza sativa', 'C  Pooled two contigs')
COMPARISON_JOBS = {
    METHODS[0]: '13449926', METHODS[1]: '13465375', METHODS[2]: '13465377',
    METHODS[3]: '13465379', METHODS[4]: '13493571',
    METHODS[8]: '13493576', METHODS[9]: '13493578',
}


def scope_cell(model, scope):
    return model['pooled'] if scope == 'pooled' else model['per_species'][scope]


def load_rows():
    reports = {}
    for key, (root, job) in SOURCES.items():
        directory = ROOT/'outputs'/root
        report = json.loads((directory/'summary.json').read_text())
        if ((directory/'STATUS').read_text().strip() != 'COMPLETED'
                or (directory/'JOBID').read_text().strip() != job
                or not report['completed'] or report['job_id'] != job
                or not report['DEV_used'] or report['sealed_data_used']
                or report['optimizer_steps_this_run'] != 0
                or report['new_predictions_or_baseline_runs_created']
                or not report['all_declared_candidates_included']
                or not report['reused_native_baselines_directly_equal_across_reports']):
            raise ValueError('Completed zero-update fixed DEV family required: '+key)
        reports[key] = report
    original, extended = reports.values()
    if (tuple(original['methods']) != METHODS[:8] or tuple(extended['methods']) != METHODS
            or extended['source_comparison_jobs'] != COMPARISON_JOBS
            or not extended['original_eight_method_family_preserved']
            or not extended['prospective_extension_not_part_of_original_eight_method_freeze']
            or extended['reference_primary_chains'] != {s:REFERENCES[s] for s in SPECIES}
            or extended['pooled_reference_primary_chains'] != REFERENCES['pooled']
            or any(original['measurements'][m] != extended['measurements'][m] for m in METHODS[:8])):
        raise ValueError('Preserve every original measurement and both preset extensions')
    rows = []
    for i, (method, label) in enumerate(zip(METHODS, LABELS)):
        model = extended['measurements'][method]
        species = [scope_cell(model, s)['exact_chain'] for s in SPECIES]
        if (set(model['per_species']) != set(SPECIES)
                or not math.isclose(model['pooled']['macro_species_exact_chain_F1'],
                                    sum(r['f1'] for r in species)/2, abs_tol=1e-12)
                or model['pooled']['worst_species_exact_chain_F1'] != min(r['f1'] for r in species)):
            raise ValueError('Fixed species and macro/worst endpoint required')
        for scope in SCOPES:
            cell = scope_cell(model, scope)
            exact = recompute(cell['exact_chain'])
            if (exact['reference'] != REFERENCES[scope]
                    or (scope == 'pooled' and any(exact[k] != sum(r[k] for r in species)
                        for k in ('tp','predicted','reference')))):
                raise ValueError('Fixed reference denominator and pooled count arithmetic required')
            rows.append(dict(method=method, label=label, scope=scope, **exact,
                original_eight_member=i < 8, native_tool=5 <= i < 8,
                predeclared_extension=i >= 8,
                source_comparison_job=COMPARISON_JOBS.get(method, COMPARISON_JOBS[METHODS[0]]),
                cumulative_optimizer_steps=None if 5 <= i < 8 else (4608 if i == 0 else 9216),
                CDS_base_F1=cell['CDS_base_unstranded_union']['f1'],
                wholly_background_chains=cell['background']['wholly_background_unique_chains'],
                background_span_FPR=cell['background']['span_FPR'],
                macro_species_exact_chain_F1=model['pooled']['macro_species_exact_chain_F1'],
                worst_species_exact_chain_F1=model['pooled']['worst_species_exact_chain_F1']))
    return rows


def style_axis(ax, scope, title):
    ax.set_ylim(9.55, -.65)
    ax.set_yticks(range(10), LABELS)
    ax.set_xlim(0, 100)
    ax.set_xticks((0, 25, 50, 75, 100))
    ax.set_title(title, loc='left', fontsize=12, fontweight='bold', pad=29)
    ax.text(0, 1.025, f'Primary references: {REFERENCES[scope]:,}',
            transform=ax.transAxes, fontsize=10, color='#5B626B')
    ax.grid(axis='x', color='#E4E7EB', linewidth=.65, zorder=0)
    ax.axhline(7.5, color='#7E8790', linewidth=.7, linestyle='--', zorder=2)
    ax.tick_params(axis='both', length=0, colors='#424850', pad=7)
    for name in ('top', 'right', 'left'):
        ax.spines[name].set_visible(False)
    ax.spines['bottom'].set_color('#7E8790')
    ax.spines['bottom'].set_linewidth(.65)


def new_figure(title):
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 9.2), sharex=True, sharey=True)
    fig.subplots_adjust(left=.205, right=.985, top=.755, bottom=.205, wspace=.17)
    fig.suptitle(title, x=.205, y=.963, ha='left', fontsize=18,
                 fontweight='bold', color='#242A31')
    fig.text(.205, .906, 'Same two raw FASTA DEV contigs; all 7,728 primary references and fixed operating points', fontsize=11)
    fig.text(.205, .122, 'Original eight methods above the dashed line; two predeclared extensions below. Order is not a ranking.', fontsize=10)
    fig.text(.205, .078, 'Pooled metrics use pooled counts. All M28 continuations: final 9,216 updates; C0: final 4,608.', fontsize=10, color='#5B626B')
    fig.text(.205, .034, 'Historically reused DEV; single fits; unmatched upstream exposure. No independent-generalization or causal claim.', fontsize=10, color='#5B626B')
    return fig, axes


def main(out):
    rows = load_rows()
    if any(out.glob('*.png')) or any(out.glob('*.svg')) or (out/'provenance.json').exists():
        raise FileExistsError('Preserve existing figure outputs')
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'svg.fonttype': 'none'})
    fig, axes = new_figure('Complete-CDS-chain accuracy on the fixed DEV panel')
    for ax, scope, title in zip(axes, SCOPES, TITLES):
        for y, method in enumerate(METHODS):
            row = next(r for r in rows if r['method'] == method and r['scope'] == scope)
            fill, edge, hatch = (('#8A949F', '#424850', None) if row['native_tool']
                                else ('white', BLUE, '///') if row['predeclared_extension']
                                else (BLUE, '#424850', None))
            value = 100*row['f1']
            ax.barh(y, value, height=.58, color=fill, edgecolor=edge,
                    linewidth=.6, hatch=hatch, zorder=3)
            ax.text(value+1.3, y, f'{value:.1f}', va='center', fontsize=10)
        style_axis(ax, scope, title)
        ax.set_xlabel('Complete-CDS-chain F1 (%)', labelpad=9)
    fig.legend(handles=[Patch(facecolor=BLUE, edgecolor='#424850', label='Original M28 family'),
                        Patch(facecolor='#8A949F', edgecolor='#424850', label='Fixed native tools'),
                        Patch(facecolor='white', edgecolor=BLUE, hatch='///', label='Predeclared extension')],
               loc='upper left', bbox_to_anchor=(.20, .864), ncol=3, frameon=False, fontsize=10)
    save_figure(fig, out, 'dev_family_complete_chain')

    fig, axes = new_figure('Complete-CDS-chain precision and recall on the fixed DEV panel')
    for ax, scope, title in zip(axes, SCOPES, TITLES):
        for y, method in enumerate(METHODS):
            row = next(r for r in rows if r['method'] == method and r['scope'] == scope)
            for metric, offset in (('precision', -.17), ('recall', .17)):
                value = 100*row[metric]
                ax.barh(y+offset, value, height=.27, color=BLUE if metric == 'precision' else 'white',
                        edgecolor=BLUE, linewidth=.7, hatch='///' if metric == 'recall' else None, zorder=3)
                ax.text(value+1.3, y+offset, f'{value:.1f}', va='center', fontsize=9)
        style_axis(ax, scope, title)
        ax.set_xlabel('Complete-CDS-chain precision / recall (%)', labelpad=9)
    fig.legend(handles=[Patch(facecolor=BLUE, edgecolor=BLUE, label='Precision: TP / unique predicted chains'),
                        Patch(facecolor='white', edgecolor=BLUE, hatch='///', label='Recall: TP / primary references')],
               loc='upper left', bbox_to_anchor=(.20, .864), ncol=2, frameon=False, fontsize=10)
    save_figure(fig, out, 'dev_family_precision_recall')

    provenance = dict(completed=True, experiment='M28-FROZEN-DEV-FAMILY-FIGURES',
        job_id=os.environ['SLURM_JOB_ID'], generated_UTC=datetime.now(timezone.utc).isoformat(),
        sources={key:dict(path=str(ROOT/'outputs'/root/'summary.json'), job_id=job)
                 for key,(root,job) in SOURCES.items()}, source_comparison_jobs=COMPARISON_JOBS,
        method_order=list(METHODS), scope_order=list(SCOPES), plotted_metric_rows=rows,
        original_eight_measurements_directly_equal=True, count_based_PRF_and_pool_macro_worst_verified=True,
        reference_ID_geometry_or_checkpoint_identity_retested=False,
        no_DNA_reference_GFF_predictions_checkpoints_or_features_opened=True,
        model_forward=False, optimizer_steps_this_run=0, DEV_used=True, DEV_results_only=True,
        sealed_data_used=False, scientific_winner_or_next_training_selected=False,
        chart_contract=dict(surface='existing private scientific manuscript PNG/SVG',
            F1_variant='three_scope_horizontal_bars_all_ten_methods',
            precision_recall_variant='three_scope_paired_horizontal_bars_all_ten_methods',
            reviewed_rows=30, numeric_marks=90, identical_percentage_axes=[0,100],
            palette_policy='single_blue_root_plus_neutral_benchmarks', blue=BLUE,
            non_color='native names, open extension hatching, recall hatching, fixed-family divider',
            provider_branding=False, uncertainty_intervals=False),
        limitation='historically reused development panel, single fits, unmatched upstream annotation exposure; not independent generalization, model selection or causality')
    with (out/'provenance.json').open('x') as handle:
        json.dump(provenance, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(dict(completed=True, job_id=provenance['job_id'], method_scope_rows=len(rows),
                         numeric_marks=90, output=str(out))), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    main(parser.parse_args().output_dir)
