#!/usr/bin/env python3
"""Local scientific figure of closed A36/A37 TRAIN decisions, not annotation F1.

Chart contract: compare selected match/nonmatch/unknown counts in both species,
and pooled matching >=5-exon pairs against the original A31_A3 baseline.
Six species/model rows and three pooled-long-chain rows; no time series/scatter.
Static Matplotlib PNG/SVG for the existing repository's scientific figures.
Two palette roots (blue/orange), neutral unknown, hatches and direct labels.
Takeaway: fewer known nonmatches coexist with fewer matching pairs; both arms
fail three retention conditions. Unmatched training inputs/capacity rule out
causal attribution. User-authored scientific figure, no provider branding.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

SPECIES = ('arabidopsis_thaliana', 'oryza_sativa')
ARMS = ('A31_A3', 'A36', 'A37')
FILES = {'A31_A3': 'a31_readout_13398521.json',
         'A36': 'a36_CPU_13408451.json', 'A37': 'a37_CPU_13410742.json'}
REFS = {'arabidopsis_thaliana': {'eligible': 2172, 'five_or_more': 854},
        'oryza_sativa': {'eligible': 744, 'five_or_more': 232}}
FIELDS = ('selected_correct', 'selected_known_nonmatch', 'selected_unknown')
COLORS = ('#2F6B9A', '#D88B40', '#BEC3CC')


def load_counts(source_dir):
    reports = {arm: json.loads((source_dir/name).read_text()) for arm, name in FILES.items()}
    counts = {s: {} for s in SPECIES}
    gates = {}
    for arm, report in reports.items():
        if (report['completed'] is not True or report['windows'] != 1491 or
                report['eligible_reference_pairs'] != REFS):
            raise ValueError('Complete original TRAIN summaries and reference-window denominators required')
        for species in SPECIES:
            old = reports['A31_A3']['selected_counts_per_species'][species]['A31_A3']
            row = report['selected_counts_per_species'][species][arm]
            if (sum(row[k] for k in FIELDS) != row['selected'] or
                    sum(row[k] for k in ('pool_correct', 'pool_known_nonmatch', 'pool_unknown')) != row['candidates'] or
                    any(row[k] != old[k] for k in ('candidates', 'pool_correct', 'pool_known_nonmatch', 'pool_unknown')) or
                    (arm != 'A31_A3' and report['selected_counts_per_species'][species]['A31_A3'] != old)):
                raise ValueError('Selected composition or unchanged stronger-baseline natural population differs')
            counts[species][arm] = row
        if arm == 'A31_A3':
            continue
        proof = ('actual_saved_labels_selections_stronger_baseline_and_five_conditions_recomputed' if arm == 'A36'
                 else 'actual_saved_labels_DP_baseline_denominators_and_five_conditions_recomputed')
        expected_sources = ({'fit': '13406144', 'readout': '13408449'} if arm == 'A36'
                            else {'fit': '13410061', 'readout': '13410739'})
        if (report['experiment'] != 'M28-'+arm+('-STRUCTURED-EMISSION-RESULT-CHECK' if arm == 'A36'
                    else '-CONTEXT-CHAIN-RESULT-CHECK') or report[proof] is not True or
                report['source_jobs'] != expected_sources or report['natural_candidates'] != 16321 or
                report['reference_injected_predictions'] != 0):
            raise ValueError('Matching completed closed-arm CPU numerical reconciliation required')
        actual = {}
        for species in SPECIES:
            old, new = counts[species]['A31_A3'], counts[species][arm]
            actual[species+'_natural_selected_correct_strictly_increased'] = new['selected_correct'] > old['selected_correct']
            actual[species+'_natural_selected_known_nonmatch_not_increased'] = new['selected_known_nonmatch'] <= old['selected_known_nonmatch']
        actual['pooled_natural_selected_five_or_more_correct_strictly_increased'] = (
            sum(counts[s][arm]['selected_five_or_more_correct'] for s in SPECIES) >
            sum(counts[s]['A31_A3']['selected_five_or_more_correct'] for s in SPECIES))
        if actual != report['predeclared_TRAIN_gate'] or report['TRAIN_preparation_supported'] is not all(actual.values()):
            raise ValueError('Closed-arm original five conditions differ from displayed counts')
        gates[arm] = actual
    if sum(counts[s]['A31_A3']['candidates'] for s in SPECIES) != 16321:
        raise ValueError('Original complete natural candidate population required')
    return reports, counts, gates


def style_axis(ax):
    ax.set_axisbelow(True); ax.grid(axis='x', color='#E8EAED', linewidth=.7)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color('#8D939C'); ax.spines[side].set_linewidth(.7)
    ax.tick_params(colors='#424850', width=.7)
    ax.set_xlabel('Window–chain pairs', color='#424850')


def main(source_dir, out):
    reports, counts, gates = load_counts(source_dir)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                        'svg.fonttype': 'none', 'axes.titlesize': 12, 'figure.facecolor': 'white'})
    fig = plt.figure(figsize=(12.8, 8.7))
    grid = fig.add_gridspec(2, 2, left=.09, right=.965, bottom=.19, top=.77,
                           height_ratios=(1, .83), hspace=.71, wspace=.30)
    max_total = max(counts[s][a]['selected'] for s in SPECIES for a in ARMS)
    titles = ('A   Arabidopsis thaliana', 'B   Oryza sativa')
    for ordinal, species in enumerate(SPECIES):
        ax = fig.add_subplot(grid[0, ordinal])
        for row_index, arm in enumerate(ARMS):
            left = 0
            for field, color, hatch in zip(FIELDS, COLORS, ('', '///', '..')):
                value = counts[species][arm][field]
                ax.barh(row_index, value, left=left, height=.58, color=color,
                        edgecolor='#424850', linewidth=.65, hatch=hatch)
                if value >= 130:
                    ax.text(left+value/2, row_index, str(value), ha='center', va='center',
                            fontsize=9, color='white' if field == 'selected_correct' else '#242A31')
                left += value
            ax.text(left+28, row_index, 'total '+str(left), va='center', fontsize=9, color='#242A31')
        ax.set_yticks(range(3), ARMS); ax.invert_yaxis()
        ax.set_xlim(0, max_total*1.14)
        ax.set_title(titles[ordinal], loc='left', pad=28, fontstyle='italic')
        ax.text(0, 1.035, f"Natural candidates {counts[species]['A31_A3']['candidates']:,} | "
                f"eligible reference pairs {REFS[species]['eligible']:,}",
                transform=ax.transAxes, fontsize=9, color='#5B626B')
        style_axis(ax)
    long_ax = fig.add_subplot(grid[1, :])
    pooled = {a: {k: sum(counts[s][a][k] for s in SPECIES) for k in FIELDS+('selected', 'selected_five_or_more_correct')}
              for a in ARMS}
    for i, arm in enumerate(ARMS):
        value = pooled[arm]['selected_five_or_more_correct']
        long_ax.barh(i, value, height=.55, color=COLORS[0], edgecolor='#424850', linewidth=.65)
        long_ax.text(value+12, i, str(value), va='center', fontsize=10, color='#242A31')
    baseline = pooled['A31_A3']['selected_five_or_more_correct']
    long_ax.axvline(baseline, color='#424850', linestyle='--', linewidth=1)
    long_ax.text(baseline+12, -.63, f'Baseline {baseline}; continuation requires strictly more',
                 fontsize=9, color='#424850')
    long_ax.set_ylim(2.5, -.7); long_ax.set_xlim(0, 1086)
    long_ax.set_yticks(range(3), ARMS)
    long_ax.set_title('C   Pooled selected matching chains with ≥5 exons', loc='left', pad=27)
    long_ax.text(0, 1.045, 'Original eligible ≥5-exon reference-window pairs: 1,086 | counts, not recall or gene F1',
                 transform=long_ax.transAxes, fontsize=9, color='#5B626B')
    style_axis(long_ax)
    fig.suptitle('Selected chain counts in closed TRAIN experiments', x=.09, y=.972,
                 ha='left', fontsize=17, fontweight='bold', color='#242A31')
    fig.text(.09, .923, '1,491 original windows · 16,321 natural candidates · same prediction pool and gain-zero selection',
             fontsize=10.5, color='#424850')
    legend = [Patch(facecolor=c, edgecolor='#424850', linewidth=.65, hatch=h, label=l)
              for c, h, l in zip(COLORS, ('', '///', '..'),
                  ('Reference-relative match (1)', 'Known nonmatch (0)', 'Protected unknown (−1)'))]
    fig.legend(handles=legend, loc='upper left', bbox_to_anchor=(.082, .899), ncol=3,
               frameon=False, fontsize=10, handlelength=2.5, columnspacing=1.8)
    fig.text(.09, .116, 'A36 and A37 each pass 2/5 original conditions: nonmatches decrease, but matching-pair and long-chain retention fail.',
             fontsize=10, color='#242A31')
    fig.text(.09, .083, 'Grain is a TRAIN window–chain pair, not a unique gene. Unknown is not an error; known nonmatch is not biological false positive.',
             fontsize=9, color='#5B626B')
    fig.text(.09, .054, 'Training differs: A36 uses positive injection / 366,191 parameters; A37 is natural-only / 612,324. Descriptive, not causal.',
             fontsize=9, color='#5B626B')
    fig.text(.09, .026, 'Sources: A31 readout 13398521; A36 CPU 13408451; A37 CPU 13410742. A38 excluded. No raw-FASTA or generalization claim.',
             fontsize=9, color='#5B626B')
    for name in ('closed_train_chain_counts.png', 'closed_train_chain_counts.svg', 'provenance.json'):
        if (out/name).exists():
            raise FileExistsError('Preserve existing scientific figure: '+str(out/name))
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out/'closed_train_chain_counts.png', dpi=300, facecolor='white')
    fig.savefig(out/'closed_train_chain_counts.svg', facecolor='white')
    plt.close(fig)
    provenance = dict(completed=True, experiment='M28-CLOSED-TRAIN-CHAIN-FIGURE',
        generated_UTC=datetime.now(timezone.utc).isoformat(), sources=FILES,
        source_experiments={a: r.get('experiment', 'A31_completed_TRAIN_readout') for a, r in reports.items()},
        windows=1491, natural_candidates=16321, eligible_reference_pairs=REFS,
        selected_counts_per_species=counts, pooled=pooled, original_five_conditions=gates,
        grain='original_TRAIN_window_natural_chain_pair_not_unique_gene',
        model_or_backbone_forward=False, model_or_feature_tensor_opened=False, optimizer_steps_this_run=0,
        DEV_used=False, test_setaria_used=False, A33_outputs_opened=False, A38_outputs_opened=False,
        palette_policy='two_roots_blue_orange_plus_neutral_unknown_with_hatches',
        limitation='Closed-arm descriptive counts, unmatched training/capacity; not causal, gene F1 or generalization.')
    with (out/'provenance.json').open('x') as handle:
        json.dump(provenance, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(dict(pooled=pooled, passed_conditions={a: sum(g.values()) for a, g in gates.items()})))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(); main(args.source_dir, args.output_dir)
