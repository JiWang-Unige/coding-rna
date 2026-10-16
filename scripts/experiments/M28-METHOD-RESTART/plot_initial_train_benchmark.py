#!/usr/bin/env python3
"""Static scientific figure of completed initial raw-FASTA TRAIN results.

Chart contract: how does final4608 C0 compare with three fixed native tools?
Twelve model/scope rows: four methods in AT, Rice and pooled six contigs.
Horizontal bars on identical 0-100% axes; no trend, uncertainty or selection.
One blue root for C0, neutral benchmark tones; names and C0 hatching also
identify methods. PNG/SVG for the existing private repository figure format.
Takeaway: initial complete-chain accuracy is below all three native tools in
both species. TRAIN diagnosis is not independent generalization or causality.
User-authored scientific figure, no provider branding or external renderer.
"""
import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SPECIES = ('arabidopsis_thaliana', 'oryza_sativa')
CONTIGS = {
    SPECIES[0]: ('NC_003070.9', 'NC_003071.7', 'NC_003076.8'),
    SPECIES[1]: ('NC_089043.1', 'NC_089044.1', 'NC_089046.1'),
}
REFERENCES = {SPECIES[0]: 17739, SPECIES[1]: 4949, 'pooled': 22688}
METHODS = ('C0 final4608', 'ANNEVO', 'Tiberius', 'Helixer')
SOURCE_FILES = {'original': 'original_13443329.json', 'native': 'native_13461307.json'}
COLORS = ('#2F6B9A', '#8A949F', '#A7AFB8', '#C4CAD0')


def recompute(row):
    tp, predicted, reference = (row[k] for k in ('tp', 'predicted', 'reference'))
    if not (0 <= tp <= min(predicted, reference) and predicted > 0 and reference > 0):
        raise ValueError('Positive denominators and feasible exact-chain counts required')
    values = dict(precision=tp/predicted, recall=tp/reference, f1=2*tp/(predicted+reference))
    if any(not math.isclose(row[k], v, rel_tol=0, abs_tol=1e-12) for k, v in values.items()):
        raise ValueError('Saved PRF differs from exact-chain counts')
    return dict(tp=tp, predicted=predicted, reference=reference, **values)


def load_rows(source_dir):
    source = {k: json.loads((source_dir/name).read_text()) for k, name in SOURCE_FILES.items()}
    for key, job in (('original', '13443329'), ('native', '13461307')):
        report = source[key]
        if (report['completed'] is not True or report['job_id'] != job
                or report['optimizer_steps_this_run'] != 0 or report['DEV_used']
                or report['sealed_data_used'] or report['reference_used_in_inference']
                or not report['unsupported_reference_strata_retained']
                or not report['training_eligibility_does_not_filter_reference']):
            raise ValueError('Matching completed original TRAIN evaluation required')
    original, native = source['original'], source['native']
    if (not original['checkpoint'].endswith('/M28-FULL-GRAMMAR-COARSE-CHANNEL-FIRST-FIT/step_004608.pt')
            or native['reference_policy'] != original['reference_policy']
            or native['reference_primary_chains'] != 22688
            or set(native['methods']) != set(METHODS[1:])):
        raise ValueError('Initial final4608 and fixed three-method comparison required')
    models = {METHODS[0]: original, **native['methods']}
    rows = []
    for method in METHODS:
        model = models[method]
        if (set(model['chromosomes']) != {c for cells in CONTIGS.values() for c in cells}
                or set(model['per_species']) != set(SPECIES)):
            raise ValueError('Exactly the original two species and six contigs required')
        for seqid, cell in model['chromosomes'].items():
            base = original['chromosomes'][seqid]
            for field in ('seqid', 'species', 'length', 'reference_primary_records',
                          'reference_assessable_isoform_records', 'retained_unsupported_reference_strata'):
                if cell[field] != base[field]:
                    raise ValueError('Reference population differs: '+seqid+' '+field)
            if (cell['CDS_base_unstranded_union']['reference'] != base['CDS_base_unstranded_union']['reference']
                    or cell['background']['bp'] != base['background']['bp']):
                raise ValueError('Reference base/background denominator differs: '+seqid)
            recompute(cell['exact_chain'])
        species_f1 = []
        for scope in (*SPECIES, 'pooled'):
            cell = model['pooled'] if scope == 'pooled' else model['per_species'][scope]
            counts = recompute(cell['exact_chain'])
            seqids = tuple(c for cells in CONTIGS.values() for c in cells) if scope == 'pooled' else CONTIGS[scope]
            if (counts['reference'] != REFERENCES[scope] or any(
                    counts[k] != sum(model['chromosomes'][c]['exact_chain'][k] for c in seqids)
                    for k in ('tp', 'predicted', 'reference'))):
                raise ValueError('Scope counts do not aggregate from the original contigs')
            if scope != 'pooled':
                species_f1.append(counts['f1'])
            rows.append(dict(method=method, scope=scope, **counts,
                f1_percent=100*counts['f1'], CDS_base_F1=cell['CDS_base_unstranded_union']['f1'],
                duplicate_chains=cell['duplicate_chain_records'],
                other_isoform_hits=cell['exact_other_assessable_isoform_not_primary_TP']))
        if (not math.isclose(model['macro_species_exact_chain_F1'], sum(species_f1)/2, abs_tol=1e-12)
                or not math.isclose(model['worst_species_exact_chain_F1'], min(species_f1), abs_tol=1e-12)):
            raise ValueError('Macro/worst-species summaries differ from the species counts')
    return source, rows


def main(source_dir, out):
    source, rows = load_rows(source_dir)
    names = ('initial_train_complete_chain.png', 'initial_train_complete_chain.svg', 'provenance.json')
    if any((out/name).exists() for name in names):
        raise FileExistsError('Preserve existing initial TRAIN figure outputs')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'svg.fonttype': 'none'})
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 5.8), sharex=True, sharey=True)
    fig.subplots_adjust(left=.145, right=.98, top=.73, bottom=.28, wspace=.18)
    for ax, scope, title in zip(axes, (*SPECIES, 'pooled'),
            ('A  Arabidopsis thaliana', 'B  Oryza sativa', 'C  Pooled six contigs')):
        selected = [next(r for r in rows if r['method'] == method and r['scope'] == scope) for method in METHODS]
        for y, (row, color) in enumerate(zip(selected, COLORS)):
            ax.barh(y, row['f1_percent'], height=.55, color=color,
                edgecolor='#424850', linewidth=.65, hatch='///' if y == 0 else None, zorder=3)
            ax.text(row['f1_percent']+1.6, y, f"{row['f1_percent']:.1f}",
                va='center', fontsize=11, color='#242A31')
        ax.set_xlim(0, 100)
        ax.set_ylim(3.65, -.65)
        ax.set_yticks(range(4), ('C0 (final 4,608)', *METHODS[1:]))
        ax.set_xticks((0, 25, 50, 75, 100))
        ax.set_xlabel('Complete-CDS-chain F1 (%)', labelpad=9)
        ax.set_title(title, loc='left', fontsize=12, fontweight='bold', pad=30)
        ax.text(0, 1.045, f"Primary references: {REFERENCES[scope]:,}",
            transform=ax.transAxes, color='#5B626B', fontsize=10)
        ax.grid(axis='x', color='#E4E7EB', linewidth=.65, zorder=0)
        ax.tick_params(axis='both', length=0, colors='#424850', pad=7)
        for name in ('top', 'right'):
            ax.spines[name].set_visible(False)
        for name in ('left', 'bottom'):
            ax.spines[name].set_color('#7E8790')
            ax.spines[name].set_linewidth(.65)
    fig.suptitle('Initial TRAIN complete-CDS-chain accuracy', x=.145, y=.965,
        ha='left', fontsize=18, fontweight='bold', color='#242A31')
    fig.text(.145, .90, 'Same six raw FASTA contigs; fixed primary-reference policy; pooled F1 from pooled counts',
        fontsize=11, color='#424850')
    fig.text(.145, .165, 'C0: fixed final 4,608 updates. Native tools: fixed existing operating points. Difficult references retained.',
        fontsize=10, color='#424850')
    fig.text(.145, .11, 'TRAIN diagnosis, not independent generalization or causal attribution. Ongoing experiment branches are not shown.',
        fontsize=10, color='#5B626B')
    fig.text(.145, .055, 'Sources: original-model evaluation 13443329; native-tool evaluation 13461307. Completed outputs only.',
        fontsize=9.5, color='#5B626B')
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out/names[0], dpi=220, facecolor='white')
    fig.savefig(out/names[1], facecolor='white')
    plt.close(fig)
    provenance = dict(completed=True, experiment='M28-INITIAL-TRAIN-BENCHMARK-FIGURE',
        job_id=os.environ['SLURM_JOB_ID'], generated_UTC=datetime.now(timezone.utc).isoformat(),
        sources=SOURCE_FILES, source_jobs={k: report['job_id'] for k, report in source.items()},
        reference_policy=source['original']['reference_policy'], plotted_rows=rows,
        counts_and_reference_metadata_reconciled=True,
        no_DNA_reference_GFF_predictions_checkpoints_or_features_opened=True,
        model_forward=False, optimizer_steps_this_run=0, DEV_used=False, sealed_data_used=False,
        limitation='Initial TRAIN diagnosis; not an independent test, matched training comparison or causal result.')
    with (out/names[2]).open('x') as handle:
        json.dump(provenance, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(completed=True, plotted_F1_percent={scope: {
        r['method']: r['f1_percent'] for r in rows if r['scope'] == scope}
        for scope in (*SPECIES, 'pooled')}), allow_nan=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    main(args.source_dir, args.output_dir)
