#!/usr/bin/env python3
"""Explain fixed A3 coverage gaps from saved DEV geometry and grammar metadata."""
import csv
import json
from bisect import bisect_left, bisect_right
from collections import Counter
from pathlib import Path
from src.m28.assembly import nearest_owner
from src.m28.candidates import Budget

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT/'outputs/M28-A2-A3-REFERENCE-TURNOVER'
OUT = ROOT/'outputs/M28-A3-REPRESENTABILITY'
ARMS = ('A2','A3')
STAGES = ('free_unique','owner_unique','selected_unique')
EXPECTED = {'arabidopsis_thaliana':5437,'oryza_sativa':2291}
WIDTH, STRIDE = 24576,12288
GRAMMAR_FLAGS = {'length_not_mod3','non_ATG_start','no_terminal_stop',
                 'internal_stop','ambiguous_CDS','split_terminal_codon','nonstandard_splice'}
ANNOTATION_FLAGS = {'same_strand_primary_overlap','longer_than_window'}
GROUPS = ('no_containing_natural_window','in_window_grammar_unsupported',
          'in_window_grammar_supported')


def geometry(chain, windows):
    """Containment uses CDS-span endpoints, valid for either genomic strand."""
    a,b = chain[0][0],chain[-1][1]
    starts,ends = [x for x,y in windows],[y for x,y in windows]
    left,right = bisect_left(ends,b),bisect_right(starts,a)
    count = max(0,right-left)
    owner = nearest_owner(chain,windows)
    x,y = windows[owner]
    return count,owner,x <= a and b <= y


def classify(row, windows, max_exons, min_intron):
    chain = tuple(tuple(x) for x in json.loads(row['CDS_intervals']))
    flags = set(filter(None,row['unsupported_strata'].split(';')))
    if flags-(GRAMMAR_FLAGS|ANNOTATION_FLAGS):
        raise ValueError('Unmapped reference metadata must not be called decoder grammar: '+str(flags))
    reasons = flags & GRAMMAR_FLAGS
    if len(chain) > max_exons:
        reasons.add('more_than_decoder_max_exons')
    if any(right[0]-left[1] < min_intron for left,right in zip(chain,chain[1:])):
        reasons.add('intron_shorter_than_decoder_minimum')
    count,owner,contains = geometry(chain,windows)
    group = GROUPS[0] if not count else GROUPS[1] if reasons else GROUPS[2]
    span = chain[-1][1]-chain[0][0]
    result = dict(row,representability_group=group,containing_natural_windows=count,
                  owner_window_index=owner,owner_window_contains_CDS_span=contains,
                  window_failure=('longer_than_window' if span > WIDTH else 'natural_grid_gap') if not count else '',
                  decoder_grammar_reasons=';'.join(sorted(reasons)),
                  annotation_overlap_flag='same_strand_primary_overlap' in flags,
                  exon_class='one' if len(chain)==1 else 'two_to_four' if len(chain)<=4 else 'five_or_more',
                  span_class='at_most_stride' if span<=STRIDE else 'stride_to_width' if span<=WIDTH else 'longer_than_width')
    for arm in ARMS:
        flags = [row[arm+'_'+stage] for stage in STAGES]
        if any(v not in ('True','False') for v in flags):
            raise ValueError('Saved membership must be explicit')
        present = [v=='True' for v in flags]
        if (present[1] and not present[0]) or (present[2] and not present[1]):
            raise ValueError('Saved reference stages must be nested')
        if group != GROUPS[2] and present[0]:
            raise ValueError('Observed exact chain contradicts representability diagnosis for '+row['id'])
    return result


def counts(rows):
    result = {'reference_chains':len(rows),'methods':{}}
    for arm in ARMS:
        values = {stage:sum(r[arm+'_'+stage]=='True' for r in rows) for stage in STAGES}
        values['missing_from_free'] = len(rows)-values['free_unique']
        values['free_coverage'] = values['free_unique']/len(rows) if rows else None
        values['states'] = dict(Counter(r[arm+'_state'] for r in rows))
        result['methods'][arm] = values
    return result


def summarize(rows):
    result = counts(rows)
    result['groups'] = {group:counts([r for r in rows if r['representability_group']==group])
                        for group in GROUPS}
    result['window_failure_counts'] = dict(Counter(r['window_failure'] for r in rows if r['window_failure']))
    result['decoder_grammar_reason_counts_nonexclusive'] = dict(Counter(
        reason for r in rows for reason in r['decoder_grammar_reasons'].split(';') if reason))
    result['annotation_overlap_flag'] = counts([r for r in rows if r['annotation_overlap_flag']])
    supported = [r for r in rows if r['representability_group']==GROUPS[2]]
    result['supported_by_exon_class'] = {key:counts([r for r in supported if r['exon_class']==key])
                                       for key in ('one','two_to_four','five_or_more')}
    result['supported_by_span_class'] = {key:counts([r for r in supported if r['span_class']==key])
                                       for key in ('at_most_stride','stride_to_width','longer_than_width')}
    result['A2_selected_lost_in_A3'] = counts([r for r in rows
        if r['A2_selected_unique']=='True' and r['A3_selected_unique']=='False'])
    return result


def main():
    if OUT.exists():
        raise FileExistsError('Use the single new representability output')
    if (SOURCE/'STATUS').read_text().strip()!='COMPLETED':
        raise ValueError('Completed full turnover source required')
    saved = json.loads((SOURCE/'summary.json').read_text())
    decoded_path = ROOT/'outputs/M28-A3-DEV-INFER'
    if (decoded_path/'STATUS').read_text().strip()!='COMPLETED':
        raise ValueError('Completed A3 inference grid required')
    decoded = json.loads((decoded_path/'summary.json').read_text())
    if (not decoded['inference_complete'] or decoded['windows']!=8690 or
        decoded['checkpoint']!=saved['source']['A3']['checkpoint'] or decoded['reference_used'] or
        decoded['test_setaria_used'] or saved['reference_primary_chains']!=7728):
        raise ValueError('Fixed source checkpoint, scope or denominator differs')
    scope = {(s['species'],s['seqid'],s['strand']):s['windows'] for s in decoded['scopes']}
    if len(scope)!=4:
        raise ValueError('Both complete DEV contigs and strands required')
    for windows in scope.values():
        if (windows!=sorted(windows) or len({tuple(w) for w in windows})!=len(windows) or
            any(b-a!=WIDTH for a,b in windows) or
            [b for a,b in windows]!=sorted(b for a,b in windows)):
            raise ValueError('Actual natural grid differs from fixed full-width DEV windows')
    with (SOURCE/'references.tsv').open() as handle:
        original = list(csv.DictReader(handle,delimiter='\t'))
    if (len(original)!=7728 or len({(r['species'],r['id']) for r in original})!=7728 or
        Counter(r['species'] for r in original)!=EXPECTED):
        raise ValueError('Every frozen primary reference must be retained exactly once')
    max_exons = decoded['configuration']['max_exons']
    min_intron = Budget().min_intron_bp
    enriched = [classify(row,scope[row['species'],row['seqid'],row['strand']],max_exons,min_intron)
                for row in original]
    result = {'experiment':OUT.name,'completed':True,
              'scope':'saved DEV representability diagnosis, not learning causality or new accuracy',
              'source':str(SOURCE.relative_to(ROOT)),'source_checkpoints':saved['source'],
              'reference_primary_chains':7728,'reference_filters_added':False,
              'reference_used_for_diagnosis':True,'new_model_inference':False,
              'optimizer_steps_this_run':0,'test_setaria_used':False,
              'rules':{'width':WIDTH,'stride':STRIDE,'max_exons':max_exons,'min_intron_bp':min_intron,
                       'same_strand_overlap_is_not_decoder_grammar':True,
                       'support_is_necessary_not_sufficient_for_greedy_query_generation':True,
                       'grammar_metadata_source':'frozen M28 primary-reference sequence checks',
                       'main_accuracy_denominator_unchanged':True},
              'pooled':summarize(enriched),
              'per_species':{species:summarize([r for r in enriched if r['species']==species])
                             for species in EXPECTED}}
    for key,group in [('pooled',result['pooled']),*result['per_species'].items()]:
        previous = saved['pooled'] if key=='pooled' else saved['per_species'][key]
        for arm in ARMS:
            for stage in STAGES:
                observed = group['methods'][arm][stage]
                expected = previous['stages'][stage]['both']+previous['stages'][stage][arm+'_only']
                if observed!=expected:
                    raise ValueError('Classification changed saved coverage accounting')
    OUT.mkdir()
    with (OUT/'references.tsv').open('x',newline='') as handle:
        writer = csv.DictWriter(handle,fieldnames=list(enriched[0]),delimiter='\t')
        writer.writeheader(); writer.writerows(enriched)
    with (OUT/'summary.json').open('x') as handle:
        json.dump(result,handle,indent=2,allow_nan=False); handle.write('\n')
    print(json.dumps(result['pooled'],indent=2),flush=True)


if __name__=='__main__':
    main()
