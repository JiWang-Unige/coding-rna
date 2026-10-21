"""Full evaluation of a c1_infer run on whole chromosomes with the project ruler (M26 CDS-assessable primary chains, exact CDS chain, CDS>=96 bp
rule as in the gate / Carbon-A scoring): pooled / per-species / per-chromosome exact-chain F1, exon-count strata (1 / 2-5 / >=6),
CDS-base F1 and background burden (a1_compare definitions).  With --compare-dev also scores the existing DEV predictions of fresh-joint
(GENERanno full grammar) and Carbon-A (filtered / unfiltered) with the SAME functions, for a same-ruler table.
Usage: c1_eval.py <run_dir> <label> <scope: DEV|HELDOUT|TRAIN4> [--compare-dev]
"""
import argparse, json, sys

import m32_common as M
import a1_common as K
sys.path.insert(0, str(K.CA_SRC))
import a1_eval as E                    # noqa: E402
import a1_compare as AC                # noqa: E402

SCOPES = {'DEV': list(M.DEV_SEQS), 'HELDOUT': list(M.HELDOUT), 'TRAIN4': [s for s in M.TRAIN_SEQS if s not in M.HELDOUT]}


def score(chains, seqs):
    r = E.summarize(chains, AC.refs, seqs)
    r['base'] = AC.base_metrics(chains, seqs)
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run')
    ap.add_argument('label')
    ap.add_argument('scope', choices=list(SCOPES))
    ap.add_argument('--compare-dev', action='store_true')
    args = ap.parse_args()
    seqs = SCOPES[args.scope]
    chains = {s: E.load_chains(args.run, s) for s in seqs}
    raw = {s: E.load_chains(args.run, s, cds96=False) for s in seqs}
    out = dict(label=args.label, run=args.run, scope=args.scope, seqs=seqs, cds96_rule=True)
    out['PC2_C1'] = score(chains, seqs)
    out['PC2_C1']['unfiltered_no_cds96_pooled'] = E.summarize(raw, AC.refs, seqs)['pooled']
    if args.compare_dev and args.scope == 'DEV':
        import torch
        from constants import CONFIDENCE_NETWORK_PATH      # noqa: E402
        out['GENERanno_fresh_joint_B'] = score({s: AC.gff_chains(AC.DEVGFF['B_freshjoint'][s], s) for s in seqs}, seqs)
        out['CarbonA_unfiltered'] = score({s: AC.carbon_dev(s, False) for s in seqs}, seqs)
        out['CarbonA_filtered'] = score({s: AC.carbon_dev(s, True) for s in seqs}, seqs)
    (M.Path(args.run) / f'eval_{args.label}_{args.scope}.json').write_text(json.dumps(out, indent=1))
    for k, v in out.items():
        if isinstance(v, dict) and 'pooled' in v:
            print(args.scope, k, round(v['pooled']['F1'], 3), {s: round(d['F1'], 2) for s, d in v['by_species'].items()},
                  {a: round(b, 2) for a, b in v['base'].items()}, {a: round(b['F1'], 2) for a, b in v['by_strata'].items()}, flush=True)


if __name__ == '__main__':
    main()
