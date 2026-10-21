"""Internal-validation labeled-F1 (protocol §4) of a run on the validation pseudo-contig.  Usage: m30_valscore.py <run_dir> <species> <k>  -> <run_dir>/valscore.json"""
import json, sys
import m30_common as M
import m30_labels as LB
import m30_score as SC

def main():
    run, species, k = sys.argv[1], sys.argv[2], float(sys.argv[3])
    contig = M.SPECIES[species]['val_contig']
    lab = LB.Labeler(species, [contig], k)
    ch = SC.load_gff_chains(f'{run}/{contig}/prediction/full_grammar.gff3', contig)
    r = LB.labeled_f1(ch, lab, contig)
    r.update(species=species, k=k, contig=contig, n_chains=len(ch), n_labeled=lab.n_labeled, n_eligible=lab.n_eligible)
    json.dump(r, open(f'{run}/valscore.json', 'w'))
    print(json.dumps(r))

if __name__ == '__main__':
    main()
