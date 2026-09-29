#!/usr/bin/env python3
"""M28 continuous-window manifests; reference metadata is not an inference input."""
import argparse
import bisect
import json
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import census as C
from eval_structure_diagnostic import parse_attrs, is_partial

WIDTH, STRIDE, FLANK = 24576, 12288, 256
STOP = {"TAA", "TAG", "TGA"}

def rc(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]

def starts(length, width=WIDTH, stride=STRIDE):
    if length <= width:
        return [0]
    values = list(range(0, length-width+1, stride))
    if values[-1] != length-width:
        values.append(length-width)
    return values

def orient_interval(a, b, start, end, strand):
    return (a-start, b-start) if strand == "+" else (end-b, end-a)

def read_allowed(path, lengths):
    seqs, buf, q = {}, [], None
    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                if q in lengths:
                    seqs[q] = "".join(buf).upper()
                q, buf = line[1:].split()[0], []
            elif q in lengths:
                buf.append(line.strip())
    if q in lengths:
        seqs[q] = "".join(buf).upper()
    if {q:len(s) for q,s in seqs.items()} != lengths:
        raise ValueError("FASTA lengths do not match allowed reference metadata")
    return seqs

def chain_checks(tx, dna):
    cds = tx["CDS"] if tx["strand"] == "+" else tx["CDS"][::-1]
    pieces = [dna[a:b] if tx["strand"] == "+" else rc(dna[a:b]) for a,b,_ in cds]
    sequence = "".join(pieces)
    failures = []
    if len(sequence)%3: failures.append("length_not_mod3")
    if sequence[:3] != "ATG": failures.append("non_ATG_start")
    if sequence[-3:] not in STOP: failures.append("no_terminal_stop")
    if any(sequence[i:i+3] in STOP for i in range(0,max(0,len(sequence)-3),3)):
        failures.append("internal_stop")
    if set(sequence)-set("ACGT"): failures.append("ambiguous_CDS")
    offset = 0
    for a,b,p in cds:
        if p not in ("0","1","2") or int(p) != (3-offset%3)%3:
            failures.append("GFF_phase_disagreement")
            break
        offset += b-a
    if len(pieces[0]) < 3 or len(pieces[-1]) < 3:
        failures.append("split_terminal_codon")
    # Native C0 supports canonical GT-AG and GC-AG; record other motifs, do not delete evaluation references.
    ordered = tx["CDS"]
    motifs = []
    for left,right in zip(ordered,ordered[1:]):
        intron = dna[left[1]:right[0]]
        if tx["strand"] == "-": intron = rc(intron)
        motifs.append(intron[:2] + "-" + intron[-2:])
    if any(m not in ("GT-AG","GC-AG") for m in motifs):
        failures.append("nonstandard_splice")
    return sorted(set(failures)), dict(Counter(motifs))

def merged(intervals):
    result = []
    for a,b in sorted(intervals):
        if b <= a: continue
        if result and a <= result[-1][1]:
            result[-1][1] = max(result[-1][1], b)
        else:
            result.append([a,b])
    return result

def probabilities(rows, enriched_fraction=.25):
    # Marginal of: pick a positive gene uniformly, then one of its existing grid windows uniformly.
    gene_windows = defaultdict(list)
    for i,r in enumerate(rows):
        for ident in r["positive_ids"]: gene_windows[ident].append(i)
    q = np.full(len(rows),1/len(rows),dtype=float)
    if gene_windows:
        q *= 1-enriched_fraction
        for indices in gene_windows.values():
            q[indices] += enriched_fraction / (len(gene_windows)*len(indices))
    return q

def sample_stratum(rows, n, rng, enriched_fraction=.25):
    q = probabilities(rows,enriched_fraction)
    return [{"window_id":rows[int(i)]["id"], "stratum_probability":float(q[i]),
             "importance_weight_uniform_within_stratum":float(1/(len(rows)*q[i]))}
            for i in rng.choice(len(rows),size=n,replace=True,p=q)]

def allocate(lengths, n=384):
    # Largest-remainder allocation, same chromosome counts on each strand.
    total = sum(lengths.values())
    exact = {q:n*length/total for q,length in lengths.items()}
    counts = {q:int(x) for q,x in exact.items()}
    for q in sorted(counts,key=lambda q:(-(exact[q]-counts[q]),q))[:n-sum(counts.values())]:
        counts[q] += 1
    return counts

def partial_cds_ids(path, allowed):
    ids = set()
    with path.open() as f:
        for line in f:
            if line.startswith("#"): continue
            fields = line.rstrip().split("\t")
            if len(fields)!=9 or fields[0] not in allowed or fields[2]!="CDS": continue
            attrs = parse_attrs(fields[8])
            if is_partial(attrs): ids.update(attrs.get("Parent","").split(","))
    return ids

def build(species, out):
    meta = C.metadata(C.yaml.safe_load(C.CONFIG.read_text()))[species]
    dna = read_allowed(meta["path"]/"genome.fa",meta["lengths"])
    coding = C.parse_annotation(meta["path"]/"reference.gff3",meta["lengths"],protein_coding_only=True)
    all_ann = C.parse_annotation(meta["path"]/"reference.gff3",meta["lengths"])
    primary = {t["id"] for t in C.primary_transcripts(coding)}
    primary_all = C.primary_transcripts(coding,complete_only=False)
    training_primary = {t["id"] for t in primary_all}
    cds_partial = partial_cds_ids(meta["path"]/"reference.gff3",meta["lengths"])
    overlaps = set()
    for q in meta["seqids"]:
        overlaps |= C.overlapping_ids([t for t in primary_all if t["seqid"]==q])
    catalog, by_qs, genes = {}, defaultdict(list), defaultdict(list)
    for g in all_ann["genes"].values():
        genes[g["seqid"]].append((g["start"],g["end"]))
    for t in coding["transcripts"].values():
        if not t["CDS"]: continue
        gene = coding["genes"][t["gene_id"]]
        partial = bool(gene["partial"] or t["partial"])
        failures, motifs = chain_checks(t,dna[t["seqid"]])
        row = {"id":t["id"],"gene_id":t["gene_id"],"seqid":t["seqid"],"strand":t["strand"],
               "CDS":t["CDS"],"span":[t["CDS"][0][0],t["CDS"][-1][1]],
               "split":meta["splits"][t["seqid"]],"primary_complete":t["id"] in primary,
               "partial":partial,"CDS_partial_attribute":t["id"] in cds_partial,
               "training_primary":t["id"] in training_primary,
               "CDS_complete_supported":not failures and t["id"] not in cds_partial,
               "sequence_structure_failures":failures,"splice_motifs":motifs,
               "overlapping_primary":t["id"] in overlaps,
               "eligible_positive":t["id"] in training_primary and t["id"] not in cds_partial and not failures and t["id"] not in overlaps,
               "exclude_as_negative":True}
        catalog[t["id"]] = row
        by_qs[t["seqid"],t["strand"]].append(row)
    out.mkdir(parents=True,exist_ok=True)
    with (out/"transcripts.jsonl").open("x") as f:
        for row in catalog.values(): f.write(json.dumps(row)+"\n")
    summaries, all_windows, exposure = {}, [], Counter()
    for split in ("train","val"):
        rows = []
        with (out/("windows_"+split+".jsonl")).open("x") as f:
            for q in meta["seqids"]:
                if meta["splits"][q] != split: continue
                L = meta["lengths"][q]
                for strand in ("+","-"):
                    txs = sorted(by_qs[q,strand],key=lambda t:t["span"][0])
                    txstarts = [t["span"][0] for t in txs]
                    for a in starts(L):
                        b = min(L,a+WIDTH)
                        local = [t for t in txs[:bisect.bisect_left(txstarts,b)] if t["span"][1] > a]
                        fully = [t for t in local if a <= t["span"][0] and t["span"][1] <= b]
                        pos = [t["id"] for t in fully if t["eligible_positive"] and a+FLANK <= t["span"][0] and t["span"][1] <= b-FLANK]
                        uncertain, partial_regions = [], []
                        for t in local:
                            x,y = max(a,t["span"][0]),min(b,t["span"][1])
                            interval = orient_interval(x,y,a,b,strand)
                            if not t["CDS_complete_supported"] or t["overlapping_primary"]:
                                partial_regions.append(interval)
                                uncertain.append(interval)
                            elif t not in fully:
                                uncertain.append(interval)
                        sequence = dna[q][a:b] if strand=="+" else rc(dna[q][a:b])
                        row = {"id":species+"|"+q+"|"+str(a)+"|"+strand,"species":species,"split":split,
                               "seqid":q,"start":a,"end":b,"width":WIDTH,"valid_bases":b-a,"strand":strand,
                               "pad_right_after_orientation":WIDTH-(b-a),"primary_complete_ids":[t["id"] for t in fully if t["primary_complete"]],
                               "positive_ids":pos,"overlapping_transcript_ids":[t["id"] for t in local],
                               "chain_unknown_intervals":merged(uncertain),"fine_label_unknown_intervals":merged(partial_regions),
                               "no_annotated_gene_either_strand":not any(x < b and y > a for x,y in genes[q]),
                               "motifs":{"start":sequence.count("ATG"),"stop":sum(sequence.count(m) for m in STOP),
                                         "donor":sequence.count("GT")+sequence.count("GC"),"acceptor":sequence.count("AG")}}
                        exposure.update(pos)
                        rows.append(row)
                        f.write(json.dumps(row)+"\n")
        all_windows.extend(rows)
        refs = [t for t in catalog.values() if t["split"]==split and t["primary_complete"]]
        observed = {i for r in rows for i in r["primary_complete_ids"]}
        eligible = {t["id"] for t in catalog.values() if t["split"]==split and t["eligible_positive"]}
        exposed = {i for r in rows for i in r["positive_ids"]}
        summaries[split] = {"windows":len(rows),"reference_complete_primary":len(refs),
                            "geometrically_covered_primary":len(observed),"eligible_positive_primary":len(eligible),
                            "positive_primary_with_256bp_flanks":len(exposed),
                            "parent_or_transcript_partial_with_complete_CDS_positives":sum(t["split"]==split and t["eligible_positive"] and t["partial"] for t in catalog.values()),
                            "unexposed_eligible_ids":sorted(eligible-exposed),
                            "reference_primary_without_training_eligibility":dict(Counter(
                                reason for t in refs for reason in (t["sequence_structure_failures"] + (["primary_overlap"] if t["overlapping_primary"] else []) + (["CDS_partial_attribute"] if t["CDS_partial_attribute"] else [])))),
                            "reference_background_windows":sum(r["no_annotated_gene_either_strand"] for r in rows),
                            "windows_with_positive":sum(bool(r["positive_ids"]) for r in rows),
                            "max_primary_per_window":max((len(r["primary_complete_ids"]) for r in rows),default=0),
                            "max_exons_per_primary":max((len(t["CDS"]) for t in refs),default=0),
                            "max_motifs_per_window":{k:max((r["motifs"][k] for r in rows),default=0) for k in ("start","stop","donor","acceptor")}}
    strata = defaultdict(list)
    for r in all_windows:
        if r["split"]=="train": strata[r["seqid"],r["strand"]].append(r)
    # 768 draws/species, half each strand; chromosomes proportional to valid sequence length.
    rng = np.random.default_rng(0)
    draws = []
    counts = allocate({q:L for q,L in meta["lengths"].items() if meta["splits"][q]=="train"})
    for key, rows in sorted(strata.items()):
        n = counts[key[0]]
        draws += [dict(x,stratum=list(key),stratum_size=len(rows),
                       probability_per_species_draw=n/768*x["stratum_probability"],
                       probability_per_paired_species_draw=n/1536*x["stratum_probability"])
                  for x in sample_stratum(rows,n,rng)]
    by_id = {r["id"]:r for r in all_windows}
    with (out/"training_draws.jsonl").open("x") as f:
        for ordinal,r in enumerate(draws):
            f.write(json.dumps(dict(r,ordinal=ordinal))+"\n")
    draw_exposure = Counter(i for r in draws for i in by_id[r["window_id"]]["positive_ids"])
    report = {"experiment":"M28-DATA-MANIFEST-R2","species":species,"width":WIDTH,"stride":STRIDE,"flank":FLANK,
              "no_model_training":True,"setaria_access":False,"test_sequence_or_labels_used":False,
              "eligible_positive_policy":"longest-CDS primary; CDS itself complete (parent/tx partial not automatic exclusion), ATG/stop, phase-consistent, no internal stop, GT/GC-AG, no same-strand primary overlap",
              "mask_policy":"unqualified/overlap fine labels and crossing-window chain targets unknown; retain certain CDS/splice local evidence; alternative-isoform conflicts masked by label builder",
              "evaluation_reference_policy":"all original complete primaries, including unsupported/long/overlapping; no denominator filtering",
              "sampling":{"draws":len(draws),"unique_windows":len({r["window_id"] for r in draws}),
                          "strata":len(strata),"seed":0,"uniform_fraction":.75,"complete_positive_fraction":.25,
                          "with_replacement":True,"unique_positive_chains_exposed":len(draw_exposure),
                          "background_draws":sum(by_id[r["window_id"]]["no_annotated_gene_either_strand"] for r in draws),
                          "positive_exposure_min_max":[min(draw_exposure.values(),default=0),max(draw_exposure.values(),default=0)],
                          "chromosome_draws_per_strand":counts,
                          "population":"equal species and strands, chromosome counts by valid length (largest remainder); importance weights undo gene-first within-stratum enrichment"},
              "splits":summaries}
    with (out/"summary.json").open("x") as f:
        json.dump(report,f,indent=2); f.write("\n")
    print(json.dumps(report,indent=2))

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--species",choices=("arabidopsis_thaliana","oryza_sativa"),required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    a=p.parse_args()
    build(a.species,a.output_dir)
