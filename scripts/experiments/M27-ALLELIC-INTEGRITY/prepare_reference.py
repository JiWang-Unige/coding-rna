#!/usr/bin/env python3
"""Reference-only feasibility preparation for the M27 allelic-integrity pilot."""
import argparse
import csv
import gzip
import json
import math
import re
import shutil
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

FASTA_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz"
GTF_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
LENGTH = 50818468
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
STOP = {"TAA", "TAG", "TGA"}


def attributes(text):
    out = defaultdict(list)
    for key, value in re.findall(r'(\w+) "([^"]*)";', text):
        out[key].append(value)
    return dict(out)


def bases(genome, positions, strand):
    seq = "".join(genome[p] for p in positions).upper()
    return seq.translate(COMP) if strand == "-" else seq


def ordered_positions(intervals, strand):
    result = []
    for start, end in sorted(intervals, reverse=(strand == "-")):
        result.extend(range(end-1, start-1, -1) if strand == "-" else range(start, end))
    return result


def merged(intervals):
    result = []
    for start, end in sorted(intervals):
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return tuple(result)


def assess_transcript(tx, genome):
    """GENCODE GTF CDS excludes stop_codon; join it explicitly for the chain."""
    tags = set(tx["attrs"].get("tag", []))
    if tags & {"cds_start_NF", "cds_end_NF"} or tx["selenocysteine"]:
        return None, "incomplete_or_selenocysteine"
    cds = sorted(tx["CDS"], reverse=(tx["strand"] == "-"))
    if not cds:
        return None, "no_CDS"
    cumulative = 0
    for start, end, phase in cds:
        if phase != (3-cumulative % 3) % 3:
            return None, "CDS_phase"
        cumulative += end-start
    cp = ordered_positions([(s,e) for s,e,_ in cds], tx["strand"])
    sp = ordered_positions(tx["stop_codon"], tx["strand"])
    ap = ordered_positions(tx["start_codon"], tx["strand"])
    if len(cp) != len(set(cp)) or len(sp) != 3 or len(set(sp)) != 3 or set(cp) & set(sp):
        return None, "CDS_stop_geometry"
    full = cp+sp
    if any(p < 0 or p >= len(genome) for p in full) or ap != cp[:3]:
        return None, "coordinates_or_start_feature"
    chain = merged([(s,e) for s,e,_ in cds] + tx["stop_codon"])
    if ordered_positions(chain, tx["strand"]) != full:
        return None, "stop_not_immediately_after_CDS_in_transcript"
    seq = bases(genome, full, tx["strand"])
    if len(seq) % 3 or set(seq)-set("ACGT") or seq[:3] != "ATG" or seq[-3:] not in STOP:
        return None, "noncanonical_or_ambiguous_ORF"
    if any(seq[i:i+3] in STOP for i in range(0,len(seq)-3,3)):
        return None, "internal_stop"
    return {"positions": full, "seq": seq, "chain": chain, "cds": cds}, "assessable"


def paired_site(tx, obj, genome):
    """Same genomic base: TAT/TAC -> synonymous Tyr or TAA, never different codons."""
    seq, pos = obj["seq"], obj["positions"]
    nc = len(seq)//3-1
    candidates = []
    for aa in range(math.ceil(0.2*nc), math.floor(0.8*nc)+1):
        j = 3*aa
        if seq[j:j+3] not in {"TAT", "TAC"}:
            continue
        q = pos[j+2]
        # The entire codon lies inside an internal CDS exon, >=30 nt from either edge.
        exon = next(((k,s,e) for k,(s,e,_) in enumerate(obj["cds"])
                     if all(s <= p < e for p in pos[j:j+3])), None)
        if exon is None:
            continue
        k, s, e = exon
        if k in {0, len(obj["cds"])-1} or min(min(pos[j:j+3])-s, e-1-max(pos[j:j+3])) < 30:
            continue
        syn = "C" if seq[j+2] == "T" else "T"
        ptc = "A"
        if tx["strand"] == "-":
            syn, ptc = syn.translate(COMP), ptc.translate(COMP)
        # Keep canonical GT/AG dinucleotides unchanged on the target gene strand.
        # This does not rule out regulatory changes or effects on the opposite strand.
        local = genome[q-1:q+2].upper()
        if tx["strand"] == "-":
            local = local.translate(COMP)[::-1]
        syn_local = syn.translate(COMP) if tx["strand"] == "-" else syn
        ptc_local = ptc.translate(COMP) if tx["strand"] == "-" else ptc
        motifs = {"GT", "AG"}
        old = tuple(local[i:i+2] if local[i:i+2] in motifs else None for i in (0,1))
        if any(tuple((local[:1]+b+local[2:])[i:i+2]
                     if (local[:1]+b+local[2:])[i:i+2] in motifs else None
                     for i in (0,1)) != old for b in (syn_local,ptc_local)):
            continue
        candidates.append((abs(aa/nc-0.5), aa, q, syn, ptc, min(q-s,e-1-q)))
    if not candidates:
        return None
    _, aa, q, syn, ptc, edge = min(candidates)
    return {"position_1based":q+1, "reference_base":genome[q].upper(),
            "synonymous_alt":syn, "stop_alt":ptc, "codon_index_1based":aa+1,
            "WT_codon":seq[3*aa:3*aa+3], "synonymous_codon":"TA"+("C" if seq[3*aa+2]=="T" else "T"),
            "stop_codon":"TAA", "distance_to_CDS_edge_nt":edge,
            "CDS_codon_fraction":aa/nc}


def choose_loci(rows, maximum=24):
    ordered = sorted(rows, key=lambda r:(r["position_1based"],r["gene_id"]))
    n = min(maximum,len(ordered))
    return [ordered[math.floor((i+0.5)*len(ordered)/n)] for i in range(n)] if n else []


def download(url, path):
    request = urllib.request.Request(url, headers={"User-Agent":"coding-rna-research/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response, path.open("xb") as out:
        shutil.copyfileobj(response,out)
        return {"url":url, "bytes":path.stat().st_size, "last_modified":response.headers.get("Last-Modified")}


def read_gtf(path, target):
    genes, txs, headers = {}, {}, []
    with gzip.open(path,"rt") as inp, target.open("x") as out:
        for line in inp:
            if line.startswith("#"):
                headers.append(line.rstrip())
                out.write(line)
                continue
            if not line.startswith("chr22\t"):
                continue
            out.write(line)
            f = line.rstrip().split("\t")
            if len(f) != 9:
                raise ValueError("GTF column count")
            attrs = attributes(f[8])
            gene = attrs["gene_id"][0]
            if f[2] == "gene":
                genes[gene] = attrs
            if "transcript_id" not in attrs:
                continue
            tid = attrs["transcript_id"][0]
            tx = txs.setdefault(tid, {"gene_id":gene, "strand":f[6], "attrs":attrs,
                                     "CDS":[], "stop_codon":[], "start_codon":[], "selenocysteine":False})
            if tx["gene_id"] != gene or tx["strand"] != f[6]:
                raise ValueError("Transcript identity/strand mismatch")
            if f[2] == "transcript":
                tx["attrs"] = attrs
            if f[2] == "CDS":
                tx["CDS"].append((int(f[3])-1,int(f[4]),int(f[7])))
            if f[2] in {"start_codon","stop_codon"}:
                tx[f[2]].append((int(f[3])-1,int(f[4])))
            if f[2] == "Selenocysteine":
                tx["selenocysteine"] = True
    return genes, txs, headers


def prepare(data, out):
    data.mkdir(parents=True, exist_ok=False)
    records = [download(FASTA_URL,data/"chr22.fa.gz"), download(GTF_URL,data/"gencode.v49.annotation.gtf.gz")]
    with gzip.open(data/"chr22.fa.gz","rt") as handle:
        if next(handle).strip() != ">chr22":
            raise ValueError("Unexpected FASTA sequence")
        genome = "".join(line.strip() for line in handle)
    if len(genome) != LENGTH or set(genome.upper())-set("ACGTN"):
        raise ValueError("Unexpected GRCh38 chr22 length/alphabet")
    # This pilot deliberately uses unmasked, sequence-only model weights. No repeat track.
    with (data/"chr22.unmasked.fa").open("x") as handle:
        handle.write(">chr22\n")
        for i in range(0,len(genome),60):
            handle.write(genome[i:i+60].upper()+"\n")
    genes, txs, headers = read_gtf(data/"gencode.v49.annotation.gtf.gz",data/"gencode.v49.chr22.gtf")
    groups, reasons, gene_reasons = defaultdict(list), Counter(), Counter()
    tx_ledger, candidates, gene_ledger = [], [], []
    for tid, tx in txs.items():
        if tx["attrs"].get("transcript_type") != ["protein_coding"]:
            continue
        obj, reason = assess_transcript(tx,genome)
        reasons[reason] += 1
        tx_ledger.append({"transcript_id":tid,"gene_id":tx["gene_id"],"status":reason})
        groups[tx["gene_id"]].append((tid,tx,obj))
    for gene in sorted(groups):
        items = groups[gene]
        valid = [(tid,tx,obj) for tid,tx,obj in items if obj is not None]
        chains = {obj["chain"] for _,_,obj in valid}
        if genes.get(gene,{}).get("gene_type") != ["protein_coding"]:
            reason = "not_protein_coding_gene"
        elif not valid:
            reason = "no_assessable_protein_coding_chain"
        elif len(chains) != 1:
            reason = "multiple_assessable_protein_coding_chains"
        elif len(valid) != len(items):
            reason = "additional_unassessable_protein_coding_transcript"
        else:
            tid, tx, obj = min(valid,key=lambda x:x[0])
            site = paired_site(tx,obj,genome)
            reason = "eligible" if site else "no_internal_matched_Tyr_site"
            if site:
                candidates.append({"gene_id":gene, "gene_name":tx["attrs"].get("gene_name",[""])[0],
                     "transcript_id":tid,"equivalent_transcripts":[x[0] for x in valid],
                     "strand":tx["strand"],"CDS_chain_0based_halfopen":obj["chain"],
                     "protein_length_aa":len(obj["seq"])//3-1, **site})
        gene_reasons[reason] += 1
        gene_ledger.append({"gene_id":gene,"status":reason,"protein_coding_transcripts":len(items),
                           "assessable_protein_coding_chains":len(chains)})
    chosen = choose_loci(candidates)
    for name,rows in (("transcript_eligibility.tsv",tx_ledger),("gene_eligibility.tsv",gene_ledger)):
        with (out/name).open("x") as handle:
            writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    result = {"stage":"reference_only_preparation_not_model_evaluation",
        "assembly":"GRCh38.p14 primary chromosome 22 / NC_000022.11 / UCSC hg38 chr22",
        "annotation":"GENCODE human release 49 comprehensive CHR, not latest release",
        "downloaded_utc":datetime.now(timezone.utc).isoformat(), "downloads":records,"GTF_headers":headers,
        "reference_length":len(genome),"original_lowercase_bases":sum(c.islower() for c in genome),
        "inference_mask_policy":"all uppercase, ANNEVO sequence-only / Tiberius no-softmask weight required",
        "gene_count_all_types":len(genes),"transcript_count_all_types":len(txs),
        "protein_coding_transcript_assessability":dict(reasons),"gene_eligibility":dict(gene_reasons),
        "eligible_loci":len(candidates),"selected_loci":len(chosen),
        "selection":"up to 24 equal-quantile loci in genomic order, reference only; no predictions read",
        "selected":chosen, "eligible":candidates,
        "no_combined_mutant_chromosome":True,"training_runs":0,"inference_runs":0,"Setaria_accessed":False,
        "scope_limits":["Tyr-to-stop same-base control is a restricted mutation class",
                        "single assessable protein-coding CDS chain is not all-isoform transcript uniqueness",
                        "unchanged canonical splice dinucleotides does not establish unchanged biological splicing",
                        "reference-selected chr22 developmental pilot, not independent clean zero-shot test"]}
    with (out/"preparation.json").open("x") as handle:
        json.dump(result,handle,indent=2)
    print(json.dumps({k:result[k] for k in ("stage","eligible_loci","selected_loci","gene_eligibility")},indent=2))


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--data-dir",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    prepare(args.data_dir,args.output_dir)
